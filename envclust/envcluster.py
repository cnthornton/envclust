#! /usr/bin/env python3
"""Cluster amplicon sequences using distribution patterns across samples to 
inform OTU generation.
"""

from __future__ import print_function, division

__author__ = 'Christopher Thornton'
__license__ = 'GPLv3'
__date__ = '2019-10-08'  #project initiation: 2014-07-23
__status__ = "Beta"
__version__ = '0.3.0'

import argparse
from math import log
import os
import random
from scipy.stats import *
import sys

def verbosity(verbose, string):
    if verbose:
        print(string)

def file_check(files, mode):
    """Check whether a file exists and/or is accessible to the user."""
    cwd = os.getcwd()
    path = cwd + '/' + files
    exists = os.path.exists(path)
    read_access = os.access(path, os.R_OK)
    access = True
    reason = None
    if mode == 'rU':
        if exists and not read_access:
            reason = "Error: you do not have the proper permissions to \
            access " + files
            access = False
        elif not exists:
            reason = "Error: " + files + " does not exist in " + cwd
            access = False
    elif mode == 'w':
        if not exists:
            try:
                fh = open(files, mode)
                fh.close()
            except IOError as e:
                reason = str(e)
                access = False
        else:
            access = False
            reason = "Error: " + files + " already exists in " + cwd + ". \
            Please move the file before continuing so that it will not be \
            overwritten" 
    return access, reason

def get_abundance(count_table):
    """Generate a list of two-tuples which contain the sequence identifier \
    and frequency that the sequence appears in the samples.
    """
    abunds = []
    abund_dict = {}
    with open(count_table, 'rU') as count:
        header = count.readline()
        for line in count:
            count_info = line.split()
            distribution = count_info[2:]
            abundance = count_info[1]
            identifier = count_info[0]
            abund_tuple = (identifier, abundance)
            abunds.append(abund_tuple)
            abund_dict[identifier] = distribution

    return abunds, abund_dict

def fasta_parse(fasta_file):
    """Return a dictionary of sequences with identifier as key and sequences \
    as values.
    """
    sequence_dict = {}
    with open(fasta_file, 'rU') as fasta:
        for line in fasta:
            if line.startswith('>'):
                identifier = line.strip('>\n')
                sequence_dict[identifier] = ''
            else:
                sequence_dict[identifier] += line.strip()

    return sequence_dict

def return_last(last):
    """Return the last item in an iterable object.
    """
    return last[-1]

def sort_by_last(tuple_list, direction):
    """Sort list of tuples by the value of the last item in tuple."""
    sorted_list = sorted(tuple_list, reverse=direction, key=return_last)
    return sorted_list
    
def distance_jc69(cand, rep):
    """Calculate the Jukes-Cantor corrected distance between a pair of 
    sequences.
    """
    num_diff = 0 #count for number of differences in bases between seqs

    for base in range(len(cand)): #each gap is penalized when calculating distance
        if cand[base] != rep[base]:
            num_diff += 1
            
    p = num_diff / len(cand)
    if p >= 3/4:
        distance = 1
    else:
        distance = -(3/4) * log(1 - (4/3) * p)

    return distance

def distance_k80(cand, rep):
    """Calculate the Kimura corrected distance between a pair of sequences.
    """ 
    ps = 0 #bases that show transitional changes
    qs = 0 #bases that show transversional changes

    pyrimidines = ['t', 'c']
    purines = ['a', 'g']
    for base in range(len(cand)):
        if cand[base] != rep[base]:
            if (cand[base] in pyrimidines and rep[base] in purines) or \
                (cand[base] in purines and rep[base] in pyrimidines):
                qs += 1
            else:
                ps += 1
                
    p = ps / len(cand)
    q = qs /len(cand)
    if q >= 1/2 or p >= (1-q)/2:
        distance = 1
    else:
        distance = -(1/2) * log(1 - 2 * p - q) - (1/4) * log(1 - 2 * q)

    return distance

def hamming_distance(cand, rep):
    """Calculate the Hamming distance between a pair of sequences.
    """
    #Return Hamming distance between equal-length sequences
    if len(cand) != len(rep):
        raise ValueError("Template and sequence index barcodes must of equal "
            "length")

    return sum(ch1 != ch2 for ch1, ch2 in zip(cand, rep))

def compare_distributions(cand, rep, iters):
    """Compare distributions to determine if they are significantly different.
    """
    obs_table = (cand, rep)

    # Remove if all cells in index contain a zero value
    zero_positions = []
    position_offset = 0
    for index in range(len(obs_table[0])):
        count = 0
        for dist in range(len(obs_table)):
            if int(obs_table[dist][index]) == 0:
                count += 1
        if count == len(obs_table):
            zero_positions.append(index - position_offset)
            position_offset += 1

    for position in zero_positions:
        for row in obs_table:
            rm_zero = row.pop(position)
 
    R, C, N, row_sums, col_sums = table_counts(obs_table)
    chisq_obs, exp_obs = chisq_test(obs_table)

    # Perform monte-carlo simulation to obtain P value if required for
    # accuracy. Otherwise, calculate p-value using chi-squared distribution
    simulate = check_cells(exp_obs)
    if simulate:
        p = simulate_chisq(obs_table, chisq_obs, iters)
    else:
        dof = (R - 1) * (C - 1)  #degrees of freedom
        cdf = chi2.cdf(chisq_obs, dof)  #cumulative distribution function
        p = 1 - cdf

    return p

def chisq_test(obs):
    """Perform a Pearson's chi-square test on a n x m contingency table."""
    R, C, N, row_sums, col_sums = table_counts(obs)

    # Compute expected values
    expected = []
    for r in row_sums:
        exp_row = [ r * c / N for c in col_sums]
        expected.append(exp_row)

    exp = tuple(expected)

    # Calculate the chi-squared test statistic
    chisq = 0
    for r in range(R):
        for c in range(C):
            chisq += (obs[r][c] - exp[r][c]) ** 2 / exp[r][c]

    return chisq, exp

def simulate_chisq(obs, chisq_obs, iters):
    """Perform a permutation chi-squared test."""
    R, C, N, row_sums, col_sums = table_counts(obs)

    chisq_sim_values = []
    row_sums.insert(0, 0)
    for iteration in range(iters):
        smpls = []
        unique_ids = []
        for c in range(C):
            identifier = 'SMPL' + str(c)
            unique_ids.append(identifier)
            composition = [identifier] * col_sums[c]
            smpls.extend(composition)

        random.shuffle(smpls)
        sim_obs = []
        for i in range(R):
            seq_dist = []
            for smpl in unique_ids:
                cell = smpls[row_sums[i]:(row_sums[i + 1] + \
                    row_sums[i])].count(smpl)
                seq_dist.append(cell)
            sim_obs.append(seq_dist)

        sim_obs = tuple(sim_obs)
        chisq_sim, exp_sim = chisq_test(sim_obs)

        if chisq_sim >= chisq_obs:
            chisq_sim_values.append(chisq_sim)

    p = len(chisq_sim_values) / iters
    return p

def table_counts(contingency_table):
    """Obtain summary information about contingency table created from \
    candidate and representative distributions.
    """
    num_rows = len(contingency_table)
    num_columns = len(contingency_table[0])

    row_sums = [sum(row) for row in contingency_table]

    column_sums = []
    for c in range(num_columns):
        col_total = 0
        for r in range(num_rows):
            col_total += contingency_table[r][c]
        column_sums.append(col_total) 

    table_total = sum(row_sums)

    return num_rows, num_columns, table_total, row_sums, column_sums

def check_cells(exp_table):
    """Verify that 80% of cells in the expected table have values of five or \
    more and that the smallest expected value is at least one. If not, \
    recommend Monte carlo simulation.
    """
    simulate = False
    count = 0
    num_cells = 0
    values = []
    for row in exp_table:
        num_cells += len(row)
        values.extend(row)
        for cell in row:
            if int(cell) < 5:
                count += 1

    percent_min_values = count / num_cells
    smallest = sorted(values)[0]

    if percent_min_values > 0.20 or smallest < 1:
        simulate = True

    return simulate

def main():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('fasta_file',
        metavar='in.fasta',
        help="input sequences in aligned FASTA format")
    parser.add_argument('count_file',
        metavar='in.count',
        help="input data table file of sequence counts. Data should be in the "
            "form species x sites, where rows represent sequence IDs and "
            "columns represent distinct samples")
    parser.add_argument('-o', '--out',
        default = sys.stdout,
        metavar='out.tsv',
        help="output tab-delimited OTU table to file [default: output to "
            "stdout]")
    parser.add_argument('-v', '--verbose', 
        action='count', 
        default=0,
        help="increase output verbosity")
    parser.add_argument('-d', '--distance',
        dest='dist_method',
        metavar='METHOD',
        choices=['jc69', 'k80', 'hamming'],
        default='hamming',
        help="model for calculating distance between sequences [default: "
            "hamming]. Options are Jukes-Cantor '69 (jc69), Kimura '80 (k80), "
            "or Hamming (hamming)")
    parser.add_argument('-i', '--iter',
        type=int,
        dest='num_iters',
        metavar='INT',
        default=10000,
        help="number of iterations to use when performing a Monte-carlo "
            "simulation if required for the chi-square test [default: 10000]")
    parser.add_argument('-m', '--max', 
        type=float, 
        metavar='DISTANCE',
        dest='max_dist', 
        default=0.10, 
        help="maximum genetic variation allowed between sequences in a "
            "population/cluster [default: 0.10]")
    parser.add_argument('-p', '--p', 
        type=float, 
        dest='p_cutoff', 
        metavar='THRESHOLD',
        default=0.05, 
        help="p-value cutoff for determining whether to reject the null "
            "hypothesis that the distribution of two sequences are "
            "statistically similar [default: 0.05]")
    parser.add_argument('--version', 
        action='version',
        version='%(prog)s ' + __version__)

    args = parser.parse_args()

    # Check program parameters
    in_files = [args.fasta_file, args.count_file]
    for in_file in in_files:
        in_access, in_reason = file_check(in_file, 'rU')
        if not in_access:
            print(in_reason)
            sys.exit(1)

    outfile = args.out
    out_access, out_reason = file_check(outfile, 'w')
    if not out_access:
        print(out_reason)
        sys.exit(1)
       
    max_dist = args.max_dist
    if max_dist <= 0 or max_dist > 1:
        print('--max-dist out of range. Value should be between 0 and 1', \
            file=sys.stderr)
        sys.exit(1)

    p_cutoff = args.p_cutoff
    if p_cutoff <= 0 or p_cutoff >= 1:
        print('--p-cutoff out of range. Value should be between 0 and 1', \
            file=sys.stderr)
        sys.exit(1)

    if args.dist_method == 'jc69':
        correction_method = distance_jc69
    elif args.dist_method == 'k80':
        correction_method = distance_k80
    else:
        correction_method = distance_hamming

    otu_dict = {}
    seq_dict = fasta_parse(args.fasta_file)
    abunds, abund_dict = get_abundance(args.count_file)
    sorted_abunds = sort_by_last(abunds, True)
    seqs_processed = 0
    num_seqs = len(sorted_abunds)
    num_otus = 0
    for candidate in sorted_abunds:
        otu_reps = otu_dict.keys()
        cand_id, cand_abund = candidate
        cand_seq = seq_dict[cand_id]

        message = 'Selecting {} as candidate sequence'.format(cand_id)
        verbosity(args.verbose, message)

        rep_distances = []
        for rep_id in otu_reps:
            rep_seq = seq_dict[rep_id]
            if len(cand_seq) != len(rep_seq):
                print('Error: sequence lengths are not the same', file=sys.stderr)
                sys.exit(1)

            distance = correction_method(cand_seq, rep_seq)
            if distance < max_dist:
                dist_tuple = (rep_id, distance)
                rep_distances.append(dist_tuple)
        
        # Compare distributions between candidate and OTU representative seqs
        sorted_reps = sort_by_last(rep_distances, False)
        added = False
        for rep in sorted_reps:
            rep_id = rep[0]
            cand_dist = [int(i) for i in abund_dict[cand_id]]
            rep_dist = [int(i) for i in abund_dict[rep_id]]
            obs_table = (cand_dist, rep_dist)
            p_value = compare_distributions(cand_dist, rep_dist, args.num_iters)

            # Add candidate to representative's OTU if distributions match
            if p_value >= p_cutoff:
                message = "Adding '{}' to cluster represented by '{}'"\
                    .format(cand_id, rep_id)
                verbosity(args.verbose, message)

                message = 'p-value: ' + str(p_value)
                verbosity(args.verbose, message)

                otu_dict[rep_id].append(cand_id)
                added = True
                break

        # Create OTU with candidate sequence as the representative sequence if
        # distribution is sufficiently unique
        if not added:
            num_otus += 1
            message = 'OTUs created: ' + str(num_otus)
            verbosity(args.verbose, message)

            otu_dict[cand_id] = [cand_id]
            message = "Creating cluster with '{}' as the representative"\
                .format(cand_id)
            verbosity(args.verbose, message)

        # Output current status
        seqs_processed += 1
        completion = int((seqs_processed / num_seqs) * 100)
        message = 'Percentage of sequence processed: {!s}%\n'.format(completion)
        verbosity(args.verbose, message)

    # Output information on run
    largest = max([len(otu_dict[i]) for i in otu_dict])
    smallest = min([len(otu_dict[j]) for j in otu_dict])
    print('Total OTUs: {!s}'.format(num_otus), file=sys.stderr)
    print('Sequences processed: {!s}'.format(seqs_processed), file=sys.stderr)
    print('Size of largest cluster: {!s}'.format(largest), file=sys.stderr)
    print('Size of smallest cluster: {!s}'.format(smallest), file=sys.stderr)
    
    # Write clusters to file
    message = "Finished processing sequences.\nWriting OTUs to '{}'."\
        .format(outfile)
    verbosity(args.verbose, message)

    fill = len(str(num_otus))
    otu_names = ['eOTU{}'.format(str(i).zfill(fill)) for i in \
        range(1, int(num_otus) + 1)]

    with open(outfile, 'w') as out_h:

        position = 0
        for cluster_rep in otu_dict:
            position += 1
            rep_name = otu_names[position]
            members = '\t'.join(otu_dict[cluster_rep])

            output = '{}\t{}\n'.format(rep_name, members)
            out.write(output)

if __name__ == '__main__':
    main()
    sys.exit(0)
