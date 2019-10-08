#! /usr/bin/env python3
"""Cluster amplicon sequences using distribution patterns across samples to 
inform OTU generation.
"""

from __future__ import print_function, division

__author__ = 'Christopher Thornton'
__date__ = '2019-10-08'  #project initiation: 2014-07-23
__version__ = '0.2.1'

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
    """\
    Generate a list of two-tuples which contain the sequence identifier and \
    frequency that the sequence appears in the samples.
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
    """\
    Return a dictionary of sequences with identifier as key and sequences as \
    values.
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
    return last[-1]

def sort_by_last(tuple_list, direction):
    """Sort list of tuples by the value of the last item in tuple."""
    sorted_list = sorted(tuple_list, reverse=direction, key=return_last)
    return sorted_list
    
def jc69_distance(cand, rep):
    """Calculate the Jukes-Cantor corrected distance between sequences."""
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

def k80_distance(cand, rep):
    """Calculate the Kimura corrected distance between sequences.""" 
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

def compare_distributions(cand, rep, iters):
    obs_table = (cand, rep)
    #remove if all cells in index contain a zero value
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
    #perform monte-carlo simulation to obtain P value if required for accuracy
    #otherwise, calculate p-value using chi-squared distribution
    simulate = check_cells(exp_obs)
    if simulate:
        p = simulate_chisq(obs_table, chisq_obs, iters)
    else:
        dof = (R - 1) * (C - 1) #degrees of freedom
        cdf = chi2.cdf(chisq_obs, dof) #the cumulative distribution function
        p = 1 - cdf
    return p

def chisq_test(obs):
    """Perform a Pearson's chi-square test on a n x m contingency table"""
    R, C, N, row_sums, col_sums = table_counts(obs)
    #compute expected values
    expected = []
    for r in row_sums:
        exp_row = [ r * c / N for c in col_sums]
        expected.append(exp_row)
    exp = tuple(expected)
    #calculate the chi-squared statistic
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
    """\
    Return information about the contigency table created from the candidate's \
    and represetative's distributions for use in other fuctions.
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
    """\
    Check if 80% of the cells in the expected table have values of five or \
    more and if the smallest expected value is at least one. If not, perform \
    monte carlo simulation.
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
        type=str,
        help="input sequences in aligned FASTA format")
    parser.add_argument('count_file',
        metavar='in.count',
        type=str,
        help="input data table file of sequence counts. Data should be in the "
            "form species x sites, where rows represent sequence IDs and "
            "columns represent distinct samples")
    parser.add_argument('-v', '--verbose', 
        action='count', 
        default=0,
        help="increase output verbosity")
    parser.add_argument('-c', '--correction',
        default='jc69',
        choices=['jc69', 'k80'],
        dest='dist_corr',
        help="distance model to use when calculating distance between "
            "sequences [default: jc69]")
    parser.add_argument('-i', '--iter',
        type=int,
        default=10000,
        dest='num_iters',
        help="number of iterations to use when performing a monte-carlo "
            "simulation if required for the chi-square test [default: 10000]")
    parser.add_argument('-m', '--max', 
        type=float, 
        default=0.10, 
        dest='max_dist', 
        help="maximum genetic variation for sequences to be allowed within "
            "the same population [default: 0.10]")
    parser.add_argument('-p', '--p-cutoff', 
        type=float, 
        default=0.05, 
        dest='p_cutoff', 
        help="p value cutoff for determining whether to reject the null "
            "hypothesis that the distribution of two sequences are "
            "statistically similar [default: 0.05]")
    parser.add_argument('-f', '--format', #add biom compatibility 
        type=str,
        default='txt',
        choices=['list', 'txt'],
        dest='out_format',
        help="output format. Options are a mothur-formatted list or classic, "
            "tab-delim otu table [default: txt]")
    parser.add_argument('--version', 
        action='version',
        version='%(prog)s ' + __version__)
    args = parser.parse_args()

    in_files = [args.fasta_file, args.count_file]
    for in_file in in_files:
        in_access, in_reason = file_check(in_file, 'rU')
        if not in_access:
            print(in_reason)
            sys.exit(1)
    if args.max_dist <= 0 or args.max_dist > 1:
        print('--max-dist out of range. Value should be between 0 and 1')
        sys.exit(1)
    elif args.p_cutoff <= 0 or args.p_cutoff >= 1:
        print('--p-cutoff out of range. Value should be between 0 and 1')
        sys.exit(1)

    otu_dict = {}
    seq_dict = fasta_parse(args.fasta_file)
    abunds, abund_dict = get_abundance(args.count_file)
    sorted_abunds = sort_by_last(abunds, True)
    seqs_processed = 0
    num_seqs = len(sorted_abunds)
    for candidate in sorted_abunds:
        otu_reps = otu_dict.keys()
        cand_id, cand_abund = candidate
        cand_seq = seq_dict[cand_id]
        output = 'Selecting ' + cand_id + ' as candidate sequence'
        verbosity(args.verbose, output)
        rep_distances = []
        for rep_id in otu_reps:
            rep_seq = seq_dict[rep_id]
            if len(cand_seq) != len(rep_seq):
                print('Error: sequence lengths are not the same')
                sys.exit(1)
            if args.dist_corr == 'jc69':
                distance = jc69_distance(cand_seq, rep_seq)
            elif args.dist_corr == 'k80':
                distance = k80_distance(cand_seq, rep_seq)
            if distance < args.max_dist:
                dist_tuple = (rep_id, distance)
                rep_distances.append(dist_tuple)
        
        #compare the distributions between candidate and OTU representative 
        #sequences
        sorted_reps = sort_by_last(rep_distances, False)
        added = False
        for rep in sorted_reps:
            rep_id = rep[0]
            cand_dist = [int(i) for i in abund_dict[cand_id]]
            rep_dist = [int(i) for i in abund_dict[rep_id]]
            obs_table = (cand_dist, rep_dist)
            p_value = compare_distributions(cand_dist, rep_dist, args.num_iters)
            #add the candidate to the rep's OTU if the distributions match
            if p_value >= args.p_cutoff:
                output = 'Adding ' + cand_id + ' to ' + rep_id + ' OTU'
                verbosity(args.verbose, output)
                output = 'P-value: ' + str(p_value)
                verbosity(args.verbose, output)
                otu_dict[rep_id].append(cand_id)
                added = True
                break

        #create OTU with candidate sequence as the representative sequence if
        #its distribution is not similar to any OTU rep's distribution
        if not added:  
            otu_dict[cand_id] = [cand_id]
            output = 'Creating OTU with ' + cand_id + ' as representative'
            verbosity(args.verbose, output)

        num_otus = len(otu_dict.keys())
        output = 'OTUs created: ' + str(num_otus)
        verbosity(args.verbose, output)
        seqs_processed += 1
        completion = int((seqs_processed / num_seqs) * 100)
        output = 'Percentage of sequence processed: ' + str(completion) + '%\n'
        verbosity(args.verbose, output)

    #write to log file
    log_file = str(args.fasta_file) + '.log'
    log_access, log_reason = file_check(log_file, 'w')
    if not log_access:
        print(log_reason)
        sys.exit(1)
    otus = 'Total OTUs: ' + str(num_otus)
    total_seqs = 'Number of sequences processed: ' + str(seqs_processed)
    largest = 'Size of largest cluster: '
    smallest = 'Size of smallest cluster: '
    with open(log_file, 'w') as log:
        log.write(total_seqs + '\n' + otus + '\n' + largest + '\n' + smallest)
    
    #write the OTUs to a file
    fasta_extensions = ('.fa', '.fna', '.fasta')
    for extension in fasta_extensions:
        if args.fasta_file.endswith(extension):
            in_name = '.'.join(args.fasta_file.split('.')[:-1])
            break
        else:
            in_name = args.fasta_file
    if args.out_format == "mothur_list":
        out_file = in_name + '.list'
    elif args.out_format == "otu_text":
        out_file = in_name + '_otu.txt'
    out_access, out_reason = file_check(out_file, 'w')
    if not out_access:
        print(out_reason)
        sys.exit(1)
    output = "Finished processing sequences\nWriting OTUs to " + out_file
    verbosity(args.verbose, output)
    fill = len(str(num_otus))
    otu_names = ['eOTU' + str(i).zfill(fill) for i in range(1, int(num_otus) \
        + 1)]
    with open(out_file, 'w') as out:
        if out_file.endswith(".txt"):
            position = 0
            for otu_rep in otu_dict.keys():
                output = otu_names[position] + '\t' + '\t'.join(otu_dict[otu_rep]) + '\n'
                out.write(output)
                position += 1
        elif out_file.endswith(".list"):
            output_list = ["dist", str(num_otus)]
            for otu_rep in otu_dict.keys():
                output_list.append(','.join(otu_dict[otu_rep]))
            output = '\t'.join(output_list)
            out.write(output)


if __name__ == '__main__':
    main()
    sys.exit(0)
