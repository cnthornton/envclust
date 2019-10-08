#! /usr/bin/env python3
"""Statistical functions for distribution-based clustering
"""

from math import log

def hamming_distance(cand, rep):
    """Calculate the Hamming distance between a pair of sequences.
    """
    #Return Hamming distance between equal-length sequences
    if len(cand) != len(rep):
        raise ValueError("Template and sequence index barcodes must of equal "
            "length")

    return sum(ch1 != ch2 for ch1, ch2 in zip(cand, rep))

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
