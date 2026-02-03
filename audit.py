import argparse
import os
import statistics
import math

import numpy as np

from collections import namedtuple

# This is a complicated way of making a string constant.
SocialChoiceFunctions = namedtuple('SocialChoiceFunctions', ['FREE_LIST_HAMILTONIAN', 'SAINTE_LAGUE'])
social_choice_fns = SocialChoiceFunctions(FREE_LIST_HAMILTONIAN='free-list-hamiltonian', SAINTE_LAGUE='sainte-lague')

# Fraction of every party's votes assumed to be unreliable. This is approximately average for Karlsruhe 2024.
UNRELIABLE = 0.6

# Maximum votes per ballot.
# FIXME votes-per-ballot should be either read in from the data, or set as a command-line param.
# This is called m in the paper, which might be confusing because there are other uses of m in this code.
VMAX = 48

# Download from https://github.com/pbstark/SHANGRLA/blob/main/shangrla/core/NonnegMean.py
# Alternatively, comment this out and use sample_size_kaplan_kolmogorov rather than
# sample_size_kaplan_kolmogorov_shangrla_update.
from shangrla.core.NonnegMean import NonnegMean

# Read summarised election data from file.
# File has the following format:
# VOTERS,Number of voters
# INVALID,Number of invalidly cast ballots
# PARTY,VOTES,SEATS
# Party name,Votes in party's tally,Seats awarded
# Party name,Votes in party's tally,Seats awarded
# ...
# Party name,Votes in party's tally,Seats awarded
#
# The total number of valid ballots cast in the election
# is computed as Number of voters - Number of invalidly
# cast ballots.
def read_data(dfile):
    data = {}
    tot_ballots = 0
    tot_voters = 0
    with open(dfile, 'r') as f:
        lines = f.readlines()

        toks1 = lines[0].split(',')
        toks2 = lines[1].split(',')
        
        tot_ballots = int(toks1[1]) - int(toks2[1])
        tot_voters = int(toks1[1])

        # Skip first three lines: party votes & seats 
        # attributions start on the 4th line.
        for i in range(3, len(lines)):
            # line = Party,Votes,Seats
            toks = lines[i].split(',')

            data[toks[0]] = (int(toks[1]), int(toks[2]))

    return data, tot_ballots, tot_voters

# This function extracts code from audit_assertion_utils.py in (a prior version of) the
# SHANGRLA repository (https://github.com/pbstark/SHANGRLA).
def sample_size_kaplan_kolmogorov(margin, prng, N, error_rate, rlimit, t=1/2, \
    g=0.1, upper_bound=1, quantile=0.5, reps=20):

    clean = 1.0/(2 - margin/upper_bound)
    # FIXME think about whether this one_vote_over is correct for multi-vote settings, or whether we
    # should be multiplying by VMAX.
    # assorter value for the single overstatement discrepancy, i.e. when the CVR contains the
    # winner and the MVR is blank/zero.
    one_vote_over = (1-0.5)/(2-margin/upper_bound)
    # but in general (VT) a single vote may give VMAX votes to the winner.
    # See Eq. 8 of the BW audit paper, and consider whether we should be introducing 'error_rate'
    # single-overstatement errors or 'error_rate' VMAX-overstatement errors.
    # I'm not sure the above expression is right except when upper_bound = 1. Based on Eq 8 of
    # the BW paper, I think it should be 0.5 / (2*upper_bound - margin), which should be right
    # for all voting types, because it's what you get when the CVR is maximal (and hence exactly
    # cancels out the upper bound) and the MVR is neutral (and hence adds 0.5 in the numerator).

    samples = [0]*reps

    for i in range(reps):
        pop = clean*np.ones(N)
        inx = (prng.random(size=N) <= error_rate)  # randomly allocate errors
        pop[inx] = one_vote_over

        sample_total = 0
        mart = (pop[0]+g)/(t+g) if t > 0 else  1
        p = min(1.0/mart,1.0)
        j = 1

        while p > rlimit and j < N:
            mart *= (pop[j]+g)*(1-j/N)/(t+g - (1/N)*sample_total)
    
            if mart < 0:
                break
            else:
                sample_total += pop[j] + g

            
            p = min(1.0/mart,1.0) if mart > 1.0 else 1.0

            j += 1;

        if p <= rlimit:
            samples[i] = j
        else:
            return np.inf 

    return np.quantile(samples, quantile)

# 'estim' is the initital estimate of the population mean, which is optional and defaults to fixed_alternative_mean.
# eta is ???
# Inputs:
# margin: the assorter margin, i.e. twice the assorter mean minus 1.
# one_over: the value of the assorter when the CVR has maximum support for the winner and the MVR (manually-observed ballot)
# is blank. In the original free-list-hamiltonian code, this is simply one_vote_over = 0.5/(2 - margin/upper_bound), however this is
# not true for more general settings with an offset and a possibility of multiple votes per ballot.
# prng: TODO do we need this?
# N: the size of the population.
# error_rate: the error rate to be simulated.
# rlimit: the risk limit.
# t: the threshold to be tested (default 1/2). Frankly unclear whether variations will get the correct answer. Strongly recommended to leave it at 1/2.
# g: a technical parameter relevant to some risk functions.
# upper_bound: the assorter upper bound.
# quantile: the quantile measured in the simulation. For example, with the default (1/2) the simulation returns the sample size
#           taken by no more than half the simulations.
# reps: the number of simulation repeats.
# Generates a population sample of length N, defaulting to 1/(2-margin/upper_bound) - this is assumed to
# be the value of the comparison assorter when there are 0 discrepancies. (See for example the BW auditing paper, Eq 7.)
#
def sample_size_comparison_assorter(margin, one_over, prng, N, error_rate, rlimit, t=1 / 2, \
                                    g=0.1, upper_bound=1, quantile=0.5, reps=20):

    # Same generation of population as original free-list-hamiltonian code.
    clean = 1.0/(2 - margin/upper_bound)
    pop = clean * np.ones(N)
    # TODO we probably don't have to insert these at random places - we can just put error_rate*N errors at the front.
    inx = (prng.random(size=N) <= error_rate)  # randomly allocate errors
    pop[inx] = one_over

    rla_worker = NonnegMean(test=NonnegMean.kaplan_kolmogorov, u = upper_bound, N = N, t = t, g=g)
    sample_size = rla_worker.sample_size(x=pop, alpha=rlimit, reps=reps, quantile=quantile)
    return sample_size

def supermajority_sample_size(hquota, seats, tot_votes, tot_ballots, \
    tot_voters, erate, rlimit, t, g, REPS, seed, rfunc, v1):

    threshold = (hquota*seats)/tot_votes

    # supermajority assertion: party p1 achieved 
    # more than 'threshold' of the vote.
    share = 1.0/(2*threshold)
    
    amean = (1.0/tot_voters)*(v1*share) - \
        tot_votes/(2*tot_voters) + 0.5

    m = 2*amean - 1

    # Estimate sample size via simulation
    sample_size = np.inf
    if rfunc == "kaplan_kolmogorov":
        prng = np.random.RandomState(seed)
        # FIXME - allow the switch more elegantly.

        one_vote_over = (1 - 0.5) / (2 - margin / upper_bound)
        sample_size =  sample_size_comparison_assorter(amean, one_over=one_vote_over, prng=prng, N=tot_ballots, \
                           error_rate=erate, rlimit=rlimit, t=t, g=g, upper_bound=share, quantile=0.5, reps=REPS)
    else:
        print("Error: function " + rfunc + " not yet incorporated. Please use Kaplan-Kolmogorov (default).")

    return sample_size, m, threshold, amean

# Sainte Lague divisors: 1, 3, 5, 7, ...
def s_l_divisor(i):
    return 2*i - 1

# unreliable is the fraction of (all parties') votes that are assumed to be unreliable.
# In general, this is probably not a fixed fraction, but varies among parties - this is just a quick estimate.
def process_sainte_lague(data, tot_votes, tot_seats, tot_ballots, tot_voters, erate, rlimit, t, g, REPS, seed, rfunc, unreliable):
    print("Doing sainte lague")
    # Check that the seat allocation is correct
    # TODO

    # TODO - note at the moment these are both just rough approximations - we should read the exact
    # values off the .xlsx spreadsheet. They'll vary for different parties, and it won't necessarily be
    # exactly true that usum = VMAX * utot.
    # Total unreliable votes
    usum = int(tot_votes * unreliable)
    # Total unreliable ballots
    utot = int(tot_ballots * unreliable)

    max_sample_size = 0
    closest_winner = ''
    closest_loser = ''

    # For each (ordered) pair of parties, compute margin for pairwise Sainte-Lague assertions
    # Assert that A's lowest winner beat B's highest loser.
    # Find the overall max estimated sample size.
    for p_A in data:
        for p_B in data:

            v_A, seats_A = data[p_A]
            v_B, seats_B = data[p_B]

            # Nothing to check if the parties are the same, or A has no winners, or B has no losers.
            if p_A == p_B or seats_A == 0 or seats_B == VMAX:
                continue

            # Sainte-lague divisors
            # for A's lowest winner
            d_W_A = s_l_divisor(seats_A)
            # for B's highest loser
            d_L_B = s_l_divisor(seats_B+1)

            # Untrusted sub-tallies.
            U_A = int(v_A * unreliable)
            U_B = int(v_B * unreliable)

            # Trusted sub-tallies.
            R_A = v_A - U_A
            R_B = v_B - U_B

            # Definition of Delta from the BW paper.
            Delta = (R_A * d_L_B / d_W_A - R_B) / usum

            # Upper bound on ballot-polling assorter.
            # This is u_A,B in the paper (Eq 4).
            upper = VMAX * (d_L_B / d_W_A + 1) / (2 * (VMAX - Delta))

            # Ballot-polling assorter mean.
            # This is \bar(h^u_A,B) in the paper (Eq 5).
            amean = (U_A * d_L_B / d_W_A - U_B) / (2 * utot * (VMAX - Delta)) \
                    + Delta / (2 * (VMAX - Delta)) \
                    + 0.5

            # Ballot-polling assorter margin
            # This is \nu_u in the BW paper (Eq 6).
            margin = 2 * (amean) - 1

            # (Apparent) mean for the comparison assorter (i.e. assuming no discrepancies).
            # Setting c=b in Eq. 7 of the BW paper:
            mean_comparison = upper / (2*upper - margin)

            # Ballot-comparison assorter margin.
            margin_comparison = 2 * mean_comparison - 1

            # Upper bound for the comparison assorter, occurring when the MVR has VMAX for the winner
            # (so h(b) = upper) and the CVR has VMAX for the loser (so h(c) = 0).
            # Using Eq 7 of the BW paper:
            upper_comparison = 2 * upper / (2*upper - margin)

            # Value of the comparison assorter when the CVR has VMAX votes and the MVR is blank.
            # This is the equivalent of a one-vote overstatement in simple plurality. See Eq 8 in the BW paper.
            vmax_over = VMAX / 2 / (VMAX - Delta) / (2*upper - margin)

            # Estimate sample size via simulation
            if rfunc == "kaplan_kolmogorov":
                prng = np.random.RandomState(seed)
                sample_size = sample_size_comparison_assorter(margin=margin_comparison, one_over=vmax_over,  prng=prng, N=tot_ballots,
                                                              error_rate=erate, rlimit=rlimit, t=t, g=g, upper_bound=upper_comparison, quantile=0.5, reps=REPS)
                print("{} lowest winner {} vs {} highest loser {}: sample size {}".format(p_A, seats_A, p_B, seats_B+1, sample_size))
                if sample_size > max_sample_size:
                    max_sample_size = sample_size
                    closest_winner = p_A
                    closest_loser = p_B
            else:
                print(
                    "Error: function " + rfunc + " not yet incorporated. Please use Kaplan-Kolmogorov (default).")

    print("Max sample size: {} lowest winner vs {} highest loser: sample size {}\n\n".format(closest_winner, closest_loser, max_sample_size))

def process_hamiltonian(data, tot_votes, tot_seats, tot_ballots, tot_voters, erate, rlimit, t, g, REPS, seed, rfunc):
    level0_max_sample = 0
    level1_max_sample = 0

    # Check that party deserved all but their last seat
    hquota = tot_votes / tot_seats
    for p1 in data:
        v1, a1 = data[p1]

        if a1 > 1:
            sample_size, m, th, am = supermajority_sample_size(hquota, a1 - 1, tot_votes, tot_ballots, tot_voters, erate,
                                                               rlimit, t, g, REPS, seed, rfunc, v1)

            level0_max_sample = max(sample_size, level0_max_sample)

            print("Level 0,{},{},{},{},{}".format(p1, a1 - 1, th, m, sample_size))

        if a1 > 0:
            v_div_q = math.floor(tot_seats * (v1 / tot_votes))

            if v_div_q > 0:
                sample_size, m, th, am = supermajority_sample_size(hquota, v_div_q, tot_votes, tot_ballots, tot_voters,
                                                                   erate, rlimit, t, g, REPS, seed, rfunc, v1)

                print("Level 1,{},{},{},{},{}".format(p1, v_div_q, th, m, sample_size))

                level1_max_sample = max(sample_size, level1_max_sample)

    # for each pair of parties, in both directions, compute
    # margin for pairwise c-diff assertion
    level2_max_sample = 0
    num_level2 = 0
    for p1 in data:
        for p2 in data:
            if p1 == p2:
                continue

            v1, a1 = data[p1]
            v2, a2 = data[p2]

            # compute 'd'
            d = (a1 - a2 - 1) / tot_seats

            # Compute mean of assorter for assertion and margin 'm'
            # formerly before change to include invalid ballots
            amean = (1.0 / tot_voters) * (((v1 - v2) - tot_votes * d + TBTS * (1 + d)) / (2 * tot_seats * (1 + d)) + iballots / 2.0)

            m = 2 * (amean) - 1

            upper = 1 / (1 + d)

            # Estimate sample size via simulation
            sample_size = np.inf
            if args.rfunc == "kaplan_kolmogorov":
                prng = np.random.RandomState(seed)
                sample_size = sample_size_kaplan_kolmogorov(m, prng, tot_ballots, erate, rlimit, t=t, g=g,
                                                            upper_bound=upper, quantile=0.5, reps=REPS)
            else:
                print(
                    "Error: function " + args.rfunc + " not yet incorporated. Please use Kaplan-Kolmogorov (default).")

            level2_max_sample = max(sample_size, level2_max_sample)

            # Print out: Party name 1, Party name 2, proportion of votes
            # in Party 1's tally, proportion of votes in Party 2's tally,
            # value of 'd', margin, and estimate of initial sample required
            # to audit the assertion.
            print("Level 2,{},{},{},{},{},{},{}".format(p1, p2, v1 / tot_votes, v2 / tot_votes, d, m, sample_size))

            num_level2 += 1

    print("Level 0, Overal ASN: {} ballots".format(level0_max_sample))
    print("Level 1, Overal ASN: {} ballots".format(level1_max_sample))
    print("Level 2, Overal ASN: {} ballots, {} assertions".format(level2_max_sample, num_level2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    # Input: data file
    parser.add_argument('-d', action='store', dest='data')

    # Input: anticipated error rate (default value maps to 2 1 vote
    # overstatements per 1000 ballots)
    parser.add_argument('-e', action='store', dest='erate', default=0.002)

    # Input: risk limit (default is 5%)
    parser.add_argument('-r', action='store', dest='rlimit', default=0.05)
    
    # Input: parameter for some risk functions
    parser.add_argument('-g', action='store', dest='g', default=0.1)

    # Input: parameter for some risk functions
    parser.add_argument('-t', action='store', dest='t', default=1/2)
    
    # Risk function to use for estimating sample size
    parser.add_argument('-rf', action='store', dest='rfunc', \
        default="kaplan_kolmogorov")

    # Input: number of repetitions to perform in simulation to determine
    # an initial sample size estimation -- the quantile of the sample
    # size is computed (default is 20 repetitions)
    parser.add_argument('-reps', action='store', dest='reps', default=20)

    # Input: seed (default is 93686630803205229070)
    parser.add_argument('-s', action='store', dest='seed', \
        default=9368663)

    # Input: social choice function (default is free-list-hamiltonian; sainte-lague also allowed)
    parser.add_argument('-scf', action='store', dest='social_choice_fn', default=social_choice_fns.FREE_LIST_HAMILTONIAN)

    args = parser.parse_args()

    # data is a mapping between party name and a (votes, seats) tuple
    data, tot_ballots, tot_voters = read_data(args.data)

    # compute total number of votes and seats
    tot_votes = 0
    tot_seats = 0

    seed = int(args.seed)
    erate = float(args.erate)
    rlimit = float(args.rlimit)

    REPS = int(args.reps)
    t = float(args.t)
    g = float(args.g)

    social_choice_fn = args.social_choice_fn

    for p,(v,s) in data.items():
        tot_votes += v
        tot_seats += s

    print("{} seats, {} voters, {} parties, {} valid ballots, "\
        "{} total votes, social choice function {}".format(tot_seats, tot_voters, len(data), \
        tot_ballots, tot_votes, social_choice_fn))

    TBTS = tot_ballots*tot_seats
    iballots = tot_voters - tot_ballots

    if social_choice_fn == social_choice_fns.FREE_LIST_HAMILTONIAN:
        process_hamiltonian(data, tot_votes, tot_seats, tot_ballots, tot_voters, erate, rlimit, t, g, REPS, seed, args.rfunc)
    elif social_choice_fn == social_choice_fns.SAINTE_LAGUE:
        print("Processing Sainte-Lague, assuming all votes are unreliable.")
        process_sainte_lague(data, tot_votes, tot_seats, tot_ballots, tot_voters, erate, rlimit, t, g, REPS, seed, args.rfunc, 1)
        print("Processing Sainte-Lague, assuming {} votes are unreliable.".format(0.8))
        process_sainte_lague(data, tot_votes, tot_seats, tot_ballots, tot_voters, erate, rlimit, t, g, REPS, seed, args.rfunc, 0.8)
        print("Processing Sainte-Lague, assuming {} votes are unreliable.".format(UNRELIABLE))
        process_sainte_lague(data, tot_votes, tot_seats, tot_ballots, tot_voters, erate, rlimit, t, g, REPS, seed, args.rfunc, UNRELIABLE)
        print("Processing Sainte-Lague, assuming {} votes are unreliable.".format(0.4))
        process_sainte_lague(data, tot_votes, tot_seats, tot_ballots, tot_voters, erate, rlimit, t, g, REPS, seed, args.rfunc, 0.4)
    else:
        print("Error: Social choice function {} not supported.".format(social_choice_fn))