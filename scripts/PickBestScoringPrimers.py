#!/usr/bin/env python3

"""
Input from stdin the file ./CandidatePrimers.blastresults.OnTargetRatioAdded.tab, and some required positional parameters

Each line in stdin input expect to be in 9 column tab delimited format:
1. chr
2. start
3. stop
4. name [PrimerName:PrimerLen:PrimerPos:PrimerTm:PrimerSeq]
5. score
6. strand
7. OffTargetCountPrediction [#]
8. OnTargetCountPrediction [#]
N. Optional extra N columns

Only the columns described in [detail] are utilized in scoring

Output to stdout the same as input with additional columns that correspond to primer scores

1. chr
2. start
3. stop
4. name
5. score
6. strand
7. OffTargetCountPrediction
8. OnTargetCountPrediction
N. Optional extra N columns
N+1. GCContentScore
N+2. TmScore
N+3. PercentOnTargetEstimate_Score
N+4. Final Score (sum(10,11,12))
N+5. TargetName

Usage:

cat ./CandidatePrimers.blastresults.OnTargetRatioAdded.tab | python PickBestScoringPrimers.py [MinPrimerLength] [MaxPrimerLength] [MinTm] [MaxTm] [MinPercentOnTargetEstimate] [Tm_weight] [PercentOnTargetEstimate_weight] [GCcontent_weight]
"""

#cat ./CandidatePrimers.blastresults.OnTargetRatioAdded.tab | python PickBestScoringPrimers.py 24 32 55 60 66 50 0.3 0.5 0.2 | head



import sys
import re
from signal import signal, SIGPIPE, SIG_DFL

signal(SIGPIPE, SIG_DFL)

def count_kmers(read, k):
    """Count kmer occurrences in a given read.

    Parameters
    ----------
    read : string
        A single DNA sequence.
    k : int
        The value of k for which to count kmers.

    Returns
    -------
    counts : dictionary, {'string': int}
        A dictionary of counts keyed by their individual kmers (strings
        of length k).

    Examples
    --------
    >>> count_kmers("GATGAT", 3)
    {'ATG': 1, 'GAT': 2, 'TGA': 1}
    """
    # Start with an empty dictionary
    counts = {}
    # Calculate how many kmers of length k there are
    num_kmers = len(read) - k + 1
    # Loop over the kmer start positions
    for i in range(num_kmers):
        # Slice the string to get the kmer
        kmer = read[i:i+k]
        # Add the kmer to the dictionary if it's not there
        if kmer not in counts:
            counts[kmer] = 0
        # Increment the count for this kmer
        counts[kmer] += 1
    # Return the final counts
    return counts

def OligoComplexity(OligoSequence):
    """
    Returns complexity (C) for an oligo of length 18 to 67 as described on wiki:
    https://en.wikipedia.org/wiki/Linguistic_sequence_complexity
    """
    U1 = len(set(OligoSequence))/float(4)
    U2 = len(count_kmers(OligoSequence, 2)) / float(len(OligoSequence)-2+1)
    U3 = len(count_kmers(OligoSequence, 3)) / float(len(OligoSequence)-3+1)
    U4 = len(count_kmers(OligoSequence, 4)) / float(len(OligoSequence)-4+1)
    return(U1*U2*U3*U4)


def PieceWiseOptimalityFunction(Min, Optimal, Max, X_in):
    """
    A piecewise function of X_in, with parameters of Min, Max, and Optimal, that will return a number between 0 and 1.
    A linear increase from 0 to 1 across the X_in range of (Min, Optimal].
    A linear decrease from 1 to 0 across the X_in range of (Optimal, Max].
    Y=0 All else.

    Used to make a GCContent score, for example, where the optimal GCContent would be 50, with a Min of 25, and a Max of 75self.

    Example:

    >>> PieceWiseOptimalityFunction(25,50,75,25)
    0

    >>> PieceWiseOptimalityFunction(25,50,75,50)
    1.0

    >>> PieceWiseOptimalityFunction(25,50,75,55)
    0.8

    """
    if Min < X_in <= Optimal:
        Y = 1/(Optimal-Min) * (X_in - Min)
    elif Optimal < X_in <= Max:
        Y = 1/(Optimal-Max) * (X_in - Max)
    else:
        Y=0
    return(Y)

# MinTm, OptimalTm, MaxTm, MinPercentOnTargetEstimate, Tm_weight, PercentOnTargetEstimate_weight, GCcontent_weight = sys.argv[1:]
MinTm, OptimalTm, MaxTm, Tm_weight, MinPercentOnTargetEstimate, PercentOnTargetEstimate_weight, GCcontent_weight, MinComplexity, Complexity_weight = sys.argv[1:]

for line in sys.stdin:
    chr, start, stop, name, score, strand, OffTargetCountPrediction, OnTargetCountPrediction, *Other = line.strip('\n').split('\t')
    TargetName,PrimerLen,PrimerPos,PrimerTm,PrimerSeq = name.split(':')
    GCContent = float((PrimerSeq.count('G') + PrimerSeq.count('C'))) / len(PrimerSeq) * 100
    GCScore = PieceWiseOptimalityFunction(20,50,80, GCContent) * float(GCcontent_weight)
    TmScore = PieceWiseOptimalityFunction(*[float(Tm) for Tm in [MinTm,OptimalTm,MaxTm,PrimerTm]]) * float(Tm_weight)
    ComplexityScore = OligoComplexity(PrimerSeq) * float(Complexity_weight)
    try:
        PercentOnTargetEstimate = float(OnTargetCountPrediction)/(float(OnTargetCountPrediction) + float(OffTargetCountPrediction))
    except ZeroDivisionError:
        PercentOnTargetEstimate = 1
    PercentOnTargetEstimateScore = PercentOnTargetEstimate * float(PercentOnTargetEstimate_weight)
    if (PercentOnTargetEstimate > float(MinPercentOnTargetEstimate)/100) and (OligoComplexity(PrimerSeq) > float(MinComplexity)):
        print(line.strip('\n') + '\t' + str(PercentOnTargetEstimate) + '\t' + str(GCScore) + '\t' + str(TmScore) + '\t' + str(PercentOnTargetEstimateScore) + '\t' + str(ComplexityScore) + '\t' + str(GCScore + TmScore + PercentOnTargetEstimateScore + ComplexityScore) + '\t' + TargetName)
