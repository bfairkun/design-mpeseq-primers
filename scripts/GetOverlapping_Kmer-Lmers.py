#!/usr/bin/env python3

"""
Input from std in the results from (bedtools getfasta -name+ -tab), and a k and l parameter

Output to std out each (k)mer, (k+1)mer, ... , (l)mer for each sequence in the input

Usage:
bedtools getfasta -bed <mybed> -fi <myfasta> -name+ -tab | python GetOverlapping_Kmer-Lmers.py <k> <l>
"""

import sys
import re
from signal import signal, SIGPIPE, SIG_DFL

signal(SIGPIPE, SIG_DFL)

k=int(sys.argv[1])
l=int(sys.argv[2])

for line in sys.stdin:
    name,seq = line.strip('\n').split('\t')
    BedName, BedChr, BedStart, BedStop, BedStrand = re.search('(^.+)::(.+?):(\d+)-(\d+)\(([+-])\)$', name).groups()
    for primerlength in range(k,l+1):
        for i in range(0,len(seq)-primerlength):
            primerseq = seq[i:i+primerlength]
            PrimerPosition = '_'.join([BedChr, str(int(BedStart)+i), str(int(BedStart)+i+primerlength),BedStrand])
            #TargetName::primerlength::PrimerPosition\tprimerseq
            print(BedName.replace(" ", "_") + "::" + str(primerlength) + "::" + PrimerPosition + "::"+ primerseq +'\t' + primerseq)
