#! /usr/bin/python

#########################################################
# Program:reversalFastaSeq.py
# Usage: Find the longest repeat segments from fasta.
# Author: michelia
# Version: 1.0
#       Editor: guoshuguang@scgene.com
#       Date: 2013-01-31
##########################################################

from __future__ import division
import argparse
from Bio import SeqIO

function = '''
Function: reverse_complement the sequence of fasta file.
'''


def main():
    parser = argparse.ArgumentParser(description=function)
    parser.add_argument("infasta", #type=str,
                        help="the input fasta file, the fastaFile can contain multi-sequence record.")
    parser.add_argument("outfile", #type=str,
                        help="the reversal fasta file out path. E.g: path to/xxx.fasta or path to/xxx.fa")
    args = parser.parse_args()
    reverseComplementFastaSeq(args.infasta, args.outfile)


def reverseComplementFastaSeq(infasta, outfile):
    seqRecordsReversal = []
    for seqRecord in SeqIO.parse(infasta, 'fasta'):
        seqRecord.id += '_reversal_complement'
        seqRecord.seq = seqRecord.seq.reverse_complement()
        seqRecordsReversal.append(seqRecord)
    SeqIO.write(seqRecordsReversal, outfile, 'fasta')

if __name__ == '__main__':
    main()
