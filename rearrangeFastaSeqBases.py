#! /usr/bin/python

#########################################################
# Program:rearrangeFastaSeqBases.py
# Usage: Find the longest repeat segments from fasta.
# Author: michelia
# Version: 1.0
#       Editor: guoshuguang@scgene.com
#       Date:
##########################################################

from __future__ import division
import argparse
from Bio import SeqIO

function = '''
Function: rearrange the fasta file Sequence from the middle .
'''


def main():
    parser = argparse.ArgumentParser(description=function)
    parser.add_argument("infile", #type=str,
                        help="the input fasta file")
    parser.add_argument("outfile", #type=str,
                        help="the rearrange fasta file out path")
    parser.add_argument("--number", type=int, default=0,
                        help="insert the number of the N, default is 10.")
    parser.add_argument("position", type=float, 
                        help="the rearrange position. When 0 < position < 1, scale rearrangement. \
                        E.g: position=0.5, from middle to rearrangement.")
    args = parser.parse_args()
    rearrange_bases(args.infile, args.outfile, args.position, args.number)

def rearrange_bases(infile, outfile, position, number=0):
    aRecord = SeqIO.read(infile, 'fasta')
    if 0 < position < 1:
        splitPoint = int(len(aRecord.seq) * position)
    else:
        splitPoint = int(position)
    aRecord.id += str(splitPoint)
    aRecord.seq = aRecord.seq[splitPoint:] + 'N'*number + aRecord.seq[:splitPoint]
    SeqIO.write([aRecord], outfile, 'fasta')

if __name__ == '__main__':
    main()
