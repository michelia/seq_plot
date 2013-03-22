#! /usr/bin/python

#########################################################
# Program:out_sam_is_paired.py
# Usage: Find the longest repeat segments from fasta.
# Author: Name
# Version: 1.0
#       Editor: name@scgene.com
#       Date: Tue Dec 25 14:45:52 CST 2012
##########################################################

from __future__ import division
import argparse
import time

function = '''
Function: From a fasta file, find the longest sequense and then save to another fasta file.
'''


def main():
    parser = argparse.ArgumentParser(version='%(prog)s 1.0',
                                    description=function,
                                    )
    parser.add_argument("infasta", #type=int
                        help="fasta file")
    parser.add_argument("outfastapath", #type=int
                        help="the outpath to saved the longest seq, format is fasta")
    args = parser.parse_args()
    print
    start_time = time.clock()
    longest_seq_from_fasta(args.infasta, args.outfastapath)
    print 'The longest sequense saved in %s file.' % args.outfastapath
    print '[INFO] Elapsed time: %.4f' % (time.clock() - start_time)
    print

def longest_seq_from_fasta(infasta, outpath):
    infasta = open(infasta)
    out = open(outpath, 'w')
    long_id = ''
    long_seq = ''
    long_lenth = 0
    for id, seq, lenth in readfasta(infasta):
        if long_lenth < lenth:
            long_lenth = lenth
            long_seq = seq
            long_id = id
    write_file(long_id, long_seq, long_lenth, out)
    out.close()
    infasta.close()


def write_file(long_id, long_seq, long_lenth, out):
    i = 0
    out.write('%s\n' % long_id)
    while True:
        if (i + 100) >= long_lenth:
            out.write(long_seq[i:])
            break
        j = i + 100
        out.write('%s\n' % long_seq[i:j])
        i = j


def readfasta(infasta):
    contig_base = ''
    id = infasta.readline().strip()
    for line in infasta:
        if line.startswith('>'):
            yield (id,contig_base,len(contig_base))
            # id = line[1:-1]
            id = line.rstrip()[1:]
            contig_base = ''
        else:
            contig_base = ''.join([contig_base, line.strip()])
    yield (id,contig_base,len(contig_base))


if __name__ == '__main__':
    main()
