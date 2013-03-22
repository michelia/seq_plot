#! /usr/bin/python

#########################################################
# Program:out_sam_is_paired.py
# Usage: Find the longest repeat segments from fasta.
# Author: Name
# Version: 1.0
#       Editor: name@scgene.com
#       Date: Friday, 14 December, 2012 11:20:15 AM CST
##########################################################

import argparse
import time

function = '''
This program's statistics the contig's N. e.g: N50, N90 ...
'''


def main():
    parser = argparse.ArgumentParser(version='%(prog)s 1.0',
                                    description=function,
                                    )
    parser.add_argument("infasta", #type=int
                        help="this is to statistical fasta file path")
    parser.add_argument("-N", type=int,
                        default=50,
                        help="If not give the N, default is N50. You can give it 60, 70, 90...")
    args = parser.parse_args()
    print
    start_time = time.clock()
    contig_N(args.infasta, args.N)
    print '[INFO] Elapsed time: %.4f' % (time.clock() - start_time)
    print

def contig_N(infasta, N=50):
    '''
    return {'longestContig:': (contigID, contigLength), 
            'contigN50': (contigID, contigLength)}
    '''
    infasta = open(infasta)
    #all of the contigs in the infasta, contigs is a list of tuble
    contigs = sorted_contigs(readfasta(infasta))
    print "The longest contig: lenth= %s; id= %s " % (contigs[0][2], contigs[0][0])
    # the sum of all bases
    base_sum = sum([lenth for id,contig_base,lenth in contigs])
    contigNstat = statistics_N(base_sum, contigs, N)  #contigNstat: (id, length)
    print "The contig's N%s: lenth= %d; id=%s" % (N, contigNstat[1], contigNstat[0])
    return {'longestContig': (contigs[0][0], contigs[0][2]), 'contigN50': (contigNstat[0], contigNstat[1])}

    infasta.close()


def statistics_N(base_sum, contigs, N=50):
    N = N / 100.0
    N_bases = 0
    for id,contig_base,lenth in contigs:
        N_bases += lenth
        if N_bases / float(base_sum) >= N:
            return id,lenth

def sorted_contigs(contigs):
    key_func = lambda contig:contig[2]
    contigs = sorted(contigs, key=key_func, reverse=True)
    return contigs


def readfasta(infasta):
    contigs = []
    contig_base = ''
    id = infasta.readline().strip()
    for line in infasta:
        if line.startswith('>'):
            contigs.append((id,contig_base,len(contig_base)))
            id = line[1:-1]
            contig_base = ''
        else:
            contig_base = ''.join([contig_base, line.strip()])
    contigs.append((id,contig_base,len(contig_base)))
    return contigs


if __name__ == '__main__':
    main()
