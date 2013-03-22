#! /usr/bin/python

#########################################################
# Program:out_sam_is_paired.py
# Usage: Find the longest repeat segments from fasta.
# Author: michelia
# Version: 1.0
#       Editor: guoshuguang@scgene.com
#       Date:  
##########################################################

from __future__ import division
import argparse
import time
from commands import getstatusoutput
import os, sys
from glob import glob
from capture_fastq_is_paired import out_sam_is_paired
from contig_N50 import contig_N

function = '''
Function: Excute SOAPdenovo and GapCloser.
'''


def main():
    parser = argparse.ArgumentParser(description=function)
    parser.add_argument("config", #type=int
                        help="SOAPdenovo configFile")
    parser.add_argument("kmerMin", type=int, default=35,
                        help="this is the infile")
    parser.add_argument("kmerMax", type=int, default=89,
                        help="this is the infile")
    parser.add_argument("output", #type=int
                        help="prefix of output file name")
    args = parser.parse_args()

# python ~/program/assembPipline.py /home/guoshuguang/node/assembPipline/refSeqFolder/pinus_massoniana_chloroplast_356999558.fasta  /home/guoshuguang/node/assembPipline/fqFolder/17_1.fq /home/guoshuguang/node/assembPipline/fqFolder/17_2.fq /home/guoshuguang/node/assembPipline/fqFolder/18_1.fq /home/guoshuguang/node/assembPipline/fqFolder/18_2.fq


    assembSOAPdenovoRecordPath = os.path.join(os.path.dirname(os.path.dirname(args.config)), 'assembSOAPdenovoRecord.txt')
    assembSOAPdenovoRecord = open(assembSOAPdenovoRecordPath, 'w')

    kmerContigValues = denovo_exec(args.config, args.kmerMin, args.kmerMax, args.output, assembSOAPdenovoRecord)
    assembSOAPdenovoRecord.write(repr(analysiskmerContigValues(kmerContigValues)))
    assembSOAPdenovoRecord.close()


def denovo_exec(config, kmerMin, kmerMax, output, assembSOAPdenovoRecord):
    SOAPdenovoPath = '/scgene/elephant/pipeline/denovo/SOAPdenovo-V1.05/bin/SOAPdenovo-127mer'
    GapCloserPath = '/scgene/elephant/pipeline/denovo/GapCloser/GapCloser'
    # denovoAssembParam = -s config -K 25  -o chloroplast_25
    kmerContigValues = {}  #contain each kmer corresponding the contigs's N50 and the longest contig
    for kmer in xrange(kmerMin, kmerMax+1, 2):
        denovoAssembCmd = SOAPdenovoPath + ' all' + ' -s ' + config + ' -K ' + str(kmer) + ' -o ' + output + '_K' + str(kmer)
        flagFile = output + '_K' + str(kmer) + '.scafSeq'
        cmd_exec(denovoAssembCmd, flagFile, assembSOAPdenovoRecord)
        kmerContigValues[kmer] = contig_N(flagFile) #statistics the contigs's N50 and the longest contig
        flagGapFile = output + '_K' + str(kmer) + '.fasta'
        GapCloserCmd = GapCloserPath + ' -b ' + config + ' -a ' + flagFile + ' -o ' + flagGapFile
        cmd_exec(GapCloserCmd, flagGapFile, assembSOAPdenovoRecord)


    for kmer in  kmerContigValues:  #record the kmerContigValues information in the assembSOAPdenovoRecord.txt
        assembSOAPdenovoRecord.write(str(kmer) + ':')        
        assembSOAPdenovoRecord.write(repr(kmerContigValues[kmer])+'\n')        
    return kmerContigValues

def analysiskmerContigValues(kmerContigValues):
    def maxN50(i):
        return kmerContigValues[i]['contigN50'][1]
    def maxLongContig(i):
        return kmerContigValues[i]['longestContig'][1]
    maxContigN50Kmer = max(kmerContigValues, key=maxN50)
    maxLongestContigKmer = max(kmerContigValues, key=maxLongContig)
    # possible Repeat assignment
    return {maxContigN50Kmer: kmerContigValues[maxContigN50Kmer], maxLongestContigKmer: kmerContigValues[maxLongestContigKmer]}

def cmd_exec(cmd, flagFile, shRecord):
    shRecord.write(cmd + '\n')
    if os.path.exists(flagFile):
        shRecord.write('    # Success 2!\n')
        return
    status, output = getstatusoutput(cmd)
    if status is not 0:
        print output
        error = ''.join(['@@@@@Error command: ', cmd])
        print error
        print 'status: ', status
        shRecord.write('    # Error!\n')
        sys.exit(1)
    shRecord.write('    # Success 1!\n')

if __name__ == '__main__':
    main()
