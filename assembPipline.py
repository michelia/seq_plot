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

function = '''
Function: Excute bwa and filter fastq.
'''


def main():
    parser = argparse.ArgumentParser(description=function)
    parser.add_argument("refSeq", #type=int
                        help="this is the infile")
    # parser.add_argument("fqDir", #type=int
    #                     help="this is the infile")
    parser.add_argument("leftReadFq1", #type=int
                        help="this is the infile")
    # parser.add_argument("assembPiplineRecordPath", #type=int
    parser.add_argument("rightReadFq1", #type=int
                        help="this is the infile")
    # parser.add_argument("assembPiplineRecordPath", #type=int
    parser.add_argument("leftReadFq2", #type=int
                        help="this is the infile")
    # parser.add_argument("assembPiplineRecordPath", #type=int
    parser.add_argument("rightReadFq2", #type=int
                        help="this is the infile")
    # parser.add_argument("assembPiplineRecordPath", #type=int
    # parser.add_argument("assembPiplineRecordPath", #type=int
    #                     help="this is the outfile")
    args = parser.parse_args()

# python ~/program/assembPipline.py /home/guoshuguang/node/assembPipline/refSeqFolder/pinus_massoniana_chloroplast_356999558.fasta  /home/guoshuguang/node/assembPipline/fqFolder/17_1.fq /home/guoshuguang/node/assembPipline/fqFolder/17_2.fq /home/guoshuguang/node/assembPipline/fqFolder/18_1.fq /home/guoshuguang/node/assembPipline/fqFolder/18_2.fq
# python ~/program/assembPipline.py /home/guoshuguang/node/assembPipline/refSeqFolder/Coniferopsida.fasta  /home/guoshuguang/node/assembPipline/fqFolder/17_1.fq /home/guoshuguang/node/assembPipline/fqFolder/17_2.fq /home/guoshuguang/node/assembPipline/fqFolder/18_1.fq /home/guoshuguang/node/assembPipline/fqFolder/18_2.fq


    assembPiplineRecordPath = os.path.join(os.path.dirname(os.path.abspath(args.leftReadFq1)), 'assembPiplineRecord.txt')
    assembPiplineRecord = open(assembPiplineRecordPath, 'w')

    bwa_exec(args.refSeq, args.leftReadFq1, args.rightReadFq1, args.leftReadFq2, args.rightReadFq2, assembPiplineRecord)

    # filter_fq(args.leftReadFq1, args.rightReadFq1, args.leftReadFq2, args.rightReadFq2, assembPiplineRecord)


    assembPiplineRecord.close()

def bwa_exec(refSeq, leftReadFq1, rightReadFq1, leftReadFq2, rightReadFq2, assembPiplineRecord):
# def bwa_exec(refSeq, fqDir, assembPiplineRecord):
    # bwaOut = outFolder + '/bwaOut'
    # if not os.path.exists(bwaOut):
    #     os.mkdir(bwaOut)
    # if os.path.exists(os.path.splitext(rightReadFq2)[0].split('_', 1)[0] + '.sam'):
    #     return
    bwaPath = '/scgene/elephant/pipeline/Align/bwa/current/bwa'
    bwaAlignParam = ' -n 8 -o 2 -q 15 -t 1 -I -L -e 10'
    bwaCmds = []
    bwaCmds.append((bwaPath + ' index ' + refSeq, refSeq+'.sa'))
    # fqFiles = [fqFile for fqFile in glob(fqDir+'/*') if '.fq' in os.path.splitext(fqFile)[1] or '.fastq' in os.path.splitext(fqFile)[1]]
    # for fqFile in fqFiles:
    #     bwaCmds.append(bwaPath + ' aln' + bwaAlignParam + ' ' + refSeq + ' ' + fqFile + ' > '  + os.path.splitext(fqFile)[0] + '.sai')
    bwaCmds.append((bwaPath + ' aln' + bwaAlignParam + ' ' + refSeq + ' ' + leftReadFq1 + ' > '  + os.path.splitext(leftReadFq1)[0] + '.sai', os.path.splitext(leftReadFq1)[0] + '.sai'))
    bwaCmds.append((bwaPath + ' aln' + bwaAlignParam + ' ' + refSeq + ' ' + rightReadFq1 + ' > '  + os.path.splitext(rightReadFq1)[0] + '.sai', os.path.splitext(rightReadFq1)[0] + '.sai'))
    bwaCmds.append((bwaPath + ' aln' + bwaAlignParam + ' ' + refSeq + ' ' + leftReadFq2 + ' > '  + os.path.splitext(leftReadFq2)[0] + '.sai', os.path.splitext(leftReadFq2)[0] + '.sai'))
    bwaCmds.append((bwaPath + ' aln' + bwaAlignParam + ' ' + refSeq + ' ' + rightReadFq2 + ' > '  + os.path.splitext(rightReadFq2)[0] + '.sai', os.path.splitext(rightReadFq2)[0] + '.sai'))
    bwaCmds.append((bwaPath + ' sampe' + ' ' + refSeq + ' ' + os.path.splitext(leftReadFq1)[0] + '.sai ' + os.path.splitext(rightReadFq1)[0] + '.sai ' + leftReadFq1 + ' ' + rightReadFq1 + ' > ' + os.path.splitext(rightReadFq1)[0].split('_', 1)[0] + '.sam', os.path.splitext(rightReadFq1)[0].split('_', 1)[0] + '.sam'))
    bwaCmds.append((bwaPath + ' sampe' + ' ' + refSeq + ' ' + os.path.splitext(leftReadFq2)[0] + '.sai ' + os.path.splitext(rightReadFq2)[0] + '.sai ' + leftReadFq2 + ' ' + rightReadFq2 + ' > ' + os.path.splitext(rightReadFq2)[0].split('_', 1)[0] + '.sam', os.path.splitext(rightReadFq2)[0].split('_', 1)[0] + '.sam'))
    for cmd, flagFile in bwaCmds:
        cmd_exec(cmd, flagFile, assembPiplineRecord)

def filter_fq(leftReadFq1, rightReadFq1, leftReadFq2, rightReadFq2, assembPiplineRecord):
    assembPiplineRecord.write('filter_fq' + '\n')
    if os.path.exists(os.path.splitext(rightReadFq2)[0]+'_filter'+os.path.splitext(rightReadFq2)[1]):
        assembPiplineRecord.write('    # Success 2!\n')
        return
    samFile1 = os.path.splitext(rightReadFq1)[0].split('_', 1)[0] + '.sam'
    samFile2 = os.path.splitext(rightReadFq2)[0].split('_', 1)[0] + '.sam'
    out_sam_is_paired(samFile1, leftReadFq1, rightReadFq1, os.path.splitext(leftReadFq1)[0]+'_filter'+os.path.splitext(leftReadFq1)[1], os.path.splitext(rightReadFq1)[0]+'_filter'+os.path.splitext(rightReadFq1)[1])
    out_sam_is_paired(samFile2, leftReadFq2, rightReadFq2, os.path.splitext(leftReadFq2)[0]+'_filter'+os.path.splitext(leftReadFq2)[1], os.path.splitext(rightReadFq2)[0]+'_filter'+os.path.splitext(rightReadFq2)[1])
    assembPiplineRecord.write('    # Success 1!\n')


def cmd_exec(cmd, flagFile, assembPiplineRecord):
    assembPiplineRecord.write(cmd + '\n')
    if os.path.exists(flagFile):
        assembPiplineRecord.write('    # Success 2!\n')
        return
    status, output = getstatusoutput(cmd)
    if status is not 0:
        print output
        error = ''.join(['@@@@@Error command: ', cmd])
        print error
        assembPiplineRecord.write('    # Error!\n')
        sys.exit(1)
    assembPiplineRecord.write('    # Success 1!\n')

if __name__ == '__main__':
    main()
