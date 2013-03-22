#encoding=utf8
#! /usr/bin/python

#########################################################
# Program:heteroSnpReadsSam.py
# Usage: 挑出包含SNP位点的sam格式的reads，此reads序列上的SNP位点的base与ref序列的碱基不一样。
# Author: michelia
# Version: 1.0
#       Editor: guoshuguang@scgene.com
#       Date:  
##########################################################

from __future__ import division
import argparse
import time
import sys
import re
from collections import namedtuple
sys.path.append('/scgene/tiger/invent/guoshuguang/program/')
sys.path.append('/scgene/tiger/invent/guoshuguang/program/moduls')
from perseSamCIGAR import read_simulate_ref, pysam_cigarstring_to_cigar
import pysam

function = '''
Function:  the sam reads contain heteroSNP.
'''

def main():
    parser = argparse.ArgumentParser(description=function,
                                    )
    parser.add_argument("snpFilePath", #type=int,
                        help="this is the infile")
    parser.add_argument("samFilePath", #type=int,
                        help="this is the outfile")
    parser.add_argument("outFilePath", #type=int,
                        help="this is the outfile")
    parser.add_argument("refName", #type=int,
                        help="this is the outfile")

    args = parser.parse_args()
    print
    start_time = time.clock()

    allReads = heter_snp_reads(args.snpFilePath, args.samFilePath, args.refName)
    writeReadToFile(allReads, arts.outFilePath)

    print args.refName
    print '[INFO] Elapsed time: %.4f' % (time.clock() - start_time)
    print 

def heter_snp_reads(snpFilePath, samFilePath, refName):
    samFile = pysam.Samfile(samFilePath, 'rb')
    patter = re.compile('\d+')
    heteroSnps = get_hetero_snp(snpFilePath)
    for heteroSnp in heteroSnps:
        for samAlign in sam_aligns(samFile, heteroSnp[0], refName):
            if samAlign.readSimulateRef[(heteroSnp[0] - samAlign.startPos)] == heteroSnp[1]:
                yield (samAlign, heteroSnp)

def writeReadToFile(allReads, outFilePath):
    outFile = open(outFilePath, 'w')
    for oneRead in allReads:
        samAlign = oneRead[0]
        heteroSnp = oneRead[1]
        outFile.write('%s\t%s\n' % (heteroSnp[0], samAlign.line))
    outFile.close()

def sam_aligns(samFile, pos, refName):
    SamLine =namedtuple('SamLine', 'qname flag rname pos mapq cigar rnext pnet tlen read qual')
    indexSam = namedtuple('indexSam', 'startPos endPos line readSimulateRef cigar')
    iter = samFile.fetch(refName, pos-1, pos)
    for align in iter:
        startPos = align.pos + 1
        readSimulateRef = read_simulate_ref(align.seq, align.cigarstring)
        endPos = startPos + len(readSimulateRef) - 1
        align_str = align_to_string(align)
        yield indexSam(*(startPos, endPos, align_str, readSimulateRef, align.cigarstring))

def align_to_string(align):
    align_str = '\t'.join(map(str,(align.qname, align.flag, align.rname, (align.pos+1), align.mapq, pysam_cigarstring_to_cigar(align.cigarstring), align.rnext, align.pnext, align.tlen, (align.seq), align.qual,tag_value_to_str(align.tags))))
    return align_str

def tag_value_to_str(tags):
    tagValueStr = ''
    for tag, value in tags:
        tagValueStr = ''.join((tagValueStr, tag, ':',str(value), '  '))
    if tagValueStr:
        tagValueStr = tagValueStr[:-1]
    return tagValueStr


def get_hetero_snp(snpFilePath):
    # heteroSnps = {}
    heteroSnps = []
    for line in (parse_SNP_file(snpFilePath)):
        if float(line[9]) != 1:
            heteroSnps.append((int(line[1]), line[7]))  #(pos, snpBase)
    return heteroSnps

def parse_SNP_file(snpFilePath):
    snpFile = open(snpFilePath)
    snpFile.readline()
    i = 0
    for line in snpFile:
        if line.strip():
            yield line.split()
    snpFile.close()
if __name__ == '__main__':
    main()
