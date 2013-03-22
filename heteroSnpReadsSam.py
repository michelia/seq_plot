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
from perseSamCIGAR import read_simulate_ref
from Bio.Seq import reverse_complement

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

    # args = parser.parse_args(['/scgene/tiger/invent/guoshuguang/dealSam/RBQ_Chr1.snp', 
    #     '/scgene/tiger/invent/guoshuguang/dealSam/NP_Chr1.sam.bam.sort.sam', 
    #     '/scgene/tiger/invent/guoshuguang/dealSam/heteroSnpsSam.txt'])
    # args = parser.parse_args(['/scgene/tiger/invent/guoshuguang/dealSam/RBQ_Chr1.snp', 
    #     '/scgene/tiger/invent/guoshuguang/dealSam/now_test.sort.sam', 
    #     '/scgene/tiger/invent/guoshuguang/dealSam/heteroSnpsSam.txt'])

    # real
    # args = parser.parse_args(['/scgene/tiger/invent/chenjiehu/ECC/SNPANGENE/SNP/RBQ/RBQ_Chr1.snp', 
    #     '/scgene/tiger/invent/guoshuguang/dealSam/NP_Chr1.sam.bam.sort.sam', 
    #     '/scgene/tiger/invent/guoshuguang/dealSam/heteroSnpsSam.txt'])

    # heter_snp_reads(args.snpFilePath, args.samFilePath, args.outFilePath)
    shs = (
            # ('/scgene/tiger/invent/chenjiehu/ECC/SNPANGENE/SNP/RBQ/RBQ_Chr2.snp',
              # '/scgene/tiger/invent/guoshuguang/dealSam/NP_Chr2.sam.bam.sort.sam',
              # '/scgene/tiger/invent/guoshuguang/dealSam/Chr2_SNP_other_relate_to_ref.txt'),
            ('/scgene/tiger/invent/chenjiehu/ECC/SNPANGENE/SNP/RBQ/RBQ_Chr3.snp',
              '/scgene/tiger/invent/guoshuguang/dealSam/NP_Chr3.sam.bam.sort.sam',
              '/scgene/tiger/invent/guoshuguang/dealSam/Chr3_SNP_other_relate_to_ref.txt'),
            # ('/scgene/tiger/invent/chenjiehu/ECC/SNPANGENE/SNP/RBQ/RBQ_Chr4.snp',
            #   '/scgene/tiger/invent/guoshuguang/dealSam/NP_Chr4.sam.bam.sort.sam',
            #   '/scgene/tiger/invent/guoshuguang/dealSam/Chr4_SNP_other_relate_to_ref.txt'),
            # ('/scgene/tiger/invent/chenjiehu/ECC/SNPANGENE/SNP/RBQ/RBQ_Chr5',
            #   '/scgene/tiger/invent/guoshuguang/dealSam/NP_Chr5.bam.sort.sam',
            #   '/scgene/tiger/invent/guoshuguang/dealSam/Chr5_other_relate_to_ref.txt'),
            # ('/scgene/tiger/invent/chenjiehu/ECC/SNPANGENE/SNP/RBQ/RBQ_Chr6.snp',
            #   '/scgene/tiger/invent/guoshuguang/dealSam/NP_Chr6.sam.bam.sort.sam',
            #   '/scgene/tiger/invent/guoshuguang/dealSam/Chr6_SNP_other_relate_to_ref.txt'),
            # ('/scgene/tiger/invent/chenjiehu/ECC/SNPANGENE/SNP/RBQ/RBQ_Chr7.snp',
            #   '/scgene/tiger/invent/guoshuguang/dealSam/NP_Chr7.sam.bam.sort.sam',
            #   '/scgene/tiger/invent/guoshuguang/dealSam/Chr7_SNP_other_relate_to_ref.txt'),
        )
    for sh in shs:
        heter_snp_reads(*sh)


def heter_snp_reads(snpFilePath, samFilePath, outFilePath):
    outFile = open(outFilePath, 'w')
    patter = re.compile('\d+')
    heteroSnps = get_hetero_snp(snpFilePath)
    for samIndex in index_sam(samFilePath):
        for key in heteroSnps:
            # if samIndex.startPos <= heteroSnps[key][0] <= samIndex.endPos:
            if samIndex.startPos <= heteroSnps[key][0] <= samIndex[1] \
                and samIndex.readSimulateRef[(heteroSnps[key][0] - samIndex.startPos)] == heteroSnps[key][1]:
                # print samIndex.startPos,heteroSnps[key][0],\
                # samIndex.readSimulateRef[(heteroSnps[key][0] - samIndex.startPos )],samIndex.endPos,\
                # samIndex.cigar
                writeRead(samIndex, heteroSnps[key], outFile)
            # else:
                # print '***', len(samIndex.readSimulateRef), (heteroSnps[key][0] - samIndex.startPos )
                # print samIndex.startPos,heteroSnps[key][0],\
                # samIndex.readSimulateRef[(heteroSnps[key][0] - samIndex.startPos )],samIndex.endPos,\
                # samIndex.cigar

    outFile.close()

def writeRead(samIndex, heteroSnp, outFile):
    outFile.write('%s\t%s\n' % (heteroSnp[0], samIndex.line))

def index_sam(samFilePath):
    indexSam = namedtuple('indexSam', 'startPos endPos line readSimulateRef cigar')
    for line in parse_sam(samFilePath):
        # if line.tlen.startswith('-'):  #reverse
        #     readSimulateRef = read_simulate_ref(reverse_complement(line.read), line.cigar)
        #     startPos = endPos - len(readSimulateRef) + 1
        #     endPos = int(line.pos)
        #     print '@@@reverse',
        #     yield indexSam(*(startPos, endPos, '\t'.join(line), readSimulateRef, line.cigar))
        # else:
        startPos = int(line.pos)
        readSimulateRef = read_simulate_ref(line.read, line.cigar)
        endPos = startPos + len(readSimulateRef) - 1
        # print '-----', startPos,endPos
        yield indexSam(*(startPos, endPos, '\t'.join(line), readSimulateRef, line.cigar))


def parse_sam(samFilePath):
    samFile = open(samFilePath)
    SamLine =namedtuple('SamLine', 'qname flag rname pos mapq cigar rnext pnet tlen read qual')
    for line in samFile:
        # print line.split()
        if line.strip():
            # yield SamLine(*(line.split()[:12]))
            yield SamLine(*(line.split()[:11]))

    samFile.close()

def get_hetero_snp(snpFilePath):
    heteroSnps = {}
    for i, line in enumerate(parse_SNP_file(snpFilePath)):
        if float(line[9]) != 1:
            heteroSnps[i] = (int(line[1]), line[7])  #(pos, snpBase)
    # print len(heteroSnps), heteroSnps
    return heteroSnps

def parse_SNP_file(snpFilePath):
    snpFile = open(snpFilePath)
    snpFile.readline()
    for line in snpFile:
        if line.strip():
            # print line.split()
            yield line.split()
    snpFile.close()
if __name__ == '__main__':
    main()
