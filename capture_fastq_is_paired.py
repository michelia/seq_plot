#! /usr/bin/python

#########################################################
# Program:capture_fastq_is_paired.py
# Usage: Find the longest repeat segments from fasta.
# Author: Name
# Version: 1.0
#       Editor: name@scgene.com
#       Date: Friday, 14 December, 2012 11:20:15 AM CST
##########################################################

import argparse
import time
import copy

function = '''
According the sam file, capture the paired fastq.
The sam file comed from bwa.
E.g: capture_fastq_is_paired.py 17.sam 17_1.fq 17_2.fq 17_1_filter.fq
17_2_filter.fq
'''


def main():
    parser = argparse.ArgumentParser(
                                    description=function,
                                    )
    parser.add_argument("samfile",
                        help="")
    parser.add_argument("fqfile1",
                        help="")
    parser.add_argument("fqfile2",
                        help="")
    parser.add_argument("outpath1",
                        help="")
    parser.add_argument("outpath2",
                        help="")
    args = parser.parse_args()
    print
    start_time = time.clock()
    out_sam_is_paired(args.samfile, args.fqfile1, args.fqfile2, args.outpath1, args.outpath2)
    print 'Result sorted in %s file.' % args.outpath1
    print 'Result sorted in %s file.' % args.outpath2
    print '[INFO] Elapsed time: %.4f' % (time.clock() - start_time)
    print


def out_sam_is_paired(samfile, fqfile1, fqfile2, outpath1, outpath2):
    out1 = open(outpath1, 'w')
    out2 = open(outpath2, 'w')
    paired_fq_ids = {}
    for read1, read2 in parse_sam(samfile):
        filter_sam(read1, read2, paired_fq_ids)
    paired_fq_ids_bak = copy.copy(paired_fq_ids)
    for onefq in parsefq(fqfile1):
        if onefq['id'] in paired_fq_ids:
            out1.write(onefq['fq'])
            del paired_fq_ids[onefq['id']]
    for onefq in parsefq(fqfile2):
        if onefq['id'] in paired_fq_ids_bak:
            out2.write(onefq['fq'])
            del paired_fq_ids_bak[onefq['id']]
    out1.close()
    out2.close()



def parsefq(fqfile):
    fqfile = open(fqfile)
    fnext = fqfile.next
    for line in fqfile:
        id = line[1:].split('/')[0]
        fq = ''.join([line, fnext(), fnext(), fnext()])
        yield {'id': id, 'fq': fq}
    fqfile.close()


def filter_sam(read1, read2, paired_fq_ids):
    read1 = read1.split('\t')
    read2 = read2.split('\t')
    if read1[2] is not '*' or read2[2] is not '*':
        paired_fq_ids.setdefault(read1[0])


def parse_sam(infile):
    infile = open(infile)
    infile_next = infile.next
    while True:
        line = infile_next()
        if not line.startswith('@'):
            break
    yield line, infile_next()
    for line in infile:
        yield line, infile_next()
    infile.close()


if __name__ == '__main__':
    main()
