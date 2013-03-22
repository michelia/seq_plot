#! /usr/bin/python

#########################################################
# Program: filter_reads.py
# Usage:Filter the correct align reads in the SAM file.
# Author: Name
# Version: 1.0
#       Editor: name@scgene.com
#       Date: Monday, 17 December, 2012 03:29:19 PM CST
##########################################################


import sys
import time
import re

USAGE = '''
Function: Filter the correct pair reads in the SAM file.
USAGE:
    %s <in.sam> > output
''' % sys.argv[0].split('/')[-1]

def main(infile, outfile):
    out = open(outfile, 'w')
    pattern = re.compile('^\d+M$')
    for read1, read2 in parse_sam(infile):
        if filter_sam(read1, read2, pattern):
            out.write(read1)
            out.write(read2)
    out.close()

def filter_sam(read1, read2, pattern):
    read1 = read1.split('\t')
    read2 = read2.split('\t')
    if read1[6] != '=' or read2[6] != '=':
        return False
    if (not pattern.search(read1[5])) or (not pattern.search(read2[5])):
        return False
    if not Q_seq(read1[10]) or not Q_seq(read2[10]):
        return False
    return True


def Q_seq(Q_values):
    for Q in Q_values:
        if (ord(Q) - 33) < 10:
            return False
    return True


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
    if len(sys.argv) != 3:
        print USAGE
    else:
        start = time.clock()
        main(sys.argv[1], sys.argv[2])
        print 'Result sorted in %s file.' % sys.argv[2]
        print '[INFO] Elapsed time: %.4f' % (time.clock() - start)
