#! /usr/bin/python

#########################################################
# Program:comparison_2_file_lines.py
# Usage: Find the longest repeat segments from fasta.
# Author: Name
# Version: 1.0
#       Editor: name@scgene.com
#       Date: Friday, 14 December, 2012 11:20:15 AM CST
##########################################################

from __future__ import division
import argparse
import time

function = '''
Function: Two files's lines number comparison value
'''


def main():
    parser = argparse.ArgumentParser(version='%(prog)s 1.0',
                                    description=function,
                                    )
    parser.add_argument("numerator_file", #type=int
                        help="this is the numerator_file")
    parser.add_argument("denominator_file", #type=int
                        help="this is the denominator_file")
    args = parser.parse_args()
    print
    compare(args.numerator_file, args.denominator_file)
    print
def compare(numerator_file, denominator_file):
    value = statistics_lines(numerator_file) / statistics_lines(denominator_file) * 100
    print '[INFO]  Comparison value: %s / %s = %.4f%%' % (numerator_file, denominator_file, value)


def statistics_lines(infile):
    with open(infile) as f:
        for count, line in enumerate(f):
            pass
        return count + 1


if __name__ == '__main__':
    main()
