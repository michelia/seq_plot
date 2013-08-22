#!/usr/bin/env python
#encoding=utf8

#########################################################
# Program: samflag.py
# Usage: Parse the flag in the samfile.
# Example: python samflag.py 163
# Flag 163:
#     1 template having multiple segments in sequencing
#     2 each segment properly aligned according to the aligner
#     32 SEQ of the next segment in the template being reversed
#     128 the last segment in the template
# Author: Name
# Version: 1.0
#       Editor: name@scgene.com
#       Modify: Add Usage
#       Date: Tuesday, 18 December, 2012 05:30:59 PM CST
##########################################################


import sys
import time

USAGE = '''Function: samflag.py
USAGE:
    %s <flag_number>
The flag_number is not greater than 2047.
1 <= flag_number <= 2047
# Example: $ python samflag.py 163
 Flag 163:
 1 template having multiple segments in sequencing
 2 each segment properly aligned according to the aligner
 32 SEQ of the next segment in the template being reversed
 128 the last segment in the template
''' % sys.argv[0].split('/')[-1]

def parse_samflag(num):
    print
    num = int(num)
    if num > 2047:
        print '$python samflag.py <flag_number>\nThe flag_number is not greater than 2047.'
        return
    print 'Flag %s:' % num
    num = bin(num)[2:][::-1]
    for i, flag in enumerate(num):
        if flag == '1':
            print_flag(i)

def print_flag(i):
    flag_id = hex(2**i)
    print '%4s\t%s' % (2**i, flag[flag_id])

flag ={'0x1': 'template having multiple segments in sequencing',
        '0x2': 'each segment properly aligned according to the aligner',
        '0x4': 'segment unmapped',
        '0x8': 'next segment in the template unmapped',
        '0x10': 'SEQ being reverse complemented',
        '0x20': 'SEQ of the next segment in the template being reversed',
        '0x40': 'the first segment in the template',
        '0x80': 'the last segment in the template',
        '0x100': 'secondary alignment',
        '0x200': 'not passing quality controls',
        '0x400': 'PCR or optical duplicate'}
# 下面是flag的另一种解释， 来源： http://picard.sourceforge.net/explain-flags.html 可以查看页面的源代码，下面的行与上面是一一对应的
        # read paired
        # read mapped in proper pair
        # read unmapped
        # mate unmapped
        # read reverse strand
        # mate reverse strand
        # first in pair
        # second in pair
        # not primary alignment
        # read fails platform/vendor quality checks
        # read is PCR or optical duplicate


if __name__ == '__main__':
    print
    if len(sys.argv) != 2:
        print USAGE
    else:
        start = time.clock()
        parse_samflag(sys.argv[1])
        print
