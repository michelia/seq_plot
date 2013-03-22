#encoding=utf8
from __future__ import division
import sys, re
sys.path.append('/scgene/tiger/invent/guoshuguang/program/')
from michelia import CsvReader, Bio
from Bio import SeqIO
import pdb
b = pdb.set_trace

def main(fa):
    '''
    所有基因的所有转录本的 fasta 序列
fa:  /scgene/tiger/invent/chenjiehu/GJXDGE/index/all.cdna.polyA.fa
    '''
    for seqRcord in SeqIO.parse(fa, 'fasta'):
        isoformID = seqRcord.id
        seq = str(seqRcord.seq)
        start = 0
        tagList = []
        while True:
            index = seq.find('CATG', start)
            if index > -1:
                tagList.append(seq[index:index+21])
                start = index + 1
            else:
                break

        if len(tagList) == 0:
            continue
        site = len(tagList)   # 3`端的第几个tag
        for tag in tagList:
            print '%s\t%s\t%s' % (isoformID, tag, str(site))
            site -= 1

if len(sys.argv) > 1:
    main(sys.argv[1])
else:
    print >>sys.stderr, "python %s <cdna>" % (sys.argv[0])
