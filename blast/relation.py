from __future__ import division
import numpy as np
import argparse
import sys
sys.setrecursionlimit(1000000)
import os
from collections import namedtuple
import matplotlib
# matplotlib.use('Agg')  # solve run on the cluster's mistakes.
import matplotlib.pyplot as plt
from itertools import count
import time
import copy
sys.path.append('/scgene/tiger/invent/guoshuguang/program/moduls')
from Bio import SeqIO
sys.path.append('/scgene/tiger/invent/guoshuguang/program/')
from michelia import vzebra, hzebra, data_dump, data_load, \
      adjust_spines, adjust_ticks_label_size, csvreader, csvwriter,\
      parse_blastout
import path

saveFilterBlastFilePath = '/scgene/tiger/invent/guoshuguang/align/blast/relative_network/head_tail_filter.blast'
outFilePath = '/scgene/tiger/invent/guoshuguang/align/blast/relative_network/relation_result.txt'
isContainOutPath = '/scgene/tiger/invent/guoshuguang/align/blast/relative_network/contain_seq.txt'
def main():
    faFilePath = '/scgene/tiger/project/130205Nb/gene/blast/index/assemble.fa'
    blastFilePath = '/scgene/tiger/project/130205Nb/gene/blast/assemble.fa.blast.filter'
    relations = seq_relation_blastout(blastFilePath, faFilePath)
    with open(outFilePath, 'w') as outFile:
        writer = csvwriter(outFile)
        for relation in relations:
            writer.writerow(relation)

def seq_relation_blastout(blastFilePath, faFilePath):
    filterblastRecorders = {}
    for blastRecorder in file_blastout(blastFilePath, faFilePath):
        filterblastRecorders[blastRecorder] = ''
    relations = []
    seqIDs = seq_ID(faFilePath)
    for seqID in seqIDs:
        relation = one_relation(seqID, filterblastRecorders)
        if len(relation) > 1:
            relations.append(relation)
    return sorted(relations, key=lambda relation: seqID_to_num(relation[0]))

def one_relation(seqID, filterblastRecorders):
    relation = {seqID:''}
    recur_relation(relation, filterblastRecorders)
    return sorted(relation.keys(), key=lambda x: seqID_to_num(x))

def recur_relation(relation, filterblastRecorders):
    '''
    recursion to find the relation
    Warning: Here the filterblastRecorders is dictionary, 
    so it's value can modify, 
    e.g: del filterblastRecorders[blastRecorder]
    '''
    for blastRecorder in filterblastRecorders:
        if blastRecorder.queryId in relation:
            relation[blastRecorder.subjectId] = ''
            del filterblastRecorders[blastRecorder]
            recur_relation(relation, filterblastRecorders)
            break
        if blastRecorder.subjectId in relation:
            relation[blastRecorder.queryId] = ''
            del filterblastRecorders[blastRecorder]
            recur_relation(relation, filterblastRecorders)
            break

def file_blastout(blastFilePath, faFilePath):
    seqLengths = seq_length(faFilePath)
    isContainOut = open(isContainOutPath, 'w')
    isContainOutWriter = csvwriter(isContainOut)    
    with open(saveFilterBlastFilePath, 'w') as saveFile:
        writer = csvwriter(saveFile)
        for blastRecorder in parse_blastout(blastFilePath):
            blastRecorderLength = blastRecorder_length(blastRecorder, seqLengths)
            isContain = is_contain(*blastRecorderLength)
            if isContain and blastRecorder.queryId != blastRecorder.subjectId:
                isContainOutWriter.writerow(blastRecorder)
            if is_head_tail(*blastRecorderLength) and not isContain:
                writer.writerow(blastRecorder)
                yield blastRecorder
    isContainOut.close()

def is_head_tail(length, rlength, queryStart, queryEnd, subjectStart, subjectEnd):
    if (queryStart == 1 and subjectEnd == rlength) or \
        (subjectStart == 1 and queryEnd == length):
        return True

def is_contain(length, rlength, queryStart, queryEnd, subjectStart, subjectEnd):
    if (queryStart == 1 and queryEnd == length) or \
        (subjectStart == 1 and subjectEnd == rlength):
        return True

def blastRecorder_length(blastRecorder, seqLengths):
    length = seqLengths[blastRecorder.queryId]
    rlength = seqLengths[blastRecorder.subjectId]
    if blastRecorder.queryStart > blastRecorder.queryEnd:
        queryStart = length + 1 - blastRecorder.queryStart
        queryEnd = length + 1 - blastRecorder.queryEnd
    else:
        queryStart = blastRecorder.queryStart
        queryEnd = blastRecorder.queryEnd
    if blastRecorder.subjectStart > blastRecorder.subjectEnd:
        subjectStart = rlength + 1 - blastRecorder.subjectStart
        subjectEnd = rlength + 1 - blastRecorder.subjectEnd
    else:
        subjectStart = blastRecorder.subjectStart
        subjectEnd = blastRecorder.subjectEnd
    return (length, rlength, queryStart, queryEnd, subjectStart, subjectEnd)

def seq_length(faFilePath):
    seqLengths = {}
    for seqRecorder in SeqIO.parse(faFilePath, 'fasta'):
        seqLengths[seqRecorder.id] = len(seqRecorder)
    # print 'seqLengths', len(seqLengths)
    return seqLengths

def seq_ID(faFilePath):
    seqIDs = {}
    for seqRecorder in SeqIO.parse(faFilePath, 'fasta'):
        seqIDs[seqRecorder.id] = ''
    # print 'seqIDs', len(seqIDs)
    return seqIDs

def seqID_to_num(seqID):
    return int(seqID.split('_')[-1])

if __name__ == '__main__':
    main()