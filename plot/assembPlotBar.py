#encoding=utf8
#! /usr/bin/python

#########################################################
# Program:assembPlotBar.py
# Usage: Plot the assemb quality and the sequencing depth relation.
# Author: michelia
# Version: 1.0
#       Editor: guoshuguang@scgene.com
#       Date:  
##########################################################

from __future__ import division
import argparse
import sys
import os
from collections import namedtuple
import matplotlib
matplotlib.use('Agg')  # solve run on the cluster's mistakes.
import matplotlib.pyplot as plt
from itertools import count
import time
sys.path.append('/scgene/tiger/invent/guoshuguang/program')
# import michelia
from michelia import vzebra, hzebra, data_dump, data_load, adjust_spines

function = '''
Plot the assemb quality and the sequencing depth relation.
主要是画横向的柱型图
'''


def main():
    parser = argparse.ArgumentParser(description=function,
                                    )
    parser.add_argument("dataPath", #type=int
                        help="this is the fasta data")
    # parser.add_argument("outfile", #type=int
    #                     help="this is the outfile")
    parser.add_argument('depthDataPath',
                        help='this is the depth data file path')
    parser.add_argument('figurePrefix', 
                        help='the figure save prefix')
    parser.add_argument('-f', '--format', default='png', 
                        choices = ['png', 'svg'],
                        dest='format',
                        help='Figure save format. If this is not, save the figure as png.')
    args = parser.parse_args()
    # args = parser.parse_args(['/home/guoshuguang/node/assembDepthPlot/mira_d_results/RF4.fa',
    #     '/home/guoshuguang/node/assembDepthPlot/mira_d_results/RF4.cns',
    #     '/home/guoshuguang/node/assembDepthPlot/RF4_assemb_depth_plot_test', '-f','png'])

    plot_fun(args.dataPath, args.depthDataPath, args.figurePrefix, args.format)

    # print
    # start_time = time.clock()
    # print 'Result sorted in %s file.' % args.outfile
    # print '[INFO] Elapsed time: %.4f' % (time.clock() - start_time)
    # print

def plot_fun(dataPath, depthDataPath, figurePrefix, format):
    data = tuple(parse_data(dataPath))
    with open(dataPath) as dataFileTemp:
        for line in dataFileTemp:
            if line.startswith('#') and 'coverage(%)' in line:
                coverage = line.split()[-1]
    ids, lengths = zip(*data)
    seqLength = sum(lengths)
    depthData = []
    depthDataPreDict = parse_depthDataPre(depthDataPath)
    for i in xrange(seqLength+1):
        depthData.append(depthDataPreDict.get(i, 0)) 
    start = 1
    leftStart = start
    bottom = 0
    compressionRatio = 1
    fig = plt.figure(figsize=(12, 6))
    ax1 = fig.add_subplot(312,position=[0.1, 0.3, 0.8, 0.03])
    # ax1.set_xlim(0, sum(lengths)*(1+0.1))
    ax1.set_xlim(0, sum(lengths)+50)
    ax2 = fig.add_subplot(311, position=[0.1, 0.34, 0.8, 0.5])  #depth plot axes
    ax2.set_ylim(ymin=-5)
    ax2.set_ylim(-1, 500)

    # ax2.set_xlim(0, sum(lengths)*(1+0.1))
    ax2.set_xlim(0, sum(lengths)+50)
    adjust_spines(ax1, ['bottom'])   
    adjust_spines(ax2, ['left', ])
    # adjust_spines(ax3, [])
    # vzebra(plt)  
    for (id, length) in data:
        if id.startswith('contig'): #plot the contig bar
            width = length/compressionRatio
            contigPlot = ax1.barh(bottom, left=leftStart, width=width, color='#0066cc',
                                height=0.5, edgecolor='none', align='center')
            leftStart += width
        else:  #plot the gap line
            width = length/compressionRatio
            gapPlot = ax1.barh(bottom, width=width, align='center',
                    left=leftStart, height=0.1, color='#00cc00', edgecolor='none')
            leftStart += width
    # depthPlot = ax2.fill_between(xrange(len(depthData)), depthData, color='#99ccff')  
    depthPlot = ax2.bar(xrange(len(depthData)), depthData, color='#99ccff', width=1, edgecolor='none')  
    ax2.set_title('Depth and Assembly information of %s (coverage: %s%%)' % (os.path.basename(dataPath).split('.')[0], coverage))
    # ax1.set_xlabel('Base position(%sbp)' % compressionRatio)
    ax1.set_xlabel('Base position(bp)')
    ax2.set_ylabel('Depth')
    # ax2.legend( [contigPlot, gapPlot], [ 'assembly', 'gap'])
    # ax3.legend( [contigPlot, gapPlot], [ 'assembly', 'gap'])
    # plt.show()
    plt.savefig(figurePrefix+'.'+format, dpi=1000)


# def compress_depth(depthData, start, end, compression=5):
#     length = end - start
#     for i in count():
#         comFrom = start + i * compression
#         comTo = comFrom + compression
#         if comTo > end:
#             yield (comFrom+comTo)/2, (sum(depthData[comFrom:])/len(depthData[comFrom:]))
#         else:
#             yield (comFrom+comTo)/2, (sum(depthData[comFrom:comTo])/compression)


def parse_data(dataPath):
    seqs = ''.join(list(parse_pre(dataPath)))
    flagID = 'gap'
    comeFrom = 0
    for i, base in enumerate(seqs):
        if flagID == 'gap':
            if base >= 'a':
                continue
            else:
                length = i - comeFrom
                yield flagID, length
                flagID = 'contig'
                comeFrom = i
        if flagID == 'contig':
            if base < 'a':
                continue
            else:
                length = i - comeFrom
                yield flagID, length
                flagID = 'gap'
                comeFrom = i
    yield flagID, (i - comeFrom)


def parse_pre(dataPath):
    with open(dataPath) as depthFile:
        for line in depthFile:
            if line.startswith('>'):
                break
        for line in depthFile:
            seq = ''.join(line.split()[1:])
            yield seq  #''.join(line.rstrip().split()[1:-1])


def parse_depthDataPre(depthDataPath):
    depthDataPreDict = {}
    with open(depthDataPath) as depthDataFile:
        for line in depthDataFile:
            line = line.split()
            depthDataPreDict[int(line[1])] = int(line[3])
    return depthDataPreDict


if __name__ == '__main__':
    main()