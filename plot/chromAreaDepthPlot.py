#! /usr/bin/python
#encoding=utf8

#########################################################
# Program: chromAreaDepth.py
# Usage: Plot the multi-chromosomes area sequnces depth 
# Author: michelia
# Version: 1.0
#       Editor: guoshuguang@scgene.com
#       Date:  
##########################################################

from __future__ import division
import sys
sys.path.append('/scgene/tiger/invent/guoshuguang/program/')
sys.path.append('/scgene/tiger/invent/guoshuguang/program/moduls')
from michelia import vzebra, hzebra, data_dump, data_load,\
         adjust_spines, adjust_ticks_label_size, CsvReader
from collections import namedtuple
import argparse
import matplotlib
import matplotlib.pyplot as plt
import time

function = '''
Plot the multi-chromosomes area sequnces depth.
'''

def main():
    parser = argparse.ArgumentParser(description=function,
        usage = 'Plot the multi-chromosomes area sequnces depth.'
                                    )
    parser.add_argument("dataPath", 
                        help="""the data file. 
                        file format:chromName start end length depth""")
    parser.add_argument('faiFilePath', 
                        help='the chrome length file. Format: chromeName  chromeLength')
    parser.add_argument('outDir', 
                        help='the Dir to save figures')
    args = parser.parse_args()
    area_depth_plot(args.faiFilePath, args.dataPath, args.outDir)

def area_depth_plot(faiFilePath, dataPath, figSaveDir):
    chromsLenDirt = chroms_len(faiFilePath)
    dataDict = parse_chrom_depth_data(dataPath)
    for chromID in chromsLenDirt:
        oneChromLen = chromsLenDirt[chromID]
        fig = plt.figure(figsize=(10, 1))
        plt.subplots_adjust(left=0.05, right=0.95, bottom=0.3, top=0.7)
        ax1 = fig.add_subplot(111)
        ax1.set_ylim(ymax=100)
        ax1.set_xlim(xmax=oneChromLen+50)
        ax1.set_xlabel('Base postion (/Mb)', fontsize=5) 
        ax1.set_ylabel('Depth', fontsize=5) 
        step = 1000000
        xticks = range(0, oneChromLen, step)
        xticks_label = map(lambda i: str(i//step), xticks)
        ax1.xaxis.set_ticks(xticks)
        ax1.xaxis.set_ticklabels(xticks_label)
        adjust_ticks_label_size(ax1, size=5)
        adjust_spines(ax1, spines=['bottom', 'left'], width=0.3)
        chromName = chromID
        ax1.set_title(chromName, fontsize=5)
        for (start, end, length, depth) in dataDict[chromName]:
            ax1.barh(0, left=start, width=length, color='#0066cc',
                    height=depth, edgecolor='none', align='center')
        savePath = figSaveDir + '/' + chromID + 'AreaDepthFigure.png'
        plt.savefig(savePath, dpi=500)

def chroms_len(faiFilePath):
    faiFile = CsvReader(faiFilePath)
    chromsLenDirt = {} # {chromID: chromLen}
    for line in faiFile.reader():
        chromsLenDirt[line[0]] = int(line[1])
    return chromsLenDirt


def parse_chrom_depth_data(dataFilePath):
    dataDict = {}
    with open(dataFilePath) as aHandle:
        aHandle.readline()
        for line in aHandle:
            line = line.split()
            if line[0] not in dataDict:
                dataDict[line[0]] = []
            dataDict[line[0]].append(map(float, line[1:]))
    return dataDict

if __name__ == '__main__':
    main()
    print 'Finish.'
