#! /usr/bin/python

#########################################################
# Program:qualityCheckPy.py
# Usage: Acording the data that jidhu'code generate, to plot the fastq quality figure.
# Author: michelia
# Version: 1.0
#       Editor: name@scgene.com
#       Date: Friday, 14 December, 2012 11:20:15 AM CST
##########################################################

from __future__ import division
import sys
import os
from collections import namedtuple
import matplotlib.pyplot as plt
from itertools import count
import argparse
import gzip
import re
import time
# from rpy2 import robjects
# from rpy2.robjects.packages import importr

# sys.path.append('/scgene/tiger/invent/guoshuguang/program')
from michelia import frange, data_dump, data_load, vzebra, hzebra


function = '''
Acording the data that jidhu'code generate, to plot the fastq quality figure.
'''

def main():
    parser = argparse.ArgumentParser(description=function)
    parser.add_argument("fq", #type=int
                        help="input fastq")
    parser.add_argument("outDir", #type=int
                        help="output folder")
    # args = parser.parse_args()
    # args = parser.parse_args(['/home/guoshuguang/node/fastqQualitPlot/PR.fq',
    #                          '/home/guoshuguang/node/fastqQualitPlot/'])
    args = parser.parse_args()

    if not os.path.exists(args.outDir):
        os.mkdir(args.outDir)
    fqname = os.path.basename(args.fq)
    dataFile = os.path.join(args.outDir, fqname + '.plot.dat')
    savePath = os.path.join(args.outDir, fqname) + 'Py'

    # qrqc_fastq_plot(args.fq, qrqcPngSavePath=(savePath+'qrqc.png')) # the qrqc plot
    read_lengths_plot(args.fq, savePath+'.lengthHist.png') # Acording the every read's length to plot histogram

    if not os.path.exists(dataFile): # if no the dataFile, perform jiehu's code to generate dataFile
        staFQ(args.fq, args.outDir) # generate dataFile
    plot(dataFile, savePath) # plot function, not contain the qrqc plot
    print

def read_lengths_plot(fqPath, pngSavePath):
    readsLenghsPicPath = os.path.join(os.path.dirname(pngSavePath), os.path.basename(fqPath)) + '.readsLenghs.piso'
    if not os.path.exists(readsLenghsPicPath):
        readsLenghs = []
        for read in parserfq(fqPath):
            readsLenghs.append(len(read[1]))
        data_dump(readsLenghsPicPath, readsLenghs)
    else:
        readsLenghs = data_load(readsLenghsPicPath)
    plot_hist(readsLenghs, pngSavePath)

def plot_hist(readsLenghs, pngSavePath):
    plt.figure(figsize=(8, 8))
    ax = plt.gca()
    plt.subplots_adjust(left=0.1, right=0.95, bottom=0.1, top=0.95)
    # plt.xlim(0, 101)
    # plt.ylim(0, 50)
    # ax.set_xticks(range(0, 101, 10))
    # ax.set_yticks(range(0, 61, 5))
    # readsLenghs = [1,2,4,5,8,9,6,2,5,4,1,2,5]
    plt.hist(readsLenghs, 50, normed=1)
    # plt.hist(readsLenghs, 100, normed=1, color='#6633AA')
    # plt.legend()
    plt.title('Sequencing Quality')
    plt.xlabel('Read length')
    plt.ylabel('Density')
    plt.savefig(pngSavePath, dpi=150)


def parserfq(fqPath):
    # fqfile = gzip.open(fqPath)
    fqfile = open(fqPath)
    next = fqfile.next
    for line in fqfile:
        yield (line.strip(), next().strip(), next().strip(), next().strip())
    fqfile.close()


def plot(dataFile, savePath):
    Data = namedtuple('Data', 
        'Cycle Base A T G C N Aperc Tperc Gperc Cperc Nperc AvgQ ModQ ErrorRatePerc')
    data = Data(*zip(*perse_date(dataFile))) # read data (jiehu's code generate dataFile)
    # the follow three function are different plot function
    quality_plot(data.AvgQ, data.ModQ, savePath)
    error_plot(data.ErrorRatePerc, savePath)
    base_perc_plot(data.Aperc, data.Tperc, data.Gperc, data.Cperc, data.Nperc, savePath)

def quality_plot(AvgQ, ModQ, savePath):
    savePath = savePath + '.quality.png'
    num = len(ModQ)
    plt.figure(figsize=(10, 8))
    ax = plt.gca()
    plt.subplots_adjust(left=0.1, right=0.95, bottom=0.1, top=0.95)
    plt.xlim(0, len(ModQ)+1)
    plt.ylim(0, 50)
    ax.set_xticks(range(0, num+1, num//10))
    ax.set_yticks(range(0, 61, 5))
    hzebra()
    plt.grid(axis='x', alpha=0.4) # just display the x axis grid line
    plt.plot(range(1, num+1), ModQ, 'g>-', color='#339999', label='Mode Quality', alpha=1, markeredgewidth=0)
    plt.plot(range(1, num+1), AvgQ, 'ro-', color='#9966CC', label='Average Quality', alpha=1, markeredgewidth=0)
    plt.legend()
    plt.title('Sequencing Quality')
    plt.xlabel('Cycles')
    plt.ylabel('Quality')
    plt.savefig(savePath, dpi=150) 

def error_plot(ErrorRatePerc, savePath):
    savePath = savePath + '.errorRate.png'
    num = len(ErrorRatePerc)
    plt.figure(figsize=(10, 8))
    ax = plt.gca()
    plt.subplots_adjust(left=0.1, right=0.95, bottom=0.1, top=0.95)
    plt.xlim(0, len(ErrorRatePerc)+1)
    plt.ylim(-(dyna_yaxis_error(ErrorRatePerc)/8/20), dyna_yaxis_error(ErrorRatePerc))
    ax.set_xticks(range(0, num+1, num//10))
    ax.set_yticks(frange(0, dyna_yaxis_error(ErrorRatePerc)+0.001, dyna_yaxis_error(ErrorRatePerc)/8))
    hzebra(plt)
    plt.grid(axis='x', alpha=0.4)
    plt.bar(range(1, num+1), ErrorRatePerc, width=0.37, label='Error Rate',
            color='red', align='center', linewidth=0, alpha=0.7)
    plt.legend()
    plt.title('Sequencing Error Rate')
    plt.xlabel('Cycles')
    plt.ylabel('Error Rate(X100)')
    plt.savefig(savePath, dpi=150)

def base_perc_plot(Aperc, Tperc, Gperc, Cperc, Nperc, savePath):
    savePath = savePath + '.GC.png'
    num = len(Aperc)  # statistics the number of sequence's cycles
    plt.figure(figsize=(12, 8))
    ax = plt.gca()
    plt.subplots_adjust(left=0.1, right=0.95, bottom=0.1, top=0.95)
    plt.xlim(0, num+1)
    plt.ylim(-1, 51)
    ax.set_xticks(range(0, num+1, num//10))
    ax.set_yticks(range(0, 51, 5))
    hzebra()
    plt.grid(axis='x', alpha=0.4) # just display the x axis grid line
    plt.plot(range(1, num+1), Aperc, 'p-', label='A', markersize=7,
            markerfacecolor='none', markeredgecolor='#000080', markeredgewidth=0.8)
    plt.plot(range(1, num+1), Tperc, '*-', label='T', markersize=8,
            markeredgewidth=0, color='#FF8C00')
    plt.plot(range(1, num+1), Gperc, 's-', label='G', markerfacecolor='none', 
                color='r', markeredgecolor='r', markeredgewidth=0.9,)
    plt.plot(range(1, num+1), Cperc, 'x-', label='C', color='g')
    plt.plot(range(1, num+1), Nperc, '2-', label='N', color='#8A2BE2')
    plt.legend()
    plt.title('GC Precentage of Sequencing')
    plt.xlabel('Cycles')
    plt.ylabel('Base Precentage(X100)')
    plt.savefig(savePath, dpi=150)  


def perse_date(dataFile):
    with open(dataFile) as data:
        for line in data:
            if line[0] != '#' and line.strip():
                yield map(float,line.split())

def qrqc_fastq_plot(filePath, qrqcPngSavePath):
    grdevices = importr('grDevices')
    rprint = robjects.globalenv.get('print')
    qrqc = importr('qrqc')
    grdevices.png(file=qrqcPngSavePath, width=512, height=512)
    fastq = qrqc.readSeqFile(filePath, type='fastq')
    rprint(qrqc.qualPlot(fastq))
    grdevices.dev_off()

def dyna_yaxis_error(ErrorRatePerc):  # get double-digit
    return ((max(ErrorRatePerc) * 100 // 10) + 1) / 10


def vzebra(plt=plt, color='#DCDCDC', alpha=0.5, align='edge', ax=None):
    """
    vertical zebra
    """
    if not ax:
        ax = plt.gca()
    xticklocs = ax.xaxis.get_ticklocs()
    xbar = xticklocs[1::2]
    # ymax = ax.get_ylim()[1] #two methods to get the max 
    ymax = ax.axis()[3]
    ax.bar(xbar, [ymax] * len(xbar), width=abs(xticklocs[1] - xticklocs[0]),
        color=color, alpha=alpha, align=align, linewidth=0)

def hzebra(plt=plt, color='#DCDCDC', alpha=0.5, align='edge', ax=None):
    """
    horizontal zebra 
    """
    if not ax:
        ax = plt.gca()
    yticklocs = ax.yaxis.get_ticklocs()
    ybar = yticklocs[1::2]
    # xmax = ax.get_xlim()[1] #two methods to get the max 
    xmax = ax.axis()[1]
    ax.barh(ybar, [xmax] * len(ybar), height=abs(yticklocs[1] - yticklocs[0]),
        color=color, alpha=alpha, align=align, linewidth=0)

# the follow is jiehu's code.
def quaStandar(min_qua, max_qua):
    qua_standar = ""
    qua_standar_diff = 0
    if min_qua == 33 and max_qua == 73:
        qua_standar = "Sanger_Phred+33"
        qua_standar_diff = 33
    elif min_qua == 59 and max_qua == 104:
        qua_standar = "Solexa_Solexa+64"
        qua_standar_diff = 64
    elif min_qua == 64 and max_qua == 104:
        qua_standar = "Illumina1.3+_Phred+64"
        qua_standar_diff = 64
    elif min_qua == 66 and max_qua == 105:
        qua_standar = "Illumina1.5+_Phred+64"
        qua_standar_diff = 64
    elif min_qua == 33 and max_qua == 74:
        qua_standar = "Illumina1.8+_Phred+33"
        qua_standar_diff == 33
    elif max_qua == 73:
        qua_standar = "Edit_Sanger_Phred+33"
        qua_standar_diff == 33
    elif max_qua == 104:
        qua_standar = "Edit_Illumina1.3+_Phred+64"
        qua_standar_diff = 64
    elif max_qua == 74:
        qua_standar = "Edit_Illumina1.8+_Phred+33"
        qua_standar_diff = 33
    elif min_qua == 33:
        qua_standar = "Edit_Phred+33"
        qua_standar_diff = 33
    elif min_qua >= 33 and min_qua < 66:
        qua_standar = "Unknown_Phred+33"
        qua_standar_diff = 33
    else:
        print >>sys.stderr, "QUALITY ERROR, unknown quality standar, raw max\
quality is %d and min quality is %d" % (max_qua, min_qua)
        sys.exit(1)
    return (qua_standar, qua_standar_diff)


def errorRate(quality, qua_diff):
    if qua_diff == 33:
        return (10 ** (quality / -10.0)) * 100
    else:
        e_rate = 10 ** (quality / -10.0)
        return e_rate / (1 + e_rate) * 100


def staFQ(fq, out_dir):
    name = os.path.basename(fq)
    dict_qua = {} #dict_qua[cycle][quality]
    dict_gc = {'A' : {}, 'T' : {}, 'G' : {}, 'C' : {}, 'N' : {}}
    list_cycle = []
    all_reads = 0
    raw_min_qua = 100
    raw_max_qua = 0
    if re.search('gz$', name):
        FQ = gzip.open(fq, 'r')
    else:
        FQ = open(fq, 'r')
    while True:
        id_1 = FQ.readline()
        if len(id_1) == 0:
            break
        seq = FQ.readline().rstrip()
        FQ.readline()
        qua = FQ.readline().rstrip()

        all_reads += 1

        length = len(seq)
        if length > len(list_cycle): #initialize cycles list
            #name a list of cycles to store information #Reversion 1.3 edit
            for x in range(len(list_cycle), length):
                list_cycle.append(0)

        for i in range(0, length):
            list_cycle[i] += 1
            base = seq[i]
            base_qua_chr = qua[i]
            base_qua = ord(base_qua_chr)
            if base_qua > raw_max_qua:
                raw_max_qua = base_qua
            if base_qua < raw_min_qua:
                raw_min_qua = base_qua

            if dict_qua.has_key(i):
                if dict_qua[i].has_key(base_qua):
                    dict_qua[i][base_qua] += 1
                else:
                    dict_qua[i][base_qua] = 1
            else:
                dict_qua[i] = {}

            if dict_gc[base].has_key(i):
                dict_gc[base][i] += 1
            else:
                dict_gc[base][i] = 1
    FQ.close()

    (qua_standar, qua_diff) = quaStandar(raw_min_qua, raw_max_qua)

    q20 = 0
    q30 = 0
    all_bases = 0
    all_quality = 0
    all_gc = 0
    all_base_noN = 0
    all_qua = 0
    all_mod_quas = {}

    outplot = open(os.path.join(out_dir, name + '.plot.dat'), 'w')
    outplot.write('#Cycle\tBase\tA\tT\tG\tC\tN\tA%\tT%\tG%\tC%\tN%\tAvgQ\tModeQ\tErrorRate%\n')
    for m in range(0, len(list_cycle)): # do analysis for each cycle
        cycle = m + 1
        bases_per_cycle = float(list_cycle[m])
        if bases_per_cycle == 0:
            print >>sys.stderr, 'error cycle', m
            sys.exit(1)
        all_bases += bases_per_cycle

        A = 0
        T = 0
        G = 0
        C = 0
        N = 0
        #not every cycle has 'ATGCN'
        if dict_gc['A'].has_key(m):
            A = dict_gc['A'][m]
        if dict_gc['T'].has_key(m):
            T = dict_gc['T'][m]
        if dict_gc['G'].has_key(m):
            G = dict_gc['G'][m]
        if dict_gc['C'].has_key(m):
            C = dict_gc['C'][m]
        if dict_gc['N'].has_key(m):
            N = dict_gc['N'][m]

        A_rate = A / bases_per_cycle * 100
        T_rate = T / bases_per_cycle * 100
        G_rate = G / bases_per_cycle * 100
        C_rate = C / bases_per_cycle * 100
        N_rate = N / bases_per_cycle * 100
        all_base_noN += (A + T + G + C)
        all_gc += (G + C)

        total_qua = 0
        mode_qua = 0
        mode_qua_bases = 0

        for quality in dict_qua[m].keys():
            cquality = quality - qua_diff
            total_qua += (cquality * dict_qua[m][quality])
            if cquality >= 20:
                q20 += dict_qua[m][quality]
            if cquality >= 30:
                q30 += dict_qua[m][quality]
            if  dict_qua[m][quality] > mode_qua_bases:
                mode_qua = cquality
                mode_qua_bases = dict_qua[m][quality]

            if cquality in all_mod_quas:  #for all mod qualitys statics
                all_mod_quas[cquality] += dict_qua[m][quality]
            else:
                all_mod_quas[cquality] = dict_qua[m][quality]

            all_qua += total_qua

        avg_qua = total_qua / bases_per_cycle

        error_rate = errorRate(avg_qua, qua_diff)

        outplot.write("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%f\n"\
                    % (cycle, bases_per_cycle, A, T, G, C, N, A_rate, T_rate, G_rate, C_rate, N_rate, avg_qua, mode_qua, error_rate))

    q20_rate = q20 / float(all_bases) * 100
    q30_rate = q30 / float(all_bases) * 100
    all_gc_rate = all_gc / float(all_base_noN) * 100
    all_avg_qua = all_qua / float(all_base_noN)
    all_mod_qua = sorted(all_mod_quas.items(), key = lambda d:d[1])[-1][0]
    avg_mod_qua_diff = abs(all_mod_qua - all_avg_qua)
    avg_error_rate = errorRate(all_avg_qua, qua_diff)
    mod_error_rate = errorRate(all_mod_qua, qua_diff)

    outplot.write("\n#File\tTotal_reads\tRead_len\tTotal_base\tGC%\t\
Quality_standar\tQ20\tQ30\tAvg_quality\tMod_quality\tabs(AvgQ-ModQ)\tAvg_error_rate%\t\
Mod_error_rate%\n")
    outplot.write("#%s\t%d\t%d\t%d\t%.2f\t%s\t%.2f\t%.2f\t%.2f\t%d\t%.2f\t%.6f\t%.6f\n"\
            % (name, all_reads, len(list_cycle), all_bases, all_gc_rate,\
                qua_standar, q20_rate, q30_rate, all_avg_qua, all_mod_qua,\
                avg_mod_qua_diff, avg_error_rate, mod_error_rate))

    outplot.close()

def parse_fq(fq):
    if re.search('gz$', name):
        FQ = gzip.open(fq, 'r')
    else:
        FQ = open(fq, 'r')   
    next = FQ.next
    for line in FQ:
        if line.strip() == None:
            break
        yield (line, next(), next(), next())


if __name__ == '__main__':
    start_time = time.clock()
    main()
    # print 'Result sorted in %s file.' % args.outfile
    print '[INFO] Elapsed time: %.4f' % (time.clock() - start_time)
    print
