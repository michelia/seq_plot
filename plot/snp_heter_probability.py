from __future__ import division
import numpy as np
import argparse
import sys
import os
from collections import namedtuple
import matplotlib
# matplotlib.use('Agg')  # solve run on the cluster's mistakes.
import matplotlib.pyplot as plt
from itertools import count
import time
sys.path.append('/scgene/tiger/invent/guoshuguang/program/')
sys.path.append('/scgene/tiger/invent/guoshuguang/program/moduls')
from michelia import vzebra, hzebra, data_dump, data_load, \
      adjust_spines, adjust_ticks_label_size, csvreader, csvwriter
import path
def main():
    # parser = argparse.ArgumentParser(description=function,
    #                                 )
    # parser.add_argument("snpFilePath", #type=int,
    #                     help="this is the infile")
    # parser.add_argument("samFilePath", #type=int,
    #                     help="this is the outfile")
    shs = (
        ('/scgene/tiger/invent/chenjiehu/ECC/SNPANGENE/SNP/RBQ/RBQ_Chr1.snp', 
           '/scgene/tiger/invent/chenjiehu/ECC/SNPANGENE/SNP/RBQA/RBQA_Chr1.snp',
           '/scgene/tiger/invent/guoshuguang/figure/snpHeterProbability/Chr1_HeterSNP_Probability.png'),
        ('/scgene/tiger/invent/chenjiehu/ECC/SNPANGENE/SNP/RBQ/RBQ_Chr2.snp', 
           '/scgene/tiger/invent/chenjiehu/ECC/SNPANGENE/SNP/RBQA/RBQA_Chr2.snp',
           '/scgene/tiger/invent/guoshuguang/figure/snpHeterProbability/Chr2_HeterSNP_Probability.png'),
        ('/scgene/tiger/invent/chenjiehu/ECC/SNPANGENE/SNP/RBQ/RBQ_Chr3.snp', 
           '/scgene/tiger/invent/chenjiehu/ECC/SNPANGENE/SNP/RBQA/RBQA_Chr3.snp',
           '/scgene/tiger/invent/guoshuguang/figure/snpHeterProbability/Chr3_HeterSNP_Probability.png'),
        ('/scgene/tiger/invent/chenjiehu/ECC/SNPANGENE/SNP/RBQ/RBQ_Chr4.snp', 
           '/scgene/tiger/invent/chenjiehu/ECC/SNPANGENE/SNP/RBQA/RBQA_Chr4.snp',
           '/scgene/tiger/invent/guoshuguang/figure/snpHeterProbability/Chr4_HeterSNP_Probability.png'),
        ('/scgene/tiger/invent/chenjiehu/ECC/SNPANGENE/SNP/RBQ/RBQ_Chr5.snp', 
           '/scgene/tiger/invent/chenjiehu/ECC/SNPANGENE/SNP/RBQA/RBQA_Chr5.snp',
           '/scgene/tiger/invent/guoshuguang/figure/snpHeterProbability/Chr5_HeterSNP_Probability.png'),
        ('/scgene/tiger/invent/chenjiehu/ECC/SNPANGENE/SNP/RBQ/RBQ_Chr6.snp', 
           '/scgene/tiger/invent/chenjiehu/ECC/SNPANGENE/SNP/RBQA/RBQA_Chr6.snp',
           '/scgene/tiger/invent/guoshuguang/figure/snpHeterProbability/Chr6_HeterSNP_Probability.png'),
        ('/scgene/tiger/invent/chenjiehu/ECC/SNPANGENE/SNP/RBQ/RBQ_Chr7.snp', 
           '/scgene/tiger/invent/chenjiehu/ECC/SNPANGENE/SNP/RBQA/RBQA_Chr7.snp',
           '/scgene/tiger/invent/guoshuguang/figure/snpHeterProbability/Chr7_HeterSNP_Probability.png'),
        ('/scgene/tiger/invent/chenjiehu/ECC/SNPANGENE/SNP/RBQ/RBQ_Chr8.snp', 
           '/scgene/tiger/invent/chenjiehu/ECC/SNPANGENE/SNP/RBQA/RBQA_Chr8.snp',
           '/scgene/tiger/invent/guoshuguang/figure/snpHeterProbability/Chr8_HeterSNP_Probability.png'),
        ('/scgene/tiger/invent/chenjiehu/ECC/SNPANGENE/SNP/RBQ/RBQ_Chr9.snp', 
           '/scgene/tiger/invent/chenjiehu/ECC/SNPANGENE/SNP/RBQA/RBQA_Chr9.snp',
           '/scgene/tiger/invent/guoshuguang/figure/snpHeterProbability/Chr9_HeterSNP_Probability.png'),
        ('/scgene/tiger/invent/chenjiehu/ECC/SNPANGENE/SNP/RBQ/RBQ_Chr10.snp', 
           '/scgene/tiger/invent/chenjiehu/ECC/SNPANGENE/SNP/RBQA/RBQA_Chr10.snp',
           '/scgene/tiger/invent/guoshuguang/figure/snpHeterProbability/Chr10_HeterSNP_Probability.png'),
        ('/scgene/tiger/invent/chenjiehu/ECC/SNPANGENE/SNP/RBQ/RBQ_Chr11.snp', 
           '/scgene/tiger/invent/chenjiehu/ECC/SNPANGENE/SNP/RBQA/RBQA_Chr11.snp',
           '/scgene/tiger/invent/guoshuguang/figure/snpHeterProbability/Chr11_HeterSNP_Probability.png'),
        ('/scgene/tiger/invent/chenjiehu/ECC/SNPANGENE/SNP/RBQ/RBQ_Chr12.snp', 
           '/scgene/tiger/invent/chenjiehu/ECC/SNPANGENE/SNP/RBQA/RBQA_Chr12.snp',
           '/scgene/tiger/invent/guoshuguang/figure/snpHeterProbability/Chr12_HeterSNP_Probability.png'),
        )
    for sh in shs:
        probability_and_plot(*sh)

def probability_and_plot(filePath, AfilePath, savePath):
    # filePath = '/scgene/tiger/invent/chenjiehu/ECC/SNPANGENE/SNP/RBQ/RBQ_Chr1.snp'
    # savePath = '/scgene/tiger/invent/guoshuguang/figure/snpHeterProbability/RBQ_Chr1_Prob.png'
    lines = []
    with open(filePath) as File:
        File.readline()
        for line in csvreader(File):
            lines.append(float(line[9]))
    lines.sort()
    # print 'lines', len(lines), lines
    # for i in lines:
    #     print 'iii:', i
    # aNum = 0
    # for one in lines:
    #     if 1 == one:#and one<= 1:
    #         aNum += 1
    # print 'aNum', aNum
    Alines = []
    with open(AfilePath)as AFile:
        AFile.readline()
        for Aline in csvreader(AFile):
#            print 'Alines', Alines
            Alines.append(float(Aline[9]))
    Alines.sort()
    probability = []
    Aprobability = []
    # probability = {}
    sPath = path.path(savePath)
    outDataFilePath = sPath.parent / sPath.namebase + '.dataFile'
    # print outDataFilePath
    outDataFile = open(outDataFilePath, 'w')
    outDataFile.write('#range\trpq\trpqa\n')
    write = csvwriter(outDataFile)
    for end in np.arange(0.05, 1.05, 0.05):
        frequence = stat_probability(lines, end)
        lines = lines[frequence:]
        probability.append(frequence)
        Afrequence = stat_probability(Alines, end)
        Alines = Alines[Afrequence:]
        Aprobability.append(frequence)
        write.writerow(('%s<=p<%s' % (end-0.05, end), frequence, Afrequence))
    plot_probability(probability, Aprobability, savePath)
    outDataFile.close()
    print 'probability',  sum(probability),probability

def stat_probability(aList, end, step=0.05):
    for i, one in enumerate(aList):
        if one >= end:
            break
    aList = aList[i:]
    return i

def plot_probability(probability, Aprobability, savePath):
    fig = plt.figure(figsize=(10, 4))
    plt.subplots_adjust(left=0.06, right=0.95, bottom=0.2, top=0.9)
    ax1 = plt.gca()
    xticks = np.arange(0, 1.05, 0.05)
    ax1.xaxis.set_ticks(xticks)
    ax1.set_xlim(xmax=1, xmin=0)
    adjust_spines(ax1, spines=['bottom', 'left'], width=0.4)
    for i, start in enumerate(np.arange(0, 1, 0.05)):
            p = plt.bar(left=start, width=0.025, color='#0066cc',
                    height=probability[i], edgecolor='none', align='edge')
            pA= plt.bar(left=start+0.025, width=0.025, color='red',
                    height=Aprobability[i], edgecolor='none', align='edge')
    plt.legend((p[0], pA[0]), ('RBQ', 'RBQA'))
    ax1.set_title(path.path(savePath).namebase)
    # ax1.set_xlabel()
    ax1.set_ylabel('Frequency')
    #ax1.set_ylabel('Probability (%)')
    plt.savefig(savePath, dpi=500)
    # plt.show()


if __name__ == '__main__':
    main()
