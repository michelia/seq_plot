from __future__ import division
import sys
sys.path.append('/scgene/tiger/invent/guoshuguang/program/moduls')
sys.path.append('/scgene/tiger/invent/guoshuguang/program/')
from michelia import blastn, makeblastdb, parse_blastout,\
    CsvReader, invert_dict
import matplotlib.pyplot as plt
from michelia import vzebra, hzebra, data_dump, data_load, \
      adjust_spines, adjust_ticks_label_size, frange
from path import path
import re

def main():
    chipFaFilePath = '/scgene/tiger/invent/guoshuguang/\
chip/compare_with_chip/seq/chip/all.fa'
    digitalChipFilePath = '/scgene/tiger/invent/guoshuguang/\
chip/compare_with_chip/seq/digital_chip/all.cdna'
    blastoutPath = '/scgene/tiger/invent/guoshuguang/\
chip/compare_with_chip/blastout_chip-digital.fmt6'
    chipExpressFilePath = '/scgene/tiger/invent/guoshuguang/\
chip/compare_with_chip/express/chip/all_expression_signal.csv'
    digitalChipExpressFilePathE = '/scgene/tiger/invent/guoshuguang/\
chip/compare_with_chip/express/digital_chip/14E0_5.GeneExpression.csv'
    digitalChipExpressFilePathL = '/scgene/tiger/invent/guoshuguang/\
chip/compare_with_chip/express/digital_chip/14L0_5.GeneExpression.csv'
    savePathE = '/scgene/tiger/invent/guoshuguang/\
chip/compare_with_chip/figure/E_fig.png'
    savePathL = '/scgene/tiger/invent/guoshuguang/\
chip/compare_with_chip/figure/L_fig.png'

    # makeblastdb(infile=digitalChipFilePath)
    # blastn(evalue=10**-100, query=chipFaFilePath,db=digitalChipFilePath,
    #         outfmt=6, out=blastoutPath, num_threads=16,)
#     compareDataPath = '/scgene/tiger/invent/guoshuguang/\
# chip/compare_with_chip/compareData.piso'
    # if path(compareDataPath).exists():
    #     compareData = data_load(compareDataPath)
    # else:
    #     compareData = compare_data(chipExpressFilePath, 
                        # digitalChipExpressFilePathE, 
    #                 idDict=id_one_to_one(blastoutPath))
    #     data_dump(compareData, compareDataPath)
    compareDataE = compare_data(chipExpressFilePath, digitalChipExpressFilePathE, 
                    blastoutPath, 40)
    compare_plot(compareDataE, savePathE)
    compareDataL = compare_data(chipExpressFilePath, digitalChipExpressFilePathL, 
                    blastoutPath, 40, flagL=1)
    compare_plot(compareDataL, savePathL)
#     dataPath = '/scgene/tiger/invent/chenjiehu/\
# GJXDGE/classify_tag/11L0_5.geneExpression_submit.xls'
#     compare_plot(parse_data(dataPath, chipExpressFilePath))
#     pass

def parse_data(dataPath, chipExpressFilePath):
    tmpData = {}
    dataFile = CsvReader(dataPath)
    for line in dataFile.reader():
        tmpData[line[0]] = float(line[-1]) * 30
    chipDataDict = chip_id_data(chipExpressFilePath)
    data = []
    for chipID in tmpData:
        if chipID in chipDataDict:
            data.append((tmpData[chipID], chipDataDict[chipID][0]))
    print 'data:', len(data)
    return data


def compare_data(chipExpressFilePath, digitalChipExpressFilePath,
                blastoutPath, multiple, flagL=False):
    idDict = id_one_to_one(blastoutPath)  #{ChipID: digitalChipID}
    compareData = []
    chipDataDict = chip_id_data(chipExpressFilePath)
    digitalChipDataDict = digital_chip_id_data(digitalChipExpressFilePath)
    # digitalChipDataDict = new_digital_chip_id_data(digitalChipExpressFilePath)
    # mapedReadsNum = len(digitalChipDataDict)# / 10**6
    for oneChipId in chipDataDict:
        if oneChipId in idDict and idDict[oneChipId] in digitalChipDataDict:
            oneDigitalId = idDict[oneChipId] 
            tmp = digitalChipDataDict[oneDigitalId] * multiple
            chipVialue = chipDataDict[oneChipId][0]
            if flagL:
                chipVialue = chipDataDict[oneChipId][1]
            compareData.append((tmp, chipVialue))
            # print 'onedata',(tmp, chip_vialue)
    # for oneDigitalId in digitalChipDataDict:
    #     if oneDigitalId in idDict and idDict[oneDigitalId] in chipDataDict:
    #         tmp = digitalChipDataDict[oneDigitalId] * 130
    #         # digitalReadsNum = mapedReadsNum * tmp
    #         chip_vialue = chipDataDict[idDict[oneDigitalId]][0]
    #         # compareData.append((digitalReadsNum, chip_vialue))
    #         compareData.append((tmp, chip_vialue))
    print len(compareData)
    return compareData

def compare_plot(compareData, savePath):
    fig = plt.figure(figsize=(4, 4))
    ax1 = plt.gca()
    max = 20000
    ax1.set_xlim(xmax=max)
    ax1.set_ylim(ymax=max)
    plt.plot((0, max),(0, max))

    # adjust_spines(ax1, spines=['bottom', 'left'], width=0.4)
    tmp, chip_vialue = zip(*compareData)
    tmp = map(float, tmp)
    chip_vialue = map(float, chip_vialue)
    print 'tmp:', len(tmp)
    print 'chip_vialue: ', len(chip_vialue)
    plt.scatter(tmp, chip_vialue, s=0.1)
    # ax1.set_xlabel()
    # ax1.set_ylabel('Frequency')
    #ax1.set_ylabel('Probability (%)')

    plt.savefig(savePath, dpi=500)
    # plt.show()

def chip_id_data(chipExpressFilePath):
    chipDataDict = {}
    ChipFile = CsvReader(chipExpressFilePath, 1)
    for line in ChipFile.reader():
        if line[2] == 'P' and line[4] == 'P':
            chipDataDict[line[0]] = (float(line[1]), float(line[3]))
    ChipFile.close()
    return chipDataDict

def digital_chip_id_data(digitalChipExpressFilePath):
    digitalChipDataDict = {}
    digitalFile = CsvReader(digitalChipExpressFilePath, 15)
    for line in digitalFile.reader():
        digitalID = line[0].partition('.')[0]
        digitalChipDataDict[digitalID] = float(line[2])
    digitalFile.close()
    return digitalChipDataDict

def new_digital_chip_id_data(digitalChipExpressFilePath):
    digitalChipDataDict = {}
    digitalFile = CsvReader(digitalChipExpressFilePath, 15)
    for line in digitalFile.reader():
        digitalChipDataDict[line[0].partition('.')[0]] = float(line[2])
    digitalFile.close()
    return digitalChipDataDict

def id_one_to_one(blastoutPath):
    ll = 'idDatall'
    if not path(ll).exists():
        idDict = {}  # {digitalChipID: ChipID}
        for blastReader in parse_blastout(blastoutPath):
            chipID = blastReader.queryId
            if chipID not in idDict:
                digitalID = blastReader.subjectId#.partition('.')[0]
                idDict[chipID] = digitalID
        # data_dump(idDict, ll)
    else:
        idDict = data_load(ll)

    print 'idDict:', len(idDict)
    print 'invert idDict:', len(invert_dict(idDict))
    return idDict

if __name__ == '__main__':
    main()