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
patter = re.compile('\d+:0.\d+')

def main():
    chipFaFilePath = '/scgene/tiger/invent/guoshuguang/\
chip/compare_with_chip/seq/chip/all.fa'
    digitalChipFilePath = '/scgene/tiger/invent/guoshuguang/\
chip/compare_with_chip/seq/digital_chip/all.cdna'
    blastoutPath = '/scgene/tiger/invent/guoshuguang/\
chip/compare_with_chip/blastout_chip-digital.fmt6'
    chipExpressFilePath = '/scgene/tiger/invent/guoshuguang/\
chip/compare_with_chip/express/chip/all_expression_signal.csv'
#-----------------------E
#     digitalChipExpressFilePath = '/scgene/tiger/invent/guoshuguang/\
# chip/compare_with_chip/express/digital_chip/14E0_5.GeneExpression.csv'
#     digitalChipExpressFilePathNew = '/scgene/tiger/invent/guoshuguang/\
# chip/compare_with_chip/express/digital_chip/14E0_5.geneExpression.xls2'

#-----------------------L
    digitalChipExpressFilePath = '/scgene/tiger/invent/guoshuguang/\
chip/compare_with_chip/express/digital_chip/14L0_5.GeneExpression.csv'
    digitalChipExpressFilePathNew = '/scgene/tiger/invent/guoshuguang/\
chip/compare_with_chip/express/digital_chip/14L0_5.geneExpression.xls2'

    saveDir = '/scgene/tiger/invent/guoshuguang/chip/\
compare_with_chip/figure/'

    compareData = compare_data(chipExpressFilePath, digitalChipExpressFilePath, 
                    blastoutPath, multiple=40, flagL=True)
    compare_plot(compareData, saveDir+'L_fig.png')
    compareDataNew = compare_data(chipExpressFilePath, digitalChipExpressFilePathNew, 
                    blastoutPath, multiple=40, flagL=True)
    compare_plot(compareDataNew, saveDir+'L_new.png','r')
    compare_mix_plot(compareData,compareDataNew ,saveDir+'L_mix.png')

# def parse_data(dataPath, chipExpressFilePath):
#     tmpData = {}
#     dataFile = CsvReader(dataPath)
#     for line in dataFile.reader():
#         tmpData[line[0]] = float(line[-1]) * 30
#     chipDataDict = chip_id_data(chipExpressFilePath)
#     data = []
#     for chipID in tmpData:
#         if chipID in chipDataDict:
#             data.append((tmpData[chipID], chipDataDict[chipID][0]))
#     print 'data:', len(data)
#     return data


def compare_data(chipExpressFilePath, digitalChipExpressFilePath,
                blastoutPath, multiple, flagL=False):
    idDict = id_one_to_one(blastoutPath)  #{ChipID: digitalChipID}
    compareData = []
    chipDataDict = chip_id_data(chipExpressFilePath)
    digitalChipDataDict = digital_chip_id_data(digitalChipExpressFilePath)
    for oneChipId in chipDataDict:
        if oneChipId in idDict and idDict[oneChipId] in digitalChipDataDict:
            oneDigitalId = idDict[oneChipId] 
            tmp = digitalChipDataDict[oneDigitalId] * multiple
            chipVialue = chipDataDict[oneChipId][0]
            if flagL:
                chipVialue = chipDataDict[oneChipId][1]
            compareData.append((tmp, chipVialue))
    print len(compareData)
    return compareData

def compare_mix_plot(compareData, dataNew, savePath, color='b'):
    fig = plt.figure(figsize=(4, 4))
    ax1 = plt.gca()
    max = 20000
    ax1.set_xlim(xmax=max)
    ax1.set_ylim(ymax=max)
    plt.plot((0, max),(0, max),c='g')

    adjust_ticks_label_size(ax1, 7)
    tmp, chip_vialue = zip(*compareData)
    tmp = map(float, tmp)
    chip_vialue = map(float, chip_vialue)
    print 'tmp:', len(tmp)
    print 'chip_vialue: ', len(chip_vialue)
    plt.scatter(tmp, chip_vialue, s=0.5, c='b', edgecolor='none')
#---------------------------New
    tmp, chip_vialue = zip(*dataNew)
    tmp = map(float, tmp)
    chip_vialue = map(float, chip_vialue)
    print 'tmp2:', len(tmp)
    print 'chip_vialue2: ', len(chip_vialue)
    plt.scatter(tmp, chip_vialue, s=0.5, c='r', edgecolor='none')
    # ax1.set_xlabel()
    # ax1.set_ylabel('Frequency')
    #ax1.set_ylabel('Probability (%)')

    plt.savefig(savePath, dpi=500)
    # plt.show()

def compare_plot(compareData, savePath, color='b'):
    fig = plt.figure(figsize=(4, 4))
    ax1 = plt.gca()
    max = 20000
    ax1.set_xlim(xmax=max)
    ax1.set_ylim(ymax=max)
    plt.plot((0, max),(0, max),c='g')

    adjust_ticks_label_size(ax1, 7)
    tmp, chip_vialue = zip(*compareData)
    tmp = map(float, tmp)
    chip_vialue = map(float, chip_vialue)
    print 'tmp:', len(tmp)
    print 'chip_vialue: ', len(chip_vialue)
    plt.scatter(tmp, chip_vialue, s=0.5, c=color, edgecolor='none')
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
    digitalFile = CsvReader(digitalChipExpressFilePath)
    for line in digitalFile.reader():
        if patter.search(line[0]):
            digitalID = line[0].split('.', 1)[0]
            for isoform in patter.findall(line[0]):
                isoformNo = isoform.split(':')[0]
                isoformID = '.'.join((digitalID, isoformNo))
                # print 'isoformID',isoformID
                digitalChipDataDict[isoformID] = float(line[2])
        else:
            digitalID = line[0]
            digitalChipDataDict[digitalID] = float(line[2])
    digitalFile.close()
    return digitalChipDataDict

# def new_digital_chip_id_data(digitalChipExpressFilePath):
#     digitalChipDataDict = {}
#     digitalFile = CsvReader(digitalChipExpressFilePath)
#     for line in digitalFile.reader():
#         digitalChipDataDict[line[0].partition('.')[0]] = float(line[2])
#     digitalFile.close()
#     return digitalChipDataDict

def id_one_to_one(blastoutPath):
    ll = 'idDatall'
    if not path(ll).exists():
        idDict = {}  # {ChipID: digitalChipID}
        for blastReader in parse_blastout(blastoutPath):
            chipID = blastReader.queryId
            if chipID not in idDict:
                digitalID = blastReader.subjectId#.split('.', 1)[0]
                idDict[chipID] = digitalID
                # print digitalID
        # data_dump(idDict, ll)
    else:
        idDict = data_load(ll)

    print 'idDict:', len(idDict)
    print 'invert idDict:', len(invert_dict(idDict))
    return idDict

if __name__ == '__main__':
    main()