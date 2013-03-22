#encoding=utf8
from __future__ import division
import sys, argparse
sys.path.append('/scgene/tiger/invent/guoshuguang/program/')
from michelia import CsvReader
import pdb
b = pdb.set_trace

def main():
    parser = argparse.ArgumentParser(
                                    )
    parser.add_argument("tag", 
                        help="""the data file. 
                        file format:chromName start end length depth""")
    parser.add_argument('gene_tag_dic', 
                        help='the chrome length file. Format: chromeName  chromeLength')
    args = parser.parse_args()
    # args = parser.parse_args(
    #     ['/scgene/tiger/invent/chenjiehu/GJXDGE/tag/11L0_5.tag', 
    #     '/scgene/tiger/invent/guoshuguang/chip/classify_tag/gene_tag_dic.txt']
    #     )

    tagSeqAmountDict = tag_seq(args.tag)
    geneIsoformDict, allTagAmountSum = isoform_class_tag(args.gene_tag_dic, tagSeqAmountDict)
    # b()
    isoformClassTagAmountDict_per = fun(geneIsoformDict, allTagAmountSum)
    tmp_print(isoformClassTagAmountDict_per, allTagAmountSum)

def tag_seq(cleanTagFilePath):
    '''
cleanTagFilePath 是 /scgene/tiger/invent/chenjiehu/GJXDGE/tag/11L0_5.tag
此数据来自于fastq， 每一行代表一次测序结果，去过接头和无用的序列等等。
    '''
    tagSeqAmountDict = {}  # {tagSeq: tagSeqAmount} 
    for reader in CsvReader(cleanTagFilePath).reader():
        tagSeq = reader[0]  #因为文件中只有一行， 但是取数据时， 还是要用 [0]
        tagSeqAmountDict[tagSeq] = tagSeqAmountDict.get(tagSeq, 0) + 1
    return tagSeqAmountDict#, allTagAmountSum

def isoform_class_tag(geneTagDictFilePath, tagSeqAmountDict):
    '''
    这个函数开始统计次数， 结果中就没有tagSeq了， 只是包含tagSeq的次数。

    classify_tag/gene_tag_dic.txt 结构
    isoformClassID    tagSeq  
    #注意：这里一个类别中可对应了多个tagSeq
    含义是把这个tag归类到哪些isoformClass上。
    所以一类转录本可能对应多个tagSeq， 即一类转录本可出现在文件gene_tag_dic.txt的多个行上。
    只是分成了多行表示。
    geneIsoformDict = {}  
    {geneID(无No): [(isoformClassNo, isoformClassTagAmountSum),]}
    一个基因所对应的 多个 转录本类别
    '''
    isoformClassTagAmountDict = {}  #{isoformClassID: 对应所有tagSeq的总次数}
    allTagAmountSum = 0 # 有用tag的总数
    for reader in CsvReader(geneTagDictFilePath).reader():
        isoformClassID = reader[0]
        tagSeq = reader[1]

        # 无论tagSeq 是否在tagSeqAmountDict中都要 初始化 为0
        isoformClassTagAmountDict.setdefault(isoformClassID, 0)
        if tagSeq in tagSeqAmountDict:
            isoformClassTagAmountDict[isoformClassID] += tagSeqAmountDict[tagSeq]
            allTagAmountSum += tagSeqAmountDict[tagSeq]

    geneIsoformDict = {}  
    # {geneID(无No): [(isoformClassNo, isoformClassTagAmountSum),]}
    # 一个基因所对应的 多个 转录本类别
    for isoformClassID in isoformClassTagAmountDict:
        geneID, isoformClassNo = isoformClassID.split('.')
        # isoformClassNo  一个转录本类别 的 所有 No (一个或多个No)
        isoformClassTagAmountSum = isoformClassTagAmountDict[isoformClassID]

        geneIsoformDict.setdefault(geneID, []).append((isoformClassNo, isoformClassTagAmountSum))
    return geneIsoformDict, allTagAmountSum

def fun(geneIsoformDict, allTagAmountSum):
    # geneIsoformDict
    # {geneID(无No): [(isoformClassNo, isoformClassTagAmountSum),]}
    # 一个基因所对应的 多个 转录本类别
    isoformClassTagAmountDict_per = {}
    for geneID in sorted(geneIsoformDict):

        # 加上了每个isoformNo 对应的的percentage
        # {isoformClassID_per: 对应所有tagSeq的总次数}

        isoformClassList = geneIsoformDict[geneID]
        # isoformClassList
        # 基因所对应的 多个 转录本类别
        # [(isoformClassNo, isoformClassTagAmountSum), ]

        if len(isoformClassList) == 1:
            # 当基因对应一个转录本类别的时候
            isoformClass = isoformClassList[0]
            one_isoform_class(geneID, isoformClass, isoformClassTagAmountDict_per)
        else:
            # 当基因对应 多个 转录本类别的时候
            many_isoform_class(geneID, isoformClassList, isoformClassTagAmountDict_per)

    return isoformClassTagAmountDict_per

def many_isoform_class(geneID, isoformClassList, isoformClassTagAmountDict_per):
    if len(isoformClassList) == 4:
        b()
    # LOC_Os01g04920   [('[1/2/3]', 3), ('3', 0), ('1', 1), ('[1/2]', 352)]
     # LOC_Os01g04920.[1:0.50/2:0.50/3:0.00]   356 107.46
     # LOC_Os01g03020    [('2', 0), ('[1/2]', 0), ('1', 95)]
     # LOC_Os01g03020.1    95  28.68
    # 当基因对应 多个 转录本类别的时候
    isoformNoAmountDict = {}  # {isoformNo: 每个 转录本 所出现的 平均次数}
    #  即 当有一个 转录本类别中 包含 多个 转录本 是 就大于 0
    gene_total_detect = 0  #  统计基因 所对应 多个转录本类别的 多个tag 所有的次数
    for isoformClass in isoformClassList:
        #  isoformClass 其中的一个转录本类别
        isoformClassNo, isoformClassTagAmountSum = isoformClass
        if isoformClassTagAmountSum == 0:
            continue
        if '/' in isoformClassNo:
            # 此转录本类别中包含 多个 转录本
            isoformClassNo = isoformClassNo.\
                replace('[', '').replace(']', '')
            isoformNoList = isoformClassNo.split('/')
            # 转录本类别中 每个 转录本 所出现的 平均次数
            percentage = isoformClassTagAmountSum / len(isoformNoList)
            gene_total_detect += isoformClassTagAmountSum
            for isoformNo in isoformNoList:
                # if geneID == 'LOC_Os01g03020':
                    # b()
                isoformNoAmountDict[isoformNo] = isoformNoAmountDict.get(isoformNo, 0) + percentage
                # isoformNoAmountDict[isoformNo] =  percentage

    if len(isoformNoAmountDict) > 0:
        # 即 此基因对应多个转录本类别， 且其中至少一个 转录本类别中包 含多个 转录本，
        # 则len(isoformNoAmountDict) > 0
        for isoformClass in isoformClassList:
            isoformClassNo, isoformClassTagAmountSum = isoformClass
            if isoformClassTagAmountSum == 0:
                continue
            if '/' not in isoformClassNo:
                # 即 此基因对应多个转录本类别， 且其中至少一个 转录本类别中包 含多个 转录本，
                # 记录 包含一个转录本的 转录本类别
                isoformNo = isoformClassNo

                if isoformNo in isoformNoAmountDict:
                    isoformNoAmountDict[isoformNo] += isoformClassTagAmountSum
                    gene_total_detect += isoformClassTagAmountSum
                else:
                    isoformClassID = '%s.%s' % (geneID, isoformClassNo)
                    isoformClassTagAmountDict_per[isoformClassID] = \
                        isoformClassTagAmountSum
    else:
        # 当len(isoformNoAmountDict)为0时，
        # 既是 这个基因有多个转录本类别，并且所有转录本类别中只包含 一个 转录本 
        # 觉得这里还是有一个bug， 假设上面的成立， 那么就只对其中一个转录本类别(且只有一个转录本)进行下面的操作
        # 如基因 LOC_Os01g03020 对应  LOC_Os01g03020.[1/2]， LOC_Os01g03020.1
        # 结果中只有  LOC_Os01g03020.1    95  28.68
        # 是否存在一个基因对应多个 转录本类别，且其中有多个转录本类别中 只有 一个 转录本
        # b()
        isoformClassID = '%s.%s' % (geneID, isoformClassNo)
        # if geneID == 'LOC_Os01g03020':
            # b()
        isoformClassTagAmountDict_per[isoformClassID] = \
                isoformClassTagAmountSum

    isoformNo_perList = [] 
    for isoformNo in sorted(isoformNoAmountDict):
        isoformNo_per = "%s:%.2f" % (isoformNo, 
                isoformNoAmountDict[isoformNo] / gene_total_detect)
        isoformNo_perList.append(isoformNo_per)

    isoformClassID_per = "%s.[%s]" % (geneID, "/".join(isoformNo_perList))
    isoformClassTagAmountDict_per[isoformClassID_per] = \
            gene_total_detect

def one_isoform_class(geneID, isoformClass, isoformClassTagAmountDict_per):
    '''
# 当基因对应一个转录本类别的时候
    '''
    isoformClassNo, isoformClassTagAmountSum = isoformClass
    if isoformClassTagAmountSum == 0:
        return
    if '/' in isoformClassNo:
        # 此转录本类别中包含 多个 转录本
        isoformClassNo = isoformClassNo.\
                replace('[', '').replace(']', '')
        isoformNoList = isoformClassNo.split('/')
        percentage = 1 / len(isoformNoList)
        isoformNoPerList = []
        for isoformNo in isoformNoList:
            isoformNoPerList.append(
                    '%s:%.2f' % (isoformNo, percentage))
        isoformClassID_per = \
                '%s.[%s]' % (geneID, '/'.join(isoformNoPerList))
        isoformClassTagAmountDict_per[isoformClassID_per] = \
                isoformClassTagAmountSum
    else:
        # 此转录本类别中包含 一个 转录本
        isoformClassID = '%s.%s' % (geneID, isoformClassNo)
        isoformClassTagAmountDict_per[isoformClassID] = \
                isoformClassTagAmountSum

def tmp_print(isoformClassTagAmountDict_per, allTagAmountSum):
    # b()
    print >>sys.stderr, allTagAmountSum
    for isoformClassID_per in sorted(isoformClassTagAmountDict_per):
        tpm = isoformClassTagAmountDict_per[isoformClassID_per] / allTagAmountSum * 1000000
        if tpm > 0:
            print "%s\t%d\t%.2f" % (isoformClassID_per, 
                        isoformClassTagAmountDict_per[isoformClassID_per], tpm)

if __name__ == '__main__':
    main()