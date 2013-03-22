#encoding=utf8
from __future__ import division
import sys, re
sys.path.append('/scgene/tiger/invent/guoshuguang/program/')
from michelia import CsvReader, Bio
from Bio import SeqIO
import pdb
b = pdb.set_trace

def main(fa, tag_count):
    '''
    所有基因的所有转录本的 fasta 序列
fa:  /scgene/tiger/invent/chenjiehu/GJXDGE/index/all.cdna.polyA.fa

    tagSeq  所测到的次数。
tag_count : /scgene/tiger/invent/chenjiehu/GJXDGE/tag/11L0_5.tag.exp
或  classify_tag/tag_count_s5.txt
    '''
    tagCountDict = {}
    for reader in CsvReader(tag_count).reader():
        tagCountDict[reader[0]] = reader[1]

    isoformTagsDict= {}
    for seqRcord in SeqIO.parse(fa, 'fasta'):
        isoformID = seqRcord.id
        seq = str(seqRcord.seq)
        isoformTagsDict[isoformID] = []   # 无论是否有tag 都建立[]列表。

        start = 0
        tagList = []
        while True:
            index = seq.find('CATG', start)
            if index > -1:
                tagList.append(seq[index:index+21])
                start = index + 1  #  从下一个碱基开始查找 CATG
            else:
                break

        if len(tagList) == 0:  # 此基因中没有tag
            continue
        site = len(tagList)  # 3`端的第几个tag
        for tag in tagList:
            if tag in tagCountDict: # 并根据测序统计文件进行判断去掉 测序文件中 没有的 tag
                isoformTagsDict[isoformID].append((tag, isoformID, str(site),\
                        tagCountDict[tag]))
                # print '%s\t%s\t%s' % ()
            site -= 1

    # isoformTagsDict  {isoformID:  [(tag, isoformID, site(3`端的第几个tag),测到的次数 )  ]}
    # 一个 基因的转录本 所对应的 多个 tag 的信息
    tagInfoDict = {}
    #  {tag : ['isoformID:位置:测的次数']} 一个 tag可能包含在 一个基因的多个转录本中
                                        # 或 多个基因中
# ('CATGAATAGAGAATCTGCAAA': ['LOC_Os08g35870.1:6:2', 'LOC_Os08g35870.2:6:2', 'LOC_Os08g35910.1:3:2'])
    for isoformID in isoformTagsDict:
        for tagInfo in isoformTagsDict[isoformID]:
            tag = tagInfo[0]
            info = ':'.join(tagInfo[1:]) # isoformID：site(3`端的第几个tag):测到的次数
                                         # LOC_Os08g35870.1:6:2
            tagInfoDict.setdefault(tag, []).append(info)
    # b()

    geneIsoNoTagsDict = {}    #{geneIsoNoID: [tag, ]}
    for tag in tagInfoDict:
        if len(tagInfoDict[tag]) > 1:
            geneIsoNoDict = {}  # 记录 一个tag 可能包含在几个基因中的 几个转录本中
                                 # {geneID: [ isoNO]}
            for info in tagInfoDict[tag]:
                geneID = info.split('.', 1)[0]
                isoNO = info.split('.', 1)[1].split(':', 1)[0]
                geneIsoNoDict.setdefault(geneID, []).append(isoNO)

            if len(geneIsoNoDict) == 1: 
                # 记录出现在一个基因中的tag, 出现在多个 基因的tag, 都去除掉
                #  geneIsoNoDict 中只有一项。
                geneID, isoNoList = geneIsoNoDict.popitem() # 弹出唯一的 项
                if geneID.startswith('LOC'):
                    geneIsoNoID = '%s.[%s]' % (geneID, '/'.join(isoNoList))
                        #  LOC_Os05g31480.[1/2]
                    geneIsoNoTagsDict.setdefault(geneIsoNoID, []).append(tag)
                    # print '%s\t%s' % (geneIsoNoID, tag)
                    #  LOC_Os05g31480.[1/2]    CATGAAGGCTAAATTTACAGT

        else:  # 记录只出现在一个 转录本的 tag
            # tag 只包含在一个转录本中 len(tagInfoDict[tag]) == 1
            # tagInfoDict[tag] 为 ['isoformID:位置:测的次数',]
            tagInfo = tagInfoDict[tag][0]  # 'isoformID:位置:测的次数'
            isoID = tagInfo.split(':')[0]  # LOC_Os05g31480.1
            geneIsoNoID = isoID
            if geneIsoNoID.startswith('LOC'):
                geneIsoNoTagsDict.setdefault(geneIsoNoID, []).append(tag)
                # print '%s\t%s' % (geneIsoNoID, tag)
                # LOC_Os01g56290.1   CATGCCTGTTCAGCTTTCTGC
    for geneIsoNoID in sorted(geneIsoNoTagsDict):
        b()
        tagList = geneIsoNoTagsDict[geneIsoNoID]
        for tag in tagList:
            print '%s\t%s' % (geneIsoNoID, tag)


if len(sys.argv) > 1:
    main(sys.argv[1], sys.argv[2])
else:
    print >>sys.stderr, "python %s <cdna> <tag_dic>" % (sys.argv[0])

