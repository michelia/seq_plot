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
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from collections import namedtuple
import argparse

function ='''
Acording the data that jidhu'code generate, to plot the fastq quality figure.
'''


def main():
    parser = argparse.ArgumentParser(description=function)
    parser.add_argument("dataFilePath",
                        help="")
    parser.add_argument("savePath",
                        help="picture savePath")
    # args = parser.parse_args()
    args = parser.parse_args()

    dataFilePath = args.dataFilePath
    Data = namedtuple('Data', 'x, y, z')
    # data = Data(*zip(*filter(f, perse_date(dataFilePath))))
    data = Data(*zip(*perse_date(dataFilePath)))

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    aScatter = ax.scatter3D(data.x, data.y, data.z, s=20, zdir='z', edgecolor='none', c='r')
    ax.set_xlabel('RBQA')
    ax.set_ylabel('RBQ')
    ax.set_zlabel('GC')
    ax.set_title('')
    ax.set_xlim(0, 3)
    ax.set_ylim(0, 3)
    ax.set_zlim(0, 1)
    plt.show()
    # plt.savefig(args.savePath, dpi=150)

def f(i):
        if i[0] > 3:
           return 
        if i[1] > 3:
           return 
        # if i[2] > 3:
        #   return 
    # if i[0] <= 3 and i[1] <= 3 and i[2] <= 3:
        return 1

def perse_date(dataFile):
    with open(dataFile) as data:
        for line in data:
            if line[0] != '#' and line.strip():
                yield map(float,line.split())

if __name__ == '__main__':
    main()