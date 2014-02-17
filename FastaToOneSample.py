#!/usr/bin/env python

"""
FastaToOneSample.py
Prepare the headers of the fasta file to agree with qiime pick_closed_reference_otus
So we need in the title of each sequence a constant string+'_'+running number
"""

__version__ = "0.9"

import argparse

from os import listdir
from os.path import isfile,join
import sys
import numpy as np
from cogent.parse.fasta import MinimalFastaParser


def FastaToOneSample(inputname,outputname,sampleid):
    ifile=open(inputname,'r')
    ofile=open(outputname,'w')
    cseq=1
    for seqid,seq in MinimalFastaParser(ifile):
        ofile.write('>'+sampleid+'_'+str(cseq)+';\n')
        ofile.write(seq+'\n')
        cseq+=1
    ofile.close()


def main(argv):
    parser=argparse.ArgumentParser(description='change the headers of the fasta so it will be one sample in qiime '+__version__)
    parser.add_argument('-i','--input',help='input fasta file')
    parser.add_argument('-o','--output',help='output file name')
    parser.add_argument('-s','--sampleid',help='the sampleid to add to each sequence',default='Sample1')
    args=parser.parse_args(argv)
    FastaToOneSample(args.input,args.output,args.sampleid)
    
if __name__ == "__main__":
    main(sys.argv[1:])                
