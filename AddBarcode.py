#!/usr/bin/env python

"""
AddBarcode.py
add the ;bacrcodelabel=XXX to a fasta file for uparse
"""

__version__ = "0.9"

import argparse

from os import listdir
from os.path import isfile,join
import sys
import numpy as np
from cogent.parse.fasta import MinimalFastaParser


def AddBarcode(inputname,outputname,barcode):
	ifile=open(inputname,'r')
	ofile=open(outputname,'w')
	for seqid,seq in MinimalFastaParser(ifile):
		ofile.write('>'+seqid+';barcodelabel='+barcode+';\n')
		ofile.write(seq+'\n')
	ofile.close()


def main(argv):
    parser=argparse.ArgumentParser(description='add the ";barcodelabel=XXX" field for uparse '+__version__)
    parser.add_argument('-i','--input',help='input fasta file (without labels)')
    parser.add_argument('-o','--output',help='output file name (with ;barcodelabel= added)')
    parser.add_argument('-b','--barcode',help='barcode label to add',default='Sample1')
    args=parser.parse_args(argv)
    AddBarcode(args.input,args.output,args.barcode)
    
if __name__ == "__main__":
    main(sys.argv[1:])                
