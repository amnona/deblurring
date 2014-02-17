#!/usr/bin/env python

"""
FastaFromBiom
Create a fasta file from a biom table
using the greengenes dataset
@author: amnon
"""

__version__ = "0.9"

import argparse

from biom.parse import parse_biom_table
from biom.util import biom_open

from os import listdir
from os.path import isfile,join
import sys
import numpy as np
from cogent.parse.fasta import MinimalFastaParser


def FastaFromBiom(biomname,ggname,colnum,outputname):
	allfreqs={}
	# load the biom table
	biom_table = parse_biom_table(biom_open(biomname,'U'))
	for values, ids, metadata in biom_table.iterObservations():
		if colnum==0:
			allfreqs[ids]=np.mean(values)
		else:
			allfreqs[ids]=values[colnum]

	print(values)
	ggfile=open(ggname,'r')
	outfile=open(outputname,'w')
	# scan greengenes for the otuids and write sequences
	for seqid,seq in MinimalFastaParser(ggfile):
		if allfreqs.has_key(seqid):
			if allfreqs[seqid]>0:
				outfile.write('>'+seqid+';size='+str(allfreqs[seqid])+'\n')
				outfile.write(seq+'\n')
			else:
				print('freq=0 for seqid '+seqid)
	outfile.close()


def main(argv):
    parser=argparse.ArgumentParser(description='Create a frequency table from fasta files version '+__version__)
    parser.add_argument('-b','--biom',help='biom table file')
    parser.add_argument('-g','--greengenes',help='greengenes fasta file containing the ggids in the biom table')
    parser.add_argument('-c','--column',help='the column to use for writing the fasta file',default=0)
    parser.add_argument('-o','--output',help='output fasta file name')
    args=parser.parse_args(argv)
    FastaFromBiom(args.biom,args.greengenes,args.column,args.output)
    
if __name__ == "__main__":
    main(sys.argv[1:])                
