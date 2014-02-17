#!/usr/bin/env python

"""
ParseUParse.py
Take from the output of uparse the otu freq. list and centroids fasta
and produce a fasta with the ;size=XXX field
@author: amnon
"""

__version__ = "0.9"

import argparse

from os import listdir
from os.path import isfile,join
import sys
import numpy as np
from cogent.parse.fasta import MinimalFastaParser


def ParseUParse(freqfilename,seqfilename,outfilename):
	allfreqs={}
	# load the otu freq file and populate the hash
	freqdat=open(freqfilename)
	# skip the header
	freqdat.readline()
	for cline in freqdat:
		cline=cline.split()
		allfreqs[cline[0]]=cline[1]

	# load the sequences
	seqfile=open(seqfilename)
	outfile=open(outfilename,'w')
	for seqid,seq in MinimalFastaParser(seqfile):
		newseqid=seqid
		if allfreqs.has_key(seqid):
			newseqid=seqid+';size='+allfreqs[seqid]
		else:
			newseqid=seqid+';size=1'
		outfile.write('>'+newseqid+'\n')
		outfile.write(seq+'\n')
	outfile.close()


def main(argv):
    parser=argparse.ArgumentParser(description='add size annotations to uparse fasta '+__version__)
    parser.add_argument('-f','--freqfile',help='uparse otu frequency file (typically otus.txt)')
    parser.add_argument('-s','--seqfile',help='uparse centroid sequence file (typically XXX.otu.fasta)')
    parser.add_argument('-o','--output',help='output fasta file name')
    args=parser.parse_args(argv)
    ParseUParse(args.freqfile,args.seqfile,args.output)
    
if __name__ == "__main__":
    main(sys.argv[1:])                
