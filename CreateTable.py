#!/usr/bin/env python

"""
CreateTable 
Create a biom like table from fasta files
Input - a list of fasta files (first one for names)
and name for the output files (.table.txt and .seq.txt)
Output - 2 files - a list of the sequences (.seq.txt) and 
a tab delimited text file of appearance of each sequence (row)
in each fasta file (column)
@author: amnon
"""

__version__ = "0.9"

import argparse

from os import listdir
from os.path import isfile,join,basename
import sys
import numpy as np
from cogent.parse.fasta import MinimalFastaParser


def CreateTable(fastanames,output,header):
	allseqs={}
	allfreqs={}
	# load all the sequences
	cseqnum=1
	for cfilename in fastanames:
		cfile=open(cfilename)
		allfreqs[cfilename]={}
		for seqid,seq in MinimalFastaParser(cfile):
			seq=seq.upper()
			# if we have a header line keep only the number os OTU name
			if header:
				allseqs[seq]=str(cseqnum)
			else:
				allseqs[seq]=str(cseqnum)+'-'+seqid
			cseqnum+=1
			# need to modify for frequency
			try:
				numseqs=float(seqid[seqid.find(';size=')+6:-1])
			except:
				numseqs=1
			if allfreqs[cfilename].has_key(seq):
				allfreqs[cfilename][seq]+=numseqs
			else:
				allfreqs[cfilename][seq]=numseqs
		cfile.close()

	# now write the table 
	outfile=open(output+'.table.txt','w')
	outfileseq=open(output+'.seq.fa','w')
	if header:
		outfile.write('OTUID')
		for cfilename in fastanames:
			cfilename=basename(cfilename)
			# need to do it more elegantly
			if cfilename[-13:]=='.fasta.ref.fa':
				outfile.write('\t'+cfilename[:-13])
			else:
				outfile.write('\t'+cfilename)
		outfile.write('\n')
	for seq,seqname in allseqs.items():
		outfileseq.write('>'+seqname+'\n')
		outfileseq.write(seq+'\n')
		outfile.write(seqname)
		for cfilename in fastanames:
			cfreq=0
			if allfreqs[cfilename].has_key(seq):
				cfreq=allfreqs[cfilename][seq]
			outfile.write('\t'+str(cfreq))
		outfile.write('\n')
	outfile.close()
	outfileseq.close()


def main(argv):
    parser=argparse.ArgumentParser(description='Create a frequency table from fasta files version'+__version__)
    parser.add_argument('-i','--input',nargs='*', dest='fasta', action='append',help='names of fasta files to use')
    parser.add_argument('-d','--indir',help='name of input directory containing .ref.fa files to use (instead of -i)',default='')
    parser.add_argument('-l','--headerline',help='include header line (for biom table compatibility',action='store_true')
    parser.add_argument('-o','--output',help='output file name (creates 2 files - .seq.fa and .table.txt)')

    args=parser.parse_args(argv)
    # if we have an input directory use it
    if len(args.indir)>0:
    	fasta=[]
    	filelist=[f for f in listdir(args.indir) if isfile(join(args.indir,f))]
        for cfile in filelist:
            if cfile[-7:]=='.ref.fa':
                fasta.append(join(args.indir,cfile))
    else:
    # otherwise use the file list
    	fasta=args.fasta[0]
    print(fasta)
    CreateTable(fasta,args.output,args.headerline)
    
if __name__ == "__main__":
    main(sys.argv[1:])                
