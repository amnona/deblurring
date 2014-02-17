#!/usr/bin/env python
"""
Created on Thu Sep 19 17:43:09 2013

@author: amnon
"""

__version__ = "1.2"

import argparse


from os import listdir
from os.path import isfile,join
from numpy import matrix,zeros
import sys
from cogent.parse.fasta import MinimalFastaParser

# prepare an ascii otu table file from all fasta files in a directory
# files are from the matlab CleanSeqs() (in AllSeqs directory)
def PrepareOTUTable(dirname,outfilename,minreads):
    filelist=[f for f in listdir(dirname) if isfile(join(dirname,f))]

    seqs=set()
    ofile=open(outfilename,'w')
    ofile.write('# from noise cleaned fasta files\n')
    
    # write the sample names
    hstr='# '
    totfiles=0;
    for cfile in filelist:
#        if cfile[-17:]=='.fasta.tuni.clean':
        if cfile[-13:]=='.fasta.ref.fa':
            hstr=hstr+'\t'+cfile[:-13]
            totfiles+=1;
    hstr+='\n'
    ofile.write(hstr)

    for cfile in filelist:
#        if cfile(-3:)=='.fa'
        if cfile[-13:]=='.fasta.ref.fa':
            print cfile
            fullname=join(dirname,cfile)
            fafile=open(fullname)
            infile=0
            for seqid,seq in MinimalFastaParser(fafile):
                numseqs=float(seqid[seqid.find(';size=')+6:-1])
                if numseqs>0:
                    seqs.add(seq)
                    infile+=1
            fafile.close()
            print 'total sequences:',len(seqs),'in this file',infile
    
    print "building matrix"
    sdict={}
    slist=[]
    cpos=0;
    for cseq in seqs:
        sdict[cseq]=cpos;
        cpos+=1
        slist.append(cseq)
        
    outmat=zeros((len(seqs),totfiles))
    cfilepos=0;
    for cfile in filelist:
#        if cfile(-3:)=='.fa'
#        if cfile[-6:]=='.clean':
        if cfile[-13:]=='.fasta.ref.fa':
            print cfile
            fullname=join(dirname,cfile)
            fafile=open(fullname)
            infile=0
            for seqid,seq in MinimalFastaParser(fafile):
                numseqs=float(seqid[seqid.find(';size=')+6:-1])
                if numseqs>0:
                    outmat[sdict[seq],cfilepos]+=numseqs
            fafile.close()
            cfilepos+=1
    print "saving"
    numok=0
    for cline in range(len(seqs)):
        if outmat[cline,:].sum()>minreads:
            numok += 1
            ofile.write(slist[cline]+'\t')
            for cnum in outmat[cline,:]:
                ofile.write(str(cnum))
                ofile.write('\t')
            ofile.write('\n')
    ofile.close()


def main(argv):
    parser=argparse.ArgumentParser(description='Make a table from the cleaned reads version '+__version__)
    parser.add_argument('dirname',help='input dir (containing .ref.fa files (from clean_indel.sh)')
    parser.add_argument('outfile',help='output file name (created in the dirname directory)')
    parser.add_argument('-r','--minreads',help='minimal number of reads in order to save to table',default=4)
    args=parser.parse_args(argv)
    PrepareOTUTable(args.dirname,args.dirname+'/'+args.outfile,float(args.minreads))

if __name__ == "__main__":
    main(sys.argv[1:])             
    
#dirname='/Users/amnon/Projects/AllSeqs/bangladesh'
#PrepareOTUTable(dirname,dirname+'/otu-bangladesh-01-hd2.txt')
