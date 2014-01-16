#!/usr/bin/env python
#!/Users/amnon/anaconda/bin/python
#!/macqiime/bin/python
#!/home/amam7564/qiime_software/python-2.7.3-release/bin/python
"""
Created on Thu Sep 19 17:43:09 2013
V1.1
@author: amnon
"""

import argparse

from os import listdir
from os.path import isfile,join
import sys
import numpy as np
#sys.path.insert(0,'/macqiime/lib/python2.7/site-packages')
from cogent.parse.fasta import MinimalFastaParser


def SeqToArray(seq):
    """ convert a string sequence to a numpy array"""
    seqa=np.zeros(len(seq),dtype=np.int8)
    for ind,base in enumerate(seq):
        if base=='A':
            seqa[ind]=0
        elif base=='a':
            seqa[ind]=0
        elif base=='C':
            seqa[ind]=1
        elif base=='c':
            seqa[ind]=1
        elif base=='G':
            seqa[ind]=2
        elif base=='g':
            seqa[ind]=2
        elif base=='T':
            seqa[ind]=3
        elif base=='t':
            seqa[ind]=3
        elif base=='-':
            seqa[ind]=4
        else:
            seqa[ind]=5
    return(seqa)


def RemoveError(seqnames,seqs,seqsnp,sfreq,readerror,meanerror):
    # we assume all sequences are of equal length
    commonlen=len(seqs[0])
    for cseq in seqs:
        if not(commonlen==len(cseq)):
            print("Not all sequences are same length!!!!")
            print(commonlen)
            print(len(cseq))
            print(cseq)
    print ("processing",len(seqs),"sequences")

    numreal=0
    for cchar in seqs[0]:
        if not (cchar=='-'):
            numreal+=1
    modfactor=pow((1-meanerror),numreal)

    # create the error profile from the read error
    # exponential independent
    fracerr=[]
    for a in range(10):
        fracerr.append(pow(readerror,a)/modfactor)
    # empirical
    fracerr=[1.0/modfactor,pow(readerror,1)/modfactor,2*pow(readerror,2)/modfactor,pow(readerror,2)/modfactor,pow(readerror,2)/modfactor,pow(readerror,2)/modfactor,pow(readerror,2)/modfactor,pow(readerror,2)/modfactor,pow(readerror,2)/modfactor,pow(readerror,2)/modfactor]
    # used for the 22 mock mixture
    fracerr=[1.0/modfactor,pow(readerror,1)/modfactor,0.01,0.01,0.01,0.005,0.005,0.005,0.005,0.005,0.005,0.001,0.001,0.001,0.001,0.001,0.001,0.0005,0.0001,0.0001]
    # used for the 44 mock mixture
    #e1=pow(readerror,1)/modfactor
    #fracerr=[1.0/modfactor,e1,e1/4,e1/5,e1/6,e1/8,e1/10,e1/15,e1/20,e1/30,e1/40,e1/50,e1/50,e1/50,e1/50,e1/50,e1/50,e1/100,e1/500,e1/500]
    print "fracerr"
    print fracerr
    print "readerror"
    print readerror
    print "modfactor"
    print modfactor        

    for idx,cseq in enumerate(seqs):
        csfreq=sfreq[cseq]
        # no need to remove neighbors if freq. is <=0
        if csfreq<=0:
            continue
        # correct for the fact that many reads are expected to be mutated
        numerr=[] 
        for a in range(len(fracerr)):
            numerr.append(fracerr[a]*csfreq)

        # if it's low level, just continue
        if numerr[1]<0.1:
            continue

        # compare to all other sequences and calculate hamming dist
        cseqnp=seqsnp[idx]
        oseqlen=len(seqs[idx].rstrip('-'))
        for idxtmp,seqnptmp in enumerate(seqsnp):
            # don't compare to ourselves (dist=0)
            if idxtmp==idx:
                continue
            # calculate the hamming distance
            hdist=np.count_nonzero(np.not_equal(seqnptmp,cseqnp))
            # if far away, don't need to correct
            if hdist>14:
                continue
            # close, so lets calculate exact distance
            numsub=0
            numindel=0
            for cpos in range(oseqlen):
                if not (cseqnp[cpos]==seqnptmp[cpos]):
                    if seqnptmp[cpos]=='-':
                        numindel+=1
                    else:
                        numsub+=1
            nerr=numerr[numsub]

            # remove errors due to (PCR?) indels (saw in 22 mock mixture)
            if numindel>0:
                nerr=nerr*0.01
            if numindel>3:
                nerr=0

            # if the effect is small - don't do anything
            if nerr<0.1:
                continue
            # met all the criteria - so correct the frequency of the neighbor
            sfreq[seqs[idxtmp]]-=nerr
    return(sfreq)




def CleanSeqs(dirname,readerror,nomod):
    print "Cleaning"
# prepare the file list in the directory
    filelist=[f for f in listdir(dirname) if isfile(join(dirname,f))]

# loop over the files
    for cfile in filelist:
        if cfile[-5:]!='.tuni':
            continue
        sfreq={}
        seqs=[]
        seqnames=[]
        seqsnp=[]
        sortfreq=[]
        print cfile
        fafile=open(join(dirname,cfile))
        for seqid,seq in MinimalFastaParser(fafile):
            # get the number or reads from the header string
            numseqs=float(seqid[seqid.find(';size=')+6:-1])
            
            # convert sequence to a numpy array
            seqa=SeqToArray(seq)
            # hash the number of reads
            sfreq[seq]=numseqs
            # and store the list of sequences (needs to be sorted)
            # we use a hash for the frequencies and a numpy array for the hamming comparisons
            seqs.append(seq)
            seqnames.append(seqid)
            seqsnp.append(seqa)
            # store the ordered list of frequencies for sorting
            sortfreq.append(numseqs)
        
        if len(seqs)==0:
            print("No sequences in file: "+cfile)
        else:
            # sort the sequences according to frequencies (descending)
            sortorder=sorted(range(len(sortfreq)),key=sortfreq.__getitem__, reverse=True)
            seqs=[seqs[i] for i in sortorder]
            seqnames=[seqnames[i] for i in sortorder]
            seqsnp=[seqsnp[i] for i in sortorder]

            # after loading the file - remove the read errors (MAIN FUNCTION)
            print "orig readerror",readerror
            cfreq=RemoveError(seqnames,seqs,seqsnp,sfreq,readerror,nomod)
            # and finally save the new fasta as a '.clean' file
            ofile=open(join(dirname,cfile+'.clean'),'w')
            for cseq in seqs:
                if round(cfreq[cseq])>0:
                    ofile.write('>aa;size='+str(int(round(cfreq[cseq])))+';\n')
                    csequp=cseq.upper()
                    for a in range(len(csequp)):
                        if not (csequp[a]=='-'):
                            ofile.write(csequp[a])
                    ofile.write('\n')
            ofile.close()
    

def main(argv):
    parser=argparse.ArgumentParser(description='Clean read errors from illumina reads')
    parser.add_argument('dirname',help='input dir (containing .tuni files unique and truncated)')
    parser.add_argument('-e','--readerror',help='readerror rate',default=0.05)
    parser.add_argument('-m','--meanerror',help='the mean error, used for original sequence estimate (default same as readerror)',default=-1)
    args=parser.parse_args(argv)
    if args.meanerror==-1:
        args.meanerror=args.readerror
    CleanSeqs(args.dirname,float(args.readerror),float(args.meanerror))
    
if __name__ == "__main__":
    main(sys.argv[1:])                