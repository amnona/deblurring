deblurring
==========
(still in alpha!)

To deblur a directory (with .fasta files), run:
clean_indel.sh dirname readlen maxreaderror meanreaderror errorprofile indelprob indelmax

recommended values:
maxreaderror - 0.02
meanreaderror - 0.005
errorprofile - 1,0.06,0.02,0.02,0.01,0.005,0.005,0.005,0.001,0.001,0.001,0.0005
indelprob - 0.01
indelmax - 3

note that the directory should contain a single fasta for each sample - so it is the result of split_fasta_on_sample_ids.py

Another option is to run in parallel on a directory using:
CleanIndelDirParallel.py dirname -l readlen -e 0.02 -n 50 -m 0.005 -d 1,0.06,0.02,0.02,0.01,0.005,0.005,0.005,0.001,0.001,0.001,0.0005
the -n gives the number of parallel processes to use

note that this curretnly works with qsub!

The main deblurring python script is:
CleanSeqsIndel.py



