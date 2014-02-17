# clean phix sequences from cleaned reads
# input:
# $1 - filename to be cleaned
# $2 - output filename
# output:
# $1.nophix.fa

# remove all matches to phix
/home/amam7564/bin/usearch7 -search_global $1 -notmatched $1.tmp.nophix.fa -id 0.95 -db /home/amam7564/data/phix/PhiX.fasta -strand both
# and reformat to 1 line
/home/amam7564/bin/fastx/fasta_formatter -i $1.tmp.nophix.fa -o $2
# remove the tmp file
rm $1.tmp.nophix.fa
