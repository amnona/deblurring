#!/bin/bash
# clean read errors from illumina fasta
# input:
# $1 - directory to clean
# $2 - the read length to trim to
# $3 - the read error upper bound
# $4 - the average read error
# $5 - the error profile (or 0 for default)
# $6 - indel probability (0.01 for illumina)
# $7 - indel maximal number (3 for illumina)

START=$(date +%s)

FILES=$1/*.fasta
echo "Trimming and de-replicating files"
echo "in $1"
for f in $FILES
do
	echo "Processing file $f..."
	# trim the sequences and throw short reads
	echo "trimming"
	/home/amam7564/bin/fastx/fastx_trimmer -i $f -o $f.trimp1.fa -t 1 -m $2
	/home/amam7564/bin/fastx/fastx_trimmer -i $f.trimp1.fa -o $f.trim.fa -l $2
	# remove duplicates and count, also remove singletons
	echo "dereplicating"
	/home/amam7564/bin/usearch7 -derep_fulllength $f.trim.fa -minuniquesize 2 -output $f.ptuni -sizeout -threads 1
	# do multiple sequence alignment
	echo "performing multiple sequence alignment"
	/home/amam7564/bin/mafft/bin/mafft --quiet --parttree --auto $f.ptuni > $f.tuni
#	/home/amam7564/bin/mafft/bin/mafft --quiet --retree 2 --maxiterate 10 $f.ptuni > $f.tuni
done

echo "Cleaning directory $1 . readlength = $2 . Errorrate=$3"
# need to add the mean error as $4
/home/amam7564/scripts/CleanSeqsIndel.py $1 -e $3 -m $4 -d $5 -i $6 --indelmax $7 --pyroseq

echo "removing chimeras and reformatting"
for f in $FILES
do
	echo "Processing $f file..."
	/home/amam7564/bin/usearch7 -uchime_denovo $f.tuni.clean -nonchimeras $f.tuni.clean.noch.fa
	/home/amam7564/bin/fastx/fasta_formatter -i $f.tuni.clean.noch.fa -o $f.ref.t.fa
	/home/amam7564/scripts/clean_phix.sh $f.ref.t.fa $f.ref.fa
done

# add the marker file for the CleanDirParallel.py so it knows we finished
echo "done" > "$1/process.finished"

echo "done"

END=$(date +%s)
DIFF=$(( $END - $START ))
echo "It took $DIFF seconds"
