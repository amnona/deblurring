# find region with V4 primers
# $1 - fasta file name (without .fasta!!!)
# $2 - length
~/bin/pprospector-1.0.1/scripts/clean_fasta.py -f $1.fasta
~/bin/pprospector-1.0.1/scripts/analyze_primers.py -f $1_filtered.fasta -p V4r -s GGACTACHVGGGTWTCTAAT
# for mock22,millions,americangut
~/bin/pprospector-1.0.1/scripts/analyze_primers.py -f $1_filtered.fasta -p V4f -s GTGCCAGCMGCCGCGGTAA
# for mock44 ds5 there is an additional A at beginning of sequence
#~/bin/pprospector-1.0.1/scripts/analyze_primers.py -f $1_filtered.fasta -p V4f -s GTGCCAGCMGCCGCGGTA
~/bin/pprospector-1.0.1/scripts/get_amplicons_and_reads.py -f $1_filtered.fasta -i V4f_$1_filtered_hits.txt:V4r_$1_filtered_hits.txt -R $2
~/bin/pprospector-1.0.1/scripts/get_amplicons_and_reads.py -f $1_filtered.fasta -i V4f_$1_filtered_hits.txt:V4r_$1_filtered_hits.txt -R $2 -d f
cp V4f_V4r_r_$2_reads.fasta $1.V4.fa
cp V4f_V4r_f_$2_reads.fasta $1.V4f.fa
~/bin/fastx/fastx_reverse_complement -i V4f_V4r_r_$2_reads.fasta -o $1.V4rc.fa
echo "done"