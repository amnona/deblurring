#run the uparse pipeline and add size annotations to fasta file
# input
# $1 - the original fasta file name
# $2 - the read length to trim to
# final output:
# $1.uparse.fa - otu sequences with size annotations (;size=XXX)

#trim
/home/amam7564/bin/fastx/fastx_trimmer -i $1 -o $1.up.trimp1.fa -t 1 -m $2
/home/amam7564/bin/fastx/fastx_trimmer -i $1.up.trimp1.fa -o $1.up.trim.fa -l $2

/home/amam7564/scripts/AddBarcode.py -i $1.up.trim.fa -o $1.up.bc.fa
/home/amam7564/bin/usearch7 -derep_fulllength $1.up.bc.fa -minuniquesize 2 -output $1.up.derep.fa -sizeout -threads 1
/home/amam7564/bin/usearch7 -sortbysize $1.up.derep.fa -output $1.up.sorted.fa -minsize 2
/home/amam7564/bin/usearch7 -cluster_otus $1.up.sorted.fa -otus $1.up.otus1.fa
/home/amam7564/bin/usearch7 -uchime_ref $1.up.otus1.fa -db /home/amam7564/data/gold/gold.fa -strand plus -nonchimeras $1.up.otus2.fa

#/home/amam7564/bin/usearch7 -uchime_denovo otus1.fa -nonchimeras otus2.fa
# Label OTU sequences OTU_1, OTU_2...
python /home/amam7564/bin/uparse/fasta_number.py $1.up.otus2.fa OTU_ > $1.up.otus.fa
# Map reads (including singletons) back to OTUs
/home/amam7564/bin/usearch7 -usearch_global $1.up.bc.fa -db $1.up.otus.fa -strand plus -id 0.97 -uc $1.up.map.uc
# Create OTU table
python /home/amam7564/bin/uparse/uc2otutab.py $1.up.map.uc > $1.up.otu_table.txt

# join to a fasta file with size=XXX added
/home/amam7564/scripts/ParseUParse.py -f $1.up.otu_table.txt -s $1.up.otus.fa -o $1.uparse.fa
