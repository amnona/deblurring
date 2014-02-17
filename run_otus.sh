# run the otu pipeline for simulations
# input:
# $1 - the fasta file name
# output:
# $1.otus-97-135.fa
/home/amam7564/scripts/FastaToOneSample.py -i $1 -o $1.forotus.fa
pick_closed_reference_otus.py -i $1.forotus.fa -r /home/amam7564/data/gg/13-5/97_otus.fasta -o otus-97-135-$1 -t /home/amam7564/data/gg/13-5/97_otu_taxonomy.txt -f -a -O 50
/home/amam7564/scripts/FastaFromBiom.py -b otus-97-135-$1/otu_table.biom -g /home/amam7564/data/gg/13-5/97_otus.fasta -o $1.otus-97-135.fa
