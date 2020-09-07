# A script to download all the transcript sequences from NCBI Refseq

for index in `seq 9`
do
wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot/human.${index}.rna.fna.gz
done
