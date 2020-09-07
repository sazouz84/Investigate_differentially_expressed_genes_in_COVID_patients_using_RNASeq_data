# Installing Bioconda if it's not available 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")


BiocManager::install("polyester")
BiocManager::install("Biostrings")

library(polyester)
library(Biostrings)

setwd("./src")
fasta_file = system.file('extdata', './fasta_ref.fa', package='polyester')

# read the big compressed fasta file (it can also read plain text fasta)
fasta =  readDNAStringSet("human.5.rna.fna.gz")

print(paste("Number of sequences in the FASTA file"))
print(length(fasta))

# We want to subset only a handful of sequences from the fasta file 
n <- 100 # here we take only the first 100 sequences from the file if you want to inlude them all just use the 'fasta' object in the further steps 

# subset the FASTA file to first 20 transcripts
small_fasta = fasta[1:n]
writeXStringSet(small_fasta, 'subset_genes.fa')

# generate random reads per transcript to get different values of FPKM
max_read <- 20 
n_groups <- 2
n_samples <- 1  # how many samples you want to simulate 
max_fold_change <- 5
coverage <- sample.int(max_read, n, replace = TRUE)
# ~20x coverage ----> reads per transcript = transcriptlength/readlength * 20
# here all transcripts will have ~equal FPKM
readspertx = round(coverage * width(small_fasta) / 100)
matrix(sample.int(max_fold_change, size = n, replace = TRUE), nrow = n, ncol = n_samples) 

# Build a matrix of random fold change for each n_samples
fold_changes <- matrix(, nrow = n, ncol = 0)
for(i in 1:n_groups) {
  fold_changes <- cbind(fold_changes, sample.int(max_fold_change, size = n, replace = TRUE))
  print(fold_changes)
}

# simulation call:
simulate_experiment('subset_genes.fa', reads_per_transcript=readspertx, 
                    num_reps=c(1,1), fold_changes=fold_changes, outdir='simulated_reads') 



