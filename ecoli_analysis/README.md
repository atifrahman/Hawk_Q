# Ampicillin resistance mapping in E. coli

## links.txt

Contains SRR IDs of short reads from E.coli strains used for mapping ampicillin resistance.

--We don't need this for q_analysis because we already have the dataset in fasta format.

## download

This will download sra files into the directory ~/ncbi/sra

--We don't need this either because we already have fasta files. However, we need to convert them to fastq format.

## countKmers

This will convert sra files into fastq files, count k-mers and 
write the names of sorted k-mer count files in 'sorted_files.txt' 
and total k-mer count in samples in 'total_kmer_counts.txt'.

## gwas_info.txt

GWAS info file for for mapping ampicillin resistance in E. coli.

-- This will now contain some values instead of case/control. This file has not been made yet.

