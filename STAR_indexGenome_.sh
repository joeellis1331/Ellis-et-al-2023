#!/bin/bash
#SBATCH --job-name=mouse_human_genome_annotate    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=j.ellis@ufl.edu     # Where to send mail
#SBATCH --ntasks=1                   # Run a single task
#SBATCH --cpus-per-task=6            # Number of CPU cores per task
#SBATCH --mem=60gb                     # Job memory request
#SBATCH --time=08:00:00               # Time limit hrs:min:sec
#SBATCH --output=serial_test_%j.log   # Standard output and error log
#SBATCH --qos=berglund-b

pwd; hostname; date

module load star/2.7.9a

#STAR --runMode genomeGenerate --runThreadN 6 --genomeDir ./mm10 --genomeFastaFiles ./FTP_GenesAndSeqFiles/Mus_musculus.GRCm38.dna.primary_assembly.fa --sjdbGTFfile ./FTP_GenesAndSeqFiles/Mus_musculus.GRCm38.101.gtf

#STAR --runMode genomeGenerate --runThreadN 6 --genomeDir ./hg38 --genomeFastaFiles ./FTP_GenesAndSeqFiles/Homo_sapiens.GRCh38.dna.primary_assembly.fa --sjdbGTFfile ./FTP_GenesAndSeqFiles/Homo_sapiens.GRCh38.101.gtf
