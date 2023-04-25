#!/bin/bash
#SBATCH --job-name=ddHEK_coreg_align    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=j.ellis@ufl.edu     # Where to send mail
#SBATCH --ntasks=1                   # Run a single task
#SBATCH --cpus-per-task=6            # Number of CPU cores per task
#SBATCH --mem=100gb                     # Job memory request
#SBATCH --time=24:00:00               # Time limit hrs:min:sec
#SBATCH --output=serial_test_%j.log   # Standard output and error log

module load star/2.7.9a

FILES="./fastqFiles/*"
counter=0

for i in $FILES; do
        if [[ $counter -eq 0 ]]; then
                f1=$i
                filename=$(basename "$i" .fastq)
                IFS='S'
                read -a splitname <<< "$filename"
                counter=1
                continue
        else
                f2=$i
                prefix="${splitname[0]::-1}"
                STAR --runThreadN 6 --readFilesCommand zcat --genomeDir ./hg38 --readFilesIn "$f1" "$f2" --outFileNamePrefix ./STAR/"$prefix" --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --alignEndsType EndToEnd
                counter=0
        fi
done
