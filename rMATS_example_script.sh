#!/bin/bash
#SBATCH --job-name=ddHEK_coreg_rMATS    # Job name
#SBATCH --mail-type=END,FAIL            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=j.ellis@ufl.edu     # Where to send mail
#SBATCH --ntasks=1                      # Run a single task
#SBATCH --cpus-per-task=1               # Number of CPU cores per task
#SBATCH --mem=10gb                      # Job memory request
#SBATCH --time=08:00:00                 # Time limit hrs:min:sec
#SBATCH --output=serial_test_%j.log     # Standard output and error log

module load mats

rmats.py --b1 ./b1_b2_files/dox0.txt --b2 ./b1_b2_files/dox100.txt --gtf Homo_sapiens.GRCh38.101.gtf\
 --od ./results/d0_vs_d100_HEK --tmp ./temp/temp_d0_d100_HEK -t paired --readLength 101 --libType fr-firststran
 
rmats.py --b1 ./b1_b2_files/dox0.txt --b2 ./b1_b2_files/dox170.txt --gtf Homo_sapiens.GRCh38.101.gtf\
 --od ./results/d0_vs_d170_HEK --tmp ./temp/temp_d0_d170_HEK -t paired --readLength 101 --libType fr-firststrand

rmats.py --b1 ./b1_b2_files/dox0.txt --b2 ./b1_b2_files/dox250.txt --gtf Homo_sapiens.GRCh38.101.gtf\
 --od ./results/d0_vs_d250_HEK --tmp ./temp/temp_d0_d250_HEK -t paired --readLength 101 --libType fr-firststrand

rmats.py --b1 ./b1_b2_files/dox0.txt --b2 ./b1_b2_files/dox280.txt --gtf Homo_sapiens.GRCh38.101.gtf\
 --od ./results/d0_vs_d280_HEK --tmp ./temp/temp_d0_d280_HEK -t paired --readLength 101 --libType fr-firststrand

rmats.py --b1 ./b1_b2_files/dox0.txt --b2 ./b1_b2_files/dox1000.txt --gtf Homo_sapiens.GRCh38.101.gtf\
 --od ./results/d0_vs_d1000_HEK --tmp ./temp/temp_d0_d1000_HEK -t paired --readLength 101 --libType fr-firststrand

