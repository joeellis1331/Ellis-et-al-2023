#!/bin/bash
#SBATCH --job-name=matt_rnamaps_cisbp    # Job name
#SBATCH --mail-type=END,FAIL            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=j.ellis@ufl.edu     # Where to send mail
#SBATCH --ntasks=1                      # Run a single task
#SBATCH --cpus-per-task=1               # Number of CPU cores per task
#SBATCH --mem=40gb                      # Job memory request
#SBATCH --time=48:00:00                 # Time limit hrs:min:sec
#SBATCH --output=serial_test_%j.log     # Standard output and error log

### $1=overall run directory, $2=se.mats.jcec, $3=fasta file, $4=maps output directory

motif_pipe() {

#rename chr column to all caps to prevent removal with sed
matt rn_cols $1/temp.SE.MATS_1.txt chr:CHR > $1/temp.SE.MATS_2.txt

#removes chr prefix for matching with fasta file
sed -i 's/chr//' $1/temp.SE.MATS_2.txt

#adds groups to table for motif analysis
matt def_cats $1/temp.SE.MATS_2.txt GROUP_PSI\
 'Inclusion=IncLevelDifference[-1.0,-0.1] FDR[0,0.05]' 'Exclusion=IncLevelDifference[0.1,1.0] FDR[0,0.05]'\
 'Background=IncLevelDifference[-0.05,0.05]' | matt add_cols $1/temp.SE.MATS_2.txt -

#actual finding motifs
matt rna_maps_cisbp $1/temp.SE.MATS_2.txt upstreamEE exonStart_0base exonEnd downstreamES CHR strand\
 GROUP_PSI[Inclusion,Exclusion,Background] 31 50 200 $3 cisbprna_regexps -d $1/$4
}

motif_pipe ./curveFit_352_run ./EventList_curvefits_352_wBkgd.txt $fasta maps_output
