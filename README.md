# Ellis-et-al-2023
-Scripts utilized in the publication Ellis-et-al-2023

-genomic GTF and FASTA files were obtained from Ensembl:
  Homo_sapiens.GRCh38.101.gtf
  Homo_sapiens.GRCh38.dna.primary_assembly.fa

-STAR alignment, rMATS, and matt were all run using UF's HiPer Gator HPCC

-AS_overview.py is the script used to generate Supplementary Figure 4

-RNAseq_MBNL1_curvefitting.py is the script used to identify, fit, and plot the 352 events skipped exon events in Supplemental File 1 and Supplemental file 2

-validate_compare_RT_RNAseq.py is the script which plotted the RT-PCR curves on top of RNAseq data for Figure 4B-E and Supplemental File 3. Additionally it conducted     all the regression analyses (Figure 4F-G; Suppl. Fig 5) and generated all swarmplots (Figure 4). See "curveValidation_inputSample_excel.xlsx" for input format


Package versions:
star              2.7.9a
rmats-turbo       4.1.0
seaborn           0.12.0
scipy             1.7.3
python            3.7.13
pandas            1.3.5
numpy             1.19.5
maxentpy          0.0.1
matplotlib        3.5.3
