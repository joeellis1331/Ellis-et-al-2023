# Scripts utilized in the publication Ellis-et-al-2023

- STAR alignment, rMATS, and matt were all run using UF's HiPer Gator HPCC

- Genomic GTF and FASTA files were obtained from Ensembl:
  - Homo_sapiens.GRCh38.101.gtf
  - Homo_sapiens.GRCh38.dna.primary_assembly.fa

- AS_overview.py
  - Used to generate Supplementary Figure 4

- RNAseq_MBNL1_curvefitting.py
  - Used to identify, fit, and plot the 352 events skipped exon events in Supplemental File 1 and Supplemental file 2

- validate_compare_RT_RNAseq.py
  - Plotted the RT-PCR curves on top of RNAseq data for Figure 4B-E and Supplemental File 3
  - Conducted all the regression analyses (Figure 4F-G; Suppl. Fig 5)
  - generated all swarmplots (Figure 4)
    - See "curveValidation_inputSample_excel.xlsx" for input format.

## Package versions:
- [star](https://github.com/alexdobin/STAR) = 2.7.9a
- [rmats-turbo](https://github.com/Xinglab/rmats-turbo) = 4.1.0
- [seaborn](https://github.com/mwaskom/seaborn) = 0.12.0
- [scipy](http://www.scipy.org/) = 1.7.3
- [python](https://www.python.org/) = 3.7.13
- [pandas](https://pandas.pydata.org/) = 1.3.5
- [numpy](https://numpy.org/) = 1.19.5
- [maxentpy](https://github.com/kepbod/maxentpy) = 0.0.1
- [matplotlib](https://matplotlib.org/) = 3.5.3
