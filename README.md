# CRISPR-Cas9 screening data analysis with mageck 

### Two processing pipelines which can quantify the read count table from raw data (before demultiplexing bam files)

#### CRISPR-UMI from Georg in IMBA

#### from Jesse in Zuber group 
some notes for this nf process pipeline:
- library id should be unique (doulbe check)
- library file: second column should be names 'group'
- barcode file: first column should be the corresponding input files (bam or fastq)
- In the default setup, the read length should be 50bp, 6bp random barcode + 4bp sample barcodes + 20bp spacer + sgRNA seq (with C padding sequence if there are not 20bp length)
- the input should be bam or fq.gz or fastq.gz, not fastq or fq
  
### Comparisons done with mageck
  