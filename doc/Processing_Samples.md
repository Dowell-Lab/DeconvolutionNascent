# Processing of Nascent Samples

Sequencing samples were processed using the Bidirectional-Flow pipeline, ([https://github.com/Dowell-Lab/Bidirectional-Flow]()) which is an extension of the Nascent-Flow pipeline ([https://github.com/Dowell-Lab/Nascent-Flow]()).
In combination these two NextFlow pipelines process nascent sequencing data from fastq files to de-novo bidirectional calls.
The documentation for these pipelines is comprehensive, and we briefly summarize the steps used here:

## Nascent Flow Steps
1. Pre-processing:
- FastQC (pre-trim) -- perform pre-trim FastQC on fastq files
1. Trimming & Mapping
- BBDuk -- trim fastq files for quality and adapters
- FastQC (post-trim) -- perform post-trim FastQC on fastq files (ensure trimming performs as expected)
- HISAT2 -- Map reads to a reference genome
1. SAMtools -- convert SAM file to BAM, index BAM, flagstat BAM, BAM --> CRAM, index CRAM
1. Quality control
- preseq -- estimate library complexity
- RSeQC -- calculate genomic coverage relative to a reference file, infer experiement (single- v. paired-end), read duplication
- Pileup.sh : BBMap Suite -- genomic coverage by chromosome, GC content, pos/neg reads, intron/exon ratio
- Picard -- Mark duplicates, GC content
- NQC -- Perform nascent-specific quality control analysis
1. Mapping Visualization
- BEDTools : non-normalized & nornmalized bedgraphs
- BEDTools and kentUtils : 5' bigwigs for dREG & normalized bigwigs for genome browser
- IGV Tools : bedGraph --> tdf

## Bidirectional Flow Steps

1. Preprocessing files for TFit and dREG
1. Run TFit and dREG
1. Count reads over genes and bidirectionals
