# Conservation Units Delineation and Predictors of Decline in Bank Swallow 

This repository contains code and workflows used in:

Robertson, E. R., Ruegg, K. C., et al. (2025). *Title of manuscript*. Journal Name.  
DOI: pending / https://doi.org/xxxx

This repository contains scripts and workflows used to process whole-genome sequencing data, analyze population structure and local adaptation, delineate conservation units, manage spatial data, contructs models explaining changes in abundance, and generate figures and tables for the manuscript. The goal is to ensure transparency and reproducibility of all analyses presented in the paper.

Paper Abstract:

## Repository structure and manuscript mapping

<pre>
BANS_adaptive_units
├── README.md          <- What you're currently reading! This has a summary of the structure of the repository and explanations for each section
├── BANS_adaptive_units.Rproj      
├── ~data/             <- The full version of this won't show up here because of space issues. Please reference the README for more information on these files.
│   ├── sample_data/   <- this is the metadata, location data, etc, for my samples
│   ├── spatial/       <- I have really large spatial files, I store these separately so I can ignore them and not have them upload to github (they're too big)
├── analysis/          <- All the workflows and code performed, summarized step by step in .Rmd files. Ordered chronologically, starting with 01.
│   ├── 01.fastq_processing/      <- Worflow for how I processed the raw sequence data
│   │   ├── 01a.fastq_processing.Rmd <- 
│   ├── 02.imputation/      <- Workflow for imputing the vcf before all downstream analysis
│   │   └──02.imputation.Rmd 
│   ├── 02.5.delineate_ESUs/      <- Workflow for delineating evolutionarily significant units (ESUs)
│   │   └──
│   ├── 03.GEA/      <- Workflow for performing the GEA analyses to identify adaptive loci
│   │   ├──
│   │   └──
│   ├── 04.delineate_AUs/      <- Workflow for delineating adaptive units (AUs)
│   │   ├──
│   ├── 05.predictors_of_decline/      <- Workflow for my spatial analysis exploring predictors of changes in abundance
│   │   ├──
│   ├── 06.genomic_offset/      <- Workflow for performing a genomic offset analysis
│   │   ├──
</pre>

## Data availability

Raw sequencing data are available from NCBI SRA under BioProject PRJNAXXXXXX.

Due to size and privacy constraints, the following are NOT included in this repository:
- Raw FASTQ files
- Intermediate BAM files
- Full VCFs (>100 GB)
- Full spatial files
- Statistical analysis output files

# Workflow Step-by-step (expanded methodology)
workflow_diagram.png




