# Chiappinelli_Lab/RNAseq/post_processing

Created: 9/15/2023

Updated: 9/15/2023

**This folder will host scripts involved in the post-processing of RNA-seq data. Specifically, this refers to the steps following QC, alignment, and annotation.**

#### _make.telescope.tetranscripts.DESeq2.input.filter.baseMean.10.py_
* This python script was originally created by James McDonald to generate a single combined count table and a DESeq2.R script from multiple telescope report files.
* This is now updated to generate a single combined count table and DESeq2.R script from telescope report files, tetranscripts count tables, or both.
* This script will take in 2 files/lists that specify sample file unique IDs:
* 1. Control Samples (Full path to each file)
  2. Treatment/Experimental Samples (Full path to each file)
  3. Annotation file Path -- Assumes GTF!!! This is the annotation file for the Telescope reports only.
  4. Output directory
  5. Optional input will be a tag for output file names and a flag for handling NA values.
  6. Optional input will be a flag for mode: combined, telescope, or tetranscripts. Default is combined.
* The input files should provide the absolute file path for each Telescope report file and TEtranscripts count table.
* The path for each Telescope report file should end in "-telescope_report.tsv".
* The path for each TEtranscripts count table should end in "-tetranscripts.cntTable".
* Please make sure the sample names match before the hyphen in the file name! This is what the script will use to match the samples.
* To run the script: python make.telescope.tetranscripts.DESeq2.input.filter.baseMean.10.py [control samples directory] "treat samples" "annotation" "output directory" -o "output name" -na "NA value" -mode "combined/telescope/tetranscripts"
