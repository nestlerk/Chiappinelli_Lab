################### Description ###################
# This script will merge TEcount and Telescope tables together. This is necessary for differential expression analysis of locus-specific (Telescope) TE expression
# using DESeq2.
# This script requires the following arguments:
# -telescope: Path to telescope tsv file (output from make.telescope.DESeq2.input.filter.baseMean.10.py script)
# -tecount: Path to a folder containing all of the tecount tables you want to merge. You do not need to specify each table individually, just the path to the folder.
# -out: Path to the output table
# Please note:
# The script will also add "ls_" to the beginning of all the Telescope transcript names to indicate they are locus specific. 
# This will help after DESeq2 analysis to distinguish the locus-specific transcripts, but be sure to remove the "ls_" before visualizing the data.
# The script also assumes the TEcount table headers have the sample name directly before the capital "A" from "Aligned.out.bam" in the header.
# You will need to remove the "transcript" column from the output table before running DESeq2.

################### Usage ###################
# Load our required packages
import argparse
import pandas as pd
import os
import glob
# Parse arguments
parser = argparse.ArgumentParser(description='Combine TEcount and Telescope tables')
parser.add_argument('-telescope', help='Path to Telescope tsv file', required=True)
parser.add_argument('-tecount', help='Path to folder containing TEcount tables', required=True)
parser.add_argument('-out', help='Output table', required=True)
args = parser.parse_args()

################### Read TEcount tables ###################
# Read in all of the tecount tables that exist in the directory specified
tecount = pd.DataFrame()
for file in glob.glob(args.tecount + "/*.cntTable"): # Loop through all the tecount tables
    df = pd.read_csv(file, sep='\t', index_col=0) # Read in the table
    tecount = pd.concat([tecount, df], axis=1) # Concatenate the tables
# Reset the index
tecount.reset_index(inplace=True)
# Rename the first column (gene/TE) to transcript
tecount.columns.values[0] = 'transcript'
# Modify the column names to include the sample name before the capital A
tecount.columns = tecount.columns.str.split("A").str[0]

################### Read Telescope tables ###################
# Read in telescope tables
telescope = pd.read_csv(args.telescope, sep='\t')
# Name the first column transcript
telescope.columns.values[0] = 'transcript'
# Set the transcript column as the index
telescope.set_index('transcript', inplace=True)
# Make all the float values integers
telescope = telescope.astype(int)
# Reset the index
telescope.reset_index(inplace=True)
# Modify the column names to include the sample name before the hyphen
telescope.columns = telescope.columns.str.split("-").str[0]
# Add "ls_" to the beginning of all the transcript names to indicate they are locus specific.
telescope['transcript'] = telescope['transcript'].apply(lambda x: 'ls_' + x)

################### Merge the two tables ###################
# Merge the two tables
combined = pd.concat([tecount, telescope], axis=0)
# Write the table to a tsv file
combined.to_csv(args.out, sep='\t', index=False)
