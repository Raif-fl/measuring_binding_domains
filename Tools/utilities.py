from plastid import BAMGenomeArray, Transcript
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import csv
import math

def ACUK_WGGA_counter(transcripts, genome):
    ACUK_count = []
    WGGA_count = []
    gene_names = []
    gene_lengths = []
    for tran in transcripts:
        # Collect the sequence for each gene
        seq = tran.get_sequence(genome)
        
        # define the coding region for each gene.
        cds_start = tran.cds_start
        cds_end = tran.cds_end

        # Change the sequence to only use the coding region
        cds_seq = seq[cds_start:cds_end]
        
        # Count the number of times the complementary sequence shows up in the DNA strand. 
        ACUK_count.append(cds_seq.count("TGAC") + cds_seq.count("TGAA"))
        WGGA_count.append(cds_seq.count("TCCT") + cds_seq.count("ACCT"))
        
        # Get the gene names and lengths for each sequence.
        gene_names.append(tran.attr["gene_name"])
        gene_lengths.append(len(cds_seq))
    
    # Save everything as a dataframe. 
    full_count = pd.DataFrame(list(zip(gene_names, gene_lengths, ACUK_count, WGGA_count)))
        
    return full_count

def find_transcript(gene,transcripts):
    '''
    A function that takes the name of a gene as input and finds 
    the corresponding transcript from a transcript list. 
    
    returns both the transcript in question and the vector of counts for that transcript.
    
    This function is still a work in progress as for now it simply gives the last 
    transcript in the list that matches the gene ID. 
    '''
    
    my_trans_list = []
    
    for i in transcripts:
            if i.attr['gene_name'] == gene:
                my_trans_list.append(i)
                index = transcripts.index(i)
                
    return my_trans_list, index



