from plastid import BAMGenomeArray, Transcript
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import csv
import math
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq

def counter(transcripts, genome, sites = ["TAC", "TAT", "TGGA", "GAC"]):
    total_cds = []
    total_utr5 = []
    total_utr3 = []
    total_all = []
    lengths_cds = []
    lengths_utr5 = []
    lengths_utr3 = []
    lengths_all = []
    gene_name = []
    gene_id = []
    for tran in transcripts:
        # Collect the sequence for each gene and transcribe to RNA sequence.
        seq = Seq(tran.get_sequence(genome), generic_dna)
        seq = seq.transcribe()
        
        # Check to ensure that it is a protein coding gene
        try:
            tran.attr["gene_name"]
        except:
            continue
        
        # define the coding region for each gene.
        cds_start = tran.cds_start
        cds_end = tran.cds_end
        cds_seq = seq[cds_start:cds_end]
        
        # define the 3' and 5' utr regions for each gene.
        utr5 = seq[0:cds_start]
        utr3 = seq[cds_end:-1]
        
        # Count the number of times the complementary sequence shows up in the DNA strand. 
        c,u5,u3,a = 0,0,0,0
        for site in sites:
            c = c + cds_seq.upper().count(site)
            u5 = u5 + utr5.upper().count(site)
            u3 = u3 + utr3.upper().count(site)
            a = a + seq.upper().count(site)
        total_cds.append(c)
        total_utr5.append(u5)
        total_utr3.append(u3)
        total_all.append(a)
        
        # Get the lengths for each gene region
        lengths_cds.append(len(cds_seq))
        lengths_utr5.append(len(utr5))
        lengths_utr3.append(len(utr3))
        lengths_all.append(len(seq))
    
        # Get the gene names for each sequence.
        gene_name.append(tran.attr["gene_name"])
        gene_id.append(tran.attr["transcript_id"])
    
    # get the count density for each region.
    density_cds = np.array(total_cds)/np.array(lengths_cds)
    density_utr5 = np.array(total_utr5)/np.array(lengths_utr5)
    density_utr3 = np.array(total_utr3)/np.array(lengths_utr3)
    density_all = np.array(total_all)/np.array(lengths_all)
    
    # Save everything as a dataframe. 
    full_count = pd.DataFrame(list(zip(gene_name, gene_id, lengths_cds, lengths_utr5, lengths_utr3, lengths_all,
                                       total_cds, total_utr5, total_utr3, total_all,
                                       density_cds, density_utr5, density_utr3, density_all)))
    full_count.columns = ["gene_name", "gene_id", "lengths_cds", "lengths_utr5", "lengths_utr3u", "lengths_all",
                                       "total_cds", "total_utr5", "total_utr3", "total_all",
                                       "density_cds", "density_utr5", "density_utr3", "density_all"]
        
    return full_count


def ACUK_WGGA_counter(transcripts, genome):
    ACUK_cds = []
    WGGA_cds = []
    ACUK_utr5 = []
    WGGA_utr5 = []
    ACUK_utr3 = []
    WGGA_utr3 = []
    ACUK_all = []
    WGGA_all = []
    lengths_cds = []
    lengths_utr5 = []
    lengths_utr3 = []
    lengths_all = []
    gene_name = []
    gene_ID = []
    for tran in transcripts:
        # Collect the sequence for each gene and transcribe to RNA sequence.
        seq = Seq(tran.get_sequence(genome), generic_dna)
        seq = seq.transcribe()
        
        # Check to ensure that it is a protein coding gene
        try:
            tran.attr["gene_name"]
        except:
            continue
        
        # define the coding region for each gene.
        cds_start = tran.cds_start
        cds_end = tran.cds_end
        cds_seq = seq[cds_start:cds_end]
        
        # define the 3' and 5' utr regions for each gene.
        utr5 = seq[0:cds_start]
        utr3 = seq[cds_end:-1]
        
        # Count the number of times the complementary sequence shows up in the DNA strand. 
        ACUK_cds.append(cds_seq.upper().count("ACUU") + cds_seq.upper().count("ACUG"))
        WGGA_cds.append(cds_seq.upper().count("AGGA") + cds_seq.upper().count("UGGA"))
        
        ACUK_utr5.append(utr5.upper().count("ACUU") + utr5.upper().count("ACUG"))
        WGGA_utr5.append(utr5.upper().count("AGGA") + utr5.upper().count("UGGA"))

        ACUK_utr3.append(utr3.upper().count("ACUU") + utr3.upper().count("ACUG"))
        WGGA_utr3.append(utr3.upper().count("AGGA") + utr3.upper().count("UGGA"))

        ACUK_all.append(seq.upper().count("ACUU") + seq.upper().count("ACUG"))
        WGGA_all.append(seq.upper().count("AGGA") + seq.upper().count("UGGA"))
    
        # Get the lengths for each gene region
        lengths_cds.append(len(cds_seq))
        
        lengths_utr5.append(len(utr5))
        
        lengths_utr3.append(len(utr3))
        
        lengths_all.append(len(seq))
    
        # Get the gene names for each sequence.
        gene_name.append(tran.attr["gene_name"])
        gene_ID.append(tran.attr["transcript_id"])
    
    # get the total counts for each region.
    total_cds = np.array(ACUK_cds) + np.array(WGGA_cds)
    total_utr5 = np.array(ACUK_utr5) + np.array(WGGA_utr5)
    total_utr3 = np.array(ACUK_utr3) + np.array(WGGA_utr3)
    total_all = np.array(ACUK_all) + np.array(WGGA_all)
    
    # get the count density for each region.
    density_cds = np.array(total_cds)/np.array(lengths_cds)
    density_utr5 = np.array(total_utr5)/np.array(lengths_utr5)
    density_utr3 = np.array(total_utr3)/np.array(lengths_utr3)
    density_all = np.array(total_all)/np.array(lengths_all)
    
    # Save everything as a dataframe. 
    full_count = pd.DataFrame(list(zip(gene_name, gene_ID, lengths_cds, lengths_utr5, lengths_utr3, lengths_all,
                                       total_cds, total_utr5, total_utr3, total_all,
                                       density_cds, density_utr5, density_utr3, density_all)))
    full_count.columns = ["gene_name", "gene_id", "lengths_cds", "lengths_utr5", "lengths_utr3u", "lengths_all",
                                       "total_cds", "total_utr5", "total_utr3", "total_all",
                                       "density_cds", "density_utr5", "density_utr3", "density_all"]
        
    return full_count

def tetra_counter(transcripts, genome, tetras = ["ACUU", "ACUG", "AGGA", "UGGA"]):
    tetra1_cds = []
    tetra2_cds = []
    tetra1_utr5 = []
    tetra2_utr5 = []
    tetra1_utr3 = []
    tetra2_utr3 = []
    tetra1_all = []
    tetra2_all = []
    lengths_cds = []
    lengths_utr5 = []
    lengths_utr3 = []
    lengths_all = []
    gene_name = []
    gene_id = []
    for tran in transcripts:
        # Collect the sequence for each gene and transcribe to RNA sequence.
        seq = Seq(tran.get_sequence(genome), generic_dna)
        seq = seq.transcribe()
        
        # Check to ensure that it is a protein coding gene
        try:
            tran.attr["gene_name"]
        except:
            continue
        
        # define the coding region for each gene.
        cds_start = tran.cds_start
        cds_end = tran.cds_end
        cds_seq = seq[cds_start:cds_end]
        
        # define the 3' and 5' utr regions for each gene.
        utr5 = seq[0:cds_start]
        utr3 = seq[cds_end:-1]
        
        # Count the number of times the complementary sequence shows up in the DNA strand. 
        tetra1_cds.append(cds_seq.upper().count(tetras[0]) + cds_seq.upper().count(tetras[1]))
        tetra2_cds.append(cds_seq.upper().count(tetras[2]) + cds_seq.upper().count(tetras[3]))
        
        tetra1_utr5.append(utr5.upper().count(tetras[0]) + utr5.upper().count(tetras[1]))
        tetra2_utr5.append(utr5.upper().count(tetras[2]) + utr5.upper().count(tetras[3]))

        tetra1_utr3.append(utr3.upper().count(tetras[0]) + utr3.upper().count(tetras[1]))
        tetra2_utr3.append(utr3.upper().count(tetras[2]) + utr3.upper().count(tetras[3]))

        tetra1_all.append(seq.upper().count(tetras[0]) + seq.upper().count(tetras[1]))
        tetra2_all.append(seq.upper().count(tetras[2]) + seq.upper().count(tetras[3]))
    
        # Get the lengths for each gene region
        lengths_cds.append(len(cds_seq))
        
        lengths_utr5.append(len(utr5))
        
        lengths_utr3.append(len(utr3))
        
        lengths_all.append(len(seq))
    
        # Get the gene names for each sequence.
        gene_name.append(tran.attr["gene_name"])
        gene_id.append(tran.attr["transcript_id"])
    
    # get the total counts for each region.
    total_cds = np.array(tetra1_cds) + np.array(tetra2_cds)
    total_utr5 = np.array(tetra1_utr5) + np.array(tetra2_utr5)
    total_utr3 = np.array(tetra1_utr3) + np.array(tetra2_utr3)
    total_all = np.array(tetra1_all) + np.array(tetra2_all)
    
    # get the count density for each region.
    density_cds = np.array(total_cds)/np.array(lengths_cds)
    density_utr5 = np.array(total_utr5)/np.array(lengths_utr5)
    density_utr3 = np.array(total_utr3)/np.array(lengths_utr3)
    density_all = np.array(total_all)/np.array(lengths_all)
    
    # Save everything as a dataframe. 
    full_count = pd.DataFrame(list(zip(gene_name, gene_id, lengths_cds, lengths_utr5, lengths_utr3, lengths_all,
                                       total_cds, total_utr5, total_utr3, total_all,
                                       density_cds, density_utr5, density_utr3, density_all)))
    full_count.columns = ["gene_name", "gene_id", "lengths_cds", "lengths_utr5", "lengths_utr3u", "lengths_all",
                                       "total_cds", "total_utr5", "total_utr3", "total_all",
                                       "density_cds", "density_utr5", "density_utr3", "density_all"]
        
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



