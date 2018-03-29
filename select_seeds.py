import pickle
import os.path
import sys
import pandas as pd
#import matplotlib.pyplot as plot
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import ClustalOmegaCommandline

def get_acc_num(str_in):
    acc_num = str_in.split('|')[1]
    return acc_num

def get_sp_name(str_in):
    split_sp_name = str(str_in).split(' ')
    sp_name = split_sp_name[0] + ' ' + split_sp_name[1]
    return sp_name

def get_shortname(str_in):
    split = str(str_in).split(' ')
    short = split[0][0] + '_' + split[1]
    return short

def select_seqs(wkbk_out):

    with open(wkbk_out, 'ab') as outfile:
        writer = pd.ExcelWriter(outfile)
    
        single_seq_srnas = []
        
        # Load the filtered BLAST results
        with open('nr_BLAST_results.p', 'rb') as infile:
            df = pickle.load(infile)

        # Load the sRNA sequences
        with open('srnas.p', 'rb') as infile:
            srna_dict = pickle.load(infile)

        # If an sRNA length dictionary exists, load it. Or make it
        if os.path.isfile('srna_len.p'):
            with open('srna_len.p', 'rb') as infile:
                srna_len = pickle.load(infile)
        else:    
            # Make an srna length dictionary
            srna_len = {name: len(sequence) for name, sequence in srna_dict.items()}
            
        # Drop all rows with a pident < 65. Resulting table will be the basis for seed selection
        df_ss = df.drop(df[df.pident < 65].index)
        print(df_ss.shape)
        
        # Group results by sRNA name 
        groupby_sRNA = df_ss.groupby(df_ss['qseqid'])
        all_seed_seqs = pd.DataFrame()
        
        # Iterate over hits one sRNA at a time
        for name, data in groupby_sRNA:
            print(name)
            
            # Drop all rows with coverage < 95
            data_good_cov = data.drop(data[data.perc_coverage < 0.95].index)
            
            # Plot a histogram of pident
            #hist = data_good_cov['pident'].plot(kind='hist')

            count, division = np.histogram(data_good_cov['pident'])
            
            # Return the numpy "count" array as a list
            bin_bound_lst = division.tolist()
            
            # Convert that list to tuples representing hist bin boundaries
            bin_bounds = list(zip(bin_bound_lst, bin_bound_lst[1:]))
            seed_seqs = pd.DataFrame()
            
            for lower, upper in bin_bounds:
            
                # Filter the results with good coverage 
                bb_filter = (data_good_cov['pident'] > lower) & (data_good_cov['pident'] < upper)
                good_cov_filtered = data_good_cov[bb_filter]
                
                # If there are samples that fall within a given histogram bin, randomly select one of them
                if not good_cov_filtered.empty:
                    selected_hit = good_cov_filtered.sample(n=1)
                    seed_seqs = seed_seqs.append(selected_hit)
            
            seq_records = []
            prealign_fn = name + '_seqs_for_seed_unaligned.fa'
            
            for i, series in seed_seqs.iterrows():
                accession = series['sseqid']
                s = series['sstart']
                e = series['send']
                description = name + ' homolog; originates from BLAST hit ' + str(i)
                seq = series['sseq']
                new_rec = SeqRecord(Seq(seq), id=accession+'/'+s+'-'+e, description=description)
                seq_records.append(new_rec)
            all_seed_seqs = all_seed_seqs.append(seed_seqs)
        
        # Write the select sequences to file.
            aligned_outfn = name + '_seqs_clusaligned.fa'
            with open(prealign_fn, 'w') as prealign_out:
                SeqIO.write(seq_records, prealign_out, 'fasta')
            # Align the sequences with Clustal and write to file
            if len(seq_records) > 1:
                clustal_align = ClustalOmegaCommandline(infile=prealign_fn, outfile=aligned_outfn, verbose=True)
                clustal_align()
            else: 
                single_seq_srnas.append(name)
        all_seed_seqs.to_excel(writer, sheet_name='seed_sequences')
        writer.save()
    if len(single_seq_srnas) > 0:
        print('These srnas have only one sequence in their alignment file:\n %s' % single_seq_srnas)
    with open('singleton_srnas.p', 'wb') as outfile:
        pickle.dump(single_seq_srnas, outfile)

if __name__ == '__main__':
    if len(sys.argv) == 3:
         select_seqs(sys.argv[1])
    else:
         print("Usage: select_seeds.py excel_workbook.xlsx")
         sys.exit(0)