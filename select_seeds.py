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

def select_seqs(wkbk_out, srnas, blast, hub_accession):
    # Load the genome count dictionary
    single_seq_srnas = []
    gen_ct_dict = pickle.load(open('sp_gen_ct.p', 'rb'))
    # Store the srna sequences of interest for later use
    with open(srnas, 'r') as srnas_in:
        srna_dict = SeqIO.to_dict(SeqIO.parse(srnas_in, 'fasta'))
    with open(wkbk_out, 'ab') as outfile:
        writer = pd.ExcelWriter(outfile)
        # Input should be filtered BLAST tabular results (top hit per query-subject pair; one HSP only)
        with open(blast, 'r') as infile:
            col_names = ['qseqid', 'sseqid', 'stitle', 'pident', 'qcovs', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'ssend', 'qframe', 'sframe', 'frames', 'evalue', 'bitscore', 'qseq', 'sseq']
            df = pd.read_csv(infile, sep='\t', header=None, names=col_names)
        # Reformat the sseqid column to show accession numbers and the stitle column to show just species names
        df['sseqid'] = df['sseqid'].apply(get_acc_num)
        df['stitle'] = df['stitle'].apply(get_sp_name)
        # Drop all rows with an evalue > 0.00001
        df = df.drop(df[df.evalue > 0.00001].index)
# PRESENCE/ABSENCE MATRIX    
###################################################################################################################################################################
        # Drop all rows with qcovs < 60. Resulting table will be the basis for presence/absence matrix
        df_pam = df.drop(df[df.qcovs < 60].index)
        # Create an empty dataframe to hold the presence/absence matrix
        pa_matrix = pd.DataFrame()
        # Group results by sRNA
        pam_sRNA_grp = df_pam.groupby(df_pam['qseqid'])
        pa_dict = {}
        for name, data in pam_sRNA_grp:
            stitle_val_cts = data['stitle'].value_counts()
            print(stitle_val_cts)
            for sp_name, num_genomes in gen_ct_dict.items():		
                if sp_name in stitle_val_cts.index:
                    num_results = stitle_val_cts.get_value(sp_name)
                    if num_results >= 0.75 * num_genomes:
                        pa_dict[sp_name] = True
                    else: pa_dict[sp_name] = False
                else: 
                    pa_dict[sp_name] = False
                    print('%s not found in BLAST results' % sp_name)
            # Add that dictionary to the empty data frame
            pa_matrix[name] = pd.Series(pa_dict)
        # Convert the Trues/Falses to ones and zeros and write this presence/absence matrix to file
        pa_matrix = pa_matrix.astype(int)
        pa_matrix.to_excel(writer, sheet_name='pa_matrix_for_GLOOME')
        # Prepare to convert pa matrix to fasta
        pa_concat = pa_matrix.apply(lambda row: ''.join(map(str, row)), axis=1)
        # Set of records for fasta
        with open('msa_for_GLOOME.fa', 'w') as msa_out:
            msa_records = []
            for i, val in pa_concat.iteritems():
                seq = Seq(str(val))
                record = SeqRecord(seq, id=i)
                msa_records.append(record)
            SeqIO.write(msa_records, msa_out, 'fasta')
    # SEED SEQUENCES
    ###################################################################################################################################################################        
        # Drop all rows with a pident < 65. Resulting table will be the basis for seed selection
        df_ss = df.drop(df[df.pident < 65].index)
        print(df_ss.shape)
        # Group results by sRNA name 
        groupby_sRNA = df_ss.groupby(df_ss['qseqid'])
        all_seed_seqs = pd.DataFrame()
        # Iterate over hits one sRNA at a time
        for name, data in groupby_sRNA:
            print(name)
            # Calculate the full length of the query 
            qlen = len(srna_dict['%s' % name].seq)
            # Calculate coverage for each hit
            data['coverage'] = data['sseq'].str.len()/qlen
            # Drop all rows with coverage < 95
            data_good_cov = data.drop(data[data.coverage < 0.95].index)
            # Plot a histogram of pident
            #hist = data_good_cov['pident'].plot(kind='hist')
            count, division = np.histogram(data_good_cov['pident'])
            # Return the numpy "count" array as a list
            bin_bound_lst = division.tolist()
            # Convert that list to tuples representing hist bin boundaries
            bin_bounds = list(zip(bin_bound_lst, bin_bound_lst[1:]))
            seed_seqs = pd.DataFrame()
            # Add the hub sRNA sequence to the seed sequence list
            hub_hit = data_good_cov.loc[data_good_cov['sseqid'] == hub_accession]
            seed_seqs = seed_seqs.append(hub_hit)
            # Before selecting hits based on bin boundaries, remove the hub hits since they'll already be included
            data_good_cov_minus_hub = data_good_cov.drop(data_good_cov[data_good_cov.sseqid == hub_accession].index) 
            for lower, upper in bin_bounds:
                # Filter the results with good coverage 
                bb_filter = (data_good_cov_minus_hub['pident'] > lower) & (data_good_cov_minus_hub['pident'] < upper)
                good_cov_filtered = data_good_cov_minus_hub[bb_filter]
                # If there are samples that fall within a given histogram bin, randomly select one of them
                if not good_cov_filtered.empty:
                    selected_hit = good_cov_filtered.sample(n=1)
                    seed_seqs = seed_seqs.append(selected_hit)
            seq_records = []
            prealign_fn = name + '_seqs_for_seed.fa'
            for i, series in seed_seqs.iterrows():
                accession = series['sseqid']
                description = name + ' homolog; originates from BLAST hit ' + str(i)
                seq = series['sseq']
                new_rec = SeqRecord(Seq(seq), id=accession, description=description)
                seq_records.append(new_rec)
            all_seed_seqs = all_seed_seqs.append(seed_seqs)
        # Write the select sequences to file.
            aligned_outfn = name + '_seqs_clusaligned.fa'
            with open(prealign_fn, 'w') as prealign_out:
                SeqIO.write(seq_records, prealign_out, 'fasta')
            # Align the sequences with Clustal and write to file
            #if len(seq_records) > 1:
                #clustal_align = ClustalOmegaCommandline(infile=prealign_fn, outfile=aligned_outfn, verbose=True)
                #clustal_align()
            #else: 
                #single_seq_srnas.append(name)
        all_seed_seqs.to_excel(writer, sheet_name='seed_sequences')
        writer.save()
        if len(single_seq_srnas) > 0:
            print('These srnas had no homologs outside of the hub:\n %s' % single_seq_srnas)
        with open('singleton_srnas.p', 'wb') as outfile:
            pickle.dump(single_seq_srnas, outfile)

select_seqs('pl_results.xlsx', 'ecoli_srnas_for_plpart01.fa', 'ecoli_sRNAS_v_rel_DB_fil.txt', 'NC_000913.3')

if __name__ == '__main__':
    if len(sys.argv) == 5:
         select_seqs(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
    else:
         print("Usage: select_seeds.py excel_workbook.xlsx srna_seqs.fa filtered_blast_results.txt hub_accession")
         sys.exit(0)