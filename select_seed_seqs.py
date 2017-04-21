import csv, sys
import pandas as pd
import matplotlib.pyplot as plot
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Align.Applications import ClustalOmegaCommandline
import ftplib
from Bio import Entrez
import os.path

hub_species = 'Streptococcus_pneumoniae'
hub_accession = 'NC_003028.3'
Entrez.email = 'hdutcher@pdx.edu'

def get_acc_num(str_in):
    acc_num = str_in.split('|')[1]
    return acc_num

def get_sp_name(str_in):
    split_sp_name = str(str_in).split(' ')
    sp_name = split_sp_name[0] + '_' + split_sp_name[1]
    return sp_name

with open('species_of_interest.txt') as sp_in:
    sp_list = list(line.rsplit('\n')[0] for line in sp_in)

# Create an excel workbook to write results to:
with open('seed_msa_seqs_and_pa_matrix.xlsx', 'wb') as outfile:
    writer = pd.ExcelWriter(outfile)
    ga_df = pd.DataFrame()
    accession_list = []
    hub_sp_accs = []
    sp_gen_ct_dict = {}
    for species in sp_list:
        print(species)
        sp_gen_ct = 0
        sp_fullname = species.split(' ')[0] + '_' + species.split(' ')[1]
        sp_shortname = species.split(' ')[0][0] + '_' + species.split(' ')[1]
        sum_fn = sp_shortname + '_asummary.txt'
        if not os.path.isfile(sum_fn):
            with open(sum_fn, 'wb') as sumfile:
                ftp = ftplib.FTP(host='ftp.ncbi.nih.gov', user='anonymous', passwd='hdutcher@pdx.edu')
                sum_dir = '/genomes/genbank/bacteria/' + sp_fullname + '/'
                ftp.cwd(sum_dir)
                ftp.retrbinary('RETR assembly_summary.txt', sumfile.write)
                ftp.quit
        with open(sum_fn, 'r') as sumfile:
            assembly_df = pd.read_csv(sumfile, sep='\t', skiprows=1, header=0)
            print(assembly_df.shape)
            # Drop all of the rows for which assembly level != Complete Genome
            assembly_df = assembly_df.drop(assembly_df[assembly_df.assembly_level != 'Complete Genome'].index)
            # Drop all of the rows for which Genbank and RefSeq assemblies are not identical
            assembly_df = assembly_df.drop(assembly_df[assembly_df.paired_asm_comp != 'identical'].index)
            print(assembly_df.shape)
            # Append retained records to assembly master dataframe
            ga_df = ga_df.append(assembly_df)
            # Get the accession.version of the GenBank assembly that is paired to the given RefSeq assembly, or vice-versa
            assembly_ids = (list(assembly_df['gbrs_paired_asm']))
            # Query NCBI databases to find accession number associated with selected chromosome assembly
            for assembly_id in assembly_ids:
                term = assembly_id + ' AND complete genome'
                refseq_id = Entrez.read(Entrez.esearch(db='nucleotide', term=term))['IdList']
                if len(refseq_id) == 0:
                    print('no Ids found for assembly: %s. Species: %s' % (assembly_id, species))
                    break
                seq_record = Entrez.efetch(db='nucleotide', id=refseq_id, retmode='xml')
                results = Entrez.read(seq_record)
                accession = results[0]['GBSeq_accession-version']
                accession_list.append(accession)
                sp_gen_ct += 1
                # if this is the hub species, save a separate list of accession numbers for the hub species only
                if sp_fullname == hub_species:
                    hub_sp_accs.append(accession)
        sp_gen_ct_dict[get_sp_name(species)] = sp_gen_ct
    # Write the full set of assembly records to file.
    ga_df = ga_df.to_excel(writer, sheet_name='selected_assemblies')
    # Download the fasta based on the genbank accession number and write to file.
    if not os.path.isfile('hub_rel_genomes.fa'):
        with open('hub_rel_genomes.fa', 'a') as rel_out:
            for acc in accession_list:
                net_handle = Entrez.efetch(db='nucleotide', id=acc, rettype='fasta', retmode='text')
                rel_out.write(net_handle.read())
    # Write hub species strains only to a separate file
    if not os.path.isfile('hub_strains.fa'):
        with open('hub_strains.fa', 'a') as strain_out:
            for acc in hub_sp_accs:
                net_handle = Entrez.efetch(db='nucleotide', id=acc, rettype='fasta', retmode='text')
                strain_out.write(net_handle.read())
    print('Coolio. Now all the genomes you need are downloaded. Make a blast database out of them and blast the sRNA sequences again them. Then input a filtered file of your blast results to continue to the next step.')

#fil_blast_name = input('Name of filtered, tabular BLAST results from sRNAs vs. all hub relatives: ')
###################################################################################################################################################################

# Input should be filtered BLAST tabular results (top hit per query-subject pair; one HSP only)
    with open('example_blast_results.csv', 'r') as infile:
        col_names = ['qseqid', 'sseqid', 'stitle', 'pident', 'qcovs', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'ssend', 'evalue', 'bitscore', 'qseq', 'sseq']
        
        # CHANGE SEP TO TAB WHEN RUNNING FOR REAL
        df = pd.read_csv(infile, sep=',', header=None, names=col_names)
        print(df.shape)
        print(df.head)
        # Reformat the sseqid column to show accession numbers and the stitle column to show just species names
        df['sseqid'] = df['sseqid'].apply(get_acc_num)
        df['stitle'] = df['stitle'].apply(get_sp_name)
        # Drop all rows with an evalue > 0.00001
        df = df.drop(df[df.evalue > 0.00001].index)
        print(df.shape)
    ###################################################################################################################################################################
        # Drop all rows with qcovs < 60. Resulting table will be the basis for presence/absence matrix
        df_pam = df.drop(df[df.qcovs < 60].index)
        print(df_pam.shape)
        # Create an empty dataframe to hold the presence/absence matrix
        pa_matrix = pd.DataFrame()
        # Group results by sRNA
        pam_sRNA_grp = df_pam.groupby(df_pam['qseqid'])
        print(sp_gen_ct_dict)
        pa_dict = {}
        for name, data in pam_sRNA_grp:
            stitle_val_cts = data['stitle'].value_counts()
            for sp_name, num_genomes in sp_gen_ct_dict.items():
                if sp_name in stitle_val_cts.index:
                    num_results = stitle_val_cts.get_value(sp_name)
                    if num_results >= 0.75 * num_genomes:
                        pa_dict[sp_name] = True
                    else: pa_dict[sp_name] = False
                else: pa_dict[sp_name] = False
            print(pa_dict)
            # Add that dictionary to the empty data frame
            pa_matrix[name] = pd.Series(pa_dict)
            print(pa_matrix.head(n=35))
        # Convert the Trues/Falses to ones and zeros
        pa_matrix = pa_matrix.astype(int)
        pa_concat = pa_matrix.apply(lambda row: ''.join(map(str, row)), axis=1)
        print(pa_concat.head)

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
            # Calculate the full length of the query (will be equivalent to qseq length for the top result)
            qlen = len(data.iloc[0]['qseq'])
            data['orig_qlen'] = qlen
            # Calculate coverage for each hit
            data['coverage'] = data['sseq'].str.len()/qlen
            print(data.shape)
            # Drop all rows with coverage < 95
            data_good_cov = data.drop(data[data.coverage < 0.95].index)
            print(data_good_cov.shape)
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
            for lower, upper in bin_bounds:
                bb_filter = (data_good_cov['pident'] > lower) & (data_good_cov['pident'] < upper)
                good_cov_filtered = data_good_cov[bb_filter]
                # If there are samples the fall within a given histogram bin, randomly select one of them
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
        all_seed_seqs.to_excel(writer, sheet_name='seed_sequences')
        writer.save()
            # Write the select sequences to file.
            #aligned_outfn = name + '_seqs_clusaligned.fa'
            #with open(prealign_fn, 'w') as prealign_out:
                #SeqIO.write(seq_records, prealign_out, 'fasta')
            # Align the sequences with Clustal and write to file
            #clustal_align = ClustalOmegaCommandline(infile=prealign_fn, outfile=aligned_outfn, verbose=True)
            #clustal_align()
            #plot.show()
# Make a pident hist of all sRNAs at once
#df_ss['pident'].hist(by=df_ss['qseqid'])
