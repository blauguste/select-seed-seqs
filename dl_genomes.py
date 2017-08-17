import sys
import os
import pickle
import ftplib
import os.path
import pandas as pd
#import matplotlib.pyplot as plot
#import numpy as np
from Bio import SeqIO
from Bio import Entrez

def get_sp_name(str_in):
    split_sp_name = str(str_in).split(' ')
    sp_name = split_sp_name[0] + '_' + split_sp_name[1]
    return sp_name

# Get the accession.version of the GenBank assembly that is paired to the given RefSeq assembly
def get_accession(assembly_id):
    term = assembly_id + ' AND complete genome NOT plasmid[Title]'
    refseq_id = Entrez.read(Entrez.esearch(db='nucleotide', term=term))['IdList']
    if len(refseq_id) == 0:
        accession = None
    else:
        seq_record = Entrez.efetch(db='nucleotide', id=refseq_id, retmode='xml')
        results = Entrez.read(seq_record)
        accession = results[0]['GBSeq_accession-version']
    return accession

def dl_seqs(email, sp_infile, wkbk_out):

    Entrez.email = email

    # Read the species names into a list
    with open(sp_infile, 'r') as sp_in:
        sp_list = list(line.rsplit('\n')[0] for line in sp_in)

    # Create an excel workbook to write results to:
    assembly_list = []
    with open(wkbk_out, 'wb') as outfile:
        writer = pd.ExcelWriter(outfile)
        ga_df = pd.DataFrame()
        accession_dict = {}
        sp_gen_ct_dict = {}
        for species in sp_list:
            print(species)
            sp_gen_ct = 0
            sp_fullname = species.split(' ')[0] + '_' + species.split(' ')[1]
            sp_shortname = species.split('_')[0][0] + '_' + species.split(' ')[1]
            sum_fn = sp_shortname + '_asummary.txt'
            # For each species, download the corresponding assembly summary from NCBI
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
                # Query NCBI databases to find accession number associated with each chromosome assembly, append to summary:
                assembly_df['paired_GB_accession'] = assembly_df['gbrs_paired_asm'].apply(get_accession)
                new_df = assembly_df[assembly_df['paired_GB_accession'].notnull()]
                accession_dict['%s' % sp_fullname] = list(new_df['paired_GB_accession'])
                assembly_list.extend(list(assembly_df['gbrs_paired_asm']))
                # Append retained records to assembly master dataframe
                ga_df = ga_df.append(assembly_df)
                # Download the fasta based on the genbank accession number and write to file.
                sp_gen_ct = 0
                with open('%s_genomes.fa' % sp_fullname, 'a') as rel_out:
                    for acc in accession_dict['%s' % sp_fullname]:
                        net_handle = Entrez.efetch(db='nucleotide', id=acc, rettype='fasta', retmode='text')
                        rel_out.write(net_handle.read())
                        sp_gen_ct += 1
            with open('%s_genomes.fa' % sp_fullname, 'r') as fastas_in:
                records = SeqIO.parse(fastas_in, 'fasta')
                sp_names = set(record.description.split(' ')[1] + ' ' + record.description.split(' ')[2] for record in records)
                print(sp_names)
                # Print a message if name in gb file doesn't match original species name
                has_alias = False
                if len(sp_names) > 1: 
                    print('This species goes by multiple names. Manually curate PA matrix.')
                if species not in sp_names and len(sp_names) >= 1:
                    print('Species alias found. Input species name: %s; Fasta species name(s): %s' % (species, sp_names))
                    alias = input('Enter species name to use as alias <Genus species>, or type <n> to proceed without alias\n')
                    if alias.lower() == 'no' or alias.lower() == 'n':
                        alias = None
                    else: has_alias = True
                else: alias = None
            sp_gen_ct_dict[species] = sp_gen_ct
            if has_alias == True:
                sp_gen_ct_dict[alias] = sp_gen_ct
        # Write the full set of assembly records to file.
        ga_df.to_excel(writer, sheet_name='selected_assemblies')
        writer.save()
    # Store genome count dictionary for later use
    with open('sp_gen_ct.p', 'wb') as sp_ct_out:
        pickle.dump(sp_gen_ct_dict, sp_ct_out)
    # Store assembly id list for later use (finding plasmids)
    with open('assembly_ids.p', 'wb') as assem_out:
        pickle.dump(assembly_list, assem_out)

if __name__ == '__main__':
    if len(sys.argv) == 4:
         dl_seqs(sys.argv[1], sys.argv[2], sys.argv[3])
    else:
         print("Usage: dl_genomes.py email@domain.com species_of_interest.txt excel_workbook_out.xlsx")
         sys.exit(0)
