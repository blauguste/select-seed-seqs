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
    unused_assemblies = []

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
            with open(sum_fn, 'rb') as sumfile:
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
                # Download the fasta based on the genbank accession number and write to file if the name on the fasta matches the species name.
                if not os.path.isfile('%s_genomes.fa' % sp_fullname):
                    with open('%s_genomes.fa' % sp_fullname, 'a') as rel_out:
                        for acc in accession_dict['%s' % sp_fullname]:
                            with Entrez.efetch(db='nucleotide', id=acc, rettype='fasta', retmode='text') as handle:
                                seq_record = SeqIO.read(handle, 'fasta')
                                org = seq_record.description.split(' ')[1] + ' ' + seq_record.description.split(' ')[2]
                                if org == species:
                                    SeqIO.write(seq_record, rel_out, 'fasta')
                                    sp_gen_ct += 1
                                else:
                                    unused_assemblies.append(seq_record.id)
                else:
                    with open('%s_genomes.fa' % sp_fullname, 'r') as rel_in:
                        # Make a dictionary of the accession numbers present in the genome FASTA
                        records = SeqIO.parse(rel_in, 'fasta')
                        accs = list(r.id for r in records)
                        for acc in accession_dict['%s' % sp_fullname]:
                            if acc not in accs:
                                unused_assemblies.append(acc)
                            else:
                                sp_gen_ct += 1
            sp_gen_ct_dict[species] = sp_gen_ct
        # Mark those assemblies with accessions that were not used because fasta name didn't match organism name
        ga_df.loc[ga_df.paired_GB_accession.isin(unused_assemblies), 'notes'] = 'dropped due to name mismatch'
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
