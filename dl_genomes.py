import sys
import os
import pickle
import ftplib
import os.path
import pandas as pd
from Bio import SeqIO
from Bio import Entrez

def get_sp_name(str_in):
    split_sp_name = str(str_in).split(' ')
    sp_name = split_sp_name[0] + '_' + split_sp_name[1]
    return sp_name

# Get the accession.version of the GenBank assembly that is paired to the given RefSeq assembly
# This looks for all associated genome assemblies, INCLUDING PLASMIDS
def get_accession(assembly_id):
    refseq_id = Entrez.read(Entrez.esearch(db='nucleotide', term=assembly_id))['IdList']
    if len(refseq_id) == 0:
        accession = None
    else:
        seq_record = Entrez.efetch(db='nucleotide', id=refseq_id, retmode='xml')
        results = Entrez.read(seq_record)
        accession = [result['GBSeq_accession-version'] for result in results]
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
            sp_shortname = species.split(' ')[0][0] + '_' + species.split(' ')[1]
            sum_fn = sp_shortname + '_asummary.txt'
            
            # For each species, download the corresponding assembly summary from NCBI
            if not os.path.isfile(sum_fn):
                with open(sum_fn, 'wb') as sumfile:
                    ftp = ftplib.FTP(host='ftp.ncbi.nih.gov', user='anonymous', passwd=email)
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
                        for accesions in accession_dict['%s' % sp_fullname]:
                            for acc in accesions:
                                with Entrez.efetch(db='nucleotide', id=acc, rettype='fasta', retmode='text') as handle:
                                    seq_record = SeqIO.read(handle, 'fasta')
                                    org = seq_record.description.split(' ')[1] + ' ' + seq_record.description.split(' ')[2]
                                    if org == species:
                                        SeqIO.write(seq_record, rel_out, 'fasta')
                                    else:
                                        unused_assemblies.append(seq_record.id)
                else:
                    with open('%s_genomes.fa' % sp_fullname, 'r') as rel_in:
                        # Make a dictionary of the accession numbers present in the genome FASTA
                        records = SeqIO.parse(rel_in, 'fasta')
                        accs = list(r.id for r in records)
                        for accesions in accession_dict['%s' % sp_fullname]:
                            for acc in accesions:
                                if acc not in accs:
                                    unused_assemblies.append(acc)
        
        print(ga_df)
        s = ga_df.apply(lambda x: pd.Series(x['paired_GB_accession']), axis=1).stack().reset_index(level=1, drop=True)
        s.name = 'accession'
        incl_df = ga_df.drop('paired_GB_accession', axis=1).join(s)
        print(incl_df)
        # Mark those assemblies with accessions that were not used because fasta name didn't match organism name
        incl_df.loc[incl_df.accession.isin(unused_assemblies), 'notes'] = 'dropped due to name mismatch'

        # Save a record of which accessions belong to which assemblies 
        identifier_dict = pd.Series(incl_df.gbrs_paired_asm.values, index=incl_df.accession).to_dict()
        print(identifier_dict)
        
        # Write the full set of assembly records to file.
        ga_df.to_excel(writer, sheet_name='selected_assemblies')
        incl_df.to_excel(writer, sheet_name='all_accessions_with_notes')
        writer.save()
    
    # Store assembly dictionary for later use
    with open('asmbly_dict.p', 'wb') as dictout:
        pickle.dump(identifier_dict, dictout)

dl_seqs('hdutcher@pdx.edu', 'test_species.txt', 'test.xlsx')

if __name__ == '__main__':
    if len(sys.argv) == 4:
         dl_seqs(sys.argv[1], sys.argv[2], sys.argv[3])
    else:
         print("Usage: dl_genomes.py email@domain.com species_of_interest.txt excel_workbook_out.xlsx")
         sys.exit(0)
