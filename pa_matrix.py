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

def overlap(start1, end1, start2, end2):
    """
    Does the range (start1, end1) overlap with (start2, end2)?
    """
    return end1 >= start2 and end2 >= start1

def compare_rows(group):
    print(group)
    winners = []
    skip = []
    if len(group) == 1:
        return group
    for i in group.index:
        if i in skip:
            continue
        for j in group.index:
            last = j == group.index[-1]
            istart = group.loc[i, 'start']
            iend = group.loc[i, 'end']
            jstart = group.loc[j, 'start']
            jend = group.loc[j, 'end']
            if overlap(istart, iend, jstart, jend):
                winner = group.loc[[i, j], 'evalue'].idxmin()
                if winner == j:
                    winners.append(winner)
                    skip.append(i)
                    break
            if last:
                winners.append(i)
    return group.loc[winners].drop_duplicates()

def select_seqs(wkbk_out, srnas, blast):

    # Load the assembly dictionary
    asmbly_dict = pickle.load(open('asmbly_dict.p', 'rb'))
    
    # Store the srna sequences of interest for later use
    with open(srnas, 'r') as srnas_in:
        srna_dict = SeqIO.to_dict(SeqIO.parse(srnas_in, 'fasta'))
    with open(wkbk_out, 'ab') as outfile:
        writer = pd.ExcelWriter(outfile)
        
        # Input should be filtered BLAST tabular results (top hit per query-subject pair; one HSP only)
        with open(blast, 'r') as infile:
            col_names = ['qseqid', 'sseqid', 'stitle', 'pident', 'qcovs', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'ssend', 'qframe', 'sframe', 'frames', 'evalue', 'bitscore', 'qseq', 'sseq']
            df = pd.read_csv(infile, sep='\t', header=None, names=col_names)
        
        # Reformat the sseqid column to show accession numbers and add a column with just the species name
        df['sseqid'] = df['sseqid'].apply(get_acc_num)
        df['sshorttitle'] = df['stitle'].apply(get_sp_name)
        
        # Add columns to reflect maximum and minimum as start/end
        df['start'] = df[['sstart','ssend']].min(axis=1)
        df['end'] = df[['sstart','ssend']].max(axis=1)

        # Create an empty dataframe to hold the presence/absence matrix
        pa_matrix = pd.DataFrame()
       
        # Sort the results by sseqid and evalue
        df.sort_values(['sseqid', 'evalue'], ascending=[True, True], inplace=True)
        
        # Identify any overlaps in BLAST results and remove them from the dataframe
        df_nr = df.groupby(['sseqid', 'sframe']).apply(compare_rows)
        df_nr.reset_index

        # Find the difference between the original results and the nr results, to be saved and examined later. Save BLAST results, too.
        overlaps = df[~df.index.isin(df_nr.index)]
        df_nr.to_excel(writer, sheet_name='BLAST_top_hit_only')
        overlaps.to_excel(writer, sheet_name='removed_due_to_overlap')

        # Group results by sRNA    
        srna_grp = df_nr.groupby(df_nr['qseqid'])

        for name, data in srna_grp:
            print(name)
            pa_dict = {}
            present = data['sseqid'].tolist()
            for key, value in asmbly_dict.items():
                if key in present:
                    pa_dict[value] = True
                else:
                    if value in pa_dict:
                        continue
                    else:
                        pa_dict[value] = False
            # Add that dictionary to the empty data frame
            pa_matrix[name] = pd.Series(pa_dict)
            print(pa_matrix)
        
        # Convert Nan values to False
        print(pa_matrix)

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
                #sp_short = get_shortname(i)
                record = SeqRecord(seq, id=i)
                msa_records.append(record)
            SeqIO.write(msa_records, msa_out, 'fasta')
        writer.save()

select_seqs('test_seeds.xls', 'srnas.fasta', 'example_blast_results.txt')