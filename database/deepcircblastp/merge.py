import pandas as pd 
from Bio import SeqIO


#%% output
def output(seqs_path, file_circblastp_csv, file_cppredict_csv, file_output):
    df_circle = pd.read_csv(file_circblastp_csv, header=0)
    df_circle.columns = ['stem', 'identity', 'evalue']

    df_predic = pd.read_csv(file_cppredict_csv, header=0)

    df_merge = pd.merge(df_circle, df_predic, how='inner', on='stem')
    df_tmp = df_merge.sort_values('evalue')
    df_tmp.reset_index(inplace=True, drop=True)

    access = []

    for stem in df_tmp['stem']:
        file_fasta = seqs_path / f'{stem}.fasta'
        seq_record = SeqIO.read(file_fasta, 'fasta') 
        seq_id = seq_record.id
        access.append(seq_id)

    df_tmp['access'] = access

    df_deepcircblastp = df_tmp[
        ['access', 'identity', 'evalue', 'isCP']
    ]
    df_deepcircblastp.to_csv(file_output, index=None)

    