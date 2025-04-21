# -*- coding: utf-8 -*-
import os, time
import threading
from pathlib import Path 
from deepcircblastp import cblast, MKEmbedding, DLPredict, merge
import argparse


#%% 主程序
def main(file_query, db_fasta, db_path, folder_results):
    # -------------> Const <-------------
    matrix_file = 'deepcircblastp/matrix/BLOSUM62.txt'
    file_circblastp_csv = Path(folder_results) / 'result_circblastp.csv'
    folder_alignments = Path(folder_results) / 'alignments'
    folder_embeddings = Path(folder_results) / 'embeddings'
    file_cppredict_csv = Path(folder_results) / 'result_DLpredict.csv'
    file_output = Path(folder_results) / 'result_deepcircblastp.csv'
    
    Path(folder_results).mkdir(parents=True, exist_ok=True)
    folder_alignments.mkdir(parents=True, exist_ok=True)
    folder_embeddings.mkdir(parents=True, exist_ok=True)
    seqs_path = Path(db_path) / 'tmp'
    seqs_path.mkdir(parents=True, exist_ok=True)
    filter_fasta = Path(db_fasta).with_stem('filtered')

    # -------------> split db fasta <-------------
    fp = os.popen(
        f"seqkit grep -s -r -v -p '[^ACDEFGHIKLMNPQRSTVWY]' {db_fasta} > {str(filter_fasta)};"
        f"seqkit split -j 10 -i -O {seqs_path} {str(filter_fasta)} "
    )
    fp.read()

    # -------------> run circblast <-------------
    now = time.strftime('%m-%d %H:%M', time.localtime())
    print(f'[{now}] Make circblast db...')
    make_circ_DB = cblast.MakeCircDB(seqs_path, matrix_file)
    make_circ_DB.db_build()
    now = time.strftime('%m-%d %H:%M', time.localtime())
    print(f'[{now}] Run circblast...')
    circblast = cblast.CircleBLASTp(file_query, seqs_path, matrix_file)
    circblast.run(folder_alignments, file_circle_csv)

    # -------------> make embeddings <-------------
    now = time.strftime('%m-%d %H:%M', time.localtime())
    print(f'[{now}] Make embeddings...')
    MKEmbedding.run(file_circblastp_csv, folder_alignments, folder_embeddings)

    # -------------> CPDetector predict <-------------
    now = time.strftime('%m-%d %H:%M', time.localtime())
    print(f'[{now}] CPDetector predicting...')
    DLPredict.run(
        folder_embeddings, 
        'deepcircblastp/models/model_deepcircblastp/model_best.pth',
        file_cppredict_csv,
    )

    # -------------> output <-------------
    now = time.strftime('%m-%d %H:%M', time.localtime())
    print(f'[{now}] DeepCircBLASTp outputting...')
    merge.output(seqs_path, file_circblastp_csv, file_cppredict_csv, file_output)

    # -------------> done <-------------
    now = time.strftime('%m-%d %H:%M', time.localtime())
    print(f'[{now}] all done.')




if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    
    parser.add_argument('file_query', type=str)
    parser.add_argument('database_fasta', type=str)
    parser.add_argument('database_path', type=str)

    parser.add_argument(
        '--folder_results', 
        type = str, 
        default = 'results'
    )


    args = parser.parse_args()
    
    main(
        args.file_query,
        args.database_fasta, 
        args.database_path, 
        args.folder_results,
    )


# python run.py query/5zoa.fasta database/test.fasta database 