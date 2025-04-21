import re, torch 
from transformers import T5Tokenizer, T5EncoderModel
import json
import argparse
import numpy as np 
import pandas as pd 
from pathlib import Path 


def make_tensor(seqs_align: list, tokenizer):
    sequence = [
        " ".join(list(re.sub(r"[-UZOB]", "X", seq)))
        for seq in seqs_align
    ]
    encoded_input = tokenizer(
        sequence, 
        add_special_tokens = True, 
        max_length = 600, 
        padding = "max_length", 
        truncation = True,
    )
    input_ids = torch.tensor(encoded_input['input_ids']).cuda()
    attention_mask = torch.tensor(encoded_input['attention_mask']).cuda()

    return input_ids, attention_mask


def run(file_circblastp_csv, folder_alignments, folder_embeddings):
    df_circle = pd.read_csv(file_circblastp_csv, header=0)
    df_circle.columns = ['stem', 'identity', 'evalue'] 
    df_tmp = df_circle.sort_values('evalue')
    df_tmp.reset_index(inplace=True, drop=True)
    df_embed = df_tmp.iloc[0:1000, :]
    
    path_model = "deepcircblastp/models/Rostlab--prot_t5_xl_half_uniref50-enc"
    tokenizer = T5Tokenizer.from_pretrained(path_model, do_lower_case=False,)
    model = T5EncoderModel.from_pretrained(path_model)
    model.half().cuda()

    # top 1k 做 embedding
    for index, row in df_embed.iterrows():
        stem = row['stem']
        file_align = folder_alignments / f'{stem}.json'
        seqs_align = json.load(file_align.open())

        input_ids, attention_mask = make_tensor(seqs_align, tokenizer)
        with torch.no_grad():
            output = model(input_ids=input_ids, attention_mask=attention_mask)
        embedding_cuda = output.last_hidden_state
        embedding = embedding_cuda.cpu().numpy()

        file_npy = Path(folder_embeddings) / f"{stem}.npy"
        np.save(file_npy, embedding)




if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--file_circblastp_csv', 
        type = str, 
        default = None,
    )
    parser.add_argument(
        '--folder_alignments', 
        type = str, 
        default = None,
    )
    parser.add_argument(
        '--folder_embeddings', 
        type = str, 
        default = None,
    )
    args = parser.parse_args()
    run(args.file_circblatp_csv, args.folder_alignments, args.folder_embeddings)

