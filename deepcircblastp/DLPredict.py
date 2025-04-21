import argparse
import numpy as np 
import pandas as pd 
from pathlib import Path 
import warnings 
warnings.filterwarnings('ignore')

import torch
from torch import nn 
import torch.nn.functional as F
from torch.utils.data import DataLoader, Dataset


class iDataset(Dataset):
    def __init__(self, folder):
        self.files = list(Path(folder).glob('*.npy'))

    def __len__(self):
        return len(self.files)

    def __getitem__(self, index):
        file_index = self.files[index]
        data = np.load(file_index)
        data_tensor = torch.from_numpy(data).float()
        stem = Path(file_index).stem
        
        return data_tensor, stem 


class MyCNN(nn.Module):
    def __init__(self, dim_in, dim_out,):
        super().__init__()
        self.conv1 = nn.Conv2d(
            in_channels=4, 
            out_channels=256,
            kernel_size=(6,dim_in),
        )
        self.conv2 = nn.Conv2d(
            in_channels=4,
            out_channels=256, 
            kernel_size=(5, dim_in), 
        )
        self.drop = nn.Dropout(0.5)
        self.fc = nn.Linear(256*2, dim_out)

    def forward(self, x):
        x_c1 = self.conv1(x)
        x_relu = F.relu(x_c1)
        x_pool = F.max_pool2d(x_relu, x_relu.size(3))
        out_c1 = x_pool.squeeze([2,3])
        x_c2 = self.conv2(x)
        x_relu = F.relu(x_c2)
        x_pool = F.max_pool2d(x_relu, x_relu.size(3))
        out_c2 = x_pool.squeeze([2,3])
        x_conv = torch.cat(
            [out_c1, out_c2,],
            axis = 1, 
        )
        x_drop = self.drop(x_conv)
        out = self.fc(x_drop)
        return out 


def run(folder_embeddings, CPDetector_path, file_cppredict_csv):
    model = MyCNN(dim_in=600, dim_out=2)
    model.load_state_dict(
        torch.load(
            CPDetector_path,
            weights_only = True,
        )
    )
    model.cuda()
    model.eval()

    batch_size = 5
    dataset = iDataset(folder_embeddings)

    iloader = DataLoader(
        dataset = dataset, 
        batch_size = batch_size, 
        shuffle = True, 
        num_workers = 0,
    )

    dict_predict = {
        'stem': [],
        'prediction': [],
    }

    for embeddings, stems in iloader:
        embeddings = embeddings.cuda()
        prediction = model(embeddings)
        prediction = torch.max(prediction.data, 1)[1].cpu().tolist()
        dict_predict['stem'].extend(stems)
        dict_predict['prediction'].extend(prediction)

    df_predict = pd.DataFrame(dict_predict)
    df_predict['isCP'] = df_predict['prediction'].map({0:"false", 1:"true"})
    df_predict.to_csv(file_cppredict_csv, index=None)




if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--folder_embeddings', 
        type = str, 
        default = None,
    )
    args = parser.parse_args()
    run(args.folder_embeddings)

