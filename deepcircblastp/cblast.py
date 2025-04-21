import json
import h5py, re
import numpy as np 
import pandas as pd 
from pathlib import Path
from Bio import SeqIO, Align
from bisect import bisect_right


 

# ========================================================================
#                            class MakeCircDB
# ========================================================================
class MakeCircDB():
    def __init__(self, db_path, matrix_file):
        self.word_len = 3
        self.db_path = Path(db_path)
        self.file_h5 = self.db_path.parent / 'indexed_db.h5'
        self.seq_id = self.db_path.parent / 'seq_id.npy'
        self.seq_name = self.db_path.parent / 'seq_name.npy'
        self.file_db_list = self.db_path.parent / 'db_list.txt'
        with open(matrix_file, 'r') as f:
            self.letters = f.readline()[:-1]

    
    def word_to_num(self):
        word_num = {}
        for i,l in enumerate(self.letters):
            word_num[l] = i
        return word_num


    def word_to_idx(self, word_num, word):
        word_index = 0
        for i,l in enumerate(word,1):
            word_index += (len(word_num))**(self.word_len-i)*word_num[l]
        return word_index
        

    def build_index(self,seq_name,seq,idx,indexed_db,seq_idx,seq_names,word_num,):
        seq_idx.append(idx)
        seq_names.append(seq_name)
        for i in range(len(seq)-self.word_len):
            word = seq[:self.word_len]
            seq = seq[1:] 
            word_index = self.word_to_idx(word_num,word)
            indexed_db[word_index].append(idx)
            idx += 1
        return idx
            

    def build_lib(self, db_path_list):
        
        indexed_db = [[] for _ in range(len(self.letters)**self.word_len)]

        word_num = self.word_to_num()
        
        seq_idx = []
        seq_names = []
        idx = 0
        
        with open(self.db_path.parent/'combined.txt', 'w') as f:
            for db_path in db_path_list:
                if len(db_path) == 0:
                    continue
                seq_name = db_path.split('/')[-1]
                with open(db_path+".fasta") as db:
                    line = db.readline()
                    seq = ""
                    for line in db:
                        if re.match(r">",line):
                            idx = self.build_index(seq_name,seq,idx,indexed_db,seq_idx,seq_names,word_num,) + self.word_len
                            seq = ""
                        else:
                            seq += line[:-1]
                    idx = self.build_index(seq_name,seq,idx,indexed_db,seq_idx,seq_names,word_num,) + self.word_len
                    f.write(seq)
        
        seq_idx.append(idx)

        # 写 seed database    
        p = 0
        with h5py.File(self.file_h5, "w") as f:
            for i in self.letters:
                for j in self.letters:
                    for k in self.letters:
                        if len(indexed_db[p]) > 0:
                            f.create_dataset(i+j+k, data = indexed_db[p])
                            p += 1
                        else:
                            indexed_db.pop(p)

        seq_idx = np.array(seq_idx)
        seq_names = np.array(seq_names)
        np.save(self.seq_id, seq_idx) # seq_id.npy has saved 
        np.save(self.seq_name, seq_names) # seq_name.npy has saved
        
        return indexed_db, seq_idx
        

    # make db
    def db_build(self):
        fasta_list = []
        for file in self.db_path.glob('*.fasta'):
            content = str(self.db_path / file.stem)
            fasta_list.append(content)
        
        self.build_lib(fasta_list)
        



# ========================================================================
#                       class GetAlignmentInfo
# ========================================================================
class GetAlignmentInfo():
    def __init__(self, str_query, str_sbjct, db_length, k=0.1, l=0.5):
        self.str_query = str_query
        self.str_sbjct = str_sbjct
        self.k = k 
        self.l = l 
        self.db_length = db_length
        self.E, self.identity, self.alignment = self.cal_i_e()

        
    def cal_i_e(self):
        aligner = Align.PairwiseAligner(
            mode = 'global', 
            open_gap_score = -11,
            extend_gap_score = -1,
        )
        seq_query = self.str_query.replace('-', '')
        seq_sbjct = self.str_sbjct.replace('-', '')
        length = len(seq_query)
        alignments = aligner.align(seq_query, seq_sbjct)
        alignment = alignments[0]
        num_same = alignment.counts().identities
        identity = num_same / alignment.shape[1]
        identity = round(identity, 4)
        self.score = alignment.score
        E = (self.k * length * self.db_length) * np.exp(-self.l * self.score)
        return E, identity, alignment




# ========================================================================
#                            class CircleBLASTp
# ========================================================================
class CircleBLASTp():
    def __init__(
            self, file_query, db_path, matrix_file, 
            gap_penalty=-11, extend_penalty=-1,
        ):
        self.word_len = 3
        self.str_query = str(SeqIO.read(file_query,'fasta').seq)
        self.db_path = db_path
        self.gap_penalty = gap_penalty
        self.extend_penalty = extend_penalty
        
        self.matrix_file = matrix_file
        self.matrix = {}
        with open(matrix_file) as f:
            letters = f.readline()[:-1]
            for line in f:
                line = line.split(',')
                self.matrix[line[0]] = int(line[1][:-1])
        self.word_num = {}
        for i,l in enumerate(letters):
            self.word_num[l] = i # each AA matched to a number (0~24)
        
        self.indexed_db = h5py.File(Path(self.db_path).parent/"indexed_db.h5", "r")

        # store the position of each sequence
        self.seq_pos = np.load(Path(self.db_path).parent/"seq_id.npy")
        # store the gene name of each sequence
        self.seq_names = np.load(Path(self.db_path).parent/"seq_name.npy")
        # store the amino acid residues of all sequence
        with open(Path(self.db_path).parent/'combined.txt') as f:
            self.seq_all = f.readline()

        self.X = -100 # during extending, lowest score accepted below best score
        self.A = 70  # max distance between two hits when they can be considered as one hit
        self.H = 0   # during denoising, min number of seeds on a hit start point that can be collected for furthur analysis
        self.N = 0   # during hit calling, min continuous seeds that can be considered a hit


    def split_query(self):
        seeds = [] 
        str_query = self.str_query
        for i in range(len(str_query) - self.word_len + 1):
            seed = str_query[:self.word_len]
            str_query = str_query[1:]
            seeds.append(seed)
        return seeds
    

    def get_seed_position(self, seeds):
        positions = []
        for seed in seeds:
            try:
                pos = self.indexed_db[seed][:]
            except KeyError:
                pos = []
            positions.append(pos)
        return positions


    def get_seq_end(self,pos):
        idx = bisect_right(self.seq_pos, pos)
        seq_end = self.seq_pos[idx] if idx < len(self.seq_pos) else float('inf')
        return seq_end
     

    def find_hit(self, positions):
        occurs = {}
        i = 0
        for position in positions:
            for pos in position:
                try:
                    occurs[pos-i].append(i)
                except KeyError:
                    occurs[pos-i] = [i]
            i += 1
            
        occurs_count = []
        occurs_clear = {}
        for key,value in occurs.items():
            if len(value) > self.H:
                occurs_count.append((key,len(value)))
                occurs_clear[key] = value
        occurs.clear()
        occurs_count.sort()
                    
        def next_gene(tmp,hits_list,last_hit_finished,hit_start,seq_end):
            if hits[j-1] - hit_start > self.N:
                tmp.append((hit_start,hit_start+occurs_count[i][0],hits[j-1]-hit_start+self.word_len))
            tmp, last_hit_finished, hits_list = add_hit(tmp, last_hit_finished, hits_list)
            tmp = []
            last_hit_finished = True
            hit_start = hits[j]
            seq_end = self.get_seq_end(hit_start+occurs_count[i][0])
            
            return tmp,hits_list,last_hit_finished,hit_start,seq_end
        
        def add_hit(tmp, last_hit_finished, hits_list):
            if len(tmp) == 0:
                return tmp, last_hit_finished, hits_list
            
            if last_hit_finished:
                hits_list.append(tmp)
                last_hit_finished = False
            else:
                last_seq_end = hits_list[-1][-1][2] + hits_list[-1][-1][1]
                gap = tmp[0][1] - last_seq_end
                if abs(gap-self.word_len) < self.A and seq_end == self.get_seq_end(last_seq_end):
                    while len(tmp) > 0:
                        last_hit_end = hits_list[-1][-1][0] + hits_list[-1][-1][2]
                        if last_hit_end >= tmp[0][0]:
                            if last_hit_end >= tmp[0][0] + tmp[0][2]:
                                tmp = tmp[1:]
                            else:
                                overlap = last_hit_end-tmp[0][0]
                                tmp[0] = (tmp[0][0]+overlap, tmp[0][1]+overlap, tmp[0][2]-overlap)
                                hits_list[-1] += tmp
                                break
                        else:
                            hits_list[-1] += tmp
                            break
                            
                else:
                    hits_list.append(tmp)
            return tmp, last_hit_finished, hits_list
            
        last_hit_finished = True
        hits_list = []
        for i in range(len(occurs_count)):
            hits = occurs_clear[occurs_count[i][0]]
            hit_start = hits[0]
            seq_end = self.get_seq_end(hit_start+occurs_count[i][0])
            tmp = []
            
            for j in range(1,len(hits)):
                if hits[j]-hits[j-1] > 1:
                    if hits[j] + occurs_count[i][0] + self.word_len > seq_end:
                        tmp,hits_list,last_hit_finished,hit_start,seq_end = next_gene(tmp,hits_list,last_hit_finished,hit_start,seq_end)
                    elif hits[j] < hits[j-1]+self.A:
                        if hits[j-1] - hit_start > self.N:
                            tmp.append((hit_start,hit_start+occurs_count[i][0],hits[j-1]-hit_start+self.word_len))
                    else:
                        tmp,hits_list,last_hit_finished,hit_start,seq_end = next_gene(tmp,hits_list,last_hit_finished,hit_start,seq_end)
                    hit_start = hits[j]
                    if hits[j] < hits[j-1]+self.A and hits[j-1] - hit_start > self.N:
                        tmp.append((hit_start, hit_start+occurs_count[i][0],hits[j-1]-hit_start+self.word_len))
            
            tmp, last_hit_finished, hits_list = add_hit(tmp, last_hit_finished, hits_list)
            
        if len(hits_list)==0:
            return -1
        else:
            return hits_list


    def extend_hit(self, hits_list, folder_alignments):
        seq_pos = self.seq_pos
        db_length = seq_pos[-1]

        hits_dict = {}

        for hit in hits_list:
            access = self.seq_names[bisect_right(seq_pos, hit[0][1])-1]
            hit_array = np.zeros([1,3])
            for tup in hit:
                hit_array = np.concatenate(
                    [np.array([tup]), hit_array],
                    axis = 0,
                )
            hit_array = np.delete(hit_array, -1, 0)
            hit = hit_array.tolist()
            if access not in hits_dict.keys():
                tmp = [] 
                for i in hit:
                    tmp.append(i)
                hits_dict[access] = tmp
            else:
                for i in hit:
                    hits_dict[access].append(i)

        df_dict = {
            'access': [],
            'iden_circle': [],
            'evalue_circle': [],
        }

        for key, value in hits_dict.items():
            alignments = []
            hits_array = np.array(value)
            tmp_array = hits_array.reshape(hits_array.shape[0], hits_array.shape[-1])
            last_col = tmp_array[:, -1]
            max_indices = np.where(last_col == np.max(last_col))[0]
            rows_max = tmp_array[max_indices, :].astype(int)
            max_identity = 0
            for hit in rows_max:
                circle_site = hit[0]
                str_query_new = self.str_query[circle_site:] + self.str_query[:circle_site]
                start = seq_pos[bisect_right(seq_pos, hit[1])-1]
                end = seq_pos[bisect_right(seq_pos, hit[1])]
                str_sbjct = self.seq_all[start:end]
                circle_site = hit[1] - start
                str_sbjct_new = str_sbjct[circle_site:] + str_sbjct[:circle_site]
                info = GetAlignmentInfo(str_query_new, str_sbjct_new, db_length)
                if info.identity > max_identity:
                    max_identity = info.identity
                    best_info = info

            # alignment 数据
            origin_info = GetAlignmentInfo(self.str_query, str_sbjct, db_length)
            alignments.extend([origin_info.alignment[0], origin_info.alignment[1]])
            alignments.extend([best_info.alignment[0], best_info.alignment[1]])
            data_json =json.dumps(alignments, indent=4)
            file_json = Path(folder_alignments) / f'{key}.json'
            file_json.write_text(data_json)

            # 绘图数据
            df_dict['access'].append(key)
            df_dict['iden_circle'].append(best_info.identity)
            df_dict['evalue_circle'].append(best_info.E)
        
        df = pd.DataFrame(df_dict)

        return df


    def run(self, folder_alignments, file_csv):
        seeds = self.split_query()
        positions = self.get_seed_position(seeds)
        hits_list = self.find_hit(positions)
        if hits_list == -1:
            print("Not found")
            return -1
        df = self.extend_hit(hits_list, folder_alignments)
        df.to_csv(file_csv, index=None)
        
            

