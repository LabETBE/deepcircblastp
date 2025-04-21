# DeepCircBLASTp

DeepCircBLASTp, an advanced framework that integrates evolutionary strategies, deep learning, and BLASTp, consists of two specialized modules CircBLASTp and CPDetector. The former optimized CP sequence searches with best-in-class performance in CP protein identification compared to four classical homologous sequence search tools. The latter demonstrates state-of-the-art capabilities in detecting and identifying CP pairs at the sequence level. Compared to existing tools, DeepCircBLASTp demonstrates exceptional efficacy in identifying CP pairs and significantly improves the retrieval of homologous sequences. Notably, it has uncovered numerous previously unreported CP pairs within the EC 3.1.1 (Carboxylic ester hydrolases) family and successfully identified a CP pair (UniProt ID: A0A382XFA1) associated with 5ZOA (PDB ID, a PETase), along with their evolutionary intermediate (UniProt ID: A0A835BI20). To extend its utility, we constructed CPDB4Enzyme, the largest CP resource to date, comprising over 638,000 CP pairs. This database enables systematic exploration of CP relationships across enzyme families, including both highly represented families such as EC 2.7.11.1 (kinases) and newly discovered ones like PETase (EC 3.1.1.101) and MHETase (EC 3.1.1.102). DeepCircBLASTp marks a significant advancement in protein sequence analysis and enzyme engineering, providing valuable insights into evolutionary biology. The tool is freely accessible at https://letbe.jiangnan.edu.cn/bioserver/deepcircblastp, serving as a powerful resource for exploring CP and homologous protein sequences. And the CPDB4Enzyme is freely accessible at https://letbe.jiangnan.edu.cn/bioserver/cpdb4enzymes. 

# os version
OS: Deepin 20.9 5.15.77-amd64-desktop

# Dependencies
- python >= 3.9
- pytorch_cuda
- transformers
- biopython>=1.81
- numpy=1
- pandas
- h5py
- requests
- seqkit (bioconda) 


# Installation steps
1. Download the required model from the following address: [models_download](https://zenodo.org/records/10.5281/zenodo.15241259.), and put it in the `deepcircblastp/models` directory 

2. Download the required `fasta` sequence as the search source, such as uniprot.fasta, and save it in the `database` directory
   


# Usage examples
```bash
python run.py query/5zoa.fasta database/test.fasta database 
```



# License
This work is licensed under the MIT license.









