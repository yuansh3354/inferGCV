import os
import math
import psutil
import argparse
import itertools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
def get_args():
    """
    Get arguments for the main scMutation script
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--gene', type=str, default='BRAF',
                        help='One gene symbols')
    parser.add_argument('-s', '--sc_expr', type=str, default='ref/GSE116237_Candidate_cancer_expr.h5',
                        help='normalize(1e6) and log1p scingle cell expression datas(cells * genes), support .h5 file the key of sc-expression store in "expr" ')
    parser.add_argument('-b', '--bulk_expr', type=str, default='ref/scMutation_ref.h5',
                        help='Bulk reference data. If user want to use own data, plz read TCGA_prepare.ipynb')
    parser.add_argument('-o', '--output', type=str, default='example',
                        help='Output filename, contain each cell\'s target gene mutation information and embedding')
    
    parser.add_argument('--batch_size', type=int, default=128,
                        help='Training loop of batch size')
    parser.add_argument('--max_epoch', type=int,  default=50,
                        help='Number of training loops')
    parser.add_argument('--seed', type=int, default=42,
                        help='Set random seed of all ')
    parser.add_argument('--lr', type=float, default=1e-4,
                        help='Set learning rate for training')
    parser.add_argument('--use_cpu', action='store_true',
                    help='Decision to use cpu(default use GPU) to training models')
    parser.add_argument('--cl_weight', type=float, default=0.9,
                    help='11111111111111')                   
    parser.add_argument('--n_neighbors_pseudo', type=int, default=30,
                    help='n_neighbors of sklearn KNN in pseudo')   

    args = parser.parse_args()
    return args

def check_mem():
    mem = psutil.virtual_memory()
    print("total memory: %.3f Gb" %(float(mem.total) / (1024*1024*1024)))
    print("used memory: %.3f Gb" %(float(mem.used) / (1024*1024*1024)))
    print("free memory: %.3f Gb" %(float(mem.free) / (1024*1024*1024)))

def p_adjust_bh(p):
    p = np.asfarray(p)
    by_descend = p.argsort()[::-1]
    by_orig = by_descend.argsort()
    steps = float(len(p)) / np.arange(len(p), 0, -1)
    q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
    return q[by_orig]

def integrate_copy_number(y, cancer_genes_df, genes, loss_df, gain_df):
    """
    Function to integrate copy number data to define gene activation or gene
    inactivation events. Copy number loss results in gene inactivation events
    and is important for tumor suppressor genes while copy number gain results
    in gene activation events and is important for oncogenes.

    Arguments:
    y - pandas dataframe samples by genes where a 1 indicates event
    cancer_genes_df - a dataframe listing bona fide cancer genes as defined by
                      the 20/20 rule in Vogelstein et al. 2013
    genes - the input list of genes to build the classifier for
    loss_df - a sample by gene dataframe listing copy number loss events
    gain_df - a sample by gene dataframe listing copy number gain events
    """

    # Find if the input genes are in this master list
    genes_sub = cancer_genes_df[cancer_genes_df['Gene Symbol'].isin(genes)]

    # Add status to the Y matrix depending on if the gene is a tumor suppressor
    # or an oncogene. An oncogene can be activated with copy number gains, but
    # a tumor suppressor is inactivated with copy number loss
    tumor_suppressor = genes_sub[genes_sub['Classification*'] == 'TSG']
    oncogene = genes_sub[genes_sub['Classification*'] == 'Oncogene']

    copy_loss_sub = loss_df[tumor_suppressor['Gene Symbol']]
    copy_gain_sub = gain_df[oncogene['Gene Symbol']]

    # Append to column names for visualization
    copy_loss_sub.columns = [col + '_loss' for col in copy_loss_sub.columns]
    copy_gain_sub.columns = [col + '_gain' for col in copy_gain_sub.columns]

    # Add columns to y matrix
    y = y.join(copy_loss_sub)
    y = y.join(copy_gain_sub)

    # Fill missing data with zero (measured mutation but not copy number)
    y = y.fillna(0)
    y = y.astype(int)
    return y

class Logger():
    def __init__(self, filename):
        self.name = filename
        self.file = open(filename, "a+", encoding='utf-8')
        self.alive = True
        self.stdout = sys.stdout
        sys.stdout = self

    def __enter__(self):
        pass

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def __del__(self):
        self.close()

    def close(self):
        if self.alive:
            sys.stdout = self.stdout
            self.file.close()
            self.alive = False

    def write(self, data):
        self.file.write(data)
        self.stdout.write(data)

    def flush(self):
        self.file.flush()
        
        
def standardization(data):
    mu = np.mean(data, axis=0)
    sigma = np.std(data, axis=0)
    return (data - mu) / sigma        
        
        
