# Create: Yuan.Sh
# Date: 2022-03-06 13:06:21
# Email: yuansh3354@163.com
# Blog: https://blog.csdn.net/qq_40966210
# Fujian Medical University

# =============================================================================
# Step.01 imports 
# =============================================================================
import os
import gc
import umap
import torch
import warnings
import argparse
import numpy as np
import pandas as pd
import torchmetrics
import scanpy as sc
import torch.nn as nn
import seaborn as sns
from pathlib import Path
from tqdm.auto import tqdm
import pytorch_lightning as pl
import matplotlib.pyplot as plt
from collections import Counter
import torch.utils.data as data
import torch.nn.functional as F
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import train_test_split
from utils.Model import *
from utils.myexptorch import *
from utils.myfunction import * 

# =============================================================================
# Step.02 Get args
# =============================================================================
gene = get_args().gene
scRNA_file = get_args().sc_expr
scMutation_ref = get_args().bulk_expr
output = get_args().output
cl_weight = get_args().cl_weight
batch_size = get_args().batch_size
max_epoch = get_args().max_epoch
seed = get_args().seed
use_cpu = get_args().use_cpu
lr = get_args().lr
n_neighbors_pseudo = get_args().n_neighbors_pseudo
if not os.path.exists(output):
    os.mkdir(output)
    
with Logger(output + '/' + gene+'_Mutation_Detail.txt'):
    print(get_args())
    print('\n')
# =============================================================================
# Step.03 Get target gene mutation information in refdata
# =============================================================================
rnaseq_full_df = pd.read_hdf(scMutation_ref, key='rnaseq_full_df')
mutation_df = pd.read_hdf(scMutation_ref,key='mutation_df')
sample_freeze = pd.read_hdf(scMutation_ref,key='sample_freeze')
mut_burden = pd.read_hdf(scMutation_ref,key='mut_burden')
copy_loss_df = pd.read_hdf(scMutation_ref,key='copy_loss_df')
copy_gain_df = pd.read_hdf(scMutation_ref,key='copy_gain_df')
cancer_genes = pd.read_table('ref/vogelstein_cancergenes.tsv')

print('Reference bulk data import aready.')
print('Memory infomation:')
check_mem()
print('Get {} mutation information.'.format(gene))
y = pd.DataFrame(mutation_df[gene])
y = integrate_copy_number(y=y, cancer_genes_df=cancer_genes,
                          genes=[gene], loss_df=copy_loss_df,
                          gain_df=copy_gain_df)
y = y.assign(total_status=y.max(axis=1))
file_names = output + '/' + gene+'_Bulk_Mutation_information.csv'
y.to_csv(file_names)
ratio = (sum(y['total_status']==1)/y.shape[0])*100
loss_weight = torch.tensor([ratio, 100 - ratio])

with Logger(output + '/' + gene+'_Mutation_Detail.txt'):
    print("{} mutation ratio in bulk tissue is {}%".format(gene, round(ratio,2)))
    print('\n')

print('free memory...')
del mutation_df, sample_freeze, mut_burden, copy_loss_df, copy_gain_df, cancer_genes
_ = gc.collect()
check_mem()

# =============================================================================
# Step.04 Data prepare for scMutation
# =============================================================================
# get Sample_id
sample_index = rnaseq_full_df.index
# read scRNA-seq
sce_full_expr = pd.read_hdf(scRNA_file, key='expr')
sc_genes = sce_full_expr.columns
sce_index = sce_full_expr.index

# get comm genes 
genes = rnaseq_full_df.columns.intersection(sc_genes)
with Logger(output + '/' + gene+'_Mutation_Detail.txt'):
    print("Total using genes {}".format(len(genes)))
    print('\n')
# filter expr data 
sce = sce_full_expr[genes]
bulk = rnaseq_full_df[genes]
groups = np.array(["bulk-seq", 'input_sce'])
groups = np.repeat(groups, [bulk.shape[0], sce.shape[0]], axis=0)
# to tensor 
bulk = DataFrame_normalize(bulk,'std')
sce = DataFrame_normalize(sce,'std')
bulk = expToTorch(bulk)
bulk = bulk.view(bulk.shape[0],1,bulk.shape[1])
label = toOneHot(y['total_status'],n_class=2).long()
sce = expToTorch(sce)
sce = sce.view(sce.shape[0],1,sce.shape[1])
pl.utilities.seed.seed_everything(seed=seed)
bulk_iter = makeDataiter(bulk, label, batch_size=batch_size, shuffle=True)
#assert len(genes) >= 3000, "The number of genes is not enough (N >= 3000) for training scMutation."

# =============================================================================
# Step.05 Training model to get sc pseduo label
# =============================================================================
use_gpu = True if not use_cpu else False
pl.utilities.seed.seed_everything(seed=seed)
model = scMutation_pseudo_label(sce=sce, gene_num=len(genes), cl_weight=cl_weight,lr=lr, loss_weight = loss_weight,use_gpu=use_gpu)
work_dir = './'
files='pseudo'
# output file
OUTPUT_DIR = './lightning_logs'
tb_logger = pl.loggers.TensorBoardLogger(save_dir='./',
                                         name='./' + files)

# set check point to choose best model
checkpoint_callback = pl.callbacks.ModelCheckpoint(
    dirpath=tb_logger.log_dir,
    filename='{epoch}-{val_acc:.4f}',
    save_top_k=1,  # keep top 5
    monitor='val_acc',  # check acc
    mode='max'
)
# train loop
if use_cpu:
    pl.utilities.seed.seed_everything(seed=seed)
    trainer = pl.Trainer(callbacks=[checkpoint_callback],max_epochs=max_epoch) 
else:
    pl.utilities.seed.seed_everything(seed=seed)
    trainer = pl.Trainer(gpus=-1,callbacks=[checkpoint_callback],
                         max_epochs=max_epoch)
pl.utilities.seed.seed_everything(seed=seed)
trainer.fit(model, bulk_iter, bulk_iter)

# =============================================================================
# Step.06 Get scMutation information
# =============================================================================
out = Path(tb_logger.log_dir)
ckpt_list = [ckpt.stem + '.ckpt' for ckpt in out.iterdir()]
ckpt = ckpt_list[np.argmax(
    [float(i.split('-')[0].split('epoch=')[1]) for i in ckpt_list])]

model = None
model = scMutation_pseudo_label(sce=sce, gene_num=len(genes), cl_weight=cl_weight, lr=lr, loss_weight = loss_weight,use_gpu=use_gpu)

ckpt_path = str(out) + '/' + ckpt
with Logger(output + '/' + gene+'_Mutation_Detail.txt'):
    print('useing pseudo_ckpt:', ckpt_path)
    print('\n')

model = model.load_from_checkpoint(ckpt_path,sce=sce, gene_num=len(genes), cl_weight=cl_weight, lr=lr, loss_weight = loss_weight,use_gpu=use_gpu)
with torch.no_grad():
    bulk_embeddings, bulk_probs,_ = model.encoder(bulk)
    sce_embeddings, sce_probs,_ = model.encoder(sce)
    embeddings =  np.concatenate((bulk_embeddings.cpu().numpy(), 
                                  sce_embeddings.cpu().numpy()))
pl.utilities.seed.seed_everything(seed=seed)
neigh = KNeighborsClassifier(n_neighbors=n_neighbors_pseudo)
neigh.fit(bulk_embeddings, y['total_status'])
pseudo_label = neigh.predict(sce_embeddings)
sc_mutation_ratio = sum(pseudo_label)/len(pseudo_label)
with Logger(output + '/' + gene+'_Mutation_Detail.txt'):
    print("single-cell {} pseudo mutation ratio: {:.2f}%".format(gene,sc_mutation_ratio*100))
    print('\n')

target = np.append(np.array(y['total_status']),pseudo_label)
probs = np.append(bulk_probs.cpu().numpy(), sce_probs.cpu().numpy(),axis=0)

df = pd.DataFrame(embeddings)
all_sample_index = sample_index.append(sce_index)
df.index = all_sample_index
data_label = np.array(["bulk-seq", "sce-seq"])
df['data'] = np.repeat(data_label, [bulk_embeddings.shape[0], sce_embeddings.shape[0]], axis=0)
df['predicted'] = np.array(pd.DataFrame(target).replace([0,1],  ['WT','MUT']))
bulk_sce_mix_embedding = output + '/ScMutation_' + gene + '_mutation_bulk_embeddings.csv'
df.to_csv(bulk_sce_mix_embedding)

# =============================================================================
# Step.09 View result 
# =============================================================================
pca = PCA(n_components=3)
pca.fit(embeddings.T)
pca_df = pd.DataFrame()
pca_df.index = all_sample_index
pca_df['PCA1'] = pca.components_[0]
pca_df['PCA2'] = pca.components_[1]
pca_df['PCA3'] = pca.components_[2]
data_label = np.array(["bulk-seq", "sce-seq"])
pca_df['data'] = np.repeat(data_label, [bulk_embeddings.shape[0], sce_embeddings.shape[0]], axis=0)
#pca_df['WT_prob'] = pd.DataFrame(probs).loc[:,0]
#pca_df['MUT_prob'] = pd.DataFrame(probs).loc[:,1]
pca_df['predict'] = np.array(pd.DataFrame(target).replace([0,1],  ['WT','MUT']))

plt.figure(figsize=(8, 6))
if pca_df[pca_df['predict'] == 'WT'].shape[0] != 0:
    plt.scatter(pca_df[pca_df['predict'] == 'WT']['PCA1'],
                pca_df[pca_df['predict'] == 'WT']['PCA2'],
                    c='#698EC3',
                    label='WT')
if pca_df[pca_df['predict'] == 'MUT'].shape[0] != 0:
    plt.scatter(pca_df[pca_df['predict'] == 'MUT']['PCA1'],
                    pca_df[pca_df['predict'] == 'MUT']['PCA2'],
                    c='#E64B35B2',
                    label='MUT')
plt.legend()
plt.savefig(output+'/'+gene +'_2D_embeddings.jpg', dpi = 300)
plt.clf()

from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure(figsize=(8,6))
ax = fig.gca(projection='3d')

if pca_df[pca_df['predict'] == 'WT'].shape[0] != 0:
	ax.scatter(pca_df[pca_df['predict'] == 'WT']['PCA1'],
	pca_df[pca_df['predict'] == 'WT']['PCA2'],
	pca_df[pca_df['predict'] == 'WT']['PCA3'],
	c='#698EC3',label='WT')
	
if pca_df[pca_df['predict'] == 'MUT'].shape[0] != 0:
	ax.scatter(pca_df[pca_df['predict'] == 'MUT']['PCA1'],
	pca_df[pca_df['predict'] == 'MUT']['PCA2'],
	pca_df[pca_df['predict'] == 'MUT']['PCA3'],
	c='#E64B35B2',label='MUT')
plt.legend()
plt.savefig(output+'/'+gene +'_3D_embeddings.jpg', dpi = 300)
pca_df.to_csv(output + '/ScMutation_' + gene + '_input_sce_embeddings.csv')

plt.figure(figsize=(8, 6))
pca_df_ = pca_df.copy()
pca_df = pca_df_[pca_df_['data']=='bulk-seq']
if pca_df[pca_df['predict'] == 'WT'].shape[0] != 0:
    plt.scatter(pca_df[pca_df['predict'] == 'WT']['PCA1'],
                pca_df[pca_df['predict'] == 'WT']['PCA2'],
                    c='#698EC3',
                    label='WT')
if pca_df[pca_df['predict'] == 'MUT'].shape[0] != 0:
    plt.scatter(pca_df[pca_df['predict'] == 'MUT']['PCA1'],
                    pca_df[pca_df['predict'] == 'MUT']['PCA2'],
                    c='#E64B35B2',
                    label='MUT')
plt.legend()
plt.savefig(output+'/'+gene +'_2D_bulk_embeddings.jpg', dpi = 300)
plt.clf()
pca_df = pca_df_[pca_df_['data']=='sce-seq']
if pca_df[pca_df['predict'] == 'WT'].shape[0] != 0:
    plt.scatter(pca_df[pca_df['predict'] == 'WT']['PCA1'],
                pca_df[pca_df['predict'] == 'WT']['PCA2'],
                    c='#698EC3',
                    label='WT')
if pca_df[pca_df['predict'] == 'MUT'].shape[0] != 0:
    plt.scatter(pca_df[pca_df['predict'] == 'MUT']['PCA1'],
                    pca_df[pca_df['predict'] == 'MUT']['PCA2'],
                    c='#E64B35B2',
                    label='MUT')
plt.legend()
plt.savefig(output+'/'+gene +'_2D_sce_embeddings.jpg', dpi = 300)
plt.clf()

