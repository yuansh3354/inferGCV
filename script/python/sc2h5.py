#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  7 00:02:53 2022

@author: yuansh
"""
import pandas as pd 
import numpy as np
import os
import re
import warnings
import argparse
import scanpy as sc
import matplotlib.pyplot as plt
import yaml
from multiprocessing import cpu_count
n_jobs = cpu_count()
warnings.filterwarnings("ignore")
warnings.simplefilter('ignore')
plt.rcParams['savefig.dpi'] = 300 #图片像素

# =============================================================================
# # my imports & settings
# =============================================================================

n_class = 2

# =============================================================================
# # get configs
# =============================================================================

parser = argparse.ArgumentParser()
parser.add_argument('--configs', type=str, default="sc2h5.yaml", help='whether to use configs')
opt = parser.parse_args()
config_path = opt.configs
seed = 42
# config_path = "configs/CaSee_Model_configs.yaml"
config_dict = yaml.safe_load(open(config_path, 'r'))
names = globals()
print('# ============================================================================= #','\n')
print("starting program...\n")
print("setting configs...\n")
# =============================================================================
# # setting configs
# =============================================================================
# data_arguments
work_dir = config_dict['data_arguments']['work_dir']
Counts_expr = config_dict['data_arguments']['Counts_expr']
use_others = config_dict['data_arguments']['use_others']
remove_genes = config_dict['data_arguments']['remove_genes']

# cell cluster
cell_annotion = config_dict['cell_annotion']

# Marker_genes
if cell_annotion:
    T_cell_marker = config_dict['Marker_genes']['T_cell']
    Fibroblast_cell_marker = config_dict['Marker_genes']['Fibroblast']
    Myeloid_cell_marker = config_dict['Marker_genes']['Myeloid']
    B_cell_marker = config_dict['Marker_genes']['B_cell']
    Endothelial_cell_marker = config_dict['Marker_genes']['Endothelial']
    Mast_cell_marker = config_dict['Marker_genes']['Mast']
    DC_cell_marker = config_dict['Marker_genes']['DC']
    Candidate_Cancer_cell_marker = config_dict['Marker_genes']['Cancer']

# save_files
files = config_dict['save_files']['files']

# ============================================================================
# create scRNA epxr
# =============================================================================
print('# ============================================================================= #','\n')
print("prepare single-cell count matrix...\n")

if not os.path.exists(work_dir+files+'_h5_file/'+files+'_Candidate_cancer_expr.h5'):
    print(work_dir+files+'_h5_file/'+files+'_Candidate_cancer_expr.h5   does not exit, starting zip h5file...\n')
    # =============================================================================
    # # create scanpy object
    # =============================================================================
    
    if not os.path.exists(work_dir + files+'_h5_file'):
        os.mkdir(work_dir + files+'_h5_file')
    if not os.path.exists(work_dir + files+'_h5_file/'+files+'_count_expr.h5'):
        print("Step.00 h5 file not exits, loading count matrix and zip to h5 file...\n")
        if Counts_expr.endswith('.txt'):
            expr = pd.read_table(work_dir + Counts_expr, index_col=0).T
        elif Counts_expr.endswith('.csv'):
            expr = pd.read_csv(work_dir + Counts_expr, index_col=0).T
        else: print('count matrix must be csv or txt format')
        
        h5 = pd.HDFStore(work_dir + files+'_h5_file/'+files+'_count_expr.h5', 'w',
                         complevel=4, complib='blosc')
        h5['expr'] = expr
        h5.close()
        #print ('creating h5 file and loading time：',run_time,' seconds')
    else: 
        print("Step.00 h5 file exits, loading h5 file...\n")
        expr = pd.read_hdf(work_dir + files+'_h5_file/'+files+'_count_expr.h5',key='expr')
        #print ('loading time：',run_time,' seconds')
        
    print("Step.01 create scanpy obj and do tsne ... \n")
    # 0.loading caount matrix
    if cell_annotion:
        marker_genes = [i for i in names if '_cell_marker' in i]
        marker_genes.sort()
        cell_types = [re.sub("_marker","",i) for i in names if '_cell_marker' in i]
        cell_types.sort()
        gene_list = [j for i in marker_genes for j in names[i]]
        
    # 1.create scanpy object
    cellinfo = pd.DataFrame(expr.index, index=expr.index, columns=['sample_index'])
    geneinfo = pd.DataFrame(expr.columns,
                            index=expr.columns,
                            columns=['genes_index'])
    sce = sc.AnnData(expr, obs=cellinfo, var=geneinfo)
    sce.var_names_make_unique()
    sce.obs_names_make_unique()
    
    # 2.clean data calculate some gene information
    mt = sce.var_names[sce.var_names.str.match(r'^MT-')]  # 线粒体DNA
    rp = sce.var_names[sce.var_names.str.match(r'^RP[SL][0-9]')]  # 核糖体DNA
    ercc = sce.var_names[sce.var_names.str.match(r'^ERCC-')]  # 外源DNA
    ncRNA = sce.var_names[sce.var_names.str.match(
        r'^[A-Z][A-Z][0-9]*\.[0-9]')]  # 匹配有小数点的基因
    LOC = sce.var_names[sce.var_names.str.match(
        r'(^LOC|LINC)[1-9]*')]  # 匹配未知的LC rna
    
    # 3.statistics
    print("Step.02 remove genes ... \n")
    print("total number of MT: {}".format(len(mt)))
    print("total number of RT: {}".format(len(rp)))
    print("total number of ERCC: {}".format(len(ercc)))
    print("total number of non-coding RNA: {}".format(len(ncRNA)))
    print("total number of LNC_RNA: {}".format(len(LOC)),'\n')
    
    # 4.remove genes 
    print('remove gene option: ',remove_genes,'\n')
    if remove_genes:
        ids = list(rp) + list(ercc) + list(ncRNA) + list(LOC) + list(mt)
        print("Total Remove Genes {}".format(len(ids)))
        use_genes = sce.var.index.values
        ids = set(use_genes) - set(ids)
        print("Number of use Genes {}".format(len(ids)),'\n')
        sce = sce[:,sce.var.genes_index.isin(ids)]

    # 5.Stander workflow
    print('Step.03 running Stander workflow...\n')
    if not os.path.exists(work_dir + 'figures/'):
        os.mkdir(work_dir + 'figures/')
    sc.pp.normalize_total(sce, target_sum=1e6)
    sc.pp.log1p(sce)
    # 6. get marker gene list expr information
    if cell_annotion:
        sce.raw = sce
        sc.pp.highly_variable_genes(sce, n_top_genes=2000)
        sc.pp.scale(sce)
        sc.tl.pca(sce, svd_solver='arpack', use_highly_variable=True, random_state=seed)
        sc.pp.neighbors(sce, random_state=seed)
        sc.tl.leiden(sce, random_state=seed)
        sc.tl.tsne(sce, n_pcs=20,n_jobs=n_jobs)
    
        ids = sc.pl.dotplot(sce,
                            gene_list,
                            groupby='leiden',
                            return_fig=True)
        ids.savefig(filename=work_dir + 'figures/' + files + '_GeneMarkers_DotPlot.pdf',show=False)
        marker_expr = ids.dot_color_df
        marker_percent = ids.dot_size_df
    
    #print ('tsne finished using：',run_time,' seconds')
    
    # =============================================================================
    # cell_annotion
    # =============================================================================
    # first annotion
    if cell_annotion:
        print("Step.04 running cell annotion...\n")
        rest_cluster_id = []
        for cell_type,marker_gene in zip(cell_types,marker_genes):
            ids = cell_type + '_cluster_id'
            names[ids] = (marker_percent[names[marker_gene]] > 0.333) & marker_expr[names[marker_gene]] > 0
            names[ids] = list(names[ids][names[ids].sum(axis=1) > 0].index)
            print(ids, ':', names[ids])
            rest_cluster_id += names[ids]
            names[cell_type] = sce.obs.leiden[sce.obs.leiden.isin(names[ids])].index
        rest_cluster_id = list(set(sce.obs.leiden) - set(rest_cluster_id))
        print('rest', ':', rest_cluster_id)
        rest = sce.obs.leiden[sce.obs.leiden.isin(rest_cluster_id)].index
        
        check_unique = []
        num_cell_types = len(cell_types)
        if not os.path.exists(work_dir + 'figures/dotplot_figures/'):
            os.mkdir(work_dir + 'figures/dotplot_figures/')
            
        for i in range(num_cell_types-1):
            for j in range(i+1, num_cell_types):
                cluster_id = set(names[cell_types[i] + '_cluster_id']
                                 ).intersection(names[cell_types[j] + '_cluster_id'])
                if len(cluster_id) != 0:
                    check_unique.append([cell_types[i],cell_types[j],cluster_id])
        
        # annotion again, adjust and remove duplicated
        unique_cells = []
        for check_list in check_unique:
            cluster_list = check_list[2]
            num_cluster = len(check_list[2])
            cells = []
            for i in range(num_cluster):
                title = check_list[0] + ' vs ' + check_list[1]
                ids = sc.pl.dotplot(sce[sce.obs.leiden == list(cluster_list)[i], ],
                                    gene_list, groupby='leiden',
                                    title=title, return_fig=True)
                cluster_expr = ids.dot_color_df.stack().idxmax()
                cluster_precent = ids.dot_size_df.stack().idxmax()
                cells += [marker_gene[:-7] for marker_gene in marker_genes if cluster_expr[1] in names[marker_gene]]
                # ids.savefig(filename='figures/dotplot_figures/'+title +
                #             "_cluster " + str(list(cluster_list)[i]), show=False)
            unique_cells.append(cells)
            
            
        rest_annotion = []
        for ids in rest_cluster_id:

            ids = sc.pl.dotplot(sce[sce.obs.leiden == ids, ],
                                gene_list,
                                groupby='leiden',
                                return_fig=True)
            cluster_expr = ids.dot_color_df.stack().idxmax()
            cluster_precent = ids.dot_size_df.stack().idxmax()
            if (ids.dot_color_df.max().max() < -0.5) | (ids.dot_size_df.max().max() < 0.05):
                rest_annotion += ['Z_Candidate']
            else:
                rest_annotion += [
                    marker_gene[:-7] for marker_gene in marker_genes
                    if cluster_expr[1] in names[marker_gene]
                ]
                
        sce.obs['cell_type'] = 'Z_Candidate'
        for cell_type in cell_types:
            sce.obs.loc[names[cell_type],'cell_type'] = cell_type
        
        for cluster_id, cell_type in zip(check_unique,unique_cells):
            num = len(cell_type)
            for ids in range(num):
                cell_annotion = cell_type[ids]
                cell_cluster = list(cluster_id[2])[ids]
                sce.obs.loc[sce.obs.leiden == cell_cluster,'cell_type'] = cell_annotion
                
        for cluster_id, cell_type in zip(rest_cluster_id,rest_annotion):
            sce.obs.loc[sce.obs.leiden == cluster_id,'cell_type'] = cell_type
        # print ('cell annotion finished using：',run_time,' seconds')
        
        # plot and save files 
        ids = sc.pl.dotplot(sce,
                      gene_list,
                      groupby='cell_type',return_fig=True,show=False)
        ids.savefig(filename=work_dir + 'figures/' + files + '_CellType_GeneMarkers_DotPlot.pdf',show=False)
        ids = sc.pl.tsne(sce, color='cell_type',return_fig=True,show=False)
        ids.savefig(fname=work_dir + 'figures/' + files + '_CellType_Tsne.pdf',show=False)
       
        sce.obs.to_csv(work_dir + files+'_h5_file/'+files+'_cell_annotion.csv')
        
        if use_others:
            use_cells = sce.obs[sce.obs['cell_type'].isin(['Candidate_Cancer_cell','Z_Candidate'])].index
        else:
            use_cells = sce.obs[sce.obs['cell_type'].isin(['Candidate_Cancer_cell'])].index
        print('num candidate cancer cell: ',len(use_cells),'\n')    
        sce = sce.raw.to_adata()
        h5 = pd.HDFStore(work_dir + files+'_h5_file/'+files+'_Candidate_cancer_expr.h5', 'w',
                         complevel=4, complib='blosc')
        h5['expr'] = sce[use_cells].to_df()
        h5.close()
    else:
        h5 = pd.HDFStore(work_dir + files+'_h5_file/'+files+'_Candidate_cancer_expr.h5', 'w',
                         complevel=4, complib='blosc')
        h5['expr'] = sce.to_df()
        h5.close()
    print('\n','# ============================================================================= #','\n')

    del sce, cellinfo, geneinfo
else:
    print(work_dir+files+'_h5_file/'+files+'_Candidate_cancer_expr.h5 exit starting BscModel...\n')    
# =============================================================================
