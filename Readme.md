# sceAPA: evaluate anomaly pathway activity status of cancer cells in scRNA-seq <img src="logo.jpeg" width="280px" align="right" />
As a core of the tumor niche, our understanding of functional programs and the nature of cancer cells remains rudimentary. One of the main reasons is strong transcriptional heterogeneity in cancer cells is much higher than that in non-cancer cells, which is kept researchers from further exploring the molecular basis and subtype of specific cancer as well as the progression trajectory of cancer cells. Although, a common approach to discriminate cancer cells from normal cells is identifying copy number variation, we cannot get more information from that to cluster hidden subclusters of cancer cells and to further understand the biological diversity of tumors. Therefore, the need is especially acute for new data analysis approach to unlock the secrets of cancer cells at the single-cell level.
> 作为肿瘤生态位的核心，我们对功能程序和癌细胞性质的理解仍处于初级阶段。其中一个主要原因是癌细胞中强烈的转录异质性远高于非癌细胞，这使得研究人员无法进一步探索特定癌症的分子基础、亚型以及癌细胞的进展轨迹。虽然区分癌细胞和正常细胞的一种常见方法是识别拷贝数变异，但我们无法从中获得更多信息来聚类隐藏的癌细胞亚群，并进一步了解肿瘤的生物多样性。因此，迫切需要一种新的数据分析方法，在单细胞水平上解开癌细胞的秘密。

Here, we present sceAPA, a transfer learning approach to evaluate anomaly pathway activity (APA) status of cancer cells in scRNA-seq. sceAPA leverages knowledge from RNA-seq, using convolutional neural network to evaluate APA status information in RNA-seq, simultaneously, to embed scRNA-seq to screen each cell’s APA status. In public scRNA-seq datasets of various human cancer types, we demonstrate that sceAPA can catch the weak signals of APA and discover that APA burden index calculated by sceAPA is highly correlated with immune-derived tumor cells in melanoma, lung cancer and colorectal cancer. Collectively, sceAPA provides an innovative approach that focuses on a more comprehensive understanding of the cluster of cancer cells in tumor research.
> 在这里，我们介绍了sceAPA，一种转移学习方法，用于评估scRNA序列中癌细胞的异常通路活性（APA）状态。sceAPA利用RNA-seq的知识，使用卷积神经网络评估RNA-seq中的APA状态信息，同时嵌入scRNA-seq以筛选每个细胞的APA状态。在各种人类癌症类型的公共scRNA-seq数据集中，我们证明sceAPA可以捕捉APA的微弱信号，并发现sceAPA计算的APA负荷指数与黑色素瘤、肺癌和结直肠癌中的免疫源性肿瘤细胞高度相关。总的来说，sceAPA提供了一种创新方法，专注于更全面地了解肿瘤研究中的癌细胞群。

<div align='center' ><b><font size='150'>Overview of sceAPA</font></b></div>
  
<img src="IMG_00001.jpg" width="3000px"/>

## Pre-requisites:

- Linux (Based on Ubuntu 20.04 LTS, Personal Computer) 
- CPU AMD Ryzen 9 3950X
- NVIDIA GeForce RTX 3090 24GB 384bit 1695MHz 19500MHz 
- Memory 128G (32GB*4) DDR4 3200MHz

### Environment and resource allocation

---

For instructions on installing anaconda on your machine (download the distribution that comes with python 3):
https://www.anaconda.com/distribution/

```
#conda env create -n CaSee -f configs/CaSee_env_info.yaml

# if there are some warnings or errors 
# you can manually install some main packages
# all pip software in all_pip_pacakges.yaml
conda create -n CaSee python==3.8.8 # python==3.8
conda activate CaSee

pip install pytorch-lightning==1.3.7 # -i https://pypi.tuna.tsinghua.edu.cn/simple
pip install scipy==1.7.0 # -i https://pypi.tuna.tsinghua.edu.cn/simple
pip install numpy==1.20.3 # -i https://pypi.tuna.tsinghua.edu.cn/simple
pip install scanpy==1.7.2 # -i https://pypi.tuna.tsinghua.edu.cn/simple
pip install scikit-learn==0.23.2 # -i https://pypi.tuna.tsinghua.edu.cn/simple
pip3 install opencv-python==4.5.2.54 # -i https://pypi.tuna.tsinghua.edu.cn/simple
pip install torchmetrics==0.3.2 # -i https://pypi.tuna.tsinghua.edu.cn/simple
pip install torchvision==0.10.0 # -i https://pypi.tuna.tsinghua.edu.cn/simple

```

And you also download pytorch https://pytorch.org/ 

Attention, if you in the Chinese mainland, plz use `pip install` instand `conda install` 

**Ubuntu**

```
pip3 install torch==1.9.0+cu111 torchvision==0.10.0+cu111 torchaudio==0.9.0 -f https://download.pytorch.org/whl/torch_stable.html


```

**MacOS**

```
pip3 install torch==1.9.0 torchvision==0.10.0 torchaudio==0.9.0


```

**Windos**

```
pip3 install torch==1.9.0+cu111 torchvision==0.10.0+cu111 torchaudio==0.9.0 -f https://download.pytorch.org/whl/cu111/torch_stable.html
```
> torch==1.9.0+cu111  
> torchvision==0.10.0+cu111  
> torchaudio==0.9.0
> ## Prepare candidate ref data

Download `ref_data.tar.xz` and unzip the file, move the whole `ref_data` into the CaSee program work_dir.


## Running Model

```
python ScMutation.py \
    -g $genes \
    -s 'ref/GSE116237_Candidate_cancer_expr.h5' \
    -b 'ref/scMutation_ref.h5' \
    -o 'Pan_cancer_Pathyway_detection/GSE116237/Pan_cancer_'$id \
    --lr 1e-4 \
    --batch_size 128 \
    --max_epoch 100 \
    --n_neighbors_pseudo 30 \
    --cl_weight 0.895
```
