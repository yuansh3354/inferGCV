# Carcinogenic Variation of Gene Confer Adaptive Evolutionary Superiority to Single Cancer Cell

It is now widely acknowledged that the adaptive evolutionary process of cancer cells is characterized by genomic and epigenetic alterations. Nevertheless, the extreme heterogeneity of cancer cells and the extreme sparsity of scRNA-seq data limited the understanding of cancer cell reconfiguration of transcriptome and adaptive evolution under selection pressure. In this study, we used a contrastive learning approach to analyze scRNA-seq data and developed a model to infer the carcinogenic variation (CGV) of the function of target gene determine the spatial mapping correlation between RNA features and genomic variation-associated transcriptome dysfunction in cancer. The dynamic transcriptional profile revealed by scRNA-seq reflects how cancer cells reconfigure their transcriptome pattern and develop various CGV patterns in response to selective pressure to gain an evolutionary advantage. We found that the plasticity of human cancer was negatively related to the transcriptome burden, and increasing transcriptome burden was associated with more chaotic CGV patterns. 
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
