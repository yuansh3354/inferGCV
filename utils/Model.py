# Model utils 
# function 
import torch.nn as nn
import torch
import pytorch_lightning as pl
import torchmetrics
import torch.nn.functional as F
from utils.myexptorch import toLabel

class sampleToMatrix(nn.Module):
    def __init__(self, n_genes,out_feature):
        self.n_genes = n_genes
        self.out_feature = out_feature
        super(sampleToMatrix,self).__init__()
        self.reshapeLayer =nn.Sequential(
            nn.Linear(self.n_genes,  self.out_feature * self.out_feature),nn.ReLU(inplace=True)
            )
    
    def forward(self, sample):
        x = self.reshapeLayer(sample)
        x = x.view(-1, 1, self.out_feature, self.out_feature)
        return x                                                             

class L1regularization(nn.Module):
    def __init__(self, weight_decay=0.2):
        super(L1regularization, self).__init__()
        self.weight_decay = weight_decay

    def forward(self, model):
        regularization_loss = 0.
        for param in model.parameters():
            regularization_loss += torch.mean(abs(param)) * self.weight_decay

        return regularization_loss
    
class Encoder(nn.Module):
    def __init__(self, gene_num):
        super(Encoder, self).__init__()
        self.input_size = gene_num
        self.to_cov = sampleToMatrix(self.input_size,20)
        self.encoder = nn.Sequential(nn.Conv2d(1, 3, kernel_size=3, stride=1, bias=False),nn.BatchNorm2d(3),
                                     nn.ReLU(inplace=True),
                                     nn.Conv2d(3, 6, kernel_size=7, stride=1, bias=False),nn.BatchNorm2d(6),
                                     nn.ReLU(inplace=True),
                                     nn.Conv2d(6, 9, kernel_size=11))
        self.decoder = nn.Sequential(nn.Linear(36 , self.input_size))
        
    def forward(self, input_x):
        input_x = self.to_cov(input_x.squeeze())
        embedding = self.encoder(input_x)
        prob = embedding.view(input_x.shape[0],-1,2).mean(dim=1)
        embedding = embedding.view(input_x.shape[0],-1)
        decoder = self.decoder(embedding)
        return embedding, prob, decoder
        
class CenterLoss(nn.Module):
    def __init__(self, w=0.95):
        super(CenterLoss, self).__init__()
        self.w = w
            
    def forward(self, embedding, label, sce):
        label = label.long()
        batch = embedding.size(0)
        dim = embedding.size(1)
        
        mask1 = torch.transpose(label, 0, 1)[0].unsqueeze(1).expand(batch, dim) 
        mask2 = torch.transpose(label, 0, 1)[1].unsqueeze(1).expand(batch, dim)
        
        embedding1 = mask1 * embedding
        embedding2 = mask2 * embedding
        
        c1 = dim *  torch.sum(embedding1, dim=0) / (mask1.sum() + 1e-5)
        c2 = dim *  torch.sum(embedding2, dim=0) / (mask2.sum() + 1e-5)
        c1_cluster = F.pairwise_distance(embedding1, c1.expand(batch,dim), p=2)
        c2_cluster = F.pairwise_distance(embedding2, c2.expand(batch,dim), p=2)
        
        c1_cluster = sum(c1_cluster * (mask1.sum(axis=1)/dim))
        c2_cluster = sum(c2_cluster * (mask2.sum(axis=1)/dim))
        c1_cluster = dim * c1_cluster / (mask1.sum() + 1e-5)
        c2_cluster = dim * c2_cluster / (mask2.sum() + 1e-5)
        
        c1_dist = F.pairwise_distance(embedding2, c1.expand(batch,dim), p=2)
        c2_dist = F.pairwise_distance(embedding1, c2.expand(batch,dim), p=2)
        
        c1_dist = sum(c1_dist * (mask2.sum(axis=1)/dim))
        c2_dist = sum(c2_dist * (mask1.sum(axis=1)/dim))
        
        c1_dist = dim * c1_dist / (mask2.sum() + 1e-5)
        c2_dist = dim * c2_dist / (mask1.sum() + 1e-5)

        c_dist = F.pairwise_distance(c2,c1) # max
        
        dist_sce_c1 = F.pairwise_distance(sce, c1.expand(sce.shape[0],dim), p=2)
        dist_sce_c2 = self.w * F.pairwise_distance(sce, c2.expand(sce.shape[0],dim), p=2)

        sce_min = torch.minimum(dist_sce_c2,dist_sce_c1).max()
        sce_max = torch.maximum(dist_sce_c2,dist_sce_c1).min()

        center_loss = 1/(c_dist + sce_max) + c1_cluster  + c2_cluster + 1/(c1_dist + c2_dist) + sce_min
        return center_loss

class scMutation_pseudo_label(pl.LightningModule):
    def __init__(self, sce=None, gene_num=None, lr=1e-4, cl_weight=None, loss_weight = None, use_gpu=True) -> None:
        super(scMutation_pseudo_label, self).__init__()
        self.gene_num = gene_num
        self.lr = lr
        self.cl_weight = cl_weight
        self.sce = sce
        self.use_gpu = use_gpu
        self.encoder = Encoder(gene_num=self.gene_num)
        self.l1_regular = L1regularization()
        self.pl_accuracy = torchmetrics.Accuracy()
        self.loss_weight = loss_weight
        self.monitor = 0
        self.center_loss = CenterLoss(w=self.cl_weight)

        if self.use_gpu:
            self.sce = self.sce.to('cuda')
            self.loss_weight = self.loss_weight.to('cuda')
    def forward(self, input_x):
        embedding, probs, decoder = self.encoder(input_x)
        return embedding, probs, decoder

    def configure_optimizers(self):
        optimizer = torch.optim.Adam(self.parameters(), lr=self.lr)
        return optimizer

    def training_step(self, train_batch, batch_ix):
        x, y = train_batch
        bulk_embeddings, prob, decoder = self.forward(x)
        sce_embeddings, _, _= self.forward(self.sce)
        loss1 = F.cross_entropy(prob, toLabel(y), weight=self.loss_weight)
        reconstruction_loss = F.mse_loss(decoder, x)
        loss2 = self.l1_regular(self.encoder)
        loss3 = self.center_loss(bulk_embeddings, y,sce_embeddings) 
        
        loss = loss1 +loss2+ loss3 + 0.01 * reconstruction_loss

        acc = self.pl_accuracy(toLabel(prob), toLabel(y))
        self.log('acc', acc, on_step=True, prog_bar=True)
        self.log('loss', loss, on_step=True, prog_bar=True)
        return loss
    
    def validation_step(self, validation_batch, batch_ix):
        x, y = validation_batch
        bulk_embeddings, prob,_ = self.forward(x)
        val_acc = self.pl_accuracy(toLabel(prob), toLabel(y))
        self.log("val_acc", val_acc, on_step=True, prog_bar=True)
        mylogdict = {'log': {'val_acc': val_acc}}
        return mylogdict
    
    def validation_epoch_end(self, output):
        val_acc = sum([out['log']['val_acc'].item() for out in output]) / len(output)
        self.monitor = round(val_acc * 100,2)





















    
    
    
   
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
