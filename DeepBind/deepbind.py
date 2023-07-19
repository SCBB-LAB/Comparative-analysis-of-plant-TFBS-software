import torch 
import torchvision
import numpy
import numpy as np
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader
#from sklearn.metrics import confusion_matrix
from torchmetrics.classification import ConfusionMatrix
from torchmetrics.classification import BinaryAccuracy
from torchmetrics.classification import BinaryAUROC
from torchmetrics.classification import BinaryRecall
from torchmetrics.classification import BinarySpecificity
from torchmetrics.classification import BinaryMatthewsCorrCoef
from torchmetrics.classification import BinaryF1Score
import sys,os
import csv
import math 
import random
from statistics import mean
import gzip
from scipy.stats import bernoulli
import torch
from sklearn import metrics
import psutil
import humanize
import os
import GPUtil as GPU
GPUs = GPU.getGPUs()

gpu = GPUs[0]


def printm():
       process = psutil.Process(os.getpid())
       print("Gen RAM Free: " + humanize.naturalsize( psutil.virtual_memory().available ), " | Proc size: " + humanize.naturalsize( process.memory_info().rss))
       print("GPU RAM Free: {0:.0f}MB | Used: {1:.0f}MB | Util {2:3.0f}% | Total {3:.0f}MB".format(gpu.memoryFree, gpu.memoryUsed, gpu.memoryUtil*100, gpu.memoryTotal))
printm()

# Device configuration
device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
nummotif=16 #number of motifs to discover
bases='ACGT' #DNA bases
batch_size=64 #fixed batch size 
dictReverse={'A':'T','C':'G','G':'C','T':'A','N':'N'} #dictionary to implement reverse-complement mode
reverse_mode=False
def seqtopad(sequence,motlen,kind='DNA'):
    rows=len(sequence)+2*motlen-2
    S=np.empty([rows,4])
    base= bases if kind=='DNA' else basesRNA
    for i in range(rows):
        for j in range(4):
            if i-motlen+1<len(sequence) and sequence[i-motlen+1]=='N' or i<motlen-1 or i>len(sequence)+motlen-2:
                S[i,j]=np.float32(0.25)
            elif sequence[i-motlen+1]==base[j]:
                S[i,j]=np.float32(1)
            else:
                S[i,j]=np.float32(0)
    return np.transpose(S)
def dinucshuffle(sequence):
    b=[sequence[i:i+2] for i in range(0, len(sequence), 2)]
    random.shuffle(b)
    d=''.join([str(x) for x in b])
    return d
def complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    complseq = [complement[base] for base in seq]
    return complseq
def reverse_complement(seq):
    seq = list(seq)
    seq.reverse()
    return ''.join(complement(seq))
class Chip():
    def __init__(self,filename,motiflen=24,reverse_complemet_mode=reverse_mode):
        self.file = filename
        self.motiflen = motiflen
        self.reverse_complemet_mode=reverse_complemet_mode            
    def openFile(self):
        train_dataset=[]
        with open(self.file, 'rt') as data:
            next(data)
            reader = csv.reader(data,delimiter='\t')
            if not self.reverse_complemet_mode:
              for row in reader:
                      train_dataset.append([seqtopad(row[2],self.motiflen),[1]])
                      train_dataset.append([seqtopad(dinucshuffle(row[2]),self.motiflen),[0]])
            else:
              for row in reader:
                      train_dataset.append([seqtopad(row[2],self.motiflen),[1]])
                      train_dataset.append([seqtopad(reverse_complement(row[2]),self.motiflen),[1]])
                      train_dataset.append([seqtopad(dinucshuffle(row[2]),self.motiflen),[0]])
                      train_dataset.append([seqtopad(dinucshuffle(reverse_complement(row[2])),self.motiflen),[0]])                           
        train_dataset_pad=train_dataset
        size=int(len(train_dataset_pad)/3)
        firstvalid=train_dataset_pad[:size]
        secondvalid=train_dataset_pad[size:size+size]
        thirdvalid=train_dataset_pad[size+size:]
        firsttrain=secondvalid+thirdvalid
        secondtrain=firstvalid+thirdvalid
        thirdtrain=firstvalid+secondvalid
        return firsttrain,firstvalid,secondtrain,secondvalid,thirdtrain,thirdvalid,train_dataset_pad
chipseq=Chip("example/"+sys.argv[1]+'_train.txt')
train1,valid1,train2,valid2,train3,valid3,alldataset=chipseq.openFile()
class chipseq_dataset(Dataset):
    def __init__(self,xy=None):
        self.x_data=np.asarray([el[0] for el in xy],dtype=np.float32)
        self.y_data =np.asarray([el[1] for el in xy ],dtype=np.float32)
        self.x_data = torch.from_numpy(self.x_data)
        self.y_data = torch.from_numpy(self.y_data)
        self.len=len(self.x_data)
    def __getitem__(self, index):
        return self.x_data[index], self.y_data[index]
    def __len__(self):
        return self.len
train1_dataset=chipseq_dataset(train1)
train2_dataset=chipseq_dataset(train2)
train3_dataset=chipseq_dataset(train3)
valid1_dataset=chipseq_dataset(valid1)
valid2_dataset=chipseq_dataset(valid2)
valid3_dataset=chipseq_dataset(valid3)
alldataset_dataset=chipseq_dataset(alldataset)
batchSize=64
if reverse_mode:
  train_loader1 = DataLoader(dataset=train1_dataset,batch_size=batchSize,shuffle=False)
  train_loader2 = DataLoader(dataset=train2_dataset,batch_size=batchSize,shuffle=False)
  train_loader3 = DataLoader(dataset=train3_dataset,batch_size=batchSize,shuffle=False)
  valid1_loader = DataLoader(dataset=valid1_dataset,batch_size=batchSize,shuffle=False)
  valid2_loader = DataLoader(dataset=valid2_dataset,batch_size=batchSize,shuffle=False)
  valid3_loader = DataLoader(dataset=valid3_dataset,batch_size=batchSize,shuffle=False)
  alldataset_loader=DataLoader(dataset=alldataset_dataset,batch_size=batchSize,shuffle=False)
else:
  train_loader1 = DataLoader(dataset=train1_dataset,batch_size=batchSize,shuffle=True)
  train_loader2 = DataLoader(dataset=train2_dataset,batch_size=batchSize,shuffle=True)
  train_loader3 = DataLoader(dataset=train3_dataset,batch_size=batchSize,shuffle=True)
  valid1_loader = DataLoader(dataset=valid1_dataset,batch_size=batchSize,shuffle=False)
  valid2_loader = DataLoader(dataset=valid2_dataset,batch_size=batchSize,shuffle=False)
  valid3_loader = DataLoader(dataset=valid3_dataset,batch_size=batchSize,shuffle=False)
  alldataset_loader=DataLoader(dataset=alldataset_dataset,batch_size=batchSize,shuffle=False)
train_dataloader=[train_loader1,train_loader2,train_loader3]
valid_dataloader=[valid1_loader,valid2_loader,valid3_loader]
# Set Hyper-parameters 
num_epochs = 1
num_classes = 10
batch_size = 100
learning_rate = 0.001
def logsampler(a,b):
        x=np.random.uniform(low=0,high=1)
        y=10**((math.log10(b)-math.log10(a))*x + math.log10(a))
        return y    
def sqrtsampler(a,b):        
        x=np.random.uniform(low=0,high=1)
        y=(b-a)*math.sqrt(x)+a
        return y      
# input of shape(batch_size,inp_chan,iW)

class ConvNet(nn.Module):
    
    def __init__(self, nummotif,motiflen,poolType,neuType,mode,dropprob,learning_rate,momentum_rate,sigmaConv,sigmaNeu,beta1,beta2,beta3,reverse_complemet_mode=reverse_mode):      
        super(ConvNet, self).__init__()
        self.poolType=poolType
        self.neuType=neuType
        self.mode=mode
        self.reverse_complemet_mode=reverse_complemet_mode
        self.dropprob=dropprob
        self.learning_rate=learning_rate
        self.momentum_rate=momentum_rate
        self.sigmaConv=sigmaConv
        self.sigmaNeu=sigmaNeu
        self.beta1=beta1
        self.beta2=beta2
        self.beta3=beta3
        self.wConv=torch.randn(nummotif,4,motiflen).to(device)
        torch.nn.init.normal_(self.wConv,mean=0,std=self.sigmaConv)
        self.wConv.requires_grad=True
        self.wRect=torch.randn(nummotif).to(device)
        torch.nn.init.normal_(self.wRect)
        self.wRect=-self.wRect
        self.wRect.requires_grad=True  
        if neuType=='nohidden':            
            if poolType=='maxavg':
                self.wNeu=torch.randn(2*nummotif,1).to(device)
            else:
                self.wNeu=torch.randn(nummotif,1).to(device)
            self.wNeuBias=torch.randn(1).to(device)
            torch.nn.init.normal_(self.wNeu,mean=0,std=self.sigmaNeu)
            torch.nn.init.normal_(self.wNeuBias,mean=0,std=self.sigmaNeu)
        else:
            if poolType=='maxavg':
                self.wHidden=torch.randn(2*nummotif,32).to(device)
            else:                
                self.wHidden=torch.randn(nummotif,32).to(device)
            self.wNeu=torch.randn(32,1).to(device)
            self.wNeuBias=torch.randn(1).to(device)
            self.wHiddenBias=torch.randn(32).to(device)
            torch.nn.init.normal_(self.wNeu,mean=0,std=self.sigmaNeu)
            torch.nn.init.normal_(self.wNeuBias,mean=0,std=self.sigmaNeu)
            torch.nn.init.normal_(self.wHidden,mean=0,std=0.3)
            torch.nn.init.normal_(self.wHiddenBias,mean=0,std=0.3)              
            self.wHidden.requires_grad=True
            self.wHiddenBias.requires_grad=True
        self.wNeu.requires_grad=True
        self.wNeuBias.requires_grad=True   
    def divide_two_tensors(self,x):
        l=torch.unbind(x)
        list1=[l[2*i] for i in range(int(x.shape[0]/2))]
        list2=[l[2*i+1] for i in range(int(x.shape[0]/2))]
        x1=torch.stack(list1,0)
        x2=torch.stack(list2,0)
        return x1,x2
    def forward_pass(self,x,mask=None,use_mask=False):        
        conv=F.conv1d(x, self.wConv, bias=self.wRect, stride=1, padding=0)
        rect=conv.clamp(min=0)
        maxPool, _ = torch.max(rect, dim=2)
        if self.poolType=='maxavg':
            avgPool= torch.mean(rect, dim=2)                          
            pool=torch.cat((maxPool, avgPool), 1)
        else:
            pool=maxPool
        if(self.neuType=='nohidden'):
            if self.mode=='training': 
                if  not use_mask:
                  mask=bernoulli.rvs(self.dropprob, size=len(pool[0]))
                  mask=torch.from_numpy(mask).float().to(device)
                pooldrop=pool*mask
                out=pooldrop @ self.wNeu
                out.add_(self.wNeuBias)
            else:
                out=self.dropprob*(pool @ self.wNeu)
                out.add_(self.wNeuBias)       
        else:
            hid=pool @ self.wHidden
            hid.add_(self.wHiddenBias)
            hid=hid.clamp(min=0)
            if self.mode=='training': 
                if  not use_mask:
                  mask=bernoulli.rvs(self.dropprob, size=len(hid[0]))
                  mask=torch.from_numpy(mask).float().to(device)
                hiddrop=hid*mask
                out=self.dropprob*(hid @ self.wNeu)
                out.add_(self.wNeuBias)
            else:
                out=self.dropprob*(hid @ self.wNeu)
                out.add_(self.wNeuBias) 
        return out,mask       
    def forward(self, x):       
        if not  self.reverse_complemet_mode:
            out,_=self.forward_pass(x)            
        else:            
            x1,x2=self.divide_two_tensors(x)
            out1,mask=self.forward_pass(x1)
            out2,_=self.forward_pass(x2,mask,True)
            out=torch.max(out1, out2)     
        return out 

         
best_SN=0
best_SP=0
best_ACC=0
best_MCC=0
best_F1=0
best_AUC=0
device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
print(device)
learning_steps_list=[4000,8000,12000,16000,20000]
for number in range(5):    
    pool_List=['max','maxavg']        
    random_pool=random.choice(pool_List)    
    neuType_list=['hidden','nohidden']
    random_neuType=random.choice(neuType_list)
    dropoutList=[0.5,0.75,1.0]         
    dropprob=random.choice(dropoutList)    
    learning_rate=logsampler(0.0005,0.05)
    momentum_rate=sqrtsampler(0.95,0.99)  
    sigmaConv=logsampler(10**-7,10**-3)   
    sigmaNeu=logsampler(10**-5,10**-2) 
    beta1=logsampler(10**-15,10**-3)
    beta2=logsampler(10**-10,10**-3)
    beta3=logsampler(10**-10,10**-3)
  
    model_auc=[[],[],[]]
    model_acc=[[],[],[]]
    model_mcc=[[],[],[]]
    model_f1=[[],[],[]]
    model_sp=[[],[],[]]
    model_sn=[[],[],[]]
    for kk in range(3):
        model = ConvNet(16,24,random_pool,random_neuType,'training',dropprob,learning_rate,momentum_rate,sigmaConv,sigmaNeu,beta1,beta2,beta3,reverse_complemet_mode=reverse_mode).to(device)
        if random_neuType=='nohidden':
            optimizer = torch.optim.SGD([model.wConv,model.wRect,model.wNeu,model.wNeuBias], lr=model.learning_rate,momentum=model.momentum_rate,nesterov=True)
        else:
            optimizer = torch.optim.SGD([model.wConv,model.wRect,model.wNeu,model.wNeuBias,model.wHidden,model.wHiddenBias], lr=model.learning_rate,momentum=model.momentum_rate,nesterov=True)
        train_loader=train_dataloader[kk]
        valid_loader=valid_dataloader[kk]
        learning_steps=0
        while learning_steps<=20000:
            model.mode='training'
            auc=[]
            for i, (data, target) in enumerate(train_loader):
                data = data.to(device)
                target = target.to(device)
                if model.reverse_complemet_mode:
                  target_2=torch.randn(int(target.shape[0]/2),1)
                  for i in range(target_2.shape[0]):
                    target_2[i]=target[2*i]
                  target=target_2.to(device)
                # Forward pass
                output = model(data)
                if model.neuType=='nohidden':
                    loss = F.binary_cross_entropy(torch.sigmoid(output),target)+model.beta1*model.wConv.norm()+model.beta3*model.wNeu.norm()
                else:
                    loss = F.binary_cross_entropy(torch.sigmoid(output),target)+model.beta1*model.wConv.norm()+model.beta2*model.wHidden.norm()+model.beta3*model.wNeu.norm()
                optimizer.zero_grad()
                loss.backward()
                optimizer.step()
                learning_steps+=1
                if learning_steps% 4000==0:       
                    with torch.no_grad():
                        model.mode='test'
                        SP=[]
                        SN=[]
                        ACC=[]
                        MCC=[]
                        F1=[]
                        AUC=[]
                        
                        for i, (data, target) in enumerate(valid_loader):
                            data = data.to(device)
                            target = target.to(device)
                            if model.reverse_complemet_mode:
                              target_2=torch.randn(int(target.shape[0]/2),1)
                              for i in range(target_2.shape[0]):
                                target_2[i]=target[2*i]
                              target=target_2.to(device)
                            # Forward pass
                            output = model(data)
                            pred_sig=torch.sigmoid(output)
                            pred=pred_sig.cpu().detach().numpy().reshape(output.shape[0])
                            labels=target.cpu().numpy().reshape(output.shape[0])
                            confmat = ConfusionMatrix(task="binary")
                            metric = BinaryAccuracy()
                            metric_AUC = BinaryAUROC()
                            y_pred=torch.tensor(pred)
                            y_label=torch.tensor(labels) 
                            sn = BinaryRecall()
                            #acc = BinaryAccuracy()        
                            #confmat = ConfusionMatrix(task="binary")
                            sp = BinarySpecificity()
                            mcc = BinaryMatthewsCorrCoef()
                            f1 = BinaryF1Score()
                                                       
                            ACC.append(metric(y_pred, y_label))                            
                            AUC.append(metrics.roc_auc_score(labels, pred))
                            MCC.append(mcc(y_pred, y_label))
                            SN.append(sn(y_pred, y_label))
                            SP.append(sp(y_pred, y_label))
                            F1.append(f1(y_pred, y_label))
                            
                         
                        a=[]
                        for i in ACC:
                            tensor_to_list=i.tolist()
                            a.append(tensor_to_list)
                        accuracy=np.mean(ACC)
                        model_acc[kk].append(np.mean(ACC)) 
                        model_sn[kk].append(np.mean(SN))   
                        model_sp[kk].append(np.mean(SP))  
                        #model_auc[kk].append(np.mean(auc))
                        model_mcc[kk].append(np.mean(MCC))
                        model_f1[kk].append(np.mean(F1))
                        model_auc[kk].append(np.mean(AUC))
                        print('ACC performance when training fold number ',kk+1, 'using ',learning_steps_list[len(model_acc[kk])-1],'learning steps = ',np.mean(ACC))
                        print('SN performance when training fold number ',kk+1, 'using ',learning_steps_list[len(model_sn[kk])-1],'learning steps = ',np.mean(SN))
                        print('SP performance when training fold number ',kk+1, 'using ',learning_steps_list[len(model_sp[kk])-1],'learning steps = ',np.mean(SP))
                        print('MCC performance when training fold number ',kk+1, 'using ',learning_steps_list[len(model_mcc[kk])-1],'learning steps = ',np.mean(MCC))
                        print('F1 performance when training fold number ',kk+1, 'using ',learning_steps_list[len(model_f1[kk])-1],'learning steps = ',np.mean(F1))
                        print('AUC performance when training fold number ',kk+1, 'using ',learning_steps_list[len(model_auc[kk])-1],'learning steps = ',np.mean(AUC))
                       

    for n in range(5):
        SN=(model_sn[0][n]+model_sn[1][n]+model_sn[2][n])/3
        SP=(model_sp[0][n]+model_sp[1][n]+model_sp[2][n])/3
        ACC=(model_acc[0][n]+model_acc[1][n]+model_acc[2][n])/3
        MCC=(model_mcc[0][n]+model_mcc[1][n]+model_mcc[2][n])/3
        F1=(model_f1[0][n]+model_f1[1][n]+model_f1[2][n])/3
        AUC=(model_auc[0][n]+model_auc[1][n]+model_auc[2][n])/3
        #print(AUC)
        if AUC>best_AUC:
            best_SN=SN
            best_SP=SP
            best_ACC=ACC
            best_MCC=MCC
            best_F1=F1
            best_AUC=AUC
            best_learning_steps=learning_steps_list[n]
            best_LearningRate=model.learning_rate
            best_LearningMomentum=model.momentum_rate
            best_neuType=model.neuType
            best_poolType=model.poolType
            best_sigmaConv=model.sigmaConv
            best_dropprob=model.dropprob
            best_sigmaNeu=model.sigmaNeu
            best_beta1=model.beta1
            best_beta2=model.beta2
            best_beta3=model.beta3
print('best_poolType=',best_poolType)
print('best_neuType=',best_neuType)
print('best_SN=',best_SN) 
print('best_SP=',best_SP)
print('best_ACC=',best_ACC)
print('best_MCC=',best_MCC)
print('best_F1=',best_F1)
print('best_AUC=',best_AUC)           
print('best_learning_steps=',best_learning_steps)      
print('best_LearningRate=',best_LearningRate)
print('best_LearningMomentum=',best_LearningMomentum)
print('best_sigmaConv=',best_sigmaConv)
print('best_dropprob=',best_dropprob)
print('best_sigmaNeu=',best_sigmaNeu)
print('best_beta1=',best_beta1)
print('best_beta2=',best_beta2)
print('best_beta3=',best_beta3)
f3 = open("output/"+sys.argv[1]+"_result.txt",'w')
f3.writelines("Train"+"\t"+"TPR"+"\t"+"TNR"+"\t"+"ACC"+"\t"+"F1"+"\t"+"MCC"+"\t"+"auc"+"\n")
f3.writelines("\t"+str(best_SN)+"\t"+str(best_SP)+"\t"+str(best_ACC)+"\t"+str(best_F1)+"\t"+str(best_MCC)+"\t"+str(best_AUC)+"\n")
f3.close()
best_hyperparameters = {'best_poolType': best_poolType,'best_neuType':best_neuType,'best_learning_steps':best_learning_steps,'best_LearningRate':best_LearningRate,
                        'best_LearningMomentum':best_LearningMomentum,'best_sigmaConv':best_sigmaConv,'best_dropprob':best_dropprob,
                        'best_sigmaNeu':best_sigmaNeu,'best_beta1':best_beta1, 'best_beta2':best_beta2,'best_beta3':best_beta3}
torch.save(best_hyperparameters,"output/"+sys.argv[1]+'_best_hyperpamarameters.pth')           
# input of shape(batch_size,inp_chan,iW)
class ConvNet_test(nn.Module):
   # import numba
   # @numba.jit
    def __init__(self,nummotif,motiflen,poolType,neuType,mode,learning_steps,learning_rate,learning_Momentum,sigmaConv,dropprob,sigmaNeu,beta1,beta2,beta3,reverse_complemet_mode):
        super(ConvNet_test, self).__init__()
        self.poolType=poolType
        self.neuType=neuType
        self.mode=mode
        self.learning_rate=learning_rate
        self.reverse_complemet_mode=reverse_complemet_mode
        self.momentum_rate=learning_Momentum
        self.sigmaConv=sigmaConv
        self.wConv=torch.randn(nummotif,4,motiflen).to(device)
        torch.nn.init.normal_(self.wConv,mean=0,std=self.sigmaConv)
        self.wConv.requires_grad=True
        self.wRect=torch.randn(nummotif).to(device)
        torch.nn.init.normal_(self.wRect)
        self.wRect=-self.wRect
        self.wRect.requires_grad=True
        self.dropprob=dropprob
        self.sigmaNeu=sigmaNeu
        self.wHidden=torch.randn(2*nummotif,32).to(device)
        self.wHiddenBias=torch.randn(32).to(device)
        if neuType=='nohidden':            
            if poolType=='maxavg':
                self.wNeu=torch.randn(2*nummotif,1).to(device)
            else:
                self.wNeu=torch.randn(nummotif,1).to(device)
            self.wNeuBias=torch.randn(1).to(device)
            torch.nn.init.normal_(self.wNeu,mean=0,std=self.sigmaNeu)
            torch.nn.init.normal_(self.wNeuBias,mean=0,std=self.sigmaNeu)
        else:
            if poolType=='maxavg':
                self.wHidden=torch.randn(2*nummotif,32).to(device)
            else:                
                self.wHidden=torch.randn(nummotif,32).to(device)
            self.wNeu=torch.randn(32,1).to(device)
            self.wNeuBias=torch.randn(1).to(device)
            self.wHiddenBias=torch.randn(32).to(device)
            torch.nn.init.normal_(self.wNeu,mean=0,std=self.sigmaNeu)
            torch.nn.init.normal_(self.wNeuBias,mean=0,std=self.sigmaNeu)
            torch.nn.init.normal_(self.wHidden,mean=0,std=0.3)
            torch.nn.init.normal_(self.wHiddenBias,mean=0,std=0.3)             
            self.wHidden.requires_grad=True
            self.wHiddenBias.requires_grad=True
        self.wNeu.requires_grad=True
        self.wNeuBias.requires_grad=True        
        self.beta1=beta1
        self.beta2=beta2
        self.beta3=beta3
    def divide_two_tensors(self,x):
        l=torch.unbind(x)
        list1=[l[2*i] for i in range(int(x.shape[0]/2))]
        list2=[l[2*i+1] for i in range(int(x.shape[0]/2))]
        x1=torch.stack(list1,0)
        x2=torch.stack(list2,0)
        return x1,x2
    def forward_pass(self,x,mask=None,use_mask=False):        
        conv=F.conv1d(x, self.wConv, bias=self.wRect, stride=1, padding=0)
        rect=conv.clamp(min=0)
        maxPool, _ = torch.max(rect, dim=2)
        if self.poolType=='maxavg':
            avgPool= torch.mean(rect, dim=2)                          
            pool=torch.cat((maxPool, avgPool), 1)
        else:
            pool=maxPool
        if(self.neuType=='nohidden'):
            if self.mode=='training': 
                if  not use_mask:
                    mask=bernoulli.rvs(self.dropprob, size=len(pool[0]))
                    mask=torch.from_numpy(mask).float().to(device)
                pooldrop=pool*mask
                out=pooldrop @ self.wNeu
                out.add_(self.wNeuBias)
            else:
                out=self.dropprob*(pool @ self.wNeu)
                out.add_(self.wNeuBias)       
        else:
            hid=pool @ self.wHidden
            hid.add_(self.wHiddenBias)
            hid=hid.clamp(min=0)
            if self.mode=='training': 
                if  not use_mask:
                    mask=bernoulli.rvs(self.dropprob, size=len(hid[0]))
                    mask=torch.from_numpy(mask).float().to(device)
                hiddrop=hid*mask
                out=self.dropprob*(hid @ self.wNeu)
                out.add_(self.wNeuBias)
            else:
                out=self.dropprob*(hid @ self.wNeu)
                out.add_(self.wNeuBias) 
        return out,mask       
    def forward(self, x):        
        if not  self.reverse_complemet_mode:
            out,_=self.forward_pass(x)
        else:            
            x1,x2=self.divide_two_tensors(x)
            out1,mask=self.forward_pass(x1)
            out2,_=self.forward_pass(x2,mask,True)
            out=torch.max(out1, out2)        
        return out    
device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
print(device)
learning_steps_list=[4000,8000,12000,16000,20000]
best_SN=0
best_SP=0
best_ACC=0
best_MCC=0
best_F1=0
best_AUC=0
best_CM=0

best_hyperparameters = torch.load('model/'+sys.argv[1]+'_best_hyperpamarameters.pth')
best_poolType=best_hyperparameters['best_poolType']
best_neuType=best_hyperparameters['best_neuType']
best_learning_steps=best_hyperparameters['best_learning_steps']
best_LearningRate=best_hyperparameters['best_LearningRate']
best_dropprob=best_hyperparameters['best_dropprob']
best_LearningMomentum=best_hyperparameters['best_LearningMomentum']
best_sigmaConv=best_hyperparameters['best_sigmaConv']
best_sigmaNeu=best_hyperparameters['best_sigmaNeu']
best_beta1=best_hyperparameters['best_beta1']
best_beta2=best_hyperparameters['best_beta2']
best_beta3=best_hyperparameters['best_beta3']
for number_models in range(6):
  model = ConvNet_test(16,24,best_poolType,best_neuType,'training',best_learning_steps,best_LearningRate,best_LearningMomentum,best_sigmaConv,best_dropprob,best_sigmaNeu,best_beta1,best_beta2,best_beta3,reverse_complemet_mode=False).to(device)
  if model.neuType=='nohidden':
      optimizer = torch.optim.SGD([model.wConv,model.wRect,model.wNeu,model.wNeuBias], lr=model.learning_rate,momentum=model.momentum_rate,nesterov=True)
  else:
      optimizer = torch.optim.SGD([model.wConv,model.wRect,model.wNeu,model.wNeuBias,model.wHidden,model.wHiddenBias], lr=model.learning_rate,momentum=model.momentum_rate,nesterov=True)
  train_loader=alldataset_loader
  valid_loader=alldataset_loader
  learning_steps=0
  while learning_steps<=best_learning_steps:  
      for i, (data, target) in enumerate(train_loader):
          data = data.to(device)
          target = target.to(device)
          if model.reverse_complemet_mode:
              target_2=torch.randn(int(target.shape[0]/2),1)
              for i in range(target_2.shape[0]):
                target_2[i]=target[2*i]
              target=target_2.to(device)
            # Forward pass
          output = model(data)          
          if model.neuType=='nohidden':
              loss = F.binary_cross_entropy(torch.sigmoid(output),target)+model.beta1*model.wConv.norm()+model.beta3*model.wNeu.norm()
          else:
              loss = F.binary_cross_entropy(torch.sigmoid(output),target)+model.beta1*model.wConv.norm()+model.beta2*model.wHidden.norm()+model.beta3*model.wNeu.norm()
          optimizer.zero_grad()
          loss.backward()
          optimizer.step()
          learning_steps+=1          
  with torch.no_grad():
      model.mode='test'
      SP=[]
      SN=[]
      ACC=[]
      MCC=[]
      F1=[]
      AUC=[]
      auc=[]
      #acc=[]
      for i, (data, target) in enumerate(valid_loader):
          data = data.to(device)
          target = target.to(device)
          if model.reverse_complemet_mode:
              target_2=torch.randn(int(target.shape[0]/2),1)
              for i in range(target_2.shape[0]):
                target_2[i]=target[2*i]
              target=target_2.to(device)
          # Forward pass
          output = model(data)
          pred_sig=torch.sigmoid(output)
          pred=pred_sig.cpu().detach().numpy().reshape(output.shape[0])
          labels=target.cpu().numpy().reshape(output.shape[0])
          auc.append(metrics.roc_auc_score(labels, pred))
          confmat = ConfusionMatrix(task="binary")
          metric = BinaryAccuracy()
          metric_AUC = BinaryAUROC()
          y_pred=torch.tensor(pred)
          y_label=torch.tensor(labels)                            
          ACC.append(metric(y_pred, y_label))                            
          AUC.append(metrics.roc_auc_score(labels, pred))
          MCC.append(mcc(y_pred, y_label))
          SN.append(sn(y_pred, y_label))
          SP.append(sp(y_pred, y_label))
          F1.append(f1(y_pred, y_label))
          sn = BinaryRecall()
          #acc = BinaryAccuracy()        
          #confmat = ConfusionMatrix(task="binary")
          sp = BinarySpecificity()
          mcc = BinaryMatthewsCorrCoef()
          f1 = BinaryF1Score()
           
      accu=(metric(y_pred, y_label))
      cm=(confmat(y_pred, y_label))
      mcc1=(mcc(y_pred, y_label))
      sn1=(sn(y_pred, y_label))
      sp1=(sp(y_pred, y_label))
      f11=(f1(y_pred, y_label))
      auc=(metric_AUC(y_pred, y_label))         
   
      model_auc[kk].append(np.mean(AUC)) 
      accuracy=np.mean(ACC)
      #model_acc[kk].append(np.mean(ACC)) 
      snn=(np.mean(SN))   
      spp=(np.mean(SP))  
      #model_auc[kk].append(np.mean(auc))
      MCCC=(np.mean(MCC))
      F11=(np.mean(F1))
      AUC=(np.mean(AUC))             
      AUC_training=np.mean(AUC)
      
      if AUC_training>best_AUC:
          state = {'conv': model.wConv,'rect':model.wRect,'wHidden':model.wHidden,'wHiddenBias':model.wHiddenBias,'wNeu':model.wNeu,'wNeuBias':model.wNeuBias }
          torch.save(state, 'output/'+sys.argv[1]+'_Model.pth')
checkpoint = torch.load('output/'+sys.argv[1]+'_Model.pth')
model = ConvNet_test(16,24,best_poolType,best_neuType,'test',best_learning_steps,best_LearningRate,best_LearningMomentum,best_sigmaConv,best_dropprob,best_sigmaNeu,best_beta1,best_beta2,best_beta3,reverse_complemet_mode=reverse_mode).to(device)
model.wConv=checkpoint['conv']
model.wRect=checkpoint['rect']
model.wHidden=checkpoint['wHidden']
model.wHiddenBias=checkpoint['wHiddenBias']
model.wNeu=checkpoint['wNeu']
model.wNeuBias=checkpoint['wNeuBias']
# Test the model 
with torch.no_grad():
      model.mode='test'
      auc=[]     
      for i, (data, target) in enumerate(valid_loader):
          data = data.to(device)
          target = target.to(device)
          if model.reverse_complemet_mode:
              target_2=torch.randn(int(target.shape[0]/2),1)
              for i in range(target_2.shape[0]):
                target_2[i]=target[2*i]
              target=target_2.to(device)
          # Forward pass
          output = model(data)
          pred_sig=torch.sigmoid(output)
          pred=pred_sig.cpu().detach().numpy().reshape(output.shape[0])
          labels=target.cpu().numpy().reshape(output.shape[0])
          #auc.append(metrics.roc_auc_score(labels, pred))                         
      AUC_training=np.mean(auc)
      print(AUC_training)
class Chip_test():
    def __init__(self,filename,motiflen=24,reverse_complemet_mode=reverse_mode):
        self.file = filename
        self.motiflen = motiflen
        self.reverse_complemet_mode=reverse_complemet_mode            
    def openFile(self):
        test_dataset=[]
        with open(self.file, 'rt') as data:
            next(data)
            reader = csv.reader(data,delimiter='\t')
            if not self.reverse_complemet_mode:
              for row in reader:
                      test_dataset.append([seqtopad(row[2],self.motiflen),[int(row[3])]])
                      test_dataset.append([seqtopad(dinucshuffle(row[2]),self.motiflen),[0]])
                      
            else:
              for row in reader:
                      test_dataset.append([seqtopad(row[2],self.motiflen),[int(row[3])]])
                      test_dataset.append([seqtopad(reverse_complement(row[2]),self.motiflen),[int(row[3])]])    
                      test_dataset.append([seqtopad(dinucshuffle(row[2]),self.motiflen),[0]])
                      test_dataset.append([seqtopad(dinucshuffle(reverse_complement(row[2])),self.motiflen),[0]])        
        return test_dataset
chipseq_test=Chip_test("example/"+sys.argv[1]+'_test.txt')
test_data=chipseq_test.openFile()
test_dataset=chipseq_dataset(test_data)
batchSize=test_dataset.__len__()
test_loader = DataLoader(dataset=test_dataset,batch_size=batchSize,shuffle=False)
with torch.no_grad():
      model.mode='test'
      #auc=[]
      SP=[]
      SN=[]
      ACC=[]
      MCC=[]
      F1=[]
      AUC=[]     
      for i, (data, target) in enumerate(test_loader):
          data = data.to(device)
          target = target.to(device)
          if model.reverse_complemet_mode:
              target_2=torch.randn(int(target.shape[0]/2),1)
              for i in range(target_2.shape[0]):
                target_2[i]=target[2*i]
              target=target_2.to(device)
          # Forward pass
          output = model(data)
          pred_sig=torch.sigmoid(output)
          pred=pred_sig.cpu().detach().numpy().reshape(output.shape[0])
          labels=target.cpu().numpy().reshape(output.shape[0])          
          confmat = ConfusionMatrix(task="binary")
          metric = BinaryAccuracy()
          metric_AUC = BinaryAUROC()
          y_pred=torch.tensor(pred)
          y_label=torch.tensor(labels)                            
          ACC.append(metric(y_pred, y_label))                            
          AUC.append(metrics.roc_auc_score(labels, pred))
          confusionmatrix= confmat(y_pred, y_label)
          MCC.append(mcc(y_pred, y_label))
          SN.append(sn(y_pred, y_label))
          SP.append(sp(y_pred, y_label))
          F1.append(f1(y_pred, y_label))
      cm=(confmat(y_pred, y_label)) 
      TN=(cm[0][0].tolist())
      FP=(cm[0][1].tolist())
      FN=(cm[1][0].tolist())
      TP=(cm[1][1].tolist())   
      accuracy_test=np.mean(ACC)
      sensitivity_test=(np.mean(SN))   
      specificity_test=(np.mean(SP))  
      mmc_score_test=(np.mean(MCC))
      f1_score_test=(np.mean(F1))                
      AUC_test=np.mean(AUC)
      
      f3 = open("output/"+sys.argv[1]+"_result.txt",'a+')       
      f3.writelines("Test"+"\t"+"TP"+"\t"+"TN"+"\t"+"FN"+"\t"+"FP"+"\t"+"TPR"+"\t"+"TNR"+"\t"+"ACC"+"\t"+"F1"+"\t"+"MCC"+"\t"+"auc"+"\n")
      f3.writelines("\t"+str(TP)+"\t"+str(TN)+"\t"+str(FN)+"\t"+str(FP)+"\t"+str(sensitivity_test)+"\t"+str(specificity_test)+"\t"+str(accuracy_test)+"\t"+str(f1_score_test)+"\t"+str(mmc_score_test)+"\t"+str(AUC_test)+"\n")
      f3.close()
