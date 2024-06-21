import tensorflow as tf
from tensorflow import keras
from tensorflow.keras.models import load_model
from tensorflow.keras.backend import clear_session
import os
import re
import numpy as np
import pandas as pd
import sys
gpus=tf.config.experimental.list_physical_devices('GPU')
if gpus:
    try:
        for gpu in gpus:
            tf.config.experimental.set_memory_growth(gpu,True)
        logical_gpus=tf.config.experimental.list_logical_devices('GPU')
        print(len(gpus),"Physical GPUs",len(logical_gpus),"Logical GPUs")
    except RuntimeError as e:
        print(e)


'''
all_models_path:the path of all models (./DenseNet_models)
fasta_path: the path of one_sequence
wpath: the path of results (./results.csv)

'''
def get_all_models(all_models_path):
    models=[]
    for root ,dirs, files in os.walk(all_models_path):
        for f in files:
            model_path=os.path.join(root,f)
            models.append(model_path)
    return models
def One_hot(dirs):
    files=open(dirs,'r')
    sample=[]
    all_bed_number=[]
    for line in files:
        if re.match('>',line) is None:
            value=np.zeros((500,4),dtype='float32')
            for index,base in enumerate(line.strip()):
                if re.match(base,'A|a'):
                    value[index,0]=1
                if re.match(base,'C|c'):
                    value[index,1]=1
                if re.match(base,'G|g'):
                    value[index,2]=1
                if re.match(base,'T|t'):
                    value[index,3]=1
            sample.append(value)
        elif re.match('>',line): 
                    all_bed_number.append(line[1:].strip())
    files.close()
    return np.array(sample),all_bed_number

def predict_one_sequence(all_models_path,fasta_path,wpath):
    sample,all_bed_number=One_hot(fasta_path)
    result={}
    count=0
    for one_model_path in get_all_models(all_models_path):
        species_name=re.findall('.*/(.*)_models',one_model_path)[-1]
        model_tf=re.findall('.*/(.*)_',one_model_path)[-1]
        count+=1
        print('Now loading the trained model of factor {} of species {} for predicting! this is the num {} model, totally 389 models!'.format(model_tf, species_name,count))
        model=load_model(one_model_path)
        pred=model.predict(sample)
        clear_session()
        result_key=species_name+'_'+model_tf
        result[result_key]=list(pred[:,0])
    result_df=pd.DataFrame(result).T.reset_index()
    result_df.columns=['TF']+all_bed_number
    print(result_df)
    result_df.to_csv(wpath,sep=',')
    return result_df
def main():
    all_models_path='./DenseNet_models'
    fasta_path=sys.argv[1]
    wpath='./results.csv'
    result_df=predict_one_sequence(all_models_path,fasta_path,wpath)

if __name__ =='__main__':
    main()

    
