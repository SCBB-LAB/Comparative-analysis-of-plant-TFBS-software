 #!/usr/bin/env python3

# -------------------------------------------------------------------------
# BIG DISCLAIMER! ONLY TESTED FOR python3
# @author: Katherine Mejia-Guerra (mm2842@cornell.edu)
# Copyright (C) 2016 Katherine Mejia-Guerra
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# -------------------------------------------------------------------------

import os
import sys
#os.system("mkdir /home/scbb/j_dap_tf/k-mer_grammar/265_dap_out/bag-of-k-mers/"+sys.argv[3]+"/")
os.system("mkdir output/")
import numpy as np
import pandas as pd
import sqlalchemy
import logging
import time
from math import log
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from itertools import product
from sklearn import metrics
from sklearn.metrics import confusion_matrix, accuracy_score, matthews_corrcoef, f1_score , average_precision_score, roc_auc_score
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.linear_model import LogisticRegression
#from sklearn.externals import joblib
import joblib
import math

def comparison(testlabel, resultslabel):
    TP = 0
    FP = 0
    TN = 0
    FN = 0
    for row1 in range(len(resultslabel)):
        if resultslabel[row1] < 0.5:
            resultslabel[row1] = 0
        else:
            resultslabel[row1] = 1
    for row2 in range(len(testlabel)):
        if testlabel[row2] == 1 and testlabel[row2] == resultslabel[row2]:
            TP = TP + 1
        if testlabel[row2] == 0 and testlabel[row2] != resultslabel[row2]:
            FP = FP + 1
        if testlabel[row2] == 0 and testlabel[row2] == resultslabel[row2]:
            TN = TN + 1
        if testlabel[row2] == 1 and testlabel[row2] != resultslabel[row2]:
            FN = FN + 1
    if TP + FN != 0:
        TPR = TP / (TP + FN)
    else:
        TPR = 0
    if TN + FP != 0:
        TNR = TN / (TN + FP)
    else:
        TNR = 0
    if TP + FP != 0:
        PPV = TP / (TP + FP)
    else:
        PPV = 0
    if TN + FN != 0:
        NPV = TN / (TN + FN)
    else:
        NPV = 0
    if FN + TP != 0:
        FNR = FN / (FN + TP)
    else:
        FNR = 0
    if FP + TN != 0:
        FPR = FP / (FP + TN)
    else:
        FPR = 0
    if FP + TP != 0:
        FDR = FP / (FP + TP)
    else:
        FDR = 0
    if FN + TN != 0:
        FOR = FN / (FN + TN)
    else:
        FOR = 0
    if TP + TN + FP + FN != 0:
        ACC = (TP + TN) / (TP + TN + FP + FN)
    else:
        ACC = 0
    if TP + FP + FN != 0:
        F1 = (2 * TP) / (2 * TP + FP + FN)
    else:
        F1 = 0
    if (TP + FP) * (TP + FN) * (TN + FP) * (TN + FN) != 0:
        MCC = (TP * TN + FP * FN) / math.sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
    else:
        MCC = 0
    if TPR != 0 and TNR != 0:
        BM = TPR + TNR - 1
    else:
        BM = 0
    if PPV != 0 and NPV != 0:
        MK = PPV + NPV - 1
    else:
        MK = 0
    return TP, FP, TN, FN, TPR, TNR, PPV, NPV, FNR, FPR, FDR, FOR, ACC, F1, MCC, BM, MK

 
def createKmerSet(kmersize):
    '''
    write all possible kmers
    :param kmersize: integer, 8
    :return uniq_kmers: list of sorted unique kmers
    '''
    kmerSet = set()
    nucleotides = ["a", "c", "g", "t"]    
    kmerall = product(nucleotides, repeat=kmersize)
    for i in kmerall:
        kmer = ''.join(i)
        kmerSet.add(kmer)
    uniq_kmers = sorted(list(kmerSet))  
    return uniq_kmers


def compute_kmer_entropy(kmer):
    '''
    compute shannon entropy for each kmer
    :param kmer: string
    :return entropy: float
    '''
    prob = [float(kmer.count(c)) / len(kmer) for c in dict.fromkeys(list(kmer))]
    entropy = - sum([ p * log(p) / log(2.0) for p in prob ])
    return round(entropy, 2)


def make_stopwords(kmersize):
    '''
    write filtered out kmers
    :param kmersize: integer, 8
    :return stopwords: list of sorted low-complexity kmers
    '''
    kmersize_filter = {5:1.3, 6:1.3, 7:1.3, 8:1.3, 9:1.3, 10:1.3}
    limit_entropy = kmersize_filter.get(kmersize)
    kmerSet = set()
    nucleotides = ["a", "c", "g", "t"]    
    kmerall = product(nucleotides, repeat=kmersize)
    for n in kmerall:
        kmer = ''.join(n)
        
        if compute_kmer_entropy(kmer) < limit_entropy:
            kmerSet.add(make_newtoken(kmer))
        else:
            continue
    stopwords = sorted(list(kmerSet))
    return stopwords

  
def createNewtokenSet(kmersize):
    '''
    write all possible newtokens
    :param kmersize: integer, 8
    :return uniq_newtokens: list of sorted unique newtokens
    ''' 
    newtokenSet = set()
    uniq_kmers = createKmerSet(kmersize)
    for kmer in uniq_kmers:
        newtoken = make_newtoken(kmer)
        newtokenSet.add(newtoken)  
    uniq_newtokens = sorted(list(newtokenSet))
    return uniq_newtokens      


def make_newtoken(kmer):
    '''
    write a collapsed kmer and kmer reverse complementary as a newtoken
    :param kmer: string e.g., "AT"
    :return newtoken: string e.g., "atnta"
    :param kmer: string e.g., "TA"
    :return newtoken: string e.g., "atnta"
    '''
    kmer = str(kmer).lower()
    newtoken = "n".join(sorted([kmer,kmer.translate(str.maketrans('tagc', 'atcg'))[::-1]]))
    return newtoken

def write_ngrams(sequence):
    '''
    write a bag of newtokens of size n
    :param sequence: string e.g., "ATCG"
    :param (intern) kmerlength e.g., 2
    :return newtoken_string: string e.g., "atnta" "gatc" "cgcg" 
    '''
    seq = str(sequence).lower()
    finalstart = (len(seq)-kmerlength)+1
    allkmers = [seq[start:(start+kmerlength)] for start in range(0,finalstart)]
    tokens = [make_newtoken(kmer) for kmer in allkmers if len(kmer) == kmerlength and "n" not in kmer]
    newtoken_string = " ".join(tokens)
    
def save_plot_prc(precision,recall, avg_prec, figure_file, name):
    '''
    make plot for precission recall
    :param precission: precission
    :param recall: recall
    :param avg_prec: avg_prec
    :param figure_file: figure_file
    :param name: name
    :return plot precission recall curve
    '''    
    plt.clf()
    title = 'Precision Recall Curve - double strand '+ name
    plt.title(title)
    plt.plot(recall, precision, label='Precission = %0.2f' % avg_prec )
    plt.legend(loc='lower right')
    plt.plot([0,1],[0,1],'r--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.savefig("output/"+sys.argv[3]+"_prc.png")
            
def save_plot_roc(false_positive_rate, true_positive_rate, roc_auc, figure_file, name):
    '''
    make plot for roc_auc
    :param false_positive_rate: false_positive_rate
    :param true_positive_rate: true_positive_rate
    :param roc_auc: roc_auc
    :param figure_file: figure_file
    :param name: name
    :return roc_auc
    '''
    plt.clf()
    title = 'Receiver Operating Characteristic - double strand '+ name
    plt.title(title)
    plt.plot(false_positive_rate, true_positive_rate, 'b', label='AUC = %0.2f'% roc_auc)
    plt.legend(loc='lower right')
    plt.plot([0,1],[0,1],'r--')
    plt.xlim([-0.1,1.2])
    plt.ylim([-0.1,1.2])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.savefig("output/"+sys.argv[3]+"_roc.png")    


if sys.argv[1] == "-help":
    print("Usage: python kgrammar_bag-of-k-mer_training_testing.py [kmersize, integer] [mode filtered, 'True', or 'False' (i.e., mode full)] [dataset_name, string]")
    print("Example: python kgrammar_bag-of-k-mer_training_testing.py 8 False FEA4")
    quit()
else:
    kmersize = sys.argv[1] #e.g 8
    if sys.argv[2] == 'True':
        filtered = True
        full = False
        mode = "_mode_filtered_"
    elif sys.argv[2] == 'False':
        filtered = False
        full = True
        mode = "mode_full_"
    dataset_name = sys.argv[3] #e.g "KN1" 
    kmerlength = int(kmersize)
    newtoken_size = 1+(kmerlength*2)
    pathname = os.path.dirname(sys.argv[0])
    WORKING_DIR = os.path.abspath(pathname)
    all_tokens = createNewtokenSet(kmerlength)
    if kmerlength > 4:
        stpwrds = make_stopwords(kmerlength)
    else:
        filtered = False
        full = True
        mode = "mode_full_"
        print("for k < 5 only full mode is available!")
    expected_tokens = len(all_tokens)
    run_id = str(int(time.time()))
    file_name = WORKING_DIR + '/output/kgrammar_bag-of-k-mers_model_' + run_id +'_' + dataset_name +'_'+str(kmerlength)+'_'+mode+'.txt'
    
    logging.basicConfig(level=logging.INFO, filename=file_name, filemode="a+",
                        format="%(asctime)-15s %(levelname)-8s %(message)s")
    logging.info("kmer_grammar_bag-of-k-mers RUN ID")
    logging.info(run_id)
    

logging.info("WORKING_DIR")
logging.info(WORKING_DIR)
logging.info("input: kmerlength")
logging.info(str(kmersize))
logging.info("input: dataset")
logging.info(str(dataset_name))
logging.info("input: filtered")
logging.info(filtered)



print('*' * 80)
print("Kgrammer run id: ", run_id)
print("-d %s -k %d -filtered %s" % (dataset_name, kmerlength, str(sys.argv[2])))

#265_new_dap
f1=open("example/"+sys.argv[3]+"_train.txt",'r')
train_labels=[]
test_labels=[]
test_tokens=[]
train_tokens=[]
for i in f1:
	z=i.split(",")
	train_tokens.append(str(z[0]))
	train_labels.append(int(z[1]))


f1.close()

dftrain = pd.DataFrame(list(zip(train_tokens, train_labels)), columns =["dna_string","bound"]) 

f1=open("example/"+sys.argv[3]+"_test.txt",'r')

for i in f1:
	z=i.split(",")
	test_tokens.append(str(z[0]))
	test_labels.append(int(z[1]))


f1.close()
dftest = pd.DataFrame(list(zip(test_tokens, test_labels)), columns =["dna_string","bound"])

print("Collecting tokens")
dftrain["tokens"] = dftrain["dna_string"].apply(write_ngrams)
dftest["tokens"] = dftest["dna_string"].apply(write_ngrams)
train_tokens = dftrain["tokens"].tolist()
test_tokens = dftest["tokens"].tolist()
print("Collecting labels")
train_labels = dftrain["bound"].tolist()
test_labels = dftest["bound"].tolist()
unique_train_labels = len(list(set(train_labels)))
unique_test_labels = len(list(set(test_labels)))


#Check that labels are as many as expected for binary classification
if unique_train_labels < 2 or unique_test_labels < 2:
    print("ERROR: Expected 2 train and test labels. Got %d train labels and %d test labels" % (unique_train_labels, unique_test_labels))
    logging.info("Unique train labels = %d" % unique_train_labels)
    logging.info("Unique test labels = %d" % unique_test_labels)
    print("log file: "+ WORKING_DIR +'/' + file_name)
    quit()

Y_DEV = np.asarray(train_labels)
Y_holdout = np.asarray(test_labels)

print("Building a vocabulary from tokens")
tmpvectorizer = TfidfVectorizer(min_df = 1 , max_df = 1.0, sublinear_tf=True,use_idf=True)
X_TFIDF_ALL =  tmpvectorizer.fit_transform(all_tokens) #newtoken sequences to numeric index.
vcblry = tmpvectorizer.get_feature_names()

if full:
    print("keeping all low-complexity k-mers")
    kmer_names = vcblry
    feature_names = np.asarray(kmer_names) #key transformation to use the fancy index into the report
else:
    print("removing %d low-complexity k-mers" % len(stpwrds))
    kmer_names = [x for x in vcblry if x not in stpwrds]
    feature_names = np.asarray(kmer_names) #key transformation to use the fancy index into the report

#Check that tokens are as many as expected math.pow(4, kmerlength)/2
if len(kmer_names) > expected_tokens:
    print("ERROR: Expected %d tokens. Obtained %d tokens" % (expected_tokens, len(kmer_names)))
    logging.info("Expecting %d tokens" % expected_tokens)
    logging.info("Feature index contains %d tokens" % len(kmer_names))
    logging.info("ERROR: expected %d tokens, got %d tokens" % (expected_tokens, len(kmer_names)))
    logging.info("ERROR: More features than expected!")
    print("log file: "+ WORKING_DIR +'/' + file_name)
    quit()
else:
    print("Expected %d tokens. Obtained %d tokens" % (expected_tokens, len(kmer_names)))
    logging.info("Feature index contains %d tokens" % len(kmer_names))

print("Extracting features from the training data using TfidfVectorizer")
vectorizer = TfidfVectorizer(min_df = 1 , max_df = 1.0, sublinear_tf=True,use_idf=True,vocabulary=kmer_names) #vectorizer for kmer frequencies
X_TFIDF_DEV =  vectorizer.fit_transform(train_tokens)

print("train_samples: %d, n_features: %d" % X_TFIDF_DEV.shape)
print("Positive n_labels: %d Negative n_labels: %d" % (train_labels.count(0),train_labels.count(1)))

logging.info("Train dataset")
logging.info("n_samples: %d, n_features: %d" % X_TFIDF_DEV.shape)
logging.info("Positive n_labels: %d Negative n_labels: %d" % (test_labels.count(0),test_labels.count(1)))

print("Extracting features from the holdout data using TfidfVectorizer")
X_TFIDF_test =  vectorizer.fit_transform(test_tokens)
print("test_samples: %d, n_features: %d" % X_TFIDF_test.shape)
print("Positive n_labels: %d Negative n_labels: %d" % (train_labels.count(0),train_labels.count(1)))

logging.info("Test dataset")
logging.info("n_samples: %d, n_features: %d" % X_TFIDF_test.shape)
logging.info("Positive n_labels: %d Negative n_labels: %d" % (test_labels.count(0),test_labels.count(1)))

print("Fiting a LogisticRegression (LR) model to the training set") 
TFIDF_LR = LogisticRegression(C=1.0, class_weight=None, dual=False, fit_intercept=True,
          intercept_scaling=1, max_iter=100, multi_class='ovr', n_jobs=1,
          penalty='l2', random_state=None, solver='liblinear', tol=0.0001,
          verbose=0, warm_start=False)
logging.info(TFIDF_LR)
TFIDF_LR.fit(X_TFIDF_DEV, Y_DEV)

LR_hold_TFIDF_pred = TFIDF_LR.predict(X_TFIDF_test) # y_pred
LR_hold_TFIDF_prob = TFIDF_LR.predict_proba(X_TFIDF_test)[:,1] # y_score

auc = metrics.roc_auc_score(Y_holdout, LR_hold_TFIDF_prob)

f3 = open("output/"+sys.argv[3]+"_grammar_result.txt",'w')
tLR_hold_TFIDF_pred = TFIDF_LR.predict(X_TFIDF_DEV) # y_pred
tLR_hold_TFIDF_prob = TFIDF_LR.predict_proba(X_TFIDF_DEV)[:,1] # y_score

auc_tr = metrics.roc_auc_score(Y_DEV, tLR_hold_TFIDF_prob)
TP, FP, TN, FN, TPR, TNR, PPV, NPV, FNR, FPR, FDR, FOR, ACC, F1, MCC, BM, MK = comparison(Y_DEV, tLR_hold_TFIDF_pred)
f3.writelines("\t"+"TP"+"\t"+"TN"+"\t"+"FN"+"\t"+"FP"+"\t"+"TPR"+"\t"+"TNR"+"\t"+"ACC"+"\t"+"F1"+"\t"+"MCC"+"\t"+"auc"+"\n")
f3.writelines("Train"+"\t"+str(TP)+"\t"+str(TN)+"\t"+str(FN)+"\t"+str(FP)+"\t"+str(TPR)+"\t"+str(TNR)+"\t"+str(ACC)+"\t"+str(F1)+"\t"+str(MCC)+"\t"+str(auc_tr)+"\n")

TP, FP, TN, FN, TPR, TNR, PPV, NPV, FNR, FPR, FDR, FOR, ACC, F1, MCC, BM, MK = comparison(Y_holdout, LR_hold_TFIDF_pred)

accuracy = metrics.accuracy_score(Y_holdout, LR_hold_TFIDF_pred)
print("accuracy", accuracy, auc)

f3.writelines("Test"+"\t"+str(TP)+"\t"+str(TN)+"\t"+str(FN)+"\t"+str(FP)+"\t"+str(TPR)+"\t"+str(TNR)+"\t"+str(ACC)+"\t"+str(F1)+"\t"+str(MCC)+"\t"+str(auc)+"\n")


f3.close()

logging.info("LR evaluation")
logging.info(metrics.classification_report(Y_holdout, LR_hold_TFIDF_pred)) #y_true, y_pred
logging.info("Accuracy")
logging.info(metrics.accuracy_score(Y_holdout, LR_hold_TFIDF_pred)) #y_true, y_pred
logging.info("ROC_AUC")
logging.info(metrics.roc_auc_score(Y_holdout, LR_hold_TFIDF_prob)) #y_true, y_score
logging.info(auc)

fpr, tpr, thresholds = metrics.roc_curve(Y_holdout, LR_hold_TFIDF_prob, pos_label=1)
roc_auc = metrics.auc(fpr, tpr)
roc_figure_file = WORKING_DIR + "/output/kgrammar_bag-of-k-mer_model_roc_" + mode + dataset_name + "_" + kmersize + "_" + run_id + ".png"
save_plot_roc(fpr, tpr, roc_auc, roc_figure_file, dataset_name)
        
precision, recall, thresholds = metrics.precision_recall_curve(Y_holdout, LR_hold_TFIDF_prob, pos_label=1)
avg_prc = metrics.average_precision_score(Y_holdout, LR_hold_TFIDF_prob)
prc_figure_file = WORKING_DIR + "/output/kgrammar_bag-of-k-mer_model_prc" + mode + dataset_name + "_" + kmersize + "_" + run_id + ".png"
save_plot_prc(precision, recall, avg_prc, prc_figure_file, dataset_name)

#Export the kmer weights from the LR classifier to a sqlite3 database
if hasattr(TFIDF_LR, 'coef_'):
    top = np.argsort(TFIDF_LR.coef_[0])[-5:] #select the top 5 index
    botton = np.argsort(TFIDF_LR.coef_[0])[:5] #select the bottom 5 index
    logging.info("database table LR_results")
    logging.info("top 5 positive kmers")
    logging.info(" ".join([ i.split('n')[0].upper() for i in feature_names[top] ]))
    logging.info(" ".join([ i.split('n')[1].upper() for i in feature_names[top] ]))
    logging.info("top 5 negative kmers")
    logging.info(" ".join([ i.split('n')[0].upper() for i in feature_names[botton] ]))
    logging.info(" ".join([ i.split('n')[1].upper() for i in feature_names[botton] ]))
    print("Saving data to database table LR_results")
    print('*' * 80)
    print("%s: %s" % ("pos kmers", " ".join([ i.split('n')[0].upper() for i in feature_names[top] ]) ))
    print("%s: %s" % ("pos kmers", " ".join([ i.split('n')[1].upper() for i in feature_names[top] ]) ))
    print() #making room
    print("%s: %s" % ("neg kmers", " ".join([ i.split('n')[0] for i in feature_names[botton] ]) ))
    print("%s: %s" % ("neg kmers", " ".join([ i.split('n')[1] for i in feature_names[botton] ]) ))
    print('*' * 80)
    print() #making room
    
    LR_weights = []
    for idx, kmer_score in enumerate(TFIDF_LR.coef_[0]):
        features = feature_names[idx].split('n')
        LR_weights.append({'kmer':features[0].upper(),'revcomp':features[1].upper(),'score':kmer_score})
    LR_weights_feature = pd.DataFrame(LR_weights)
    outengine = 'sqlite:///'+ WORKING_DIR + "/output/kgrammar_bag-of-k-mers_weights" + mode + dataset_name + "_" + kmersize + "_" + run_id + ".db"
    disk_engine = sqlalchemy.create_engine(outengine) 
    LR_weights_feature.to_sql('LR_results', disk_engine, if_exists='append')
    df = pd.read_sql_query('SELECT * FROM LR_results LIMIT 3', disk_engine)
    print(df.head())

#Saving models to disk for model persistency
LR_model_file = WORKING_DIR + "/output/kgrammar_bag-of-k-mers_LR_" + mode + dataset_name + "_" + kmersize + "_" + run_id +".pkl"
joblib.dump(TFIDF_LR, LR_model_file) 
logging.info("Saving model to disk")
logging.info(LR_model_file)
#Later load back the pickled model (possibly in another Python process) with:
#clf = joblib.load('filename.pkl') 
    
print("log file: "+ file_name)
print("model file LR: "+LR_model_file)
print("plot ROC for LR model: "+roc_figure_file)
print("plot PRC for LR model: "+prc_figure_file)
print("result database: "+outengine)
print('*' * 80)
print()
