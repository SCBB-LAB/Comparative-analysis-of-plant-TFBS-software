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
import numpy as np
import pandas as pd
import sqlalchemy
import logging
import multiprocessing
import gensim
import time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from copy import deepcopy
from itertools import product
from sklearn import metrics
from sklearn.metrics import roc_curve 
from sklearn.metrics import auc

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
def make_newtoken(kmer):
    '''make_newtoken'''
    s = set('acgt') #if odd characters return token
    kmer = str(kmer).lower()
    if not set(kmer).difference(s):
        newtoken = "n".join(sorted([kmer,kmer.translate(str.maketrans('tagc', 'atcg'))[::-1]]))
        return newtoken
    else:
        return "token"
    

def split_len(seq, length):
    '''To slide througth a sequence and extract the k-mers'''
    return [''.join(x) for x in zip(*[list(seq[z::length]) for z in range(length)])]


def merge_two_dicts(x, y):
    '''Given two dicts, merge them into a new dict as a shallow copy'''
    z = x.copy()
    z.update(y)
    return z


def write_sentences(sequence):
    '''write_sentences with newtokens'''
    kmersize = kmerlength
    seq = str(sequence).lower()
    sentences = []
    if len(seq) > 0:
        first_sentence_kmers = split_len(seq, kmersize)
        alltokens = [make_newtoken(kmer) for kmer in first_sentence_kmers if len(kmer) == kmersize]
        first_sentence_newtokens = [newtoken for newtoken in alltokens if newtoken != "token"]
        sentences.append(first_sentence_newtokens) #each sentence is a list
        n = kmersize-1
        while n >= 1:
            next_sentence_kmers = split_len(seq[n:], kmersize)
            alltokens = [make_newtoken(kmer) for kmer in next_sentence_kmers if len(kmer) == kmersize]
            next_sentence_newtokens = [newtoken for newtoken in alltokens if newtoken != "token"]
            sentences.append(next_sentence_newtokens)
            n=n-1
    return sentences
    

def write_kmer_sentences(sequence):
    '''write_sentences with kmers'''
    kmersize = kmerlength
    seq = str(sequence).lower()
    sentences = []
    if len(seq) > 0:
        first_sentence_kmers = split_len(seq, kmersize)
        alltokens = [make_newtoken(kmer) for kmer in first_sentence_kmers if len(kmer) == kmersize]
        first_sentence_newtokens = [newtoken for newtoken in alltokens if newtoken != "token"]
        sentences.append(first_sentence_newtokens) #each sentence is a list
        n = kmersize-1
        while n >= 1:
            next_sentence_kmers = split_len(seq[n:], kmersize)
            alltokens = [make_newtoken(kmer) for kmer in next_sentence_kmers if len(kmer) == kmersize]
            next_sentence_newtokens = [newtoken for newtoken in alltokens if newtoken != "token"]
            sentences.append(next_sentence_newtokens)
            n=n-1
    return sentences


def createKmerSet(kmersize):
    '''write all possible kmers'''
    kmerSet = set()
    nucleotides = ["a", "c", "g", "t"]    
    kmerall = product(nucleotides, repeat=kmersize)
    for i in kmerall:
        kmer = ''.join(i)
        kmerSet.add(kmer)
    uniq_kmers = sorted(list(kmerSet))  
    return uniq_kmers


def createNewtokenSet(kmersize):
    '''write all possible newtokens'''   
    newtokenSet = set()
    uniq_kmers = createKmerSet(kmersize)
    for kmer in uniq_kmers:
        newtoken = make_newtoken(kmer)
        newtokenSet.add(newtoken)  
    uniq_newtokens = sorted(list(newtokenSet))
    return uniq_newtokens      
    
            
def sentences_in_category(list_dictionaries, category):
    '''Given a dictionary per document (sequence) return a generator that yield the rated sentences'''
    for dictionary in list_dictionaries:
        if dictionary['bound'] in category:
        #if dictionary['range_category'] in category:
            for sentence in dictionary['range_sentences']:
                yield sentence
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

def sentences_for_vocab(uniq_tokens):
    '''Given a dictionary per document (sequence) return a generator that yield the rated sentences'''            
    for token in uniq_tokens:
        sentence = [token, token, token]
        yield sentence


#TODO: parallelize me! Needs to be improved, so far everything is done in-memory!
def docprob(seqtest, models):
    '''
    docprob is the actual classifier, as uses the w2v output and bayes inversion
    it takes two arguments
    * A list of documents, each of which is a list of sentences
    * The candidate word2vec models (each potential class)
    '''
    docs = [r['range_sentences'] for r in seqtest]
    docs_cats = pd.Series([r['bound'] for r in seqtest])
    #docs_cats = pd.Series([r['range_category'] for r in seqtest])
    sentlist = [s for d in docs for s in d]
    llhd = np.array( [ m.score(sentlist, len(sentlist)) for m in models ] )
    lhd = np.exp(llhd - llhd.max(axis=0)) # subtract row max to avoid numeric overload
    prob = pd.DataFrame( (lhd/lhd.sum(axis=0)).transpose() )
    prob["seq"] = [i for i,d in enumerate(docs) for s in d]
    prob = prob.groupby("seq").mean()
    prob['true_category'] = docs_cats.values
    prob['predict'] = np.where(prob[1] <= prob[0], 0, 1)
    prob['predict_proba'] = np.where(prob[1] <= prob[0], prob[0], prob[1])
    prob.columns = ["category_0","category_1","true_category", "predict","predict_proba"]
    return prob

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
    plt.savefig(figure_file)
            
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
    plt.savefig(figure_file)    

if sys.argv[1] == "-help" : #or 
    print("Usage: python kgrammar_vector-k-mers_training_testing.py [kmersize, integer] [windowsize, integer] [kmer_parsing, 'True' for kmer, 'False' for newtokens] [dataset_name, string]")
    print("Example: python kgrammar_vector-k-mers_training_testing.py 8 5 False FEA4")
    quit()
else:
    kmerlength = int(sys.argv[1]) #e.g 8
    windowsize = int(sys.argv[2]) #e.g 5
    if sys.argv[3] == 'True':
        kmer_parsing = True #in case that vectors for each k-mer are the desired output
    elif sys.argv[3] == 'False':
        kmer_parsing = False
    dataset_name = sys.argv[4] #e.g "FEA4" 
    run_id = str(int(time.time()))
    pathname = os.path.dirname(sys.argv[0])
    WORKING_DIR = os.path.abspath(pathname)
    file_name = WORKING_DIR + '/output/kgrammar_vector-k-mers_model_' + run_id +'_' + dataset_name +'_'+str(kmerlength)+'_'+str(windowsize)+'.txt'
    logging.basicConfig(level=logging.INFO, filename=file_name, filemode="a+",
                        format="%(asctime)-15s %(levelname)-8s %(message)s")
    logging.info("kmer_grammar_vector-k-mers RUN ID")
    logging.info(run_id)

logging.info("WORKING_DIR")
logging.info(WORKING_DIR)
logging.info("input: kmerlength")
logging.info(kmerlength)
logging.info("input: newtoken size")
logging.info(kmerlength*2)
logging.info("input: dataset")
logging.info(str(dataset_name))  
   
#collect all data in the three required formats
print('*' * 80)
print("Dataset: ",dataset_name)


#265_new_dap
f1=open("example/"+sys.argv[4]+"_train.txt",'r')
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

f1=open("example/"+sys.argv[4]+"_test.txt",'r')

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
categories = list(set(train_labels))
#print(unique_train_labels)
unique_test_labels = len(list(set(test_labels)))


if len(categories) != 2:
    logging.info("Unexpected number of categories for binary classifier")
    quit()

if kmer_parsing:
    dftrain["range_sentences"] = dftrain["dna_string"].apply(write_kmer_sentences).astype('object')
    dftest["range_sentences"] = dftest["dna_string"].apply(write_kmer_sentences).astype('object')
else:
    dftrain["range_sentences"] = dftrain["dna_string"].apply(write_sentences).astype('object')
    dftest["range_sentences"] = dftest["dna_string"].apply(write_sentences).astype('object')

logging.info("Train dataset")
logging.info(dftrain.shape)
logging.info("Holdout dataset")
logging.info(dftest.shape)

holdout_set = dftest.to_dict('records')
DEV_set = dftrain.to_dict('records')
alldata = holdout_set + DEV_set #required to build complete vocabularies (word2vec and vectorizer)

##Preparing the embedding models for raw data
assert gensim.models.word2vec.FAST_VERSION > -1 #This will be painfully slow otherwise
num_features = 300 # Word vector dimensionality
min_word_count = 0 # Minimum word count
num_workers = multiprocessing.cpu_count() #use all my cores
context = windowsize  #context word window
hs=3
negative = 0
iterations=30 # Sweeps of Stochastic Gradient Descending through the data

#model with word2vec
print("Preparing the base model ...")
basemodel = gensim.models.Word2Vec(workers=num_workers, iter=iterations, hs=hs, negative=negative,
                              size=num_features, min_count = min_word_count, window = context)

if kmer_parsing:
    all_tokens = createKmerSet(kmerlength)
    expected_newtokens = len(all_tokens)
else:
    all_tokens = createNewtokenSet(kmerlength)
    expected_newtokens = len(all_tokens)

logging.info("Building vocabulary")
basemodel.build_vocab(sentences_for_vocab(all_tokens))
logging.info(basemodel)


if len(basemodel.wv.vocab) > expected_newtokens:
    print("ERROR: Expected %d tokens. Obtained %d tokens" % (expected_newtokens, len(basemodel.wv.vocab)))
    logging.info("Expecting %d tokens" % expected_newtokens)
    logging.info("Feature index contains %d tokens" % len(basemodel.vocab))
    logging.info("ERROR: expected %d tokens, got %d tokens" % (expected_newtokens, len(basemodel.wv.vocab)))
    logging.info("ERROR: More features than expected!")
    print("log file: "+ WORKING_DIR +'/' + file_name)
    quit()
else:
    print(basemodel)
    print("Vocabulary info: Expected %d tokens. Obtained %d tokens" % (expected_newtokens, len(basemodel.wv.vocab)))
    print()
    logging.info("Expecting %d tokens" % expected_newtokens)
    logging.info("Feature index contains %d tokens" % len(basemodel.wv.vocab))
    

catmodels = [deepcopy(basemodel) for each in categories]
print("Iterating through categories to build the model ...")
model_file_base = WORKING_DIR + "/output/kgrammar_vector-k-mer_model_"+ dataset_name + "_" + run_id + "_kmersize_" +  str(kmerlength) + '_windowsize_'+  str(windowsize) 
#model_file_base = WORKING_DIR + "/output/vector-k-mers/" + dataset_name + "/kgrammar_vector-k-mer_model_"+ dataset_name + "_" + run_id + "_kmersize_" +  str(kmerlength) + '_windowsize_'+  str(windowsize) 
for category in categories:

    print(category)
    t0 = time.time()
    logging.info("building model for category: "+str(category))
    slist = list(sentences_in_category(DEV_set, [category]))
    print(category, "category (", len(slist), ")")
    catmodels[category].train(slist, total_examples=len(slist), epochs=basemodel.iter)
    duration = time.time() - t0
    logging.info("done in %fs" % (duration))
    print("done in %fs" % (duration))
    logging.info("saving model for category: "+str(category)) #saving models 
    model_file =model_file_base + '_category_' + str(category)
    catmodels[category].save(model_file)
    logging.info(model_file)
    print()

#predict with word2vec and bayesinversion for the holdout dataset.
np.random.shuffle( holdout_set ) 
print("Predicted labels for holdout set")
w2v_hold_df = docprob(holdout_set, catmodels )
model_test = WORKING_DIR + "/output/kgrammar_vector-k-mer_test_results_"+ dataset_name + "_" + run_id + "_kmersize_" +  str(kmerlength) + '_windowsize_'+  str(windowsize) + '_category_' + str(category)
w2v_hold_df.to_pickle(model_test+'.pkl')
w2v_hold_df.to_csv(model_test+'.csv', sep='\t', encoding='utf-8', index=False)
logging.info(w2v_hold_df.head(3))
w2v_hold_true = w2v_hold_df["true_category"].tolist()
w2v_hold_pred = w2v_hold_df["predict"].tolist()
w2v_hold_prob = w2v_hold_df["category_1"].tolist()

print("Model Evaluation:")
print(metrics.classification_report(w2v_hold_true, w2v_hold_pred))
print("Accuracy score")
print(metrics.accuracy_score(w2v_hold_true, w2v_hold_pred))
print("ROC_AUC")
print(metrics.roc_auc_score(w2v_hold_true, w2v_hold_prob))

auc = metrics.roc_auc_score(w2v_hold_true, w2v_hold_prob)
TP, FP, TN, FN, TPR, TNR, PPV, NPV, FNR, FPR, FDR, FOR, ACC, F1, MCC, BM, MK = comparison(w2v_hold_true, w2v_hold_pred)
f3 = open("output/"+sys.argv[4]+"_vector_result.txt",'w')
f3.writelines("\t"+"TP"+"\t"+"TN"+"\t"+"FN"+"\t"+"FP"+"\t"+"TPR"+"\t"+"TNR"+"\t"+"ACC"+"\t"+"F1"+"\t"+"MCC"+"\t"+"auc"+"\n")
f3.writelines("Test"+"\t"+str(TP)+"\t"+str(TN)+"\t"+str(FN)+"\t"+str(FP)+"\t"+str(TPR)+"\t"+str(TNR)+"\t"+str(ACC)+"\t"+str(F1)+"\t"+str(MCC)+"\t"+str(auc)+"\n")
f3.close()
logging.info("Evaluation report")
logging.info(metrics.classification_report(w2v_hold_true, w2v_hold_pred))
logging.info("ROC_AUC")
logging.info(metrics.roc_auc_score(w2v_hold_true, w2v_hold_prob))
logging.info("Accuracy score")
logging.info(metrics.accuracy_score(w2v_hold_true, w2v_hold_pred))

fpr, tpr, thresholds = roc_curve(w2v_hold_true, w2v_hold_prob, pos_label=1)
#roc_auc = auc(fpr, tpr)
roc_figure_file = WORKING_DIR + "/output/kgrammar_vector-k-mer_model_roc_" + dataset_name + "_" + str(kmerlength) + "_" + run_id + ".png"
save_plot_roc(fpr, tpr, auc, roc_figure_file, dataset_name)

precision, recall, thresholds = metrics.precision_recall_curve(w2v_hold_true, w2v_hold_prob, pos_label=1)
avg_prc = metrics.average_precision_score(w2v_hold_true, w2v_hold_prob)
prc_figure_file = WORKING_DIR + "/output/kgrammar_vector-k-mer_model_prc_" + dataset_name + "_" + str(kmerlength) + "_" + run_id + ".png"
save_plot_prc(precision, recall, avg_prc, prc_figure_file, dataset_name)

print("log file: "+ file_name)
for category in categories:
    print("model files: "+ model_file_base + str(category))
print("plot ROC : "+roc_figure_file)
print("plot PRC : "+prc_figure_file)
print("results for testing: "+model_test+'.csv')
print('*' * 80)
print()
