#!/usr/bin/env python
# encoding: utf-8

import numpy as np
import os
import time
from optparse import OptionParser
from keras.preprocessing.text import Tokenizer
#from keras.preprocessing.sequence import pad_sequences
from tensorflow.keras.preprocessing.sequence import pad_sequences
from keras.models import Sequential
from keras.layers import Dense, LSTM, Dropout, Bidirectional, Flatten, GRU
#from keras.layers.embeddings import Embedding
from tensorflow.keras.layers import Embedding
from keras.callbacks import ModelCheckpoint, EarlyStopping
from keras.optimizers import RMSprop,SGD,Adam,Adagrad
from sklearn import metrics
from sklearn.metrics import roc_auc_score ,accuracy_score
from gensim.models import Word2Vec
from keras.models import load_model
os.environ["CUDA_VISIBLE_DEVICES"] = "-1"

np.random.seed(12345)

################################################################################
# dbrmodel.py
#
#
# DBRNN model run.
################################################################################

################################################################################
# main
################################################################################
def main():
    usage = 'usage'
    parser = OptionParser(usage)
    #parser.add_option('-p', dest='file_path', default='/home/szhen/szhen/DBLSTM/test2', help='train data file path')
    parser.add_option('-p', dest='file_path', default='/DATA/kegru_code/sampledata', help='train data file path')
    parser.add_option('-g', dest='gpu', type=int, default=0, help='using which gpu')
    parser.add_option('-n', dest='name', default='init', type=str, help='first word about file name')
    parser.add_option('-k', dest='k', default=5, type=int, help='kmer length')
    parser.add_option('-s', dest='s', default=2, type=int, help='stride when slicing k-mers')
    parser.add_option('-b', dest='batchsize', type=int, default=200, help='size of one batch')
    parser.add_option('-i', dest='init', action='store_true', default=True, help='initialize vector')
    parser.add_option('-t', dest='trainable', action='store_true', default=True, help='embedding vectors trainable')
    parser.add_option('-T', dest='test', action='store_true', default=False, help='only test step')
    parser.add_option('-l', dest='layer', default=1, type=int, help='the number of layers')
    parser.add_option('-B', dest='bidirectional', default=True, help='stride when slicing k-mers')
    parser.add_option('-u', dest='u', default=50, type=int, help='the number of rnn unit')
    parser.add_option('-r', dest='out', type=str, help='result out file')
    parser.add_option('-O', dest='optim', type=str, help='hyper parm')
    parser.add_option('-N', dest='moname', type=str, help='hyper parm')
    parser.add_option('-U', dest='un_num', type=str, help='hyper parm')
    (options,args) = parser.parse_args()

    #hyper patm
    moname = options.moname
    un_num = options.un_num
    optim = options.optim

    #load train seq data
    print('Loading seq data...')
    pos_name = options.file_path + '/' + options.name + '_pos.fa' 
    neg_name = options.file_path + '/' + options.name + '_neg.fa'
    #pos_name = options.file_path + '/' + options.name + '_pos_' + str(options.k) + 'gram_' + str(options.s) + 'stride'
    #neg_name = options.file_path + '/' + options.name + '_neg_' + str(options.k) + 'gram_' + str(options.s) + 'stride'
    pos_seqs = [line[:-2] for line in open(pos_name) if len(line.split()) > 15]
    neg_seqs = [line[:-2] for line in open(neg_name) if len(line.split()) > 15]
    seqs1=[]
    seqs = pos_seqs + neg_seqs
    
    for i in seqs:
        z=i.split(" ")
        seqs1.append(z)
    length=len(z[1])
    lens = [len(line.split()) for line in seqs]
    n_seqs = len(lens)
    y = np.array([1] * len(pos_seqs) + [0] * len(neg_seqs))

    print('model load')

    wmodel = Word2Vec(seqs1, size=100)

    print('model ready')

    print('Tokenizing seqs...')
    MAX_LEN = 1000
    NB_WORDS = 20000
    tokenizer = Tokenizer(nb_words=NB_WORDS)
    tokenizer.fit_on_texts(seqs)
    sequences = tokenizer.texts_to_sequences(seqs)
    X = pad_sequences(sequences, maxlen=MAX_LEN)

    print(f'Train padded shape: {X.shape}')
    kmer_index = tokenizer.word_index
    print('Found %s unique tokens.' % len(kmer_index))

    print('Spliting train, valid, test parts...')
    indices = np.arange(n_seqs)
    np.random.shuffle(indices)
    X = X[indices]
    y = y[indices]
    #n_tr = int(n_seqs * 0.85)
    #n_va = int(n_seqs * 0.05)
    n_te = n_seqs
    
    #X_train = X[:n_tr]
    #y_train = y[:n_tr]
    #X_valid = X[n_tr:n_tr + n_va]
    #y_valid = y[n_tr:n_tr + n_va]
    X_test = X[-n_te:]
    y_test = y[-n_te:]
    #print(X_test)
    #print(y_test)
    embedding_vector_length = 100
    nb_words = min(NB_WORDS, len(kmer_index))  # kmer_index starting from 1
    print('Building model...')
    model = Sequential()
    
    if options.init:
        print('initialize embedding layer with w2v vectors')
        embedding_matrix = np.zeros((nb_words + 1, embedding_vector_length))      
        for kmer, i in kmer_index.items():
            if i > NB_WORDS:
                continue
            elif (len(kmer)==length):                
                kmer=kmer.upper()
                embedding_matrix[i] = wmodel[kmer]
                
                
        print('embedding layers trainable %s' % options.trainable)
        model.add(Embedding(nb_words + 1,
                            embedding_vector_length,
                            weights=[embedding_matrix],
                            input_length=MAX_LEN,
                            trainable=options.trainable))
    else:
    
        model.add(Embedding(nb_words + 1,
                            embedding_vector_length,
                            input_length=MAX_LEN))
    model.add(Dropout(0.2))
    #model.add(Bidirectional(RCG.RevComGRU(50, return_sequences=True)))
    #model.add(Dropout(0.2))
    model.add(Bidirectional(GRU(int(80), return_sequences=True)))
    model.add(Dropout(0.2))
    model.add(Dense(20, activation='relu'))
    model.add(Dropout(0.5))
    model.add(Flatten())
    model.add(Dense(1, activation='sigmoid'))
    model.compile(loss='binary_crossentropy', optimizer="Adam", metrics=['accuracy'])
    print((model.summary()))

    model_path = options.file_path + '/model/' + options.name + '_bestmodel_' + str(options.k) + '_withlstm.hdf5'

    start = 0
    end = 0
    if not options.test:
        start = time.time()
        checkpointer = ModelCheckpoint(filepath=model_path, verbose=1,
                                       save_best_only=True)
        earlystopper = EarlyStopping(monitor='val_loss', patience=6, verbose=1)

        print('Training model...')
        #history_callback = model.fit(X_train, y_train, epochs=100, batch_size=options.batchsize, shuffle=True,
       #
        #          validation_data=(X_valid, y_valid),
         #         callbacks=[checkpointer, earlystopper],
          #        verbose=1)
        #end = time.time()
    
    print('Testing model...')
    model = load_model('data/model/ABF2_bestmodel_5_withlstm.hdf5')    
    #prob = model.predict([x_test_1,x_test_2])
    #model.load_weights(model_path)
    #tresults = model.evaluate(X_test, y_test)
    y_pred = model.predict(X_test, batch_size=options.batchsize, verbose=1)    
    #print(y_pred)
    a_file = open(options.name+"_test.txt", "a+")
    for row in y_pred:
        np.savetxt(a_file, row)   
    y = y_test
    
    print('Calculating AUC...')
    auroc = metrics.roc_auc_score(y, y_pred)
    auprc = metrics.average_precision_score(y, y_pred)
    #print(auroc, auprc)
    
    labels=np.array(y_pred)
    
    labels[labels>=0.5]=1
    labels[labels<0.5]=0
    acc = accuracy_score(y,labels)
    tn, fp, fn, tp = metrics.confusion_matrix(y,labels).ravel()
    # Specificity or true negative rate
    specificity = tn / (tn+fp)
    # Sensitivity, recall, or true positive rate
    sensitivity = tp/(tp+fn)
    mcc=metrics.matthews_corrcoef(y,labels)
    f1=metrics.f1_score(y,labels)
    #print(tn, fp, fn, tp)
    #print(sensitivity)
    #print(specificity)
    #print(acc)
    #print(f1)
    #print(mcc)
    #print(auroc)
   
    print('save results to local file')
    runtime = end - start
    curout = options.file_path.split('/')
    curoutn = curout[-1]
    #myhist = history_callback.history
    f3 = open("result/"+options.name+"_stats",'w')
    
    f3.writelines("TP"+"\t"+"TN"+"\t"+"FN"+"\t"+"FP"+"\t"+"TPR"+"\t"+"TNR"+"\t"+"ACC"+"\t"+"F1"+"\t"+"MCC"+"\t"+"auc"+"\n")
    f3.writelines(str(tp)+"\t"+str(tn)+"\t"+str(fp)+"\t"+str(fn)+"\t"+str(sensitivity)+"\t"+str(specificity)+"\t"+str(acc)+"\t"+str(f1)+"\t"+str(mcc)+"\t"+str(auroc)+"\n")

    #all_hist = np.asarray([myhist["loss"], myhist["accuracy"], myhist["val_loss"], myhist["val_accuracy"]]).transpose()
    #np.savetxt(os.path.join("result/"+options.name+ "_training_history.txt"), all_hist, delimiter="\t",
             #  header='loss\taccuracy\tval_loss\tval_accuracy')
    comstr = 'echo ' +  '\t' + str(auroc) + '\t' + str(auprc) + '\t' + str(curoutn) + '\t' + str(moname) + ' >> ' + str(options.out)
    os.system(comstr)

if __name__ == '__main__':
    main()
