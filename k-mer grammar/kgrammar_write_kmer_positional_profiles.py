#!/usr/bin/env python3

# -------------------------------------------------------------------------
## #!/workdir/mm2842/anaconda/bin/python
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
import re
import pandas as pd
import numpy as np
import pickle
import time
import logging
import sqlalchemy

from scipy.stats import gaussian_kde
from itertools import product


def make_newtoken(kmer):
    kmer = kmer.lower()
    newtoken = "n".join(sorted([kmer,kmer.translate(str.maketrans('tagc', 'atcg'))[::-1]]))
    return newtoken

def createTokenSet(kmersize, token=True):
    kmerSet = set()
    nucleotides = ["a", "c", "g", "t"]    
    kmerall = product(nucleotides, repeat=kmersize)
    for i in kmerall:
        kmer = ''.join(i)
        kmerSet.add(kmer)
    uniq_kmers = sorted(list(kmerSet))
    if token:
        print("returning %d kmers" % (len(uniq_kmers)))
        return uniq_kmers
    else:
        newtokenSet = set()
        for kmer in uniq_kmers:
            newtoken = make_newtoken(kmer)
            newtokenSet.add(newtoken)  
        uniq_newtokens = sorted(list(newtokenSet))
        print("returning %d newtokens" % (len(uniq_newtokens)))
        return uniq_newtokens   

def newtoken_kmer(lst_nwtk, knwtkn):
    return [(nwtk.split("n")[0].upper(), nwtk.split("n")[1].upper()) for nwtk in lst_nwtk if nwtk != knwtkn]

def load_sequences(dataset_name, positive=True):
    inengine = 'sqlite:////Users/mm2842/notebooks/CNNs_NLP/data/model_db_all_datasets/' + dataset_name + '/data_model.db'
    dbcon = sqlalchemy.create_engine(inengine)
    if positive:
        query = "SELECT * FROM train WHERE bound = 1 UNION ALL SELECT * FROM test WHERE bound = 1"
        dfquery = pd.read_sql_query(query, dbcon)
        dfquery.columns = ["chr_num","left_idx","right_idx","dna_string","bound"]
    else:
        query = "SELECT * FROM train WHERE bound = 0 UNION ALL SELECT * FROM test WHERE bound = 0"
        dfquery = pd.read_sql_query(query, dbcon)
        dfquery.columns = ["chr_num","left_idx","right_idx","dna_string","bound"]
    return dfquery


if sys.argv[1] == "-help" : #or 
    print("Usage: python kgrammar_write_kmer_positional_profiles.py [kmersize, integer] [dataset_name, string]")
    print("Example: python kgrammar_write_kmer_positional_profiles.py 8 FEA4")
    quit()
else:
    kmerlength = int(sys.argv[1]) #e.g 8
    dataset_name = sys.argv[2] #e.g "FEA4"
    pathname = os.path.dirname(sys.argv[0])
    run_id = str(int(time.time()))
    WORKING_DIR = os.path.abspath(pathname)
    file_name = WORKING_DIR + '/output/'+ dataset_name+'/kmer_positional_profiles_' + run_id +'_' + dataset_name +'_'+str(kmerlength)+'.txt'
    logging.basicConfig(level=logging.INFO, filename=file_name, filemode="a+",
                        format="%(asctime)-15s %(levelname)-8s %(message)s")
    logging.info("kmer_positional_profiles RUN ID")
    logging.info(run_id)
    

logging.info("WORKING_DIR")
logging.info(WORKING_DIR)
logging.info("input: kmerlength")
logging.info(str(kmerlength))
logging.info("input: dataset")
logging.info(str(dataset_name))
logging.info("input dataset")
logging.info(dataset_name)     
uniq_newtokens = createTokenSet(kmerlength, token=False)
logging.info("total newtokens")
logging.info(len(uniq_newtokens))

dfpositive = load_sequences(dataset_name)
dfnegative = load_sequences(dataset_name, positive=False)

logging.info("load positive")
logging.info(dfpositive.shape)
logging.info("load negative")
logging.info(dfnegative.shape)

sequences = dfpositive['dna_string'].tolist()
nullsequences = dfnegative['dna_string'].tolist()

density_tokens = dict()
nulldensity_tokens = dict()
xs = np.linspace(0.0, 1.0, num=300)
notfound = [0 for i in xs]
lambda_val = 0.1

logging.info("lambda")
logging.info(lambda_val)
logging.info("linspace")
logging.info(len(xs))

for token in uniq_newtokens:
    combo = [i.upper() for i in token.split("n")]
    label = " ".join(combo)
    matcher = re.compile('|'.join(combo))
    matches = []
    nullmatches = []
    for idx, seq in enumerate(sequences):
        nseq = nullsequences[idx]
        for match in matcher.finditer(seq):
            if len(match.group()) >= 1:
                mmid_val = match.start()+int((match.end()-match.start())/2)
                md = {"idseq":idx,
                      "start":match.start(), 
                      "end":match.end(),
                      "middle":mmid_val,
                      "length":len(match.group()),
                      "repeat":match.group()}
                if md:
                      matches.append(md)
        for nullmatch in matcher.finditer(nseq):
            if len(nullmatch.group()) >= 1:
                nullmmid_val = nullmatch.start()+int((nullmatch.end()-nullmatch.start())/2)
                nullmd = {"idseq":idx,
                      "start":nullmatch.start(), 
                      "end":nullmatch.end(),
                      "middle":nullmmid_val,
                      "length":len(nullmatch.group()),
                      "repeat":nullmatch.group()}
                if nullmd:
                      nullmatches.append(nullmd)
            
    if matches:
        stats_semkmer = pd.DataFrame(matches)
        skcoord = [i/300 for i in stats_semkmer['middle'].tolist()]
        sdensity = gaussian_kde(skcoord)
        sdensity.covariance_factor = lambda : lambda_val
        sdensity._compute_covariance()
        density_tokens[label] = sdensity(xs)
    else:
        density_tokens[label] = notfound
        
    if nullmatches:
        nullstats_semkmer = pd.DataFrame(nullmatches)
        nullskcoord = [i/300 for i in nullstats_semkmer['middle'].tolist()]
        nullsdensity = gaussian_kde(nullskcoord)
        nullsdensity.covariance_factor = lambda : lambda_val
        nullsdensity._compute_covariance()
        nulldensity_tokens[label] = nullsdensity(xs)
    else:
        nulldensity_tokens[label] = notfound    

picklefile = 'kgrammar_write_kmer_positional_profiles'+ dataset_name +'_'+str(kmerlength)+'.pkl'
picklepath = WORKING_DIR + '/output/'+ dataset_name +picklefile

with open(picklepath, 'wb') as f:
    pickle.dump(density_tokens, f)

nullpicklefile = '/kgrammar_write_kmer_positional_profiles'+ dataset_name +'_'+str(kmerlength)+'.pkl'
nullpicklepath = WORKING_DIR + '/output/'+ dataset_name + nullpicklefile

with open(nullpicklepath, 'wb') as f:
    pickle.dump(nulldensity_tokens, f)

print("density files are ready")
print("*" * 80)    

   