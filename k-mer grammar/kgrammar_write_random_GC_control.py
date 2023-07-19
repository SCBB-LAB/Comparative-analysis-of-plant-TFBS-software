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

import sys
import os
import csv
import sqlalchemy
import pandas as pd
import numpy as np
import logging
import time
from pathlib import Path

if sys.argv[1] == "-help":
    print("Usage: python kgrammar_write_random_GC_control.py [data_set_name] [genome_version, string] [input_file, string] [output_file, string] [radious, int]")
    print("Example: python kgrammar_write_random_GC_control.py ZmB73_AGPv3 train.csv control_train.csv 1000000")
    quit()
else:
    pathname = os.path.dirname(sys.argv[0])
    WORKING_DIR = os.path.abspath(pathname)
    genome_version = str(sys.argv[1])
    data_set_name = str(sys.argv[2])
    input_file = str(sys.argv[3]) #a comma separated values file, a fasta file is not acceptable
    output_file = str(sys.argv[4])
    radious = int(sys.argv[5])
    run_id = str(int(time.time()))
    file_name = WORKING_DIR + '/output/kgrammar_write_random_GC_control_' + run_id +'.txt'
    logging.basicConfig(level=logging.INFO, filename=file_name, filemode="a+",
                        format="%(asctime)-15s %(levelname)-8s %(message)s")
    logging.info("kgrammar_write_random_GC_control RUN ID")
    logging.info(run_id)

logging.info("WORKING_DIR")
logging.info(WORKING_DIR)
logging.info("log_file")
logging.info(file_name)
logging.info("genome_version")
logging.info(str(genome_version))
logging.info("input: data_set_name")
logging.info(str(data_set_name))
logging.info("input: input_file")
logging.info(str(input_file))
logging.info("input: output_file")
logging.info(str(output_file))
logging.info("input: radious") #value to add at the boundaries of each region to search for controls 
logging.info(str(radious))

#gemome file is a tab separated file with chromosome name and chromosome length for the AGPv3 version of the reference maize genome
genome_dict = {}
genome_file = WORKING_DIR+'/input_databases/'+genome_version+'_genome_file.tsv'
path_to_genome_file = Path(genome_file)

#check if gemome file doesn't exist or if empty before continue
if path_to_genome_file.is_file():
    logging.info("genome file")
    logging.info(genome_file)
    with open(genome_file) as genomefile:
        for line in genomefile:
            fields = line.strip().split()
            genome_dict[fields[0]] = int(fields[1])
else:
    raise ValueError("file is missing! something is wrong with path to: "+genome_file);

#db with index to optimize the query
input_db_file = '/input_databases/'+genome_version+'_control_regions.sql'
engine = 'sqlite:///'+WORKING_DIR+input_db_file

#check if db file doesn't exist or if empty before continue
if os.path.exists(os.path.join(WORKING_DIR,input_db_file)):
    if os.stat(os.path.join(WORKING_DIR,input_db_file)).st_size == 0:
        logging.info("empty database file")
        raise ValueError("file is empty! something is wrong with path to: "+engine);
    else:
        logging.info("database for control regions")
        logging.info(engine)
        dbcon = sqlalchemy.create_engine(engine)
else:
    logging.info("input file doesn't exist")
    raise ValueError("file doesn't exist! something is wrong with path to: "+engine);

#check if input file doesn't exist or if empty before continue
if os.path.exists(os.path.join(WORKING_DIR,input_file)):
    if os.stat(os.path.join(WORKING_DIR,input_file)).st_size == 0:
        logging.info("empty input file")
        raise ValueError("input file is empty! something is wrong with file: "+input_file);
    else:
        print("Input file: ", os.path.join(WORKING_DIR,input_file))
else:
    logging.info("input file doesn't exist")
    raise ValueError("input file is missing! something is wrong with path to: "+input_file);

id_dict = {}
sample_dicts = {}

#input file is expected to be a csv file without header
#[chromosome],[start],[end],[dna_string]
#fields after dna_string will be ignored
with open(os.path.join(WORKING_DIR,input_file)) as csvfile:
    filereader = csv.reader(csvfile, delimiter=',')
    for row in filereader:
        query_region = str(row[3])+" "+str(row[0])+":"+str(row[1])+"-"+str(row[2]) 

        try_close = True
        try_again = False
        
        up_boundary = genome_dict.get(row[0])
        
        if try_close:
            values=[]
            
            #chromosome
            values.append(int(row[0]))
            
            #if using the values in the example the lower boundary would be 1 Mb from the start
            if int(row[1]) - radious < 0:
                values.append(0)
            else:
                values.append(int(row[1]) - radious)
            
            #if using the values in the example the upper boundary would be 1 Mb from the end
            if int(row[2]) + radious > up_boundary:
                values.append(up_boundary)
            else:
                values.append(int(row[2]) + radious)
            #GC%
            values.append(round(float(row[5])-0.01,3))
            values.append(round(float(row[5])+0.01,3))
            
            #length query region
            queryLength = int(row[2])-int(row[1])
            
            #if queryLength = 300, table2query have no overlapping regions of length 300 that doesn't include the queries (i.e., norep_regions_300)
            #in the database there are only tables availables for lengths 100 and 300
            table2query = "regions_" + str(queryLength) 
            query = ("SELECT * from "+table2query+" WHERE chr = %d AND (START BETWEEN %d AND %d) AND (GCperc BETWEEN %f AND %f) AND N = 0" % tuple(values))
            df_query = pd.read_sql_query(query, dbcon)
            if df_query.shape[0] >= 1:
                random_sample = df_query.take(np.random.permutation(len(df_query))[:1])
                random_id = "chr"+random_sample['chr'].map(str)+":"+random_sample['start'].map(str)+"-"+random_sample['end'].map(str)
                sample_id = random_id.to_string(index=False)
                if sample_id in sample_dicts:
                    try_again = True
                else:
                    random_sample.insert(3, 'id', sample_id, allow_duplicates=False)
                    random_sample.to_csv(os.path.join(WORKING_DIR,output_file), mode='a',
                                 header=False,index=False, sep=',')
                    id_dict[str(row[3])] = sample_id
                    sample_dicts[sample_id] = str(row[3])
                    try_again = False
            else:
                values=[]
                values.append(int(row[0])) 
                values.append(0) #try as hard as possible
                values.append(up_boundary) #try as hard as possible
                values.append(round(float(row[5])-0.01,3))
                values.append(round(float(row[5])+0.01,3))
                query = ("SELECT * from "+table2query+" WHERE chr = %d AND (START BETWEEN %d AND %d) AND (GCperc BETWEEN %f AND %f) AND N = 0" % tuple(values))
                df_rquery = pd.read_sql_query(query, dbcon)
                if df_rquery.shape[0] >= 1:
                    random_sample = df_rquery.take(np.random.permutation(len(df_rquery))[:1])
                    random_id = "chr"+random_sample['chr'].map(str)+":"+random_sample['start'].map(str)+"-"+random_sample['end'].map(str)
                    sample_id = random_id.to_string(index=False)
                    if sample_id in sample_dicts:
                        if df_rquery.shape[0] > 1:
                            try_again = True
                        else:
                            logging.info("No control available for: ",query_region)
                            print("No more regions available in the same chromosome: ", str(row[3]))
                            id_dict[str(row[3])] = "NA"
                            try_again = False
                    else: 
                        random_sample.insert(3, 'id', sample_id, allow_duplicates=False)
                        random_sample.to_csv(os.path.join(WORKING_DIR,output_file), mode='a',
                                 header=False,index=False, sep=',')
                        id_dict[str(row[3])] = sample_id
                        sample_dicts[sample_id] = str(row[3])
                        try_again = False
                else:
                    logging.info("No control available for: ",query_region)
                    print("Region without a control in chromosome: ", str(row[3]))
                    id_dict[str(row[3])] = "NA"
                    try_again = False
                    
        while try_again:
            values=[]
            values.append(int(row[0])) 
            values.append(0) #try as hard as possible
            values.append(up_boundary) #try as hard as possible
            values.append(round(float(row[5])-0.01,3))
            values.append(round(float(row[5])+0.01,3))
            query = ("SELECT * from "+table2query+" WHERE chr = %d AND (START BETWEEN %d AND %d) AND (GCperc BETWEEN %f AND %f) AND N = 0" % tuple(values))
            df_rquery = pd.read_sql_query(query, dbcon)
            if df_rquery.shape[0] >= 1:
                random_sample = df_rquery.take(np.random.permutation(len(df_rquery))[:1])
                random_id = "chr"+random_sample['chr'].map(str)+":"+random_sample['start'].map(str)+"-"+random_sample['end'].map(str)
                sample_id = random_id.to_string(index=False)
                if sample_id in sample_dicts:
                    try_again = True
                else:
                    random_sample.insert(3, 'id', sample_id, allow_duplicates=False)
                    random_sample.to_csv(os.path.join(WORKING_DIR,output_file), mode='a',
                                 header=False,index=False, sep=',')
                    id_dict[str(row[3])] = sample_id
                    sample_dicts[sample_id] = str(row[3])
                    try_again = False
            else:
                logging.info("No control available for: ",query_region)
                print("Unlikely, but got out without a control in chromosome: ", str(row[3]), " for query_region: ", query_region )
                id_dict[str(row[3])] = "NA"
                try_again = False


#make a log file with the pairs of the ids for the input file and the query results
df_ids = pd.DataFrame([[key,value] for key,value in id_dict.items()],columns=["provided_id","query_res"])
file_pairs = str(data_set_name) + "_ids_sample_ids_nulls.csv"
df_ids.to_csv(os.path.join(WORKING_DIR,file_pairs), mode='a', header=True,index=False, sep=',')
print("IDs pairs writen to file: ", os.path.join(WORKING_DIR,file_pairs))
print("Output file is ready: ", os.path.join(WORKING_DIR,output_file))