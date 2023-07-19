import numpy as np

from sklearn.model_selection import train_test_split
import sys, os
import sys,os
infile_pos=[x.strip() for x in open(sys.argv[1]).readlines()]
infile_neg=[x.strip() for x in open(sys.argv[2]).readlines()]

p_train, p_test, n_train, n_test = train_test_split(infile_pos, infile_neg, train_size=0.70, shuffle=False)

f1= open(sys.argv[3],"w+")
for i in p_train:
	f1.write(i+'\n')
f1.close()

f1= open(sys.argv[3],"a+")
for i in n_train:
	f1.write(i+'\n')
f1.close()

f1= open(sys.argv[4],"w+")
for i in p_test:
	f1.write(i+'\n')
f1.close()

f1= open(sys.argv[4],"a+")
for i in n_test:
	f1.write(i+'\n')
f1.close()
