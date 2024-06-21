import pandas as pd, os, sys
filename=str(sys.argv[1])
seq=[i.strip() for i in os.popen("grep -v '>' temp.fa").readlines()]

master=[]
for i in range(len(seq[0])):
	l1=[]	
	a = {  "A": 0, "C": 0,  "G": 0,  "T": 0}
	for j in seq:
		k=j[i].upper()
		if k=='A' or k=='T' or k=='C' or k=='G':
			a[k]+=1
		else:
			a['A']+=1
	for j in a:
		l1.append(a[j]/len(seq))
	master.append(l1)


df=pd.DataFrame(master,columns =['A', 'C', 'G','T'])
df=df.T
df.to_csv(filename+".csv")
os.system("sed -i '1d' "+filename+".csv")
os.system("Rscript 1.R "+filename+".csv "+filename+".pdf")

master=[]
for i in range(len(seq[0])):
	l1=[]	
	a = {  "A": 0, "C": 0,  "G": 0,  "T": 0}
	for j in seq:
		k=j[i].upper()
		if k=='A' or k=='T' or k=='C' or k=='G':
			a[k]+=1
		else:
			a['A']+=1
	for j in a:
		l1.append(a[j])
	master.append(l1)


df=pd.DataFrame(master,columns =['A', 'C', 'G','T'])
df=df.T
df.to_csv(filename+"_pfm.csv")

