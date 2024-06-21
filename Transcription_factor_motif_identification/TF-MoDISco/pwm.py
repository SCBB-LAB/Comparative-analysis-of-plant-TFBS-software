import os, sys, numpy as np
sfile=[x.strip() for x in open(sys.argv[1]).readlines()] #fasta file
fastadict={};
for line in sfile:
    if not line:
        continue
    if line.startswith('>'):
        sname = line
        if line not in fastadict:
            fastadict[line] = ''
        continue
    fastadict[sname] += line


infile=[x.split() for x in open(sys.argv[2]).readlines()] # pwm file A, C, T, G only not A C G T
f1=np.array(infile,dtype="float32")
m=[]
zz=[]
for i in open(sys.argv[2]):
	z=i.split("\t")
y=len(z)
for j in range(0,y):
	z=[]
	for i in range(0,4):
		z.append(float(infile[i][j]))
	zz.append(z)
	m.append(max(z))

monk=np.array(m).sum()
def mix(string):
	mz=[]
	for i, j in enumerate(string):
		if j=='A':
			mz.append(zz[i][0])
		elif j=='C':
			mz.append(zz[i][1])
		elif j=='T':
			mz.append(zz[i][2])
		elif j=='G':
			mz.append(zz[i][3])
	return (np.array(mz).sum()*100)/monk

ids = list(fastadict.keys())
seq = list(fastadict.values())
#f3=open(sys.argv[3],'w')
f3=open(sys.argv[3]+".fa",'w')
one=[]
for k, j in enumerate(seq):
	sequence=[j[i:i+y] for i in range(len(j)-(y-1))]
	for m, n in enumerate(sequence):
		app = mix(n.upper())
		if app>=75:
#			f3.writelines(str(ids[k])+"\t"+str(m)+"\t"+str(m+y)+"\t"+str(n)+"\t"+str(app)+"\n")
			f3.writelines(">seq\n"+str(n)+"\n")
			one.append(str(ids[k]))

#f3.close()
f3.close()

one = set(one)
print((len(one)/len(seq))*100)

