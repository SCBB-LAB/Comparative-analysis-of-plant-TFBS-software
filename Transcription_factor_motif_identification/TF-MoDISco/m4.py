import os, sys, pandas as pd
file1=str(sys.argv[1])
f1=[i.split() for i in os.popen("sed '1d' jaspar/"+file1+".jaspar|sed 's+\[++'|sed 's+]++'").readlines()]
f2=[i.split() for i in os.popen("awk '{if(NR==1) print $2}' jaspar/"+file1+".jaspar").readlines()][0][0]
ran = len(f1[0])
def minimum(a, n):
	maxpos = a.index(max(a))
	if maxpos==0:
		return("A")
	elif maxpos==1:
		return("C")
	elif maxpos==2:
		return("G")
	else:
		return("T")

ss=[]
for i in range(1,ran):
	a=int(f1[0][i])
	b=int(f1[1][i])
	c=int(f1[2][i])
	d=int(f1[3][i])
	s = [a, b, c, d]
	ss.append(minimum(s, len(s)))


print(f2+"\t"+"".join(ss))
