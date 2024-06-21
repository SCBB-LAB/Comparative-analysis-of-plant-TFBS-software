import os, sys, pandas as pd
import pandas as pd
#file1=str(sys.argv[1])
file1=str(sys.argv[1])

f1=[i.split() for i in os.popen("sed '1,11d' "+file1+".txt|sed '/^$/d'").readlines()]
df = pd.DataFrame(f1)
df =df.T
df.to_csv(file1+'.csv', index=False)

f1=[i.split() for i in os.popen("sed '1d' "+ sys.argv[3]+"_modisco" + file1+".csv| tr ',' '\t'").readlines()]
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
for i in range(0,ran):
	a=float(f1[0][i])
	b=float(f1[1][i])
	c=float(f1[2][i])
	d=float(f1[3][i])
	s = [a, b, c, d]
	ss.append(minimum(s, len(s)))


f2 = []
f2.append(f1[0])
f2.append(f1[1])
f2.append(f1[3])
f2.append(f1[2])

df = pd.DataFrame(f2)
df.to_csv(file1+'.csv', index=False)
os.system("sed '1d' "+file1+".csv| tr ',' '\t' >"+file1+"t.csv")
os.system("mv "+file1+"t.csv "+file1+".csv")


f3=[i.split() for i in os.popen("python3.8 pwm.py "+sys.argv[2]+" "+file1+".csv temp").readlines()][0][0]

if float(f3)>=75:
	print("".join(ss), f3)
	os.system("python3.8 m6.py "+file1)
#	os.system("weblogo -X NO -Y NO --errorbars NO --format PDF --fineprint \"\"  -C \"#CB2026\" A A -C \"#34459C\" C C -C \"#FBB116\" G G -C \"#0C8040\" T T < temp.fa > "+file1+".pdf")
