#module load system/Python-3.6.3
#python
import csv
import sys
import os
from functools import reduce
import glob
import numpy as np
import pandas as pd

if(sys.argv[3]=='remove'):
	w=open(sys.argv[2],'w')
	#w=open('triallelic_remove.txt','w')
	with open(sys.argv[1]) as f: 
	#with open('alt2.txt') as f: 
		for a in f:
			ai=a.strip('\n').split(' ')
			a=ai[5:len(ai)]
			a1=[]
			a2=[]
			for i in range(0,len(a)):
				a1.append(int(a[i].strip('[').strip(']').split(',')[0]))
				a2.append(int(a[i].strip('[').strip(']').split(',')[1]))
			n_a1=len(a1)-a1.count(0)
			n_a2=len(a2)-a2.count(0)
			if(n_a1==n_a2):
				s=ai[0]+' '+ai[1]
				w.write(s+'\n')
			if(n_a1>n_a2 and n_a2<0.1*len(a)):
				s=ai[0]+' '+ai[1]+' '+ai[2]+' '+ai[3].split(',')[0]
			elif(n_a1<n_a2 and n_a1<0.1*len(a)):
				s=ai[0]+' '+ai[1]+' '+ai[2]+' '+ai[3].split(',')[1]
			elif(n_a1>0.1*len(a) and n_a2>0.1*len(a)):
				s=ai[0]+' '+ai[1]
				w.write(s+'\n')

elif(sys.argv[3]=='recode'):	
	w=open(sys.argv[2],'w')
	#w=open('triallelic_recode.txt','w')
	with open(sys.argv[1]) as f: 
	#with open('alt2.txt') as f: 
		for a in f:
			ai=a.strip('\n').split(' ')
			a=ai[5:len(ai)]
			a1=[]
			a2=[]
			for i in range(0,len(a)):
				a1.append(int(a[i].strip('[').strip(']').split(',')[0]))
				a2.append(int(a[i].strip('[').strip(']').split(',')[1]))
			n_a1=len(a1)-a1.count(0)
			n_a2=len(a2)-a2.count(0)
			if(n_a1==n_a2):
				s=ai[0]+' '+ai[1]
			if(n_a1>n_a2 and n_a2<0.1*len(a)):
				s=ai[0]+' '+ai[1]+' '+ai[2]+' '+ai[3].split(',')[0]
				w.write(s+'\n')
			elif(n_a1<n_a2 and n_a1<0.1*len(a)):
				s=ai[0]+' '+ai[1]+' '+ai[2]+' '+ai[3].split(',')[1]
				w.write(s+'\n')
			elif(n_a1>0.1*len(a) and n_a2>0.1*len(a)):
				s=ai[0]+' '+ai[1]
