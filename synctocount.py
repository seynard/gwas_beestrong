#module load system/Python-3.6.3
#python
import csv
import sys
import os
from functools import reduce
import glob
import numpy as np
import pandas as pd

f=open(sys.argv[1],"r")
#f=open('BS17-0561.count','r')
lines=f.read().splitlines()
d=[]
c_ref=[]
c_alt=[]
for N in range(1,len(lines)):
	Line=lines[N]
	info_ref=Line.split(' ')[2]
	info_alt=Line.split(' ')[3]
	count_A=int(Line.split(' ')[4].split(':')[0])
	count_T=int(Line.split(' ')[4].split(':')[1])
	count_C=int(Line.split(' ')[4].split(':')[2])	
	count_G=int(Line.split(' ')[4].split(':')[3])
	count_N=int(Line.split(' ')[4].split(':')[4])
	count_del=int(Line.split(' ')[4].split(':')[5])
	l_dict={'A':count_A,'T':count_T,'C':count_C,'G':count_G,'N':count_N,'del':count_del}
	if(sum(l_dict.values())>0 and count_del/sum(l_dict.values())>0.1):
		l_dict={k:0 for k in l_dict}
	elif(count_del>0 and count_del/sum(l_dict.values())=<0.1):
		l_dict['del']=0
	if(sum(l_dict.values())>0 and count_N/sum(l_dict.values())>0.1):
		l_dict={k:0 for k in l_dict}
	elif(count_N>0 and count_N/sum(l_dict.values())=<0.1):
		l_dict['N']=0
	if(len(info_alt)==1):
		allele=[info_ref,info_alt]
		nt_count=[k for k, v in l_dict.items() if v != 0]
		if(len(nt_count)==len(allele)):
			if(allele!=nt_count):
				not_in=[x for x in nt_count if x not in allele]
				for x in not_in:
					if(sum(l_dict.values())>0 and l_dict[x]/sum(l_dict.values())>0.1):
						l_dict={k:0 for k in l_dict}
					else:
						l_dict[x]=0
		elif(len(nt_count)>len(allele)):
			not_allele=list(set(nt_count)-set(allele))
			count_not_allele=sum(dict((k,l_dict[k]) for k in not_allele if k in l_dict).values())
			if(sum(l_dict.values())>0 and count_not_allele/sum(l_dict.values())>0.1):
				l_dict={k:0 for k in l_dict}
			else:	
				for x in not_allele:
					l_dict[x]=0
		elif(len(nt_count)<len(allele)):	
			not_in=[x for x in nt_count if x not in allele]	
			for x in not_in:
				if(sum(l_dict.values())>0 and l_dict[x]/sum(l_dict.values())>0.1):
					l_dict={k:0 for k in l_dict}
				else:
					l_dict[x]=0								
	elif(len(info_alt)>1):	
		info_alt1=info_alt.split(',')[0]
		info_alt2=info_alt.split(',')[1]
		info_alt=info_alt1+','+info_alt2
		allele=[info_ref,info_alt1,info_alt2]
		nt_count=[k for k, v in l_dict.items() if v != 0]
		if(len(nt_count)==len(allele)):
			if(allele!=nt_count):
				not_in=[x for x in nt_count if x not in allele]
				for x in not_in:
					if(sum(l_dict.values())>0 and l_dict[x]/sum(l_dict.values())>0.1):
						l_dict={k:0 for k in l_dict}
					else:
						l_dict[x]=0			
		elif(len(nt_count)>len(allele)):
			not_allele=list(set(nt_count)-set(allele))
			count_not_allele=sum(dict((k,l_dict[k]) for k in not_allele if k in l_dict).values())
			if(sum(l_dict.values())>0 and count_not_allele/sum(l_dict.values())>0.1):
				l_dict={k:0 for k in l_dict}
			else:	
				for x in not_allele:
					l_dict[x]=0
		elif(len(nt_count)<len(allele)):	
			not_in=[x for x in nt_count if x not in allele]	
			for x in not_in:
				if(sum(l_dict.values())>0 and l_dict[x]/sum(l_dict.values())>0.1):
					l_dict={k:0 for k in l_dict}
				else:
					l_dict[x]=0			
	depth=sum(l_dict.values())
	ref=l_dict[info_ref]
	alt=sum(list(dict((k, l_dict[k]) for k in info_alt if k in l_dict).values()))
	if(ref+alt!=depth):
		print('error in ',Line.split(' ')[0],Line.split(' ')[1])
	d.append(sum(l_dict.values()))
	c_ref.append(l_dict[info_ref])
	c_alt.append(list(dict((k, l_dict[k]) for k in info_alt if k in l_dict).values()))

np.savetxt(sys.argv[2]+'/tmp_depth_'+lines[0].split(' ')[4]+'.txt', np.c_[d], delimiter=' ', header=lines[0].split(' ')[4],fmt='%s',comments='')
np.savetxt(sys.argv[2]+'/tmp_count_ref_'+lines[0].split(' ')[4]+'.txt', np.c_[c_ref], delimiter=' ', header=lines[0].split(' ')[4],fmt='%s',comments='')
np.savetxt(sys.argv[2]+'/tmp_count_alt_'+lines[0].split(' ')[4]+'.txt', np.c_[c_alt], delimiter=' ', header=lines[0].split(' ')[4],fmt='%s',comments='')
