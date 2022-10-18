#cd /work/genphyse/dynagen/seynard/GWAS/infestation/
#module load system/Python-3.7.4
#python
import csv
import sys
import os
from functools import reduce
import glob
import numpy as np
import pandas as pd
import statistics
import dask.dataframe as dd
import collections

w=open('average_depth_colony.txt','w')
with open('data/data_extra/depth_all_nus.txt') as f: 
	for df in f:
		df=df.strip('\n').split(' ')
		if(df[0]=='CHROM'):
			print('header')
			line='col nb_col nb_col_depth0 nb_col_depth_non0 ave_depth min_depth max_depth sd_depth'
			w.write(line+'\n')
			break

head=df[4:len(df)]		
for x in range(0,len(head)):
	df=pd.read_csv('data/data_extra/depth_all_nus.txt',usecols=head[x:(x+1)],header=0,sep=' ')
	df=df[head[x]].tolist()
	df_len=len(df)
	col_id=head[x]
	nb_0=df.count(0)
	df_n0=df
	if nb_0 > 0 :
		while 0 in df_n0: df_n0.remove(0)
	df_n0_len=len(df_n0)
	line=str(col_id)+' '+str(df_len)+' '+str(nb_0)+' '+str(df_n0_len)+' '+str(round(statistics.mean(df_n0),3))+' '+str(min(df_n0))+' '+str(max(df_n0))+' '+' '+str(round(statistics.stdev(df_n0),3))
	w.write(line+'\n')

w.close()

head=df[4:len(df)]		
x=1
df=dd.read_csv('data/data_extra/depth_all_nus.txt',usecols=head[x:(x+1)],header=0,sep=' ')
#df=df[head[x]].tolist()
df_len=len(df)
col_id=head[x]
nb_0=df.count(0)
df_n0=df
if nb_0 > 0 :
	while 0 in df_n0: df_n0.remove(0)
df_n0_len=len(df_n0)
line=str(col_id)+' '+str(df_len)+' '+str(nb_0)+' '+str(df_n0_len)+' '+str(round(statistics.mean(df_n0),3))+' '+str(min(df_n0))+' '+str(max(df_n0))+' '+' '+str(round(statistics.stdev(df_n0),3))
w.write(line+'\n')

