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

w=open('average_depth_snp.txt','w')
with open('data/depth.txt') as f: 
	for df in f:
		df=df.strip('\n').split(' ')
		if(df[0]=='CHROM'):
			print('header')
			line='rs nb_col nb_col_depth0 nb_col_depth_non0 ave_depth min_depth max_depth sd_depth'
			w.write(line+'\n')
		else: 
			chr=df[0]
			pos=df[1]
			df=df[4:len(df)]
			df=[ int(x) for x in df ]
			df_len=len(df)
			nb_0=df.count(0)
			df_n0=df
			if nb_0 > 0 :
				while 0 in df_n0: df_n0.remove(0)
			df_n0_len=len(df_n0)
			line=str(chr)+':'+str(pos)+' '+str(df_len)+' '+str(nb_0)+' '+str(df_n0_len)+' '+str(round(statistics.mean(df_n0),3))+' '+str(min(df_n0))+' '+str(max(df_n0))+' '+' '+str(round(statistics.stdev(df_n0),3))
			w.write(line+'\n')
w.close()
