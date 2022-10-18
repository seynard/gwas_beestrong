#module load system/Python-3.6.3
#python
import csv
import sys
import os
from functools import reduce
import glob
import numpy as np
import pandas as pd
import statistics

w=open(sys.argv[3],'w')
#w=open('work/GWAS/bam_BeeStrongMOSAR/depth_control.txt','w')
with open(sys.argv[1]) as f1, open(sys.argv[2]) as f2: 
#with open('work/GWAS/bam_BeeStrongMOSAR/depth3.txt') as f1, open('work/GWAS/bam_BeeStrongMOSAR/count_ref3.txt') as f2: 
	for d,cr in zip(f1,f2):
		d=d.strip('\n').split(' ')
		cr=cr.strip('\n').split(' ')
		if ('BS' in d[4]):
			head=str(d).strip('[').strip(']').replace(',','').replace("'",'')
			s=head.split(' ')[0]+' '+head.split(' ')[1]+' '+'mean'+' '+'median'
			w.write(s+'\n')	
		else:
			dd=str(d).strip('[').strip(']').replace(',','').replace("'",'').split(' ')
			crr=str(cr).strip('[').strip(']').replace(',','').replace("'",'').split(' ')
			d_n=list(map(int,d[4:len(d)]))
			mean_d=round(sum(d_n)/len(d_n),3)
			median_d=statistics.median(d_n)
			s=d[0]+' '+d[1]+' '+str(mean_d)+' '+str(median_d)
			w.write(str(s)+'\n')
	w.close()

