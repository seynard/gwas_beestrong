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
import bz2
from bz2 import BZ2File as bzopen
from decimal import *
setcontext(ExtendedContext)
    
w=open(sys.argv[3],'w')
#w=open('work/GWAS/bam_BeeStrongMOSAR/freq.txt','w')
with open(sys.argv[1]) as f1, open(sys.argv[2]) as f2: 
#with open('depth_final.txt') as f1, open('count_ref_final.txt') as f2: 
	for d,cr in zip(f1,f2):
#		d=d.decode("utf-8")
#		cr=cr.decode("utf-8")		
		d=d.strip('\n').split(' ')
		cr=cr.strip('\n').split(' ')
		if (len(d[3])>1 and 'BS' in d[6]):
			head=str(d).strip('[').strip(']').replace(',','').replace("'",'')
			w.write(head+'\n')	
		else:
			d_n=list(map(float,d[4:len(d)]))
			cr_n=list(map(float,cr[4:len(cr)]))					
#			d_n2=['NA' if x==0 else x for x in d_n]
#			cr_n2=[-999 if x==0 else x for x in cr_n]
#			freq=['NA' if d_n2 == 'NA' else x/y for x,y in zip(cr_n,d_n2)]
			freq=[Decimal(x)/Decimal(y) for x,y in zip(cr_n,d_n)]
			freq_round=[ '%.3f' % elem for elem in freq]
#			freq_uniq=set(freq_round)
#			if (len(freq_uniq)==1):
#				print(d[0]+' '+d[1])
#			else:
			s=d[0]+' '+d[1]+' '+d[2]+' '+d[3]+' '+str(freq_round).strip('[').strip(']').replace(',','').replace("'",'')
			w.write(str(s)+'\n')
	w.close()
	
