#### missing genotypes are 'imputed' as mean genotype value for the rest of the population #####
#### missing genotypes are frequency for which depth = count = 0 ####  
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

w=open(sys.argv[2],'w')
#w=open('work/GWAS/BeeStrongHAV3_1/freq_imputed.txt','w')
#for f in open('work/GWAS/BeeStrongHAV3_1/freq.txt'):
for f in open(sys.argv[1]):
	freq=f.strip('\n').split(' ')
	if ('BS' in freq[6]):
		head=str(freq).strip('[').strip(']').replace(',','').replace("'",'')
		w.write(str(head)+'\n')	
	else:
		Freq=np.array(freq)
		Freq=Freq[4:len(Freq)]
		p_miss=(Freq=='nan').sum()/len(Freq)
		if (p_miss<0.05):
			f_cor=[x for x in Freq if x != 'nan' ]
			f_cor=[float(item) for item in f_cor]
			m_f_cor=sum(f_cor)/len(f_cor)
			Freq2=[m_f_cor if x == 'nan' else x for x in Freq]
			freq_uniq=set(Freq2)
#			if (len(freq_uniq)==1):
#				print(freq[0]+' '+freq[1])
#			elif (max(Freq2)>1):
			Freq2=[float(item) for item in Freq2]
			if (max(Freq2)>1):
				print(freq[0]+' '+freq[1]+' still freq error')
			else:
				s=freq[0]+' '+freq[1]+' '+freq[2]+' '+freq[3]+' '+str(Freq2).strip('[').strip(']').replace(',','').replace("'",'')
				w.write(str(s)+'\n')

w.close()
