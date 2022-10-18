import csv
import sys
import os
import glob
import bz2
import statistics
import numpy as np

dir_in=sys.argv[1]
inp=sys.argv[2]
outp=sys.argv[3]
filter_maf=float(sys.argv[4])
filter_missing=float(sys.argv[5])

w=open(dir_in+'/'+outp,'w')
if ('bz2' in inp):
	for f in bz2.open(dir_in+'/'+inp, 'rt'):
		freq=f.strip('\n').split(' ')
		if (freq[0]==''):
			freq=freq[1:len(freq)]
		if ('BS' in freq[4]):
			head=str(freq).strip('[').strip(']').replace(',','').replace("'",'')
			w.write(head+'\n')
		elif (all(i is '' for i in freq[4:len(freq)])):
			print('empty')	
		else:
			freq[4:len(freq)]=[float(item) for item in freq[4:len(freq)]]
			F=np.nanmean(freq[4:len(freq)])
			F=round(F,3)
			if (F > 0.5 ):
				MAF=1-F
			else :
				MAF=F
			M=np.count_nonzero(np.isnan(freq[4:len(freq)]))
			M=M/len(freq[4:len(freq)])	
			if ( MAF > filter_maf and M < filter_missing):
				freq[4:len(freq)]=[element * 2 for element in freq[4:len(freq)]]
				f=str(freq).strip('[').strip(']').replace(',','').replace("'",'')
				w.write(f+'\n')	
else :
	for f in open(dir_in+'/'+inp, 'rt'):
		freq=f.strip('\n').split(' ')
		if (freq[0]==''):
			freq=freq[1:len(freq)]
		if ('BS' in freq[4]):
			head=str(freq).strip('[').strip(']').replace(',','').replace("'",'')
			w.write(head+'\n')	
		elif (all(i is '' for i in freq[4:len(freq)])):
			print('empty')	
		else:
			freq[4:len(freq)]=[float(item) for item in freq[4:len(freq)]]
			F=np.nanmean(freq[4:len(freq)])
			F=round(F,3)
			if (F > 0.5 ):
				MAF=1-F
			else :
				MAF=F
			M=np.count_nonzero(np.isnan(freq[4:len(freq)]))
			M=M/len(freq[4:len(freq)])	
			if ( MAF > filter_maf and M < filter_missing):
				freq[4:len(freq)]=[element * 2 for element in freq[4:len(freq)]]
				f=str(freq).strip('[').strip(']').replace(',','').replace("'",'')
				w.write(f+'\n')	
w.close()
