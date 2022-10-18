#module load system/Python-3.6.3
#python
import csv
import sys
import os
from functools import reduce
import glob
import numpy as np
import pandas as pd

f=open(sys.argv[3],'w')
#f=open('depth3.txt','w')

if('depth' in sys.argv[3]):
	with open(sys.argv[1]) as f1, open(sys.argv[2]) as f2: 
	#with open('depth2.txt') as f1, open('count_ref2.txt') as f2: 
		for d,cr in zip(f1,f2):
			d=d.strip('\n').split(' ')
			cr=cr.strip('\n').split(' ')
			if ('BS' in d[4]):
				head=str(d).strip('[').strip(']').replace(',','').replace("'",'')
				f.write(head+'\n')
			elif ('BS' not in d):
				dd=' '.join(str(u) for u in d)
				crr=' '.join(str(u) for u in cr)
				ds=list(map(int,d[4:len(d)]))
				crs=list(map(int,cr[4:len(cr)]))
				if(ds==crs):
					print(d[0],' ',d[1],' monoallelic')
				else:
					f.write(str(dd)+'\n')
elif('ref' in sys.argv[3]):				
	with open(sys.argv[1]) as f1, open(sys.argv[2]) as f2: 
	#with open('depth2.txt') as f1, open('count_ref2.txt') as f2: 
		for d,cr in zip(f1,f2):
			d=d.strip('\n').split(' ')
			cr=cr.strip('\n').split(' ')
			if ('BS' in d[4]):
				head=str(d).strip('[').strip(']').replace(',','').replace("'",'')
				f.write(head+'\n')
			elif ('BS' not in d):
				dd=' '.join(str(u) for u in d)
				crr=' '.join(str(u) for u in cr)
				ds=list(map(int,d[4:len(d)]))
				crs=list(map(int,cr[4:len(cr)]))
				if(ds==crs):
					print(cr[0],' ',cr[1],' monoallelic')
				else:
					f.write(str(crr)+'\n')
