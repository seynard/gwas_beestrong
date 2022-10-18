import numpy as np
import pandas as pd
import sys
prfx=sys.argv[1]
wped=open(sys.argv[2]+'.tped','w')

with open(sys.argv[1],'r') as f: 
	for l in f:
		l=l.strip('\n').split(',')
		if(l[0]=='CHROM:POS'):
			id=l[3:len(l)]
			id=[i.replace(':','_') for i in id]
			d={'fid':1,'iid':id,'father':0,'mother':0,'sex':0,'pheno':0}
			df=pd.DataFrame(data=d)
			df.to_csv(sys.argv[2]+'.tfam',header=False,index=False,sep=' ',mode='w')
		else:		
			CHROM=l[0].split(':')[0]
			name=l[0]
			cm=0
			POS=l[0].split(':')[1]
			REF=l[1]
			ALT=l[2]
			l=l[3:len(l)]
			G=[]
			for x in l:
				if(x=='2'):
					G.append(ALT+' '+ALT)
				elif(x=='1'):
					G.append(REF+' '+ALT)
				elif(x=='0'):
					G.append(REF+' '+REF)
				else:
					G.append(0+' '+0)
			G=' '.join(G)
			s=str(CHROM)+' '+name+' '+str(cm)+' '+str(POS)+' '+str(G)
			wped.write(str(s)+'\n')		

wped.close()

