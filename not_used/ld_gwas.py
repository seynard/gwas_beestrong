#module load system/Python-3.6.3
#python
import os
import sys
from ld_estimator import pairwise_ld
dir_in=sys.argv[1]
dir_out=sys.argv[2]
pop=sys.argv[3]
CHR=sys.argv[4]

w=open(dir_out+'/ld_'+pop+'_'+CHR+'.txt','w')
for f in open(dir_in+'/egs2_'+pop+'_'+CHR+'.txt'):
	x=f.strip('\n').split(' ')
	if ('BS' in x[4]):
		print('BS')
	else :	
		snpx=str(x[0])+':'+str(x[1])
		X=x[4:len(x)]
		X=[float(i) for i in X]
		for f in open(dir_in+'/egs2_'+pop+'_'+CHR+'.txt'):
			y=f.strip('\n').split(' ')
			if ('BS' in y[4]):
				print('BS')
			else :	
				snpy=str(y[0])+':'+str(y[1])
				Y=y[4:len(y)]
				Y=[float(i) for i in Y]
				is_haploid=[False]*len(X)
				Xround=[round(num) for num in X]
				Xround=[(0,0) if item == 0 else item for item in Xround]
				Xround=[(0,1) if item == 1 else item for item in Xround]
				Xround=[(1,1) if item == 2 else item for item in Xround]
				Yround=[round(num) for num in Y]
				Yround=[(0,0) if item == 0 else item for item in Yround]
				Yround=[(0,1) if item == 1 else item for item in Yround]
				Yround=[(1,1) if item == 2 else item for item in Yround]
				ld=pairwise_ld(Xround,Yround,is_haploid)
				if ((ld!=None and ld.dprime>0) or (ld!=None and ld.r_squared>0)):
					str_ld=snpx+' '+snpy+' '+str(ld.dprime)+' '+str(ld.r_squared)
					w.write(str_ld+'\n')
										

