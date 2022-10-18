#module load system/Python-3.6.3
#python
import statistics
import sys
import bz2

maf_min=float(sys.argv[2])
#w=open('egs1.filter','w')
w=open(sys.argv[1]+'_'+str(maf_min)+'.txt.filter','w')
#maf_min=0.01

#with open('queen_geno1.egs') as f:
#with open(sys.argv[1]) as f:
with bz2.open(sys.argv[1]+'.txt.bz2','rt') as f:
	header_line=next(f)
	if ( 'freq' in sys.argv[1] ):
		s=header_line.strip('\n').replace(' ',',')
	else:
		s=header_line.strip('\n')	
	w.write(s+'\n')	
	for egs in f:
		if ( 'freq' in sys.argv[1] ):
			egs=egs.strip('\n').split(' ')
			chr=egs[0]
			pos=egs[1]
			ref=egs[2]
			alt=egs[3]
			egs_num=[float(i) for i in egs[4:len(egs)]]
		else:
			egs=egs.strip('\n').split(',')
			ind=egs[0]
			chr=egs[1]
			pos=egs[2]
			ref=egs[3]
			alt=egs[4]
			egs_num=[float(i) for i in egs[5:len(egs)]]
		F=statistics.mean(egs_num)
		F=round(F,3)
#		if (min(egs_num)==0):
#			egs_num=[F if i==0 else i for i in egs_num]
#			F=statistics.mean(egs_num)
#			F=round(F,3)		
		if (F > 0.5 ):
			MAF=1-F
		else:
			MAF=F
		if ( MAF > maf_min ):
			if ( 'freq' in sys.argv[1] ):		
				s=chr+','+pos+','+ref+','+alt+','+str(egs_num).strip('[').strip(']').replace("'",'')
			else:
				s=ind+','+chr+','+pos+','+ref+','+alt+','+str(egs_num).strip('[').strip(']').replace("'",'')
			w.write(s+'\n')
