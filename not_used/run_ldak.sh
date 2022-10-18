#! /bin/bash
############################################ parameters ################################################## 
j=${1}
geno=${2}
########################################################################################################## 
./ldak5.1.linux --cut-weights sections$j --bfile $geno --chr $j
./ldak5.1.linux --calc-weights-all sections$j --bfile $geno --chr $j
