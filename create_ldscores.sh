#!/bin/bash
programs_dir=$1
annot_cell=$2
LDSC_PATH="/working/lab_miguelr/xochitlY/scRNAxGWAS/LDSCORE"
bfile_path=$LDSC_PATH/Phase3_1000G
hapmap_path=$LDSC_PATH/hapmap3_snps

module load ldsc/1.0.0

IFS="
"

for line in `ls $annot_cell | awk '{print $1}' | sort | uniq`;
do
   annot_module=`echo $line | awk '{print $1}'`
   echo $annot_cell $annot_module
   for enhancer_method in `ls $annot_cell/$annot_module | awk '{print $1}' | sort | uniq`;
   do
      annot_dir=`echo $enhancer_method | awk '{print $1}'`
      echo $annot_dir
      if [ ! -d $annot_cell/$annot_module/$annot_dir ]
      then
				mkdir $annot_cell/$annot_module/$annot_dir
			fi
			for chrom in {1..22}
      do
			  if [ ! -f $annot_cell/$annot_module/$annot_dir/$annot_dir.$chrom.l2.ldscore.gz ]
	      then
	        ldsc.py --bfile $bfile_path/1000G.EUR.QC.$chrom --l2 --ld-wind-cm 1 --yes-really --annot $annot_cell/$annot_module/$annot_dir/$annot_dir.$chrom.annot.gz --print-snps $hapmap_path/hm.$chrom.snp --out $annot_cell/$annot_module/$annot_dir/$annot_dir.$chrom
	      fi
			done
   done
done
