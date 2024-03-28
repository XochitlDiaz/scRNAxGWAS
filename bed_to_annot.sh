#!/bin/bash
module load bedtools/2.29.0
module load gcc/6.4.0
module load python/3.11.4

bed_dir=$1
bimfile_path=$2
geneset=$3
annot_path="$4/${geneset}_annots"
FOLDER="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"


IFS="
"
	
if [ ! -d $annot_path/$geneset ];then
	if [ ! -d $annot_path ];then
		mkdir $annot_path
	fi
	mkdir $annot_path/$geneset
fi

for bedline in `ls $bed_dir/$geneset/ | cat | sort | uniq | cut -f 1 -d '.'`;
do
	bedname=`echo $bedline | awk '{print $1}'`
	if [ ! -d $annot_path/$geneset/$bedname ];then
	    mkdir $annot_path/$geneset/$bedname
	fi
	if [ ! -f $annot_path/$geneset/$bedname/$bedname.22.annot.gz ];then
	    python  $FOLDER/bedgraph_to_annot.py --bedname $bedname --bedfile_path $bed_dir/$geneset --bimfile_path $bimfile_path --annot_path $annot_path/$geneset/$bedname
	fi
done
