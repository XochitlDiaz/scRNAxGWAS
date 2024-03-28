annot_cell=$1
LDSC_PATH="/working/lab_miguelr/xochitlY/scRNAxGWAS/LDSCORE"
baseline_cell=$LDSC_PATH/baselineLD_v2.1_annots
weights_path=$LDSC_PATH/weights_hm3_no_hla
freq_path=$LDSC_PATH/1000G_Phase3_frq
sumstats_cell=$LDSC_PATH/SUMSTATS
output_cell_pre=$2

module load ldsc/1.0.0

IFS="
"

sumstats_path=$LDSC_PATH/SUMSTATS


if [ ! -d $output_cell_pre ]
then
    mkdir $output_cell_pre
fi

output_cell=$output_cell_pre

if [ ! -d $output_cell ]
then
    mkdir $output_cell
fi

echo $output_cell
for line in `ls $annot_cell | awk '{print $1}' | sort | uniq`;
do
    annot_module=`echo $line | awk '{print $1}'`
    echo $annot_cell $annot_module
    if [ ! -d $annot_cell/$annot_module ]
    then
        echo "Error: annotation module directory not found" > ldsc_logfile.log
        exit 100
    fi
    if [ ! -d $output_cell/$annot_module ]
    then
        mkdir $output_cell/$annot_module
    fi
    for ll in `ls $annot_cell/$annot_module | awk '{print $1}' | sort | uniq`;
    do
	annot_dir=`echo $ll | awk '{print $1}'`
	echo $annot_dir
	if [ ! -d $annot_cell/$annot_module/$annot_dir ]
	then
            echo "Error: annotation module directory not found" > ldsc_logfile.log
            exit 101
	fi
	if [ ! -d $output_cell/$annot_module/$annot_dir ];then
            mkdir $output_cell/$annot_module/$annot_dir
	fi
	for step in `ls $sumstats_cell | awk '{print $1}' | sort | uniq`;
	do
            sumstats_file=`echo $step | awk '{print $1}'`
            echo $sumstats_cell $sumstats_file
            if [ ! -f $sumstats_cell/$sumstats_file ];then
		echo "Error: sumstats file not found" > ldsc_logfile.log
		exit 102
            fi
            if [ ! -f $output_cell/$annot_module/$annot_dir/$sumstats_file.results ];then
            cmd="ldsc.py  --h2 $sumstats_cell/$sumstats_file --ref-ld-chr $annot_cell/$annot_module/$annot_dir/$annot_dir.,$baseline_cell/baselineLD.  --frqfile-chr $freq_path/1000G.EUR.QC. --w-ld-chr $weights_path/weights.hm3_noMHC. --overlap-annot --print-coefficients --print-delete-vals --out $output_cell/$annot_module/$annot_dir/$sumstats_file"
            sbatch --time=120:00 --mem=20000 --output=reg_max.out --error=reg_max.err -p short -c 1 --wrap="$cmd"
            fi
				done
     done
done
