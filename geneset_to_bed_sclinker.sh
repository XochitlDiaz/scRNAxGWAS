geneset_dir=$1
bed_dir=$2


module load gcc/6.4.0
module load R/4.2.0

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$DIR"



IFS="
"

for ll in `ls -1 $geneset_dir | sed 's/\.txt//g' | awk '{print $1}' | sort | uniq`;
do
      annot_name=`echo $ll | awk '{print $1}'`
      echo $annot_name
      Rscript geneset_to_bed_sclinker.R  $geneset_dir $bed_dir $annot_name
      #cmd="Rscript geneset_to_bed_sclinker.R  $geneset_dir $bed_dir $annot_name"
      #sbatch --time=90:00 --mem=20000 --output=geneset_to_bed.out --error=geneset_to_bed.err -p short -c 1 --wrap="$cmd"
done
