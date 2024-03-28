#bedops_cell=[BEDOPS_DIRECTORY]/BEDOPS/bin
#bedtools_cell=[BEDTOOLS_DIRECTORY]/BEDTOOLS/bedtools2/bin
input_cell=$1  #sclinker_beds/{Celltype} [directory where bed files are located]
names=`ls $input_cell | cut -f 1 -d '.'`

module load bedtools/2.29.0
module load bedops/2.4.41

for name in $names
do
bedtools sort -i $input_cell/$name.bed > $input_cell/$name.2.bed
bedtools merge -i $input_cell/$name.2.bed -c 4 -o max > $input_cell/$name.3.bed
mv $input_cell/$name.3.bed $input_cell/$name.bed
rm $input_cell/$name.2.bed
done
