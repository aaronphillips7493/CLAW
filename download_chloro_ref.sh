declare -a arr=("NC_041421.1")

for i in "${arr[@]}"
do
	ncbi-acc-download --format fasta $i
done

mkdir -p chloro_assembly/reference

touch chloro_table_3.fasta
for j in *.fa
do
	echo ${j}
	grep '^>' ${j} >> chloro_table_3.fasta
	grep -v '^>' ${j} > temp
	tr -d '\n' < temp > temp2
	cat temp2 >> chloro_table_3.fasta
	cat temp2 >> chloro_table_3.fasta
	sed -i -e '$a\' chloro_table_3.fasta
done

mv chloro_table_3.fasta chloro_assembly/reference/
