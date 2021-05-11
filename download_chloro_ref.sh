declare -a arr=("NC_041421.1")

for i in "${arr[@]}"
do
	ncbi-acc-download --format fasta $i
done

touch chloro_table.fasta
for j in *.fa
do
	echo ${j}
	grep '^>' ${j} >> chloro_table.fasta
	grep -v '^>' ${j} > temp
	tr -d '\n' < temp > temp2
	cat temp2 >> chloro_table.fasta
	cat temp2 >> chloro_table.fasta
	sed -i -e '$a\' chloro_table.fasta
done
