region=$1
offset=$(expr $2 - 1)
file=$3
genome=$4
partref=$5
chr=$6
pathDB=$7
pathData=$8

predb=$(echo $file | sed 's/.fa$//')

db=$pathDB/$predb
makeblastdb -in $partref -parse_seqids -dbtype nucl -out $db

prefix=$(echo $pathData/$file | sed 's/.fa$//')
prefixSM=$(echo $prefix | awk -F'/' '{ print $5 }' )
prefixID=$(echo $prefixSM | awk -F'.' '{ print $1 }')
regionID=$(echo $region | tr ":" "-")
prefixname=$(echo $file | sed 's/.fa$//')

blastn -query $pathData/$file -db $db -out $prefix.blast.xml -outfmt 5 -dust no -evalue 1e-100 -word_size 5
blastn -query $pathData/$file -db $db -out $prefix.blast.sam -outfmt '17 std SQ SR' -dust no -evalue 1e-100 -word_size 5

samtools view -h $prefix.goodname.blast.sam | samtools addreplacerg -r '@RG\tID:'$prefixID"-"$regionID'\tSM:'$prefixSM'' - | sed "s/$region/$chr/" | \
awk '{ if (substr($1,0,1)=="@") { print $0 } else { print $4+('$offset')"\t"$0 } }' | \
awk '{ if (substr($1,0,1)=="@") { print $0 } else { for (i=2;i<=NF;i++) { if (i==5) { printf "\t"$1 } else { printf "\t"$i } } print "" } }' | \
sed 's/^\t//' | sed 's/\(\t[0-9]\+\t\)[0-9]\+H/\1/' | \
sed 's/\([MDIN]\)[0-9]\+H\t/\1\t/' | samtools calmd -b - $genome | samtools sort - > $prefix.blast.bam
 
samtools index $prefix.blast.bam

rm $pathDB/*

