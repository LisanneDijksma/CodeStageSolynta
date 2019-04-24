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
#to do gapopen en gapextension gebruiken 
#python module vinden die hsp om kan zetten naar goede sam regel
#Toevoegen van RG tag zodat je kunt onderscheiden in de gemergde bam file welke read bij welk genoom hoort

#Probleem is SR zorgt wel voor omdraaien subject en query, maar hierdoor wordt de refseq gebruikt ipv de subject sequentie
#misschien is het mogelijk beter handmatig omdraaien query en reference

blastn -query $pathData/$file -db $db -out $prefix.blast.xml -outfmt 5 -dust no -evalue 1e-100 -word_size 5
blastn -query $pathData/$file -db $db -out $prefix.blast.sam -outfmt '17 std SQ SR' -dust no -evalue 1e-100 -word_size 5


#samtools view -h $prefix.blast.sam | sed "s/$name/$prefix/" > $prefix.help.blast.sam

#name=$(samtools view -h $pathData/$prefix.blast.sam | awk ' { if (substr($1,0,2)!="@") { print $1 } } ')
./setQueryBamName.sh $prefixname $region $pathData $prefix.blast.sam

samtools view -h $prefix.goodname.blast.sam | samtools addreplacerg -r '@RG\tID:'$prefixID"-"$regionID'\tSM:'$prefixSM'' - | sed "s/$region/$chr/" | \
awk '{ if (substr($1,0,1)=="@") { print $0 } else { print $4+('$offset')"\t"$0 } }' | \
awk '{ if (substr($1,0,1)=="@") { print $0 } else { for (i=2;i<=NF;i++) { if (i==5) { printf "\t"$1 } else { printf "\t"$i } } print "" } }' | \
sed 's/^\t//' | sed 's/\(\t[0-9]\+\t\)[0-9]\+H/\1/' | \
sed 's/\([MDIN]\)[0-9]\+H\t/\1\t/' | samtools calmd -b - $genome | samtools sort - > $prefix.blast.bam
 
samtools index $prefix.blast.bam

rm $pathDB/*

