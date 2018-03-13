
head -n1 Filtered_Fusions.tsv > header

sed '1d' Filtered_Fusions.tsv | sort -k 4,4 -k 1,1 | awk '{print >> $4; close($4)}' -

for i in `ls CCRC*`
do
	cat header $i > Fusions_in_${i}.txt
done
rm -f CCRC*

for i in `ls UCEC*`
do
        cat header $i > Fusions_in_${i}.txt
done
rm -f UCEC*

