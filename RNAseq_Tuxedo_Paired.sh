#!/bin/bash

root=$1
samples=$2
reference=$3
annotation=$4
index=$5
proc=$6

#if [ "$7" != "" ]; then
#	uniqueness=$7
#else 
#	uniqueness
#fi

start=$PWD

#mkdir 1-HiSat_alignment
#mkdir 1-HiSat_alignment/NonConcordant
#cd 1-HiSat_alignment
#for i in $(cat ../$samples);do
#	echo "1-HiSat alignment:" $i
#	hisat2 -p $proc --no-discordant --no-mixed --met-file $i_metrics.txt --un-conc-gz NonConcordant/$i --dta -x ../$index -1 ../$root/*$i*_1.fq.gz -2 ../$root/*$i*_2.fq.gz -S $i.sam
#	samtools sort -@ $proc -o $i.bam $i.sam
#	rm $i.sam
#done
#cd ..

#--met-file $i"_metrics"

mkdir 2-Stringtie
cd 2-Stringtie
for i in $(cat ../$samples);do
	echo "2-Stringtie quantification:" $i
	#echo ../1-HiSat_alignment/$i.bam
	stringtie -p $proc -G ../$annotation -o $i.gtf -l $i ../1-HiSat_alignment/$i.bam
done

#echo "Collecting..."
#for i in $(cat ../$samples);do
#	echo $i.gtf >> mergelist.txt
#done
#ls *.gtf > mergelist.txt
#for i in $(cat ../../Against_BristolN_EN-WB235/2-Stringtie/mergelist.txt);do
#	echo ../../Against_BristolN_EN-WB235/2-Stringtie/$i >> mergelist_Metformin.txt
#done
#cp mergelist.txt mergelist_Combined.txt
#cat mergelist_Metformin.txt >> mergelist_Combined.txt

echo "Merging..."
#stringtie --merge -p $proc -G ../$annotation -o stringtie_merged.gtf mergelist.txt

stringtie --merge -p $proc -G ../$annotation -o stringtie_merged_DR_only.gtf mergelist.txt


echo "Comparing..."
gffcompare -r ../$annotation -G -o Comparison stringtie_merged_DR_only.gtf
cd ..

mkdir 3-Ballgown-DRonly
cd 3-Ballgown-DRonly
for i in $(cat ../../Heinz_DR/Samples.txt);do
	echo "3-Ballgown parsing:" $i
	stringtie ../1-HiSat_alignment/$i.bam -G ../2-Stringtie/stringtie_merged_DR_only.gtf -l $i -o $i/$i.gtf -p $proc -e -B
done
cd ..

#mkdir 3-Ballgown-Met
#cd 3-Ballgown-Met
#for i in $(cat ../../Samples.txt);do
#	echo "3-Ballgown parsing:" $i
#	stringtie ../../Against_BristolN_EN-WB235/1-HiSat_alignment/$i.bam -G ../2-Stringtie/stringtie_merged.gtf -l $i -o $i/$i.gtf -p $proc -e -B
#done
#cd ..

#prepDE.py -i 3-Ballgown-Met -l 51 -g Gene_count_Met.csv -t Transcript_count_Met.csv

prepDE.py -i 3-Ballgown-DRonly -l 100 -g Gene_count_DRonly.csv -t Transcript_count_DRonly.csv
