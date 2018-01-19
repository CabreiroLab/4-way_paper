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

mkdir 1-HiSat_alignment
cd 1-HiSat_alignment
for i in $(cat ../$samples);do
	echo "1-HiSat alignment:" $i
	if [ "$7" == "" ]; then 
		hisat2 -p $proc --dta -x ../$index -U ../$root/*$i*.fastq.gz -S $i.sam
		samtools sort -@ $proc -o $i.bam $i.sam
		rm $i.sam
	elif [ "$7" == "ST" ]; then 
		hisat2 -p $proc --dta -x ../$index -U ../$root/*$i*.fastq.gz -S $i.sam
		samtools sort -@ $proc -o $i_unfiltered.bam $i.sam
		samtools view -bq 4 $i_unfiltered.bam > $i.bam
		rm $i_unfiltered.bam
		rm $i.sam
	elif [ "$7" == "HS" ]; then 
		hisat2 -p $proc -k 1 --dta -x ../$index -U ../$root/*$i*.fastq.gz -S $i.sam
		samtools sort -@ $proc -o $i.bam $i.sam
		rm $i.sam
	else
		echo "Unknown alignment filtering"
		exit 1
	fi
done
cd ..

#--met-file $i"_metrics"

mkdir 2-Stringtie
cd 2-Stringtie
for i in $(cat ../$samples);do
	echo "2-Stringtie quantification:" $i
	#echo ../1-HiSat_alignment/$i.bam
	stringtie -p $proc -G ../$annotation -o $i.gtf -l $i ../1-HiSat_alignment/$i.bam
	
done

echo "Collecting..."
ls *.gtf > mergelist.txt
echo "Merging..."
stringtie --merge -p $proc -G ../$annotation -o stringtie_merged.gtf mergelist.txt
echo "Comparing..."
gffcompare -r ../$annotation -G -o Comparison stringtie_merged.gtf
cd ..

mkdir 3-Ballgown
cd 3-Ballgown
for i in $(cat ../$samples);do
	echo "3-Ballgown parsing:" $i
	stringtie ../1-HiSat_alignment/$i.bam -G ../2-Stringtie/stringtie_merged.gtf -l $i -o $i/$i.gtf -p $proc -e -B
done
cd ..

prepDE.py -i 3-Ballgown -l 51
