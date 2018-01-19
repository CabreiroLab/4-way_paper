#!/bin/bash

proc=16

mkdir 3-Ballgown-Met
cd 3-Ballgown-Met
echo "Metformin dataset...."
for i in $(cat ../../Samples.txt);do
	echo "3-Ballgown parsing:" $i
	stringtie ../../Against_BristolN_EN-WB235/1-HiSat_alignment/$i.bam -G ../2-Stringtie/stringtie_merged.gtf -l $i -o $i/$i.gtf -p $proc -e -B
done
cd ..

mkdir 3-Ballgown-DR
cd 3-Ballgown-DR

echo "DR dataset...."
for i in $(cat ../../Heinz_DR/Samples.txt);do
	echo "3-Ballgown parsing:" $i
	stringtie ../1-HiSat_alignment/$i.bam -G ../2-Stringtie/stringtie_merged.gtf -l $i -o $i/$i.gtf -p $proc -e -B
done

cd ..


prepDE.py -i 3-Ballgown-Met -l 51 -g Gene_count_Met.csv -t Transcript_count_Met.csv

prepDE.py -i 3-Ballgown-DR -l 100 -g Gene_count_DR.csv -t Transcript_count_DR.csv


#stringtie ../1-HiSat_alignment/$i.bam -G ../2-Stringtie/stringtie_merged.gtf -l $i -o $i/$i.gtf -p 16 -e -B
