#! /usr/bin/env bash

# create corresponding table for GRCh name and UCSC name
rm -rf GCF_000001405.13.assembly.txt
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/GCF_000001405.13.assembly.txt
python make_ucsc_grch.py GCF_000001405.13.assembly.txt > grch2ucsc.txt

# for GRCh38 (hg38)
# ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/GCF_000001405.31.assembly.txt

# for GRCm38 (mm10)
# wget http://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/refGene.txt.gz
# wget http://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/ensGene.txt.gz

rm -rf refGene.txt.gz
rm -rf ensGene.txt.gz
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/ensGene.txt.gz

# branchpoint (Mercer et al., Genome Research, 2015)
wget http://genome.cshlp.org/content/suppl/2014/12/16/gr.182899.114.DC1/Supplemental_DataS2.bed.gz

echo "python listGene.py refGene.txt.gz ref | sort -k1,1 -k2,2n -k3,3n - > refGene.bed"
python listGene.py refGene.txt.gz ref | sort -k1,1 -k2,2n -k3,3n - > refGene.bed

echo "python listGene.py ensGene.txt.gz ens | sort -k1,1 -k2,2n -k3,3n - > ensGene.bed"
python listGene.py ensGene.txt.gz ens | sort -k1,1 -k2,2n -k3,3n - > ensGene.bed

echo "python listExon.py refGene.txt.gz ref | sort -k1,1 -k2,2n -k3,3n - > refExon.bed"
python listExon.py refGene.txt.gz ref | sort -k1,1 -k2,2n -k3,3n - > refExon.bed

echo "python listExon.py ensGene.txt.gz ens | sort -k1,1 -k2,2n -k3,3n - > ensExon.bed"
python listExon.py ensGene.txt.gz ens | sort -k1,1 -k2,2n -k3,3n - > ensExon.bed

echo "python listCoding.py refGene.txt.gz | sort -k1,1 -k2,2n -k3,3n > refCoding.bed"
python listCoding.py refGene.txt.gz | sort -k1,1 -k2,2n -k3,3n > refCoding.bed

# compressing  and indexing
echo "bgzip -f refGene.bed"
bgzip -f refGene.bed

echo "bgzip -f ensGene.bed"
bgzip -f ensGene.bed

echo "bgzip -f refExon.bed"
bgzip -f refExon.bed

echo "bgzip -f ensExon.bed"
bgzip -f ensExon.bed

echo "bgzip -f refCoding.bed"
bgzip -f refCoding.bed

echo "tabix -f -p bed refGene.bed.gz"
tabix -f -p bed refGene.bed.gz 

echo "tabix -f -p bed ensGene.bed.gz"
tabix -f -p bed ensGene.bed.gz 

echo "tabix -f -p bed refExon.bed.gz"
tabix -f -p bed refExon.bed.gz 

echo "tabix -f -p bed ensExon.bed.gz"
tabix -f -p bed ensExon.bed.gz 

echo "tabix -f -p bed refCoding.bed.gz"
tabix -f -p bed refCoding.bed.gz

# branchpoint
sort -k1,1 -k2,2n Supplemental_DataS2.bed | bgzip -c > branchpoint_mercer.bed.gz
tabix -p bed branchpoint_mercer.bed.gz

