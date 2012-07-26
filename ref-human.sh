#!/bin/sh
##
echo `date` ">> " "setting the settings.."
#Please set reference genome path here.
#refPATH=/Please/set/reference/genome/path/here/
refPATH=`pwd`
export PATH=$PATH:$refPATH
cd $refPATH
###############################################################################
# for human
###############################################################################
genome_ucsc="http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz"
echo `date` ">> " "downloading human reference genome from UCSC"
echo `date` ">> " "curl -O $genome_ucsc"
curl -O $genome_ucsc
#Gunzip the tarball
tar -zxvf chromFa.tar.gz
#Remove the haplotype, unmapped, and random chromosomes
rm *random*
rm *Un*
rm *hap*
#Concatenate the different chromosomes into one chromosome.
cat chr*.fa >hg19.fa
# build the indexes for bwa
bwa index -a bwtsw hg19.fa
echo `date` ">> " "done"