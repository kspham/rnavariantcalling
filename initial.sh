#!/bin/bash
pushd `dirname $0` > /dev/null
SCRIPTPATH=`pwd -P`
popd > /dev/null

STAR=$SCRIPTPATH/STARref
HISAT2=$SCRIPTPATH/HISAT2ref
BIT2FA=$SCRIPTPATH/bin/./twoBitToFa

if [ ! -d $STAR ]; then
mkdir $STAR
mkdir $HISAT2
fi

#if [ ! -d $STAR/hg19 ] || [ ! -f $STAR/hg19STARindex.tar.gz ]; then
if [ ! -d $STAR/hg19 ]; then
cd $STAR && wget -O hg19STARindex.tar.gz --no-check-certificate --no-proxy --timestamping 'https://www.encodeproject.org/files/ENCFF069ZCO/@@download/ENCFF069ZCO.tar.gz'
tar -xvzf hg19STARindex.tar.gz
rm STARindex.tar.gz
mv out hg19
fi

#if [ ! -d $STAR/mm10 ] || [ ! -f $STAR/mm10STARindex.tar.gz ]; then
if [ ! -d $STAR/mm10 ]; then
cd $STAR && wget -O mm10STARindex.tar.gz --no-check-certificate --no-proxy --timestamping 'https://www.encodeproject.org/files/ENCFF518RJA/@@download/ENCFF518RJA.tar.gz'
tar -xvzf mm10STARindex.tar.gz
rm mm10STARindex.tar.gz
mv out mm10
fi

if [ ! -d $HISAT2/hg19 ]; then
cd $HISAT2 && wget --timestamping 'ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/hg19.tar.gz' --no-check-certificate
tar -xvzf hg19.tar.gz
rm hg19.tar.gz
fi

if [ ! -d $HISAT2/mm10 ]; then
cd $HISAT2 && wget --timestamping 'ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/mm10.tar.gz' --no-check-certificate
tar -xvzf mm10.tar.gz
rm mm10.tar.gz
fi

if [ ! -f $SCRIPTPATH/lib/hg19.2bit ];then
cd $SCRIPTPATH/lib/
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit
chmod +x $BIT2FA
$BIT2FA hg19.2bit hg19.fa
fi

if [ ! -f "$SCRIPTPATH/lib/mm10.2bit" ];then
cd $SCRIPTPATH/lib/
wget http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.2bit
chmod +x $BIT2FA
$BIT2FA mm10.2bit mm10.fa
fi

if [ ! -f $SCRIPTPATH/lib/All_20151104.vcf.gz ];then
cd $SCRIPTPATH/lib/
wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b146_GRCh37p13/VCF/All_20151104.vcf.gz
tabix -p vcf All_20151104.vcf.gz
fi

if [ ! -f $SCRIPTPATH/lib/Mus_musculus.vcf.gz ];then
cd $SCRIPTPATH/lib/
wget ftp://ftp.ensembl.org/pub/release-84/variation/vcf/mus_musculus/Mus_musculus.vcf.gz
tabix -p vcf Mus_musculus.vcf.gz 
fi

if [ ! -d $SCRIPTPATH/bin/snpEff/data ]; then
cd $SCRIPTPATH/bin/snpEff
wget http://downloads.sourceforge.net/project/snpeff/databases/v4_2/snpEff_v4_2_GRCm38.82.zip
unzip snpEff_v4_2_GRCm38.82.zip
rm snpEff_v4_2_GRCm38.82.zip
wget http://downloads.sourceforge.net/project/snpeff/databases/v4_2/snpEff_v4_2_GRCh37.75.zip
unzip snpEff_v4_2_GRCh37.75.zip
rm snpEff_v4_2_GRCh37.75.zip
fi
