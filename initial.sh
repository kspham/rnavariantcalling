#!/bin/bash
pushd `dirname $0` > /dev/null
SCRIPTPATH=`pwd -P`
popd > /dev/null

STAR=$SCRIPTPATH/STARref
HISAT2=$SCRIPTPATH/HISAT2ref
BIT2FA=$SCRIPTPATH/bin/./twoBitToFa

if [ ! -d "$STAR" ]; then
mkdir $STAR
mkdir $HISAT2
fi

if [ ! -d "$STAR/out" ]; then
cd $STAR && wget -O STARindex.tar.gz --no-check-certificate --no-proxy --timestamping 'https://www.encodeproject.org/files/ENCFF069ZCO/@@download/ENCFF069ZCO.tar.gz'
tar -xvzf STARindex.tar.gz
rm STARindex.tar.gz
fi

if [ ! -d "$HISAT2/grch37" ]; then
cd $HISAT2 && wget --timestamping 'ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grch37.tar.gz' --no-check-certificate
tar -xvzf grch37.tar.gz
rm grch37.tar.gz
fi

if [ ! -f "$SCRIPTPATH/lib/hg19.fa" ];then
cd $SCRIPTPATH/lib/
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit
chmod +x $BIT2FA
$BIT2FA hg19.2bit hg19.fa
fi

if [ ! -f "$SCRIPTPATH/lib/All_20151104.vcf.gz" ]; then
cd $SCRIPTPATH/lib
if wget 'ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/All_20151104.vcf.gz'; then
echo "VCF database downloaded";
else
rm All_20151104.vcf.gz
echo "Download error! Please run initial again";
fi
fi

if [ ! -f "$SCRIPTPATH/lib/All_20151104.vcf.gz.tbi" ]; then
cd $SCRIPTPATH/lib
wget 'ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/All_20151104.vcf.gz.tbi';
fi

if [ ! -f "$SCRIPTPATH/bin/java64.tar.gz" ]; then
cd $SCRIPTPATH/bin
wget -O java64.tar.gz 'http://javadl.sun.com/webapps/download/AutoDL?BundleId=111741'
tar -xvzf java64.tar.gz
fi
