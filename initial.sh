pushd `dirname $0` > /dev/null
SCRIPTPATH=`pwd -P`
popd > /dev/null

STAR=$SCRIPTPATH/STARref
HISAT2=$SCRIPTPATH/HISAT2ref
BIT2FA=$SCRIPTPATH/bin/./twoBitToFa
mkdir $STAR
mkdir $HISAT2

cd $STAR && wget -O STARindex.tar.gz --no-check-certificate --no-proxy --timestamping 'https://www.encodeproject.org/files/ENCFF069ZCO/@@download/ENCFF069ZCO.tar.gz'
tar -xvzf STARindex.tar.gz
rm STARindex.tar.gz

cd $HISAT2 && wget --timestamping 'ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grch37.tar.gz' --no-check-certificate
tar -xvzf grch37.tar.gz
rm grch37.tar.gz

cd $SCRIPTPATH/lib/
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit
#wget ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz

gunzip Homo_sapiens.GRCh37.75.gtf.gz

chmod +x $BIT2FA
$BIT2FA hg19.2bit hg19.fa


