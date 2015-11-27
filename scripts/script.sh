VCFTOOLS=$1
HISAT2=$2
STAR=$3
EDITINGSITES=$4
OUTDIR=$5

for f in $HISAT2/*.vcf; do
	bgzip -c $f > $f.gz
	tabix -p vcf $f.gz
done

for f in $STAR/*.vcf; do
	bgzip -c $f > $f.gz
	tabix -p vcf $f.gz
done

for f in $STAR/*.vcf.gz; do
	$VCFTOOLS --gzvcf $f --gzdiff $HISAT2/"${f##*/}" --diff-site --out $f
	awk '{for (i=1; i<=NF; i++) {if ($i==".") next}; if (NF) print}' $f.diff.sites_in_files | awk '{print $1"\t"$2}' | sed '1d' > $f.common
done

for f in $STAR/*.vcf.gz; do
	
	$VCFTOOLS --gzvcf $f --positions $f.common --recode --recode-INFO-all --out "${f%.*}"
done

mv $STAR/*recode* $OUTDIR

rm $STAR/*.vcf.*
rm $HISAT2/*.vcf.*

