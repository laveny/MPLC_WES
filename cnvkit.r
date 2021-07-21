
cd ~/miniconda2/bin/
source ./activate
source activate cnvkit



Tumor_bam=~/MPLC/bam/tumor.bam
Normal_bam=~/MPLC/bam/normal.bam
outdir=~/MPLC/cnvkit/${patient}/
snv_file=~/MPLC



if [ ! -d $outdir ]
then mkdir -p $outdir
fi

if [ ! -d $outdir/call ]
then mkdir -p $outdir/call
fi

cd $outdir

for i in $Tumor_bam
do
fn=`basename $i`
sample=${fn%%.bam}


cnvkit.py batch $i --normal $Normal_bam \
    -m hybrid \
    --targets ex_region.sort.bed --annotate hg38.refFlat.txt \
    --fasta Homo_sapiens_assembly38.fasta  \
    --output-reference $outdir/${sample}_reference.cnn --output-dir $outdir \
   --diagram --scatter
   
cnvkit.py call ${sample}.cns -v ${snv_file}/${sample}.vcftools_selected.vcf -m threshold -t=-1.1,-0.4,0.3,0.7 -o $outdir/call/${sample}.call.cns
