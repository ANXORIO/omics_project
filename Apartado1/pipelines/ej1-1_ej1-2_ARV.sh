alias fastqc='perl /home/vant/FastQC/fastqc'

#################### CONTROL DE CALIDAD CON FASTQC #################### 
cd /home/vant/transcriptomics-project/Apartado1
mkdir -p ./fastQC

for sid in $(ls input/*.fastq | cut -d "." -f1 | sed "s:input/::" | sort | uniq)
do
    fastqc -o ./fastQC --noextract input/$sid.chr21.fastq
done
#################### ALINEAMIENTO DE LAS LECTURAS CON HISAT2 ####################
export PATH=/home/vant/miniconda3/pkgs/hisat2-2.2.1-h87f3376_4:$PATH

# si no se ha hecho el indexado del genoma de referencia:
hisat2-build --seed 123 -p 2 input/REF/chr21.fa input/REF/chr21

cd input
mkdir -p hisat2 

for sid in $(ls *.fastq | cut -d "." -f1)
do
	hisat2 --new-summary --summary-file hisat2/$sid.hisat2.summary  --rna-strandness R --seed 123 --phred33 -p 2 -k 1 -x REF/chr21 -U $sid.chr21.fastq -S hisat2/$sid.sam
done

#################### SAMTOOLS: VISUALIZAR, ORDENAR E INDEXAR ####################
export PATH=/home/vant/miniconda3/pkgs/samtools-1.7-1:$PATH
cd hisat2
for sid in $(ls *.sam  | cut -d "." -f1)
do
	samtools view -bh -S $sid.sam > $sid.bam
	samtools sort $sid.bam -o $sid.sorted.bam
	samtools index $sid.sorted.bam
done


#################### CONTAR LECTURAS CON HTseq-count  ####################
cd ..
mkdir -p htseq

for sid in $(ls *.fastq | cut -d "." -f1)
do
	htseq-count --format=bam --stranded=reverse --mode=intersection-nonempty --minaqual=10 --type=exon --idattr=gene_id --additional-attr=gene_name /home/vant/hisat2/$sid.sorted.bam /home/vant/transcriptomics-project/Apartado1/input/REF/GRCh38.gencode.v38.annotation.for.chr21.gtf > htseq/$sid.htseq
done

