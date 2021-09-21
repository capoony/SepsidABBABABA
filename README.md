# Bioinformatic pipeline for ABBA-BABA analyses of genomic allele frequency data of _Sepsis cynipsea_ and _S. neocynipsea_ in Europe

## (1) Trimming of raw reads, which are deposited on the SRA archive under [PRJNA612154](https://www.ncbi.nlm.nih.gov/bioproject?LinkName=biosample_bioproject&from_uid=14363816) was carried out with [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic).

```bash
java -jar trimmomatic-0.36.jar \
PE -threads 20 \
input_R1.fastq.gz \
input_R2.fastq.gz \
output_R1.fastq.gz \
output_R1_un.fastq.gz \
output_R2.fastq.gz \
output_R2_un.fastq.gz \
ILLUMINACLIP:AdapterSeq_new.fa:2:30:7:5:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:40
```

## (2) Mapping of the trimmed reads was carried out with bwa mem following the pipeline in [Kapun _et al._ 2020](https://academic.oup.com/mbe/article/37/9/2661/5837682) as described in [here](https://github.com/capoony/DrosEU_pipeline) using the reference genome of _S. thoracica_, which can be found [here](www.cgae.de/seto_01_genome_masked.fasta).

## (3) After merging the BAM files with _samtools mpileup_ as described [here](https://github.com/capoony/DrosEU_pipeline), we called SNPs using the heuristic SNP caller [PoolSNP](https://github.com/capoony/PoolSNP) using the following parameters.

```bash
sh PoolSNP.sh  \
mpileup=/ABBA_BABA.mpileup.gz \
reference=seto_01_genome.fasta \
names=PhC,PtC,SoC,ZuC,MoC,GeN,HoN,SoN,MoN,IZuC,IZuN,Sor  \
max-cov=0.9 \
min-cov=10 \
min-count=20 \
min-freq=0.01 \
miss-frac=0 \
jobs=10 \
base-quality=15 \
output=ABBA_BABA \
badsites=0
```
