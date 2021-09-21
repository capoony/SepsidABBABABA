# Bioinformatic pipeline for ABBA-BABA analyses of genomic allele frequency data of _Sepsis cynipsea_ and _S. neocynipsea_ in Europe

## (1) Trimming of raw reads was carried out with [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)

```bash
java -jar trimmomatic-0.36.jar \
PE -threads 20 \
input_R1.fastq.gz \
input_R2.fastq.gz \
output_R1.fastq.gz \
output_R1_un.fastq.gz \
output_R2.fastq.gz \
output_R2_un.fastq.gz \
ILLUMINACLIP:/Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/AdapterSeq_new.fa:2:30:7:5:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:40
```

## (2) Mapping of the trimmed reads was carried out with bwa mem following the pipeline in [Kapun _et al._ 2020](https://academic.oup.com/mbe/article/37/9/2661/5837682) as described in [here](https://github.com/capoony/DrosEU_pipeline) using the reference genome of _S. thoracica_, which can be found [here](www.cgae.de/seto_01_genome_masked.fasta).

## (3) SNPs were called using the heuristic SNP caller [PoolSNP](https://github.com/capoony/PoolSNP)

```bash
sh /Volumes/MartinResearch2/Wolf2019/scripts/PoolSNP-v1.07/PoolSNP.sh  \
mpileup=/Volumes/MartinResearch2/Wolf2019/mapping/thoracica/ABBA_BABA.mpileup.gz \
reference=/Volumes/MartinResearch2/Wolf2019/references/seto_01_genome.fasta \
names=ABBASorenbergCyn,ABBAZurichCyn,CevennesCyn,EstoniaCyn,LenzerheideCyn,LudwigshafenCyn,MaggiaCyn,MonteCeneriCyn,PetroiaCyn,SorenbergCyn,WesterwaldCyn,ZurichCyn,ABBASorenbergNeo,ABBAZurichNeo,BassettNeo,CevennesNeo,GeschinenNeo,HospentalNeo,KentuckyNeo,MontanaNeo,SierravilleNeo,SorenbergNeo,SorSim,Sor \
max-cov=/Volumes/MartinResearch2/Wolf2019/mapping/thoracica/ABBA_BABA-cov-0.9.txt \
min-cov=10 \
min-count=20 \
min-freq=0.01 \
miss-frac=0 \
jobs=10 \
base-quality=15 \
output=/Volumes/MartinResearch2/Wolf2019/mapping/thoracica/ABBA_BABA \
badsites=0
```
