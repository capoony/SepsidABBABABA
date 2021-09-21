#!/bin/sh

#\ \ AtheneNew.sh
#\ \
#
#\ \ Created\ by\ Martin\ Kapun\ on\ 16.07.19.
#\ \
##\ trim\ on\ server

Folder=(Cyn Neo)
Old=(CYN_FRA NEO_TUS)
New=(CYN_FRA NEO_FRA)

for i in ${!Folder[*]}

do

java -jar /Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/Trimmomatic/trimmomatic-0.36.jar \
PE -threads 20 \
/Volumes/MartinResearch3/Wolf2019/FASTQ/${Folder[i]}/20180601.B-${Old[i]}_R1.fastq.gz \
/Volumes/MartinResearch3/Wolf2019/FASTQ/${Folder[i]}/20180601.B-${Old[i]}_R2.fastq.gz \
/Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/trimmed/${New[i]}_R1.fastq.gz \
/Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/trimmed/${New[i]}_R1_un.fastq.gz \
/Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/trimmed/${New[i]}_R2.fastq.gz \
/Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/trimmed/${New[i]}_R2_un.fastq.gz \
ILLUMINACLIP:/Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/AdapterSeq_new.fa:2:30:7:5:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:40

done

Folder=(Cyn Neo)
Old=(CYN_FRA NEO_TUS)
New=(CYN_FRA NEO_FRA)

for i in ${!Folder[*]}

do

java -jar /Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/Trimmomatic/trimmomatic-0.36.jar \
PE -threads 20 \
/Volumes/MartinResearch3/Wolf2019/FASTQ/${Folder[i]}/20180601.B-${Old[i]}_R1.fastq.gz \
/Volumes/MartinResearch3/Wolf2019/FASTQ/${Folder[i]}/20180601.B-${Old[i]}_R2.fastq.gz \
/Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/trimmed/${New[i]}_R1.fastq.gz \
/Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/trimmed/${New[i]}_R1_un.fastq.gz \
/Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/trimmed/${New[i]}_R2.fastq.gz \
/Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/trimmed/${New[i]}_R2_un.fastq.gz \
ILLUMINACLIP:/Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/AdapterSeq_new.fa:2:30:7:5:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:40

done



mkdir /Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/map

# align your reads to the reference genome, introgression

bwa index /Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/ref/SEOR_allpaths.fasta

for i in CYN NEO

do

#bwa mem -t 5 -k 20 /Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/ref/SEOR_allpaths.fasta \
#/Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/trimmed/${i}_FRA_R1.fastq.gz \
#/Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/trimmed/${i}_FRA_R2.fastq.gz | \
#samtools view -Sbh -q 10 > /Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/map/${i}_FRA.bam

sambamba sort \
-o /Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/map/${i}_FRA-sort.bam \
-t 20 \
-p \
/Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/map/${i}_FRA.bam

done

for i in Bassett Geschinen Hospental Kentucky Montana Sierraville Sorenberg

do

bwa mem -t 22 -k 20 /Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/ref/SEOR_allpaths.fasta \
/Volumes/MartinResearch2/Wolf2019/FASTQ/3_Trimmo/Neo/${i}Neo/${i}Neo_R1_trimPair.fastq.gz \
/Volumes/MartinResearch2/Wolf2019/FASTQ/3_Trimmo/Neo/${i}Neo/${i}Neo_R2_trimPair.fastq.gz | \
samtools view -Sbh -q 10 > /Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/map/${i}Neo_ORTHO.bam

sambamba sort \
-o /Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/map/${i}Neo_ORTHO_sorted.bam \
-t 22 \
-p \
/Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/map/${i}Neo_ORTHO.bam

done


### OK and now something completely different

for Pop in Westerwald Petroia Maggia Ludwigshafen Lenzerheide MonteCeneri Zurich #Estonia Sorenberg Cevennes

do

echo /Volumes/MartinResearch2/Wolf2019/FASTQ/3_Trimmo/Cyn/${Pop}Cyn/${Pop}Cyn_R1_trimPair.fastq.gz

sh /Volumes/MartinResearch2/Wolf2019/scripts/pipeline/mapping_noTrim_noRealign.sh \
/Volumes/MartinResearch2/Wolf2019/FASTQ/3_Trimmo/Cyn/${Pop}Cyn/${Pop}Cyn_R1_trimPair.fastq.gz \
/Volumes/MartinResearch2/Wolf2019/FASTQ/3_Trimmo/Cyn/${Pop}Cyn/${Pop}Cyn_R2_trimPair.fastq.gz \
${Pop} \
/Volumes/MartinResearch2/Wolf2019/mapping/thoracica/${Pop}Cyn \
bwa \
/Volumes/MartinResearch2/Wolf2019/references/seto_01_genome

done

for Pop in Bassett Geschinen Hospental Kentucky Sierraville #Montana Sorenberg Cevennes

do

echo /Volumes/MartinResearch2/Wolf2019/FASTQ/3_Trimmo/Neo/${Pop}Neo/${Pop}Neo_R1_trimPair.fastq.gz

sh /Volumes/MartinResearch2/Wolf2019/scripts/pipeline/mapping_noTrim_noRealign.sh \
/Volumes/MartinResearch2/Wolf2019/FASTQ/3_Trimmo/Neo/${Pop}Neo/${Pop}Neo_R1_trimPair.fastq.gz \
/Volumes/MartinResearch2/Wolf2019/FASTQ/3_Trimmo/Neo/${Pop}Neo/${Pop}Neo_R2_trimPair.fastq.gz \
${Pop} \
/Volumes/MartinResearch2/Wolf2019/mapping/thoracica/${Pop}Neo \
bwa \
/Volumes/MartinResearch2/Wolf2019/references/seto_01_genome

done

for Pop in ABBASorenbergCyn ABBASorenbergNeo ABBAZurichCyn ABBAZurichNeo

do

echo /Volumes/MartinResearch2/Wolf2019/FASTQ/3_Trimmo/ABBA/${Pop}/${Pop}_R1_trimPair.fastq.gz

sh /Volumes/MartinResearch2/Wolf2019/scripts/pipeline/mapping_noTrim_noRealign.sh \
/Volumes/MartinResearch2/Wolf2019/FASTQ/3_Trimmo/ABBA/${Pop}/${Pop}_R1_trimPair.fastq.gz \
/Volumes/MartinResearch2/Wolf2019/FASTQ/3_Trimmo/ABBA/${Pop}/${Pop}_R2_trimPair.fastq.gz \
${Pop} \
/Volumes/MartinResearch2/Wolf2019/mapping/thoracica/${Pop} \
bwa \
/Volumes/MartinResearch2/Wolf2019/references/seto_01_genome

done


NAME=(Estonia Sorenberg Cevennes Montana Sorenberg Cevennes)
Species=(Cyn Cyn Cyn Neo Neo Neo)

NAME=(Westerwald Petroia Maggia Ludwigshafen Lenzerheide MonteCeneri Zurich)
Species=(Cyn Cyn Cyn Cyn Cyn Cyn Cyn)

NAME=(Bassett Geschinen Hospental Kentucky Sierraville)
Species=(Neo Neo Neo Neo Neo)

for i in ${!NAME[*]} #ABBASorenbergCyn ABBASorenbergNeo ABBAZurichCyn ABBAZurichNeo

do


## f) re-align around InDels

java -Xmx20g -jar /Volumes/MartinResearch2/Wolf2019/scripts/pipeline/scripts/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar \
-T IndelRealigner \
-R /Volumes/MartinResearch2/Wolf2019/references/seto_01_genome.fasta \
-I /Volumes/MartinResearch2/Wolf2019/mapping/thoracica/${NAME[i]}${Species[i]}/mapping/${NAME[i]}-dedup.bam \
-targetIntervals /Volumes/MartinResearch2/Wolf2019/mapping/thoracica/${NAME[i]}${Species[i]}/mapping//realign_list/${NAME[i]}.list \
-o /Volumes/MartinResearch2/Wolf2019/mapping/thoracica/${NAME[i]}${Species[i]}/mapping/${NAME[i]}-dedup_InDel.bam &

done

for i in ABBASorenbergCyn ABBASorenbergNeo ABBAZurichCyn ABBAZurichNeo

do


## f) re-align around InDels

java -Xmx20g -jar /Volumes/MartinResearch2/Wolf2019/scripts/pipeline/scripts/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar \
-T IndelRealigner \
-R /Volumes/MartinResearch2/Wolf2019/references/seto_01_genome.fasta \
-I /Volumes/MartinResearch2/Wolf2019/mapping/thoracica/${i}/mapping/${i}-dedup.bam \
-targetIntervals /Volumes/MartinResearch2/Wolf2019/mapping/thoracica/${i}/mapping//realign_list/${i}.list \
-o /Volumes/MartinResearch2/Wolf2019/mapping/thoracica/${i}/mapping/${i}-dedup_InDel.bam &

done

### now merge and call SNPs

echo '''
/Volumes/MartinResearch2/Wolf2019/mapping/thoracica/ABBASorenbergCyn/mapping/ABBASorenbergCyn-dedup_InDel.bam
/Volumes/MartinResearch2/Wolf2019/mapping/thoracica/ABBAZurichCyn/mapping/ABBAZurichCyn-dedup_InDel.bam
/Volumes/MartinResearch2/Wolf2019/mapping/thoracica/CevennesCyn/mapping/Cevennes-dedup_InDel.bam
/Volumes/MartinResearch2/Wolf2019/mapping/thoracica/EstoniaCyn/mapping/Estonia-dedup_InDel.bam
/Volumes/MartinResearch2/Wolf2019/mapping/thoracica/LenzerheideCyn/mapping/Lenzerheide-dedup_InDel.bam
/Volumes/MartinResearch2/Wolf2019/mapping/thoracica/LudwigshafenCyn/mapping/Ludwigshafen-dedup_InDel.bam
/Volumes/MartinResearch2/Wolf2019/mapping/thoracica/MaggiaCyn/mapping/Maggia-dedup_InDel.bam
/Volumes/MartinResearch2/Wolf2019/mapping/thoracica/MonteCeneriCyn/mapping/MonteCeneri-dedup_InDel.bam
/Volumes/MartinResearch2/Wolf2019/mapping/thoracica/PetroiaCyn/mapping/Petroia-dedup_InDel.bam
/Volumes/MartinResearch2/Wolf2019/mapping/thoracica/SorenbergCyn/mapping/Sorenberg-dedup_InDel.bam
/Volumes/MartinResearch2/Wolf2019/mapping/thoracica/WesterwaldCyn/mapping/Westerwald-dedup_InDel.bam
/Volumes/MartinResearch2/Wolf2019/mapping/thoracica/ZurichCyn/mapping/Zurich-dedup_InDel.bam
/Volumes/MartinResearch2/Wolf2019/mapping/thoracica/ABBASorenbergNeo/mapping/ABBASorenbergNeo-dedup_InDel.bam
/Volumes/MartinResearch2/Wolf2019/mapping/thoracica/ABBAZurichNeo/mapping/ABBAZurichNeo-dedup_InDel.bam
/Volumes/MartinResearch2/Wolf2019/mapping/thoracica/BassettNeo/mapping/Bassett-dedup_InDel.bam
/Volumes/MartinResearch2/Wolf2019/mapping/thoracica/CevennesNeo/mapping/Cevennes-dedup_InDel.bam
/Volumes/MartinResearch2/Wolf2019/mapping/thoracica/GeschinenNeo/mapping/Geschinen-dedup_InDel.bam
/Volumes/MartinResearch2/Wolf2019/mapping/thoracica/HospentalNeo/mapping/Hospental-dedup_InDel.bam
/Volumes/MartinResearch2/Wolf2019/mapping/thoracica/KentuckyNeo/mapping/Kentucky-dedup_InDel.bam
/Volumes/MartinResearch2/Wolf2019/mapping/thoracica/MontanaNeo/mapping/Montana-dedup_InDel.bam
/Volumes/MartinResearch2/Wolf2019/mapping/thoracica/SierravilleNeo/mapping/Sierraville-dedup_InDel.bam
/Volumes/MartinResearch2/Wolf2019/mapping/thoracica/SorenbergNeo/mapping/Sorenberg-dedup_InDel.bam
/Volumes/MartinResearch2/Wolf2019/mapping/thoracica/SorSim/mapping/SorSim-dedup_InDel.bam
/Volumes/MartinResearch2/Wolf2019/mapping/thoracica/Sor/mapping/SorSim-dedup_InDel.bam
''' > /Volumes/MartinResearch2/Wolf2019/mapping/thoracica/bamlist.txt

## merge populations with samtools

samtools mpileup -B \
-f /Volumes/MartinResearch2/Wolf2019/references/seto_01_genome.fasta \
-b /Volumes/MartinResearch2/Wolf2019/mapping/thoracica/bamlist.txt \
-q 20 \
-Q 20 \
| gzip > /Volumes/MartinResearch2/Wolf2019/mapping/thoracica/ABBA_BABA.mpileup.gz

## call SNPs with PoolSNP

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

sh /Volumes/MartinResearch2/Wolf2019/scripts/PoolSNP-v1.07/PoolSNP.sh  \
mpileup=/Volumes/MartinResearch2/Wolf2019/mapping/thoracica/ABBA_BABA.mpileup.gz \
reference=/Volumes/MartinResearch2/Wolf2019/references/seto_01_genome.fasta \
names=ABBASorenbergCyn,ABBAZurichCyn,CevennesCyn,EstoniaCyn,LenzerheideCyn,LudwigshafenCyn,MaggiaCyn,MonteCeneriCyn,PetroiaCyn,SorenbergCyn,WesterwaldCyn,ZurichCyn,ABBASorenbergNeo,ABBAZurichNeo,BassettNeo,CevennesNeo,GeschinenNeo,HospentalNeo,KentuckyNeo,MontanaNeo,SierravilleNeo,SorenbergNeo,SorSim,Sor \
max-cov=/Volumes/MartinResearch2/Wolf2019/mapping/thoracica/ABBA_BABA-cov-0.9.txt \
min-cov=10 \
min-count=20 \
min-freq=0.01 \
miss-frac=0.1 \
jobs=10 \
base-quality=15 \
output=/Volumes/MartinResearch2/Wolf2019/mapping/thoracica/ABBA_BABA_missfrac0.1 \
badsites=0


python2.7 /Volumes/MartinResearch3/Tad2019/Patrick/scripts/filter-pos-from-vcf.py \
NA \
/Volumes/MartinResearch2/Wolf2019/references/seto_repeat_annotation.gff3 \
/Volumes/MartinResearch2/Wolf2019/mapping/thoracica/ABBA_BABA_missfrac0.1.vcf.gz \
| gzip > /Volumes/MartinResearch2/Wolf2019/mapping/thoracica/ABBA_BABA-filtered.vcf.gz

## convert to sync file
gunzip -c /Volumes/MartinResearch2/Wolf2019/mapping/thoracica/ABBA_BABA-filtered.vcf.gz  | parallel \
-k \
--pipe \
-j20 \
--no-notice \
--cat python /Volumes/MartinResearch3/Tad2019/Patrick/scripts/vcf2sync.py \
--input {} \
| gzip > /Volumes/MartinResearch2/Wolf2019/mapping/thoracica/ABBA_BABA-filtered.sync.gz


## calculate MAF
gunzip -c /Volumes/MartinResearch2/Wolf2019/mapping/thoracica/ABBA_BABA-filtered.sync.gz \
| parallel \
-k \
--pipe \
-j 20 \
--no-notice \
--cat python2.7 /Volumes/MartinResearch3/Tad2019/Patrick/scripts/AFbyAllele.py \
{} \
| gzip > /Volumes/MartinResearch2/Wolf2019/mapping/thoracica/ABBA_BABA-filtered.af.gz

## calculate MAF
gunzip -c /Volumes/MartinResearch2/Wolf2019/mapping/thoracica/ABBA_BABA-filtered.sync.gz \
| parallel \
-k \
--pipe \
-j 20 \
--no-notice \
--cat python /Users/mkapun/Documents/GitHub/ABBA-BABA-4-AF/SYNC2AF.py \
--sync {} \
--MinCov 10 \
| gzip > /Volumes/MartinResearch2/Wolf2019/mapping/thoracica/ABBA_BABA-filtered_mc10.af.gz


## first test the script

python /Volumes/MartinResearch2/Wolf2019/scripts/ABBA_BABA.py --AlleleFrequencies /Volumes/MartinResearch2/Wolf2019/scripts/TestABBA_BABA/testinput.txt --output /Volumes/MartinResearch2/Wolf2019/scripts/TestABBA_BABA/testoutput.txt --SNPs 5 --order 2,3,4,5

python3 /Volumes/MartinResearch2/Wolf2019/scripts/ABBA_BABA.py --AlleleFrequencies /Volumes/MartinResearch2/Wolf2019/mapping/thoracica/ABBA_BABA-filtered.af.gz --output /Volumes/MartinResearch2/Wolf2019/scripts/TestABBA_BABA/PeC_CeC_CeN_Sor.txt --SNPs 100 --order -16,-22,-9,-1


## now calculate ABBA-BABA in multiple different combinations

# at first test the three sympatric populations using Montana neocynipsea as H1

mkdir /Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/AbbaBaba_MoN

python /Volumes/MartinResearch2/Wolf2019/scripts/ABBA_BABA_driver.py \
--input /Volumes/MartinResearch2/Wolf2019/mapping/thoracica/ABBA_BABA-filtered.af.gz \
--names ISoC,IZuC,CeC,EsC,LeC,LuC,MaC,MoC,PeC,SoC,WeC,ZuC,ISoN,IZuN,BaN,CeN,GeN,HoN,KeN,MoN,SiN,SoN,SorS,Sor \
--combinations MoN,SoN,SoC+MoN,ISoN,ISoC+MoN,SoN,ISoC+MoN,ISoN,SoC+MoN,IZuN,IZuC+MoN,IZuN,ZuC+MoN,CeN,CeC \
--outgroup Sor \
--SNPs 100,500,1000 \
--output /Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/AbbaBaba_MoN

# Then test the three sympatric populations using Sierraville neocynipsea as H1

mkdir /Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/AbbaBaba_SiN

python /Volumes/MartinResearch2/Wolf2019/scripts/ABBA_BABA_driver.py \
--input /Volumes/MartinResearch2/Wolf2019/mapping/thoracica/ABBA_BABA-filtered.af.gz \
--names ISoC,IZuC,CeC,EsC,LeC,LuC,MaC,MoC,PeC,SoC,WeC,ZuC,ISoN,IZuN,BaN,CeN,GeN,HoN,KeN,MoN,SiN,SoN,SorS,Sor \
--combinations SiN,SoN,SoC+SiN,ISoN,ISoC+SiN,SoN,ISoC+SiN,ISoN,SoC+SiN,IZuN,IZuC+SiN,IZuN,ZuC+SiN,CeN,CeC \
--outgroup Sor \
--SNPs 100,500,1000 \
--output /Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/AbbaBaba_SiN

# Then test the three sympatric populations using Geschinen neocynipsea as H1

mkdir /Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/AbbaBaba_GeN

python2.7 /Volumes/MartinResearch2/Wolf2019/scripts/ABBA_BABA_driver.py \
--input /Volumes/MartinResearch2/Wolf2019/mapping/thoracica/ABBA_BABA-filtered.af.gz \
--names ISoC,IZuC,CeC,EsC,LeC,LuC,MaC,MoC,PeC,SoC,WeC,ZuC,ISoN,IZuN,BaN,CeN,GeN,HoN,KeN,MoN,SiN,SoN,SorS,Sor \
--combinations GeN,SoN,SoC+GeN,ISoN,ISoC+GeN,SoN,ISoC+GeN,ISoN,SoC+GeN,IZuN,IZuC+GeN,IZuN,ZuC+GeN,CeN,CeC \
--outgroup Sor \
--SNPs 100,500,1000 \
--output /Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/AbbaBaba_GeN

# Finally test the three sympatric populations using Hospental neocynipsea as H1

mkdir /Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/AbbaBaba_HoN

python2.7 /Volumes/MartinResearch2/Wolf2019/scripts/ABBA_BABA_driver.py \
--input /Volumes/MartinResearch2/Wolf2019/mapping/thoracica/ABBA_BABA-filtered.af.gz \
--names ISoC,IZuC,CeC,EsC,LeC,LuC,MaC,MoC,PeC,SoC,WeC,ZuC,ISoN,IZuN,BaN,CeN,GeN,HoN,KeN,MoN,SiN,SoN,SorS,Sor \
--combinations HoN,SoN,SoC+HoN,ISoN,ISoC+HoN,SoN,ISoC+HoN,ISoN,SoC+HoN,IZuN,IZuC+HoN,IZuN,ZuC+HoN,CeN,CeC \
--outgroup Sor \
--SNPs 100,500,1000 \
--output /Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/AbbaBaba_HoN

# now flip the tree topology and use Estonia cynipase as H1

mkdir /Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/AbbaBaba_EsC

python2.7 /Volumes/MartinResearch2/Wolf2019/scripts/ABBA_BABA_driver.py \
--input /Volumes/MartinResearch2/Wolf2019/mapping/thoracica/ABBA_BABA-filtered.af.gz \
--names ISoC,IZuC,CeC,EsC,LeC,LuC,MaC,MoC,PeC,SoC,WeC,ZuC,ISoN,IZuN,BaN,CeN,GeN,HoN,KeN,MoN,SiN,SoN,SorS,Sor \
--combinations EsC,SoC,SoN+EsC,ISoC,ISoN+EsC,ISoC,SoN+EsC,SoC,ISoN+EsC,IZuC,IZuN+EsC,ZuC,IZuN+EsC,CeC,CeN \
--outgroup Sor \
--SNPs 100,500,1000 \
--output /Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/AbbaBaba_EsC

mkdir /Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/AbbaBaba_EsC_mc10

python2.7 /Volumes/MartinResearch2/Wolf2019/scripts/ABBA_BABA_driver.py \
--input /Volumes/MartinResearch2/Wolf2019/mapping/thoracica/ABBA_BABA-filtered_mc10.af.gz \
--names ISoC,IZuC,CeC,EsC,LeC,LuC,MaC,MoC,PeC,SoC,WeC,ZuC,ISoN,IZuN,BaN,CeN,GeN,HoN,KeN,MoN,SiN,SoN,SorS,Sor \
--combinations EsC,SoC,SoN+EsC,ISoC,ISoN+EsC,ISoC,SoN+EsC,SoC,ISoN+EsC,IZuC,IZuN+EsC,ZuC,IZuN+EsC,CeC,CeN \
--outgroup Sor \
--SNPs 100,500,1000 \
--window NA \
--output /Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/AbbaBaba_EsC_mc10


# Then test the three sympatric populations using Petroia cynipsea as H1

mkdir /Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/AbbaBaba_PeC

python2.7 /Volumes/MartinResearch2/Wolf2019/scripts/ABBA_BABA_driver.py \
--input /Volumes/MartinResearch2/Wolf2019/mapping/thoracica/ABBA_BABA-filtered.af.gz \
--names ISoC,IZuC,CeC,EsC,LeC,LuC,MaC,MoC,PeC,SoC,WeC,ZuC,ISoN,IZuN,BaN,CeN,GeN,HoN,KeN,MoN,SiN,SoN,SorS,Sor \
--combinations PeC,SoC,SoN+PeC,ISoC,ISoN+PeC,ISoC,SoN+PeC,SoC,ISoN+PeC,IZuC,IZuN+PeC,ZuC,IZuN+PeC,CeC,CeN \
--outgroup Sor \
--SNPs 100,500,1000 \
--output /Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/AbbaBaba_PeC

# Then test the three sympatric populations using Westerwald cynipsea as H1

mkdir /Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/AbbaBaba_WeC

python /Volumes/MartinResearch2/Wolf2019/scripts/ABBA_BABA_driver.py \
--input /Volumes/MartinResearch2/Wolf2019/mapping/thoracica/ABBA_BABA-filtered.af.gz \
--names ISoC,IZuC,CeC,EsC,LeC,LuC,MaC,MoC,PeC,SoC,WeC,ZuC,ISoN,IZuN,BaN,CeN,GeN,HoN,KeN,MoN,SiN,SoN,SorS,Sor \
--combinations WeC,SoC,SoN+WeC,ISoC,ISoN+WeC,ISoC,SoN+WeC,SoC,ISoN+WeC,IZuC,IZuN+WeC,ZuC,IZuN+WeC,CeC,CeN \
--outgroup Sor \
--SNPs 100,500,1000 \
--output /Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/AbbaBaba_WeC

for i in /Users/martinkapun/Desktop/ABBA/500/*.stat ; do awk '{print FILENAME"\t"$0}' $i >> /Users/martinkapun/Desktop/ABBA/stats.txt ; done



### plot AFS of Isofemale lines

echo '''

library(plotrix)

AF=read.table(gzfile("/Volumes/MartinResearch2/Wolf2019/mapping/thoracica/ABBA_BABA-filtered.af.gz", open = "r"),header=F,na.strings = "NA")

MinAl <- function(x) { if(x>0.5){return (1-x)}else{return(x)}}

H=list()
X=sapply(na.omit(AF[,4]),MinAl)
H$"ISoC"<-X[X!=0]
X=sapply(na.omit(AF[,5]),MinAl)
H$"IZuC"<-X[X!=0]
X=sapply(na.omit(AF[,6]),MinAl)
#H$"CeC"<-X[X!=0]
X=sapply(na.omit(AF[,16]),MinAl)
H$"ISoN"<-X[X!=0]
X=sapply(na.omit(AF[,17]),MinAl)
H$"IZuN"<-X[X!=0]

pdf("/Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/Isofemale_AFS.pdf",width=8,height=6)
par(las=1,cex=1,yaxt="n")
M=multhist(H, breaks=10,freq=F,col=grey.colors(4),legend=T)
legend("topright",legend=names(H),col=grey.colors(4),pch=15,cex=0.75)
t<-axTicks(2)
par(yaxt="s")
axis(2,at=t,labels=round(t/rowSums(M[[2]])[1],2))
box()
dev.off()

''' > /Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/Isofemale_AFS.r

Rscript /Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/Isofemale_AFS.r


## count mapped reads

for i in /Volumes/MartinResearch2/Wolf2019/FASTQ/3_Trimmo/*/*/*R1*

do

echo ${i##*/} `gunzip -c "$i" | wc -l | awk '{print $1}'` >> /Volumes/MartinResearch2/Wolf2019/mapping/thoracica/RawReads.txt

done


## count raw reads

for i in /Volumes/MartinResearch2/Wolf2019/mapping/thoracica/*/mapping/*InDel.bam

do

echo ${i##*/} `samtools view -c "$i" | awk '{print $1}'` >> /Volumes/MartinResearch2/Wolf2019/mapping/thoracica/mappedReads.txt

done

## count raw mapped reads

for i in /Volumes/MartinResearch2/Wolf2019/mapping/thoracica/*/mapping/*raw.bam

do

echo ${i##*/} `samtools view -c -F 260 "$i" | awk '{print $1}'` >> /Volumes/MartinResearch2/Wolf2019/mapping/thoracica/RawMappedReads.txt

done

## count mapped reads

for i in /Volumes/MartinResearch2/Wolf2019/FASTQ/Seor/SEOR_250_trimmed_paired_1.fq.gz

do

echo ${i##*/} `gunzip -c "$i" | wc -l | awk '{print $1}'` >> /Volumes/MartinResearch2/Wolf2019/mapping/thoracica/RawReads.txt

done


### calculate pi
sed 's/ /:/g' /Volumes/MartinResearch2/Wolf2019/analyses/checkFASTA/seto.txt > /Volumes/MartinResearch2/Wolf2019/analyses/checkFASTA/seto.cont

mkdir /Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/pi

python2.7 /Volumes/MartinResearch2/Wolf2019/scripts/pipeline4Publication/scripts/TrueWindows.py \
--badcov /Volumes/MartinResearch2/Wolf2019/mapping/thoracica/ABBA_BABA_BS.txt.gz \
--te /Volumes/MartinResearch2/Wolf2019/references/seto_repeat_annotation.gff3 \
--window 200000 \
--step 200000 \
--chromosomes /Volumes/MartinResearch2/Wolf2019/analyses/checkFASTA/seto.cont \
--output /Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/pi/200k_windows.txt

python2.7 /Volumes/MartinResearch2/Wolf2019/scripts/pipeline4Publication/scripts/PopGen-var.py \
--input /Volumes/MartinResearch2/Wolf2019/mapping/thoracica/ABBA_BABA-filtered.sync.gz \
--pool-size 50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50 \
--min-count 2 \
--window 200000 \
--step 200000 \
--sitecount /Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/pi/200k_windows.txt-200000-200000.txt \
--min-sites-frac 0.75 \
--output /Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/pi/popgen



### blast regions in 5% tails against Drosophila AA sequences

mkdir /Volumes/MartinResearch2/Wolf2019/references/DmelProt

# make BLAST DB

makeblastdb -in /Volumes/MartinResearch2/Wolf2019/references/DmelProt/dmel-all-translation-r6.31.fasta -parse_seqids -dbtype prot -parse_seqids

for i in /Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/*/*/*_500SNPs.abba

do

echo $i

sh /Volumes/MartinResearch2/Wolf2019/shell/doBlast.sh $i &

NPROC=$(($NPROC+1))
if [ "$NPROC" -ge 15 ]; then
wait
NPROC=0
fi

done


mkdir /Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/GO

for i in /Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/*/*/*_500SNPs.abba_Top.go

do

echo $i

python /Volumes/MartinResearch2/Wolf2019/scripts/FixGO.py $i >> /Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/GO/All_Go2.txt

done

## calculate coverages

mkdir /Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/cov

python2.7 /Volumes/MartinResearch3/Tad2019/Patrick/scripts/mpileup2cov.py \
/Volumes/MartinResearch2/Wolf2019/mapping/thoracica/ABBA_BABA.mpileup.gz \
/Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/cov/ABBA_BABA.cov \
ISoC,IZuC,CeC,EsC,LeC,LuC,MaC,MoC,PeC,SoC,WeC,ZuC,ISoN,IZuN,BaN,CeN,GeN,HoN,KeN,MoN,SiN,SoN,SorS,Sor

python2.7 /Volumes/MartinResearch2/Wolf2019/scripts/MergeCovNStuff.py \
/Volumes/MartinResearch2/Wolf2019/analyses/checkFASTA/seto.txt \
/Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/cov/ABBA_BABA.cov_mean \
> /Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/cov/ABBA_BABA_cov_mean_len.txt

rscript /Volumes/MartinResearch2/Wolf2019/scripts/findCoveragePeaks.r \
/Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/cov/ABBA_BABA_cov_mean_len


## test new window function

python2.7 /Volumes/MartinResearch2/Wolf2019/scripts/ABBA_BABA_driver.py \
--input /Volumes/MartinResearch2/Wolf2019/mapping/thoracica/ABBA_BABA-filtered.af.gz \
--names ISoC,IZuC,CeC,EsC,LeC,LuC,MaC,MoC,PeC,SoC,WeC,ZuC,ISoN,IZuN,BaN,CeN,GeN,HoN,KeN,MoN,SiN,SoN,SorS,Sor \
--combinations EsC,SoC,SoN+EsC,ISoC,ISoN+EsC,ISoC,SoN+EsC,SoC,ISoN+EsC,IZuC,IZuN+EsC,ZuC,IZuN+EsC,CeC,CeN \
--outgroup Sor \
--SNPs NA \
--window 100000 \
--output /Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/AbbaBaba_EsC_window

mkdir /Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/AbbaBaba_EsC_window_200k

python2.7 /Volumes/MartinResearch2/Wolf2019/scripts/ABBA_BABA_driver.py \
--input /Volumes/MartinResearch2/Wolf2019/mapping/thoracica/ABBA_BABA-filtered.af.gz \
--names ISoC,IZuC,CeC,EsC,LeC,LuC,MaC,MoC,PeC,SoC,WeC,ZuC,ISoN,IZuN,BaN,CeN,GeN,HoN,KeN,MoN,SiN,SoN,SorS,Sor \
--combinations EsC,SoC,SoN+EsC,ISoC,ISoN+EsC,ISoC,SoN+EsC,SoC,ISoN+EsC,IZuC,IZuN+EsC,ZuC,IZuN+EsC,CeC,CeN \
--outgroup Sor \
--SNPs NA \
--window 200000 \
--output /Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/AbbaBaba_EsC_window_200k

mkdir  /Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/AbbaBaba_EsC_window_200k_mc10

python2.7 /Volumes/MartinResearch2/Wolf2019/scripts/ABBA_BABA_driver.py \
--input /Volumes/MartinResearch2/Wolf2019/mapping/thoracica/ABBA_BABA-filtered_mc10.af.gz \
--names ISoC,IZuC,CeC,EsC,LeC,LuC,MaC,MoC,PeC,SoC,WeC,ZuC,ISoN,IZuN,BaN,CeN,GeN,HoN,KeN,MoN,SiN,SoN,SorS,Sor \
--combinations EsC,SoC,SoN+EsC,ISoC,ISoN+EsC,ISoC,SoN+EsC,SoC,ISoN+EsC,IZuC,IZuN+EsC,ZuC,IZuN+EsC,CeC,CeN \
--outgroup Sor \
--SNPs NA \
--window 200000 \
--output /Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/AbbaBaba_EsC_window_200k_mc10

mkdir /Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/AbbaBaba_PeC_100k

python2.7 /Volumes/MartinResearch2/Wolf2019/scripts/ABBA_BABA_driver.py \
--input /Volumes/MartinResearch2/Wolf2019/mapping/thoracica/ABBA_BABA-filtered.af.gz \
--names ISoC,IZuC,CeC,EsC,LeC,LuC,MaC,MoC,PeC,SoC,WeC,ZuC,ISoN,IZuN,BaN,CeN,GeN,HoN,KeN,MoN,SiN,SoN,SorS,Sor \
--combinations PeC,SoC,SoN+PeC,ISoC,ISoN+PeC,ISoC,SoN+PeC,SoC,ISoN+PeC,IZuC,IZuN+PeC,ZuC,IZuN+PeC,CeC,CeN \
--outgroup Sor \
--SNPs NA \
--window 100000 \
--output /Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/AbbaBaba_PeC_100k


mkdir /Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/AbbaBaba_PeC_200k

python2.7 /Volumes/MartinResearch2/Wolf2019/scripts/ABBA_BABA_driver.py \
--input /Volumes/MartinResearch2/Wolf2019/mapping/thoracica/ABBA_BABA-filtered.af.gz \
--names ISoC,IZuC,CeC,EsC,LeC,LuC,MaC,MoC,PeC,SoC,WeC,ZuC,ISoN,IZuN,BaN,CeN,GeN,HoN,KeN,MoN,SiN,SoN,SorS,Sor \
--combinations PeC,SoC,SoN+PeC,ISoC,ISoN+PeC,ISoC,SoN+PeC,SoC,ISoN+PeC,IZuC,IZuN+PeC,ZuC,IZuN+PeC,CeC,CeN \
--outgroup Sor \
--SNPs NA \
--window 200000 \
--output /Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/AbbaBaba_PeC_200k

mkdir /Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/AbbaBaba_GeN_100k

python2.7 /Volumes/MartinResearch2/Wolf2019/scripts/ABBA_BABA_driver.py \
--input /Volumes/MartinResearch2/Wolf2019/mapping/thoracica/ABBA_BABA-filtered.af.gz \
--names ISoC,IZuC,CeC,EsC,LeC,LuC,MaC,MoC,PeC,SoC,WeC,ZuC,ISoN,IZuN,BaN,CeN,GeN,HoN,KeN,MoN,SiN,SoN,SorS,Sor \
--combinations GeN,SoN,SoC+GeN,ISoN,ISoC+GeN,SoN,ISoC+GeN,ISoN,SoC+GeN,IZuN,IZuC+GeN,IZuN,ZuC+GeN,CeN,CeC \
--outgroup Sor \
--SNPs NA \
--window 100000 \
--output /Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/AbbaBaba_GeN_100k

mkdir /Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/AbbaBaba_GeN_200k

python2.7 /Volumes/MartinResearch2/Wolf2019/scripts/ABBA_BABA_driver.py \
--input /Volumes/MartinResearch2/Wolf2019/mapping/thoracica/ABBA_BABA-filtered.af.gz \
--names ISoC,IZuC,CeC,EsC,LeC,LuC,MaC,MoC,PeC,SoC,WeC,ZuC,ISoN,IZuN,BaN,CeN,GeN,HoN,KeN,MoN,SiN,SoN,SorS,Sor \
--combinations GeN,SoN,SoC+GeN,ISoN,ISoC+GeN,SoN,ISoC+GeN,ISoN,SoC+GeN,IZuN,IZuC+GeN,IZuN,ZuC+GeN,CeN,CeC \
--outgroup Sor \
--SNPs NA \
--window 200000 \
--output /Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/AbbaBaba_GeN_200k


# Finally test the three sympatric populations using Hospental neocynipsea as H1

mkdir /Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/AbbaBaba_HoN_100k

python2.7 /Volumes/MartinResearch2/Wolf2019/scripts/ABBA_BABA_driver.py \
--input /Volumes/MartinResearch2/Wolf2019/mapping/thoracica/ABBA_BABA-filtered.af.gz \
--names ISoC,IZuC,CeC,EsC,LeC,LuC,MaC,MoC,PeC,SoC,WeC,ZuC,ISoN,IZuN,BaN,CeN,GeN,HoN,KeN,MoN,SiN,SoN,SorS,Sor \
--combinations HoN,SoN,SoC+HoN,ISoN,ISoC+HoN,SoN,ISoC+HoN,ISoN,SoC+HoN,IZuN,IZuC+HoN,IZuN,ZuC+HoN,CeN,CeC \
--outgroup Sor \
--SNPs NA \
--window 100000 \
--output /Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/AbbaBaba_HoN_100k

mkdir /Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/AbbaBaba_HoN_200k

python2.7 /Volumes/MartinResearch2/Wolf2019/scripts/ABBA_BABA_driver.py \
--input /Volumes/MartinResearch2/Wolf2019/mapping/thoracica/ABBA_BABA-filtered.af.gz \
--names ISoC,IZuC,CeC,EsC,LeC,LuC,MaC,MoC,PeC,SoC,WeC,ZuC,ISoN,IZuN,BaN,CeN,GeN,HoN,KeN,MoN,SiN,SoN,SorS,Sor \
--combinations HoN,SoN,SoC+HoN,ISoN,ISoC+HoN,SoN,ISoC+HoN,ISoN,SoC+HoN,IZuN,IZuC+HoN,IZuN,ZuC+HoN,CeN,CeC \
--outgroup Sor \
--SNPs NA \
--window 200000 \
--output /Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/AbbaBaba_HoN_200k


mkdir  /Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/AbbaBaba_Zu_window_200k_mc10

python2.7 /Volumes/MartinResearch2/Wolf2019/scripts/ABBA_BABA_driver.py \
--input /Volumes/MartinResearch2/Wolf2019/mapping/thoracica/ABBA_BABA-filtered_mc10.af.gz \
--names ISoC,IZuC,CeC,EsC,LeC,LuC,MaC,MoC,PeC,SoC,WeC,ZuC,ISoN,IZuN,BaN,CeN,GeN,HoN,KeN,MoN,SiN,SoN,SorS,Sor \
--combinations ZuC,SoC,SoN+IZuC,SoC,SoN+ZuC,ISoC,ISoN+IZuC,ISoC,ISoN+ZuC,ISoC,SoN+IZuN,SoN,SoC+IZuN,ISoN,ISoc \
--outgroup Sor \
--SNPs NA \
--window 200000 \
--output /Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/AbbaBaba_Zu_window_200k_mc10


### climate

mkdir /Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/climate

### now convert coordinates in decimal format and store as coordinates.txt


echo '''

# load installed package to the R project
require(raster)
# first load WC bio variables at the resolution of 2.5 deg
biod <- getData("worldclim2", var="bio", res=2.5)
tmind <- getData("worldclim", var="tmin", res=2.5)
tmaxd <- getData("worldclim", var="tmax", res=2.5)
precd <- getData("worldclim", var="prec", res=2.5)

# read csv file with geographic coordinates
geod<-read.table("/Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/climate/coordinates.txt", header=T, stringsAsFactors=F)

# extact for each coordinate bio clim variables
bio<-extract(biod, geod[,c(3,2)])
tmin<-extract(tmind, geod[,c(3,2)])
tmax<-extract(tmaxd, geod[,c(3,2)])
precd<-extract(precd, geod[,c(3,2)])

# create a full dataset
bio.data<-cbind(geod,bio,tmin,tmax,precd)

# save into external file
write.table(bio.data,file="/Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/climate/coordinates-GIS.txt",sep="\t", row.names=FALSE ,quote=FALSE)


xx<- c(c(1:12,12:1),c(1:12,12:1),c(1:12,12:1))
yy <- c(c(tmin[1,]/10,rev(tmax[1,]/10)),c(tmin[2,]/10,rev(tmax[2,]/10)),c(tmin[3,]/10,rev(tmax[3,]/10)))

pdf("/Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/climate/climate.pdf",width=12,height=6)

plot   (xx, yy, type = "n", xlab = "Months", ylab = "Temperature")

polygon(c(1:12,12:1),c(tmin[1,]/10,rev(tmax[1,]/10)),col = rgb(0,0,1,0.2),border="blue")
polygon(c(1:12,12:1),c(tmin[2,]/10,rev(tmax[2,]/10)),col = rgb(1,0,0,0.2),border="red")
#polygon(c(1:12,12:1),c(tmin[3,]/10,rev(tmax[3,]/10)),col = rgb(0.2,0.2,0.2,0.2), border = "green")

par(new = TRUE)

plot(1:12, precd[1,], type = "l", col = "blue",axes = FALSE, xlab = "", ylab = "",ylim=c(min(precd),max(precd)),lty=2,lwd=4)
points(1:12, precd[2,], type = "l", col = "red",lty=2,lwd=4)
#points(1:12, precd[3,], type = "l", col = "green",lty=2,lwd=2)

axis(side = 4, at = pretty(precd))
mtext("precipitation", side = 4, line = 3)

dev.off()

pdf("/Volumes/MartinResearch2/Wolf2019/analyses/AtheneNew/climate/climate_SWISS.pdf",width=12,height=6)

plot   (xx, yy, type = "n", xlab = "Months", ylab = "Temperature")

polygon(c(1:12,12:1),c(tmin[1,]/10,rev(tmax[1,]/10)),col = rgb(0,0,1,0.2),border="blue")
#polygon(c(1:12,12:1),c(tmin[2,]/10,rev(tmax[2,]/10)),col = rgb(1,0,0,0.2),border="red")
polygon(c(1:12,12:1),c(tmin[3,]/10,rev(tmax[3,]/10)),col = rgb(0.2,0.2,0.2,0.2), border = "green")

par(new = TRUE)

plot(1:12, precd[1,], type = "l", col = "blue",axes = FALSE, xlab = "", ylab = "",ylim=c(min(precd),max(precd)),lty=2,lwd=4)
#points(1:12, precd[2,], type = "l", col = "red",lty=2,lwd=4)
points(1:12, precd[3,], type = "l", col = "green",lty=2,lwd=2)

axis(side = 4, at = pretty(precd))
mtext("precipitation", side = 4, line = 3)


dev.off()
