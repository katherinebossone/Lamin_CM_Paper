#!/bin/bash1

module load bedtools/ samtools/ bwa/

conda activate pairtools-1.1.0

## LB1 het rep 1 ##
bash docker-4dn-hic/scripts/run-bwa-mem.sh rawdata/sample4/8450_2_P0_5_LB1Fd_R1.fastq rawdata/sample4/8450_2_P0_5_LB1Fd_R2.fastq /mnt/sequence/home/bossone/4DNFI823LSI8.bwaIndex.tgz . 8450_2 16

bash docker-4dn-hic/scripts/run-pairsam-parse-sort.sh /mnt/sequence/home/bossone/HiC/8450_2.bam /mnt/sequence/home/bossone/HiC/mm10.sizes . 8450_2 16 lz4c

bash docker-4dn-hic/scripts/run-pairsam-merge.sh 8450_2 16 8450_2.sam.pairs.gz

bash docker-4dn-hic/scripts/run-pairsam-markasdup.sh 8450_2.merged.sam.pairs.gz 8450_2

bash docker-4dn-hic/scripts/run-pairsam-filter.sh 8450_2.marked.sam.pairs.gz 8450_2 mm10.sizes

bash docker-4dn-hic/scripts/run-merge-pairs.sh 8450_2 8450_2.marked.sam.pairs.gz
bash docker-4dn-hic/scripts/run-addfrag2pairs.sh -i 8450_2.pairs.gz -r mm10_DpnII.txt -o 8450_2

bash docker-4dn-hic/scripts/run-cooler.sh 8450_2.pairs.gz mm10.sizes 1000 1 8450_2 2 -u 1000,2000,5000,10000,25000,50000,100000,250000,500000,1000000,2500000,5000000,10000000

java -jar /mnt/sequence/home/bossone/HiC/docker-4dn-hic/juicer_tools_1.22.01.jar pre -n 8450_2.ff.pairs.gz 8450_2_prefix.hic mm10.sizes -r 1000,2000,5000,10000,25000,50000,100000,250000,500000,1000000,2500000,5000000,10000000 -q 0

bash docker-4dn-hic/scripts/run-cool2multirescool.sh -i 8450_2.cool -p 1 -o 8450_2 -c 10000000 -u 1000,2000,5000,10000,25000,50000,100000,250000,500000,1000000,2500000,5000000,10000000

bash docker-4dn-hic/scripts/run-add-hicnormvector-to-mcool.sh 8450_2_prefix.hic 8450_2.multires.cool .

## LB1 het rep 2 ##
bash docker-4dn-hic/scripts/run-bwa-mem.sh rawdata/sample4/8450_4_P0_5_LB1Fd_R1.fastq rawdata/sample4/8450_4_P0_5_LB1Fd_R2.fastq /mnt/sequence/home/bossone/4DNFI823LSI8.bwaIndex.tgz . 8450_4 16

bash docker-4dn-hic/scripts/run-pairsam-parse-sort.sh /mnt/sequence/home/bossone/HiC/8450_4.bam /mnt/sequence/home/bossone/HiC/mm10.sizes . 8450_4 16 lz4c

bash docker-4dn-hic/scripts/run-pairsam-merge.sh 8450_4 16 8450_4.sam.pairs.gz

bash docker-4dn-hic/scripts/run-pairsam-markasdup.sh 8450_4.merged.sam.pairs.gz 8450_4

bash docker-4dn-hic/scripts/run-pairsam-filter.sh 8450_4.marked.sam.pairs.gz 8450_4 mm10.sizes

bash docker-4dn-hic/scripts/run-merge-pairs.sh 8450_4 8450_4.marked.sam.pairs.gz
bash docker-4dn-hic/scripts/run-addfrag2pairs.sh -i 8450_4.pairs.gz -r mm10_DpnII.txt -o 8450_4

bash docker-4dn-hic/scripts/run-cooler.sh 8450_4.pairs.gz mm10.sizes 1000 1 8450_4 2 -u 1000,2000,5000,10000,25000,50000,100000,250000,500000,1000000,2500000,5000000,10000000

java -jar /mnt/sequence/home/bossone/HiC/docker-4dn-hic/juicer_tools_1.22.01.jar pre -n 8450_4.ff.pairs.gz 8450_4_prefix.hic mm10.sizes -r 1000,2000,5000,10000,25000,50000,100000,250000,500000,1000000,2500000,5000000,10000000 -q 0

bash docker-4dn-hic/scripts/run-cool2multirescool.sh -i 8450_4.cool -p 1 -o 8450_4 -c 10000000 -u 1000,2000,5000,10000,25000,50000,100000,250000,500000,1000000,2500000,5000000,10000000

bash docker-4dn-hic/scripts/run-add-hicnormvector-to-mcool.sh 8450_4_prefix.hic 8450_4.multires.cool .

## LB1 mutant rep 1 ##
bash docker-4dn-hic/scripts/run-bwa-mem.sh rawdata/sample4/8450_2_P0_5_LB1Fd_R1.fastq rawdata/sample4/8450_2_P0_5_LB1Fd_R2.fastq /mnt/sequence/home/bossone/4DNFI823LSI8.bwaIndex.tgz . 8450_2 16

bash docker-4dn-hic/scripts/run-pairsam-parse-sort.sh /mnt/sequence/home/bossone/HiC/8450_2.bam /mnt/sequence/home/bossone/HiC/mm10.sizes . 8450_2 16 lz4c

bash docker-4dn-hic/scripts/run-pairsam-merge.sh 8450_2 16 8450_2.sam.pairs.gz

bash docker-4dn-hic/scripts/run-pairsam-markasdup.sh 8450_2.merged.sam.pairs.gz 8450_2

bash docker-4dn-hic/scripts/run-pairsam-filter.sh 8450_2.marked.sam.pairs.gz 8450_2 mm10.sizes

bash docker-4dn-hic/scripts/run-merge-pairs.sh 8450_2 8450_2.marked.sam.pairs.gz
bash docker-4dn-hic/scripts/run-addfrag2pairs.sh -i 8450_2.pairs.gz -r mm10_DpnII.txt -o 8450_2

bash docker-4dn-hic/scripts/run-cooler.sh 8450_2.pairs.gz mm10.sizes 1000 1 8450_2 2 -u 1000,2000,5000,10000,25000,50000,100000,250000,500000,1000000,2500000,5000000,10000000

java -jar /mnt/sequence/home/bossone/HiC/docker-4dn-hic/juicer_tools_1.22.01.jar pre -n 8450_2.ff.pairs.gz 8450_2_prefix.hic mm10.sizes -r 1000,2000,5000,10000,25000,50000,100000,250000,500000,1000000,2500000,5000000,10000000 -q 0

bash docker-4dn-hic/scripts/run-cool2multirescool.sh -i 8450_2.cool -p 1 -o 8450_2 -c 10000000 -u 1000,2000,5000,10000,25000,50000,100000,250000,500000,1000000,2500000,5000000,10000000

bash docker-4dn-hic/scripts/run-add-hicnormvector-to-mcool.sh 8450_2_prefix.hic 8450_2.multires.cool .

## LB1 mutant rep 2 ##
bash docker-4dn-hic/scripts/run-bwa-mem.sh rawdata/sample4/8450_4_P0_5_LB1Fd_R1.fastq rawdata/sample4/8450_4_P0_5_LB1Fd_R2.fastq /mnt/sequence/home/bossone/4DNFI823LSI8.bwaIndex.tgz . 8450_4 16

bash docker-4dn-hic/scripts/run-pairsam-parse-sort.sh /mnt/sequence/home/bossone/HiC/8450_4.bam /mnt/sequence/home/bossone/HiC/mm10.sizes . 8450_4 16 lz4c

bash docker-4dn-hic/scripts/run-pairsam-merge.sh 8450_4 16 8450_4.sam.pairs.gz

bash docker-4dn-hic/scripts/run-pairsam-markasdup.sh 8450_4.merged.sam.pairs.gz 8450_4

bash docker-4dn-hic/scripts/run-pairsam-filter.sh 8450_4.marked.sam.pairs.gz 8450_4 mm10.sizes

bash docker-4dn-hic/scripts/run-merge-pairs.sh 8450_4 8450_4.marked.sam.pairs.gz
bash docker-4dn-hic/scripts/run-addfrag2pairs.sh -i 8450_4.pairs.gz -r mm10_DpnII.txt -o 8450_4

bash docker-4dn-hic/scripts/run-cooler.sh 8450_4.pairs.gz mm10.sizes 1000 1 8450_4 2 -u 1000,2000,5000,10000,25000,50000,100000,250000,500000,1000000,2500000,5000000,10000000

java -jar /mnt/sequence/home/bossone/HiC/docker-4dn-hic/juicer_tools_1.22.01.jar pre -n 8450_4.ff.pairs.gz 8450_4_prefix.hic mm10.sizes -r 1000,2000,5000,10000,25000,50000,100000,250000,500000,1000000,2500000,5000000,10000000 -q 0

bash docker-4dn-hic/scripts/run-cool2multirescool.sh -i 8450_4.cool -p 1 -o 8450_4 -c 10000000 -u 1000,2000,5000,10000,25000,50000,100000,250000,500000,1000000,2500000,5000000,10000000

bash docker-4dn-hic/scripts/run-add-hicnormvector-to-mcool.sh 8450_4_prefix.hic 8450_4.multires.cool .


# deactivate environment
conda deactivate
