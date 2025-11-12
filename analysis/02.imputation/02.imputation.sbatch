#!/bin/bash
#
#SBATCH --job-name=impute
#SBATCH --output=impute.%j.out
#SBATCH --error=impute.%j.err
#SBATCH -t 24:00:00
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node 24
#SBATCH --mem=90G
#SBATCH --mail-type=END,FAIL
#SBATCH  --mail-user=ericacnr@colostate.edu

#echo commands to stdout

set -x
source ~/.bashrc

BEAGLE_JAR="/projects/ericacnr@colostate.edu/mambaforge/envs/bioinf/java/beagle.27Jan18.7e1.jar"
vcf="/scratch/alpine/ericacnr@colostate.edu/BANS/01.fastq_processing/mega-non-model-wgs-snakeflow/results/bqsr-round-0/vcf_final/BANS.all.ds6x.pass-maf-0.05.SNP.miss.8.vcf.gz"

java -Xss5m -Xmx90g -jar $BEAGLE_JAR \
  gt=$vcf \
  out=BANS.all.ds6x.pass-maf-0.05.SNP.miss.8_imputed4.1 \
  dr2=TRUE \
  nthreads=48
