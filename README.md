# whole_exome_squence_snakemake

## Conda Enviornment
```
conda install -c bioconda -c conda-forge bowtie sra-tools samtools bcftools vcftools snpeff multiqc qualimap fastp bowtie2 bwa bowtie bedtools snakemake
```

## GATK Download
```
wget  https://github.com/broadinstitute/gatk/releases/download/4.0.6.0/gatk-4.0.6.0.zip
unzip gatk-4.0.6.0.zip
cd /GATK/gatk-4.0.6.0
./gatk --help
```

## Files needed for analysis
### Homo_sapiens_assembly38
```
wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.fasta.gz  
wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.fasta.fai  
unzip Homo_sapiens_assembly38.fasta.gz
```

### dbsnp
```
wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/dbsnp_146.hg38.vcf.gz
wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/dbsnp_146.hg38.vcf.gz.tbi
```
### indels
```
wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz 
wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi
```

# How To Execute
## File Tree (fastq files are saved in rawfastq file)
>WES_snakemake.py<br>
>rawfastq
>>SRR13018671_1.fastq.gz<br>
>>SRR13018671_2.fastq.gz<br>
>>SRR10720394_1.fastq.gz<br>
>>SRR10720394_2.fastq.gz<br>

```
genomes=path-to-bwa-index GATK=path-to-gatk \
ref=path-to-Homo_sapiens_assembly38.fasta \
snp=path-to-dbsnp_146.hg38.vcf.gz \
indel=path-to-path-to-Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
```
# example
```
genomes=/mnt/d/WES/GATK/hg38/bwa_index/ GATK=/mnt/d/WES/GATK/gatk-4.1.7.0/gatk \
ref=/mnt/d/WES/GATK/hg38/Homo_sapiens_assembly38.fasta \
snp=/mnt/d/WES/GATK/hg38/dbsnp_146.hg38.vcf.gz \
indel=/mnt/d/WES/GATK/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
```
# result
![image](https://github.com/XiaoqiLuo/whole_exome_squence_snakemake/blob/main/%E5%BE%AE%E4%BF%A1%E6%88%AA%E5%9B%BE_20210806201251.png)
