############################################################################################
# 2021-07-28
# Xiaoqi Luo
# WES pipeline
############################################################################################
# Terminal:
# genomes=/mnt/d/lxq/Training/WES/GATK/hg38/bwa_index/ GATK=/mnt/d/lxq/Training/WES/GATK/gatk-4.1.7.0/gatk \
# ref=/mnt/d/lxq/Training/WES/GATK/hg38/Homo_sapiens_assembly38.fasta \
# snp=/mnt/d/lxq/Training/WES/GATK/hg38/dbsnp_146.hg38.vcf.gz \
# indel=/mnt/d/lxq/Training/WES/GATK/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz


files=os.listdir(config['workspace']+'/rawfastq/')
SRR = []
for f in files[0:3:2]:
    SRR.append(re.split(r'_', f)[0])
#print(SRR)
rule all:
    input:
        expand(
            config['workspace'] + '/mapped/{sample}_raw.vcf',
            sample=SRR
		)
  
rule quality:
	output:
		fq1=config['workspace'] + '/qc/{sample}_1_trim.fastq.gz',
		fq2=config['workspace'] + '/qc/{sample}_2_trim.fastq.gz',
		json=config['workspace'] + '/qc/{sample}_fastp.json',
		html=config['workspace'] + '/qc/{sample}_fastp.html'
	input:
		fq1=config['workspace'] + '/rawfastq/{sample}_1.fastq.gz',
		fq2=config['workspace'] + '/rawfastq/{sample}_2.fastq.gz'
	log:
		config['workspace'] + '/logs/{sample}_trim.log'
	shell:
		'fastp -i {input.fq1} -I {input.fq2} '
		'-o {output.fq1} -O {output.fq2} '
		'-g -q 5 -u 50 -n 15 '
		'-j {output.json} -h {output.html} 2>{log}'

rule mapping:
	output:
		config['workspace'] + "/mapped/{sample}.bam"
	input:
		fq1=config['workspace'] + '/qc/{sample}_1_trim.fastq.gz',
		fq2=config['workspace'] + '/qc/{sample}_2_trim.fastq.gz'
	log:
		config['workspace'] + '/logs/{sample}_bwa_mem.log'
	params:
		index=config['genomes']+'gatk_hg38', #genomes: D:\lxq\Training\WES\GATK\hg38\bwa_index
		extra=r"-R '@RG\tID:lane1\tPL:illumina\tLB:WGS\tSM:Illumina'",
	shell:
		'bwa mem {params.extra} {params.index} '
		'{input.fq1} {input.fq2} '
		'| samtools sort -o {output} -'

rule markdup:
	output:
		bam=config['workspace'] + '/mapped/{sample}_marked.bam',
		metrics=config['workspace'] + '/mapped/{sample}_marked.metrics.txt'
	input:
		config['workspace'] + '/mapped/{sample}.bam'
	log:
		config['workspace'] + '/logs/{sample}_markdup.log'
	params:
		extra=r' "-Xmx2G -Djava.io.tmpdir=./" ',
		GATK=config['GATK'] 
	shell:
		'{params.GATK} --java-options {params.extra} MarkDuplicates '
		'-I {input} -O {output.bam} -M {output.metrics}'

rule fixmat:
    output:
        config['workspace'] + '/mapped/{sample}_marked_fixed.bam'
    input:
        config['workspace'] + '/mapped/{sample}_marked.bam'
    log:
        config['workspace'] + '/logs/{sample}_fix.log'
    params:
        extra=r' "-Xmx2G -Djava.io.tmpdir=./" ',
        GATK=config['GATK'] 
    shell:
        '{params.GATK} --java-options {params.extra} FixMateInformation '
        '-I {input} -O {output} -SO coordinate'

rule baserecal:
    output:
        config['workspace'] + '/mapped/{sample}_recal.table'
    input:
        config['workspace'] + '/mapped/{sample}_marked_fixed.bam'
    log:
        config['workspace'] + '/logs/{sample}_baserecal.log'
    params:
        extra=r' "-Xmx2G -Djava.io.tmpdir=./" ',
        GATK=config['GATK'],
        ref=config['ref'],
        snp=config['snp'],
        indel=config['indel']
    shell:
        '{params.GATK} --java-options {params.extra} BaseRecalibrator ' 
        '-R {params.ref} --input {input} '
        '--known-sites {params.snp} --known-sites {params.indel} '
        '-O {output} '

rule bqsr:
    output:
        config['workspace'] + '/mapped/{sample}_bqsr.bam'
    input:
        bam=config['workspace'] + '/mapped/{sample}_marked_fixed.bam',
        table=config['workspace'] + '/mapped/{sample}_recal.table'
    log:
        config['workspace'] + '/logs/{sample}_fqsr.log'
    params:
        extra=r' "-Xmx2G -Djava.io.tmpdir=./" ',
        GATK=config['GATK'],
        ref=config['ref']
    shell:
        '{params.GATK} --java-options {params.extra} ApplyBQSR '
        '-R {params.ref} --input {input.bam} '
        '-bqsr {input.table} -O {output}'
    
rule call:
    output:
        config['workspace'] + '/mapped/{sample}_raw.vcf'
    input:
        config['workspace'] + '/mapped/{sample}_bqsr.bam'
    log:
        config['workspace'] + '/logs/{sample}_call.log'
    params:
        extra=r' "-Xmx2G -Djava.io.tmpdir=./" ',
        GATK=config['GATK'],
        ref=config['ref'],
        snp=config['snp']
    shell:
        '{params.GATK} --java-options {params.extra} HaplotypeCaller '
        '-R {params.ref} --input {input} '
        '--dbsnp {params.snp} -O {output}'
