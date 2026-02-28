#! /usr/bin/env nextflow
nextflow.enable.dsl=2

/* 
###################
Parameters
###################
*/

// All parameters are specified in the json file

/* 
###################
Processes
###################
*/

process Nr_raw_reads {

	input:
	path fastq

	output:
	tuple env(sample_id), val('Raw_reads'), env(nr_reads), emit: nr_raw_ch

	script:
	"""
	sample_id=\$(basename $fastq .fastq.gz)
	nr_lines=\$(gunzip -c $fastq | wc -l)
	nr_reads=\$((nr_lines/4))
	"""

}

process FastQC {

container "biocontainers/fastqc:v0.11.9_cv8"
publishDir "${params.outdir}/01_FastQC", mode: 'copy'

input:
path fastq

output:
path "*_fastqc.{zip,html}"

script:
"""
fastqc --extract -t ${params.threads} $fastq
"""

}

process MultiQC {

container "quay.io/biocontainers/multiqc:1.33--pyhdfd78af_0"
publishDir "${params.outdir}/02_MultiQC", mode: 'copy'

input:
path fastqc_reports

output:
path "multiqc_report.html"

script:
"""
multiqc $fastqc_reports
"""

}

process Trimmomatic {

container "biocontainers/trimmomatic:v0.38dfsg-1-deb_cv1"
publishDir "${params.outdir}/03_Trimmed", mode: 'copy'

input:
path fastq
path adapters

output:
path "*.fastq.gz", emit: trimming_ch

script:
"""
TrimmomaticSE -phred${params.phred} -threads ${params.threads} \
${fastq} \
"trimmed_${fastq}" \
ILLUMINACLIP:${adapters}:${params.adapters_params[0]}:\
${params.adapters_params[1]}:${params.adapters_params[2]} \
LEADING:${params.leading} TRAILING:${params.trailing} \
SLIDINGWINDOW:${params.slidingwindow[0]}:${params.slidingwindow[1]} \
MINLEN:${params.minlen}
"""

}

process Nr_trimmed_reads {

	input:
	path trimmed_fastq

	output:
	tuple env(sample_id), val('Trimmed_reads'), env(nr_reads), emit: nr_trimmed_ch

	script:
	"""
	sample_id=\$(basename $trimmed_fastq .fastq.gz)
	nr_lines=\$(gunzip -c $trimmed_fastq | wc -l)
	nr_reads=\$((nr_lines/4))
	"""
}

process FastQC2 {

container "biocontainers/fastqc:v0.11.9_cv8"
publishDir "${params.outdir}/04_FastQC_Trimmed", mode: 'copy'

input:
path trimmed_fastq

output:
path "*_fastqc.{zip,html}"

script:
"""
fastqc --extract -t ${params.threads} $trimmed_fastq
"""

}

process MultiQC2 {

container "quay.io/biocontainers/multiqc:1.33--pyhdfd78af_0"
publishDir "${params.outdir}/05_MultiQC_Trimmed", mode: 'copy'

input:
path fastqc_reports

output:
path "multiqc_report.html"

script:
"""
multiqc $fastqc_reports
"""

}

process Bismark_Align {

container "quay.io/biocontainers/bismark:0.25.1--hdfd78af_0"
publishDir "${params.outdir}/06_Bismark_Aligned", mode: 'copy'

input:
path trimmed_fastq
path genome

output:
path "*.bam", emit: Aligned_ch
path "*.txt"

script:
"""
bismark --genome $genome $trimmed_fastq \
--multicore ${params.threads} --phred${params.phred}-quals \
--bowtie2 -N ${params.seedmss} -L ${params.seedlen} \
--score_min L,${params.intercept_slope[0]},${params.intercept_slope[1]}
"""

}

process Samtools {

	container "quay.io/biocontainers/samtools:1.23--h96c455f_0"
	publishDir "${params.outdir}/07_Samtools", mode: 'copy'

	input:
	path bam_file
 
	output:
	path "*_filtered.bam", emit: Samtools_ch
	path "*_filtered_sorted.bam"
	path "*_filtered_sorted.bam.bai"

	script:
	"""
	samtools view -bq ${params.mapq} $bam_file > ${bam_file.simpleName}_filtered.bam
	samtools sort ${bam_file.simpleName}_filtered.bam -o ${bam_file.simpleName}_filtered_sorted.bam
	samtools index ${bam_file.simpleName}_filtered_sorted.bam
	"""
}

process Bismark_Meth_Extractor {

container "quay.io/biocontainers/bismark:0.25.1--hdfd78af_0"
publishDir "${params.outdir}/08_Bismark_Methylation_Extractor", mode: 'copy'

input:
path sorted_bam

output:
path "*.txt"

script:
"""
bismark_methylation_extractor $sorted_bam --comprehensive --multicore ${params.threads}
"""

}

/*
###################
Workflow
###################
*/

workflow {

// Define channels
fastq_ch = Channel.fromPath("${params.input}/*.fastq.gz")

// Calculate nr of raw reads and view in terminal
Nr_raw_reads(fastq_ch)
Nr_raw_reads.out.nr_raw_ch.view()

// Run FastQC
FastQC(fastq_ch)

// Collect output from FastQC and run MultiQC
MultiQC(FastQC.out.collect())

// Run trimmomatic
Trimmomatic(fastq_ch, params.adapters)

// Count nr of reads
Nr_trimmed_reads(Trimmomatic.out.trimming_ch)

// Reformat channel and view in terminal
Nr_trimmed_reads.out.nr_trimmed_ch.view()

// Rerun FastQC
FastQC2(Trimmomatic.out.trimming_ch)

// Rerun MultiQC post-trim
MultiQC2(FastQC2.out.collect())

// Run Bismark alignment
Bismark_Align(Trimmomatic.out.trimming_ch,params.genome)

// Run Samtools 
Samtools(Bismark_Align.out.Aligned_ch)

// Run Methylation Extractor
Bismark_Meth_Extractor(Samtools.out.Samtools_ch)

}