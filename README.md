# BSS-Pipeline-nf
Reproducible Nextflow workflow that takes raw bisulfite sequencing FASTQ files and outputs genome-wide methylation calls.

# ⚙️ Workflow
<img width="721" height="436" alt="BSS-pipeline drawio (1)" src="https://github.com/user-attachments/assets/64dee5ea-07c8-4bd9-900d-2b34fe91f926" />

# ❗Dependencies 
The pipeline runs in Docker containers, requiring only Nextflow and Docker to be installed.

# ▶️ Running the script
Before running the pipeline, a few input parameters must be specified.

**User-specific parameters**    
These parameters should be updated in the JSON configuration file according to the user's environment:
- Input: Path to the folder containing only raw .fastq.gz files.
- Outdir: Path to the directory where output files will be written.
- Genome: Path to the folder containing the Bismark genome index (this folder should include the Bisulfite_Genome directory and the reference .fa file(s) (as created by `bismark_genome_preparation`).
- Adapters: Path to the file containing adapter sequences to be trimmed from the reads.

**Optional Settings**
- Threads: Number of threads to use.
- Phred: Specify the encoding of the user's FASTQ files - phred-33 or phred-64.
- Leading: Quality score threshold for trimming bases from the leading end of the reads.
- Trailing: Quality score threshold for trimming bases from the trailing end of the reads.
- Slidingwindow: [window length, average quality cut-off] for Trimmomatic.
- Minlen: Minimum read length to retain after trimming.
- Adapters_params: Adapter clip options [Seed Mismatches, Palindrome Clip Threshold, Simple Clip Threshold].
- Seedmms: Maximum number of mismatches allowed in the Bowtie seed (for Bismark).
- Seedlen: Seed length used for Bowtie alignment (for Bismark).
- Intercept_slope: Intercept and slope for the Bowtie minimum score function [intercept, slope].
- Mapq: Minimum MAPQ to retain after Samtools filtering.

**Running the script**    
To run the pipeline in Bash:
1. Make sure the .json and .config files are in the same directory as the .nf file.
2. There are separate Nextflow scripts for single-end (SE) and paired-end (PE) data — use the appropriate one.

Example:
`nextflow run BSS-pipelineSE.nf -params-file BSS-pipeline-params.json`

# ✅ Versioning
This pipeline uses Docker containers that include specific versions of the tools:
- FastQC v0.11.9
- MultiQC v1.33
- Trimmomatic v0.38
- Bismark v0.25.1
- Samtools v1.23

The pipeline itself was developed and tested with:
- Nextflow v25.10.2
- Docker v29.2.1
