# Maurano Lab bwa NGS processing pipeline
## Usage
submit.sh genomesToMap analysisType samplePrefix BS
1. genomesToMap - comma-separated list of reference genomes
2. analysisType - broken into two verbs, separated by a comma. Must be chosen from this list:
  * processingCommand: aggregate, aggregateRemarkDups, mapBwaAln, mapBwaMem, bamintersect, none
  * sampleType: atac, dnase, chipseq, dna, capture, none
3. samplePrefix - prefix to be appended to output, possible including directory tree and filename prefix
4. BS - sample identifier used to identify fastq files and name output
5. sampleAnnotation - semicolon-separated list of key-value pairs. Enclose the whole thing in double quotes if metadata includes spaces.

requires inputs.txt in current directory to contain a list of fastq file locations

## Dependencies:
1. Major module dependencies listed in submit.sh
2. SGE-based grid submission system (we use with SLURM using [sge2slurm](https://github.com/mauranolab/sge2slurm))
3. genomeinfo.sh containes locations of references, parameters, etc. per supported genome
