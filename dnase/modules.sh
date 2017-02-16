#!/bin/bash
set -e -o pipefail

echo "Loading modules"

module load picard/1.140
module load FastQC/0.11.4
module load bedops/2.4.19
module load bwa/0.7.15
module load htslib/1.2.1
module load samtools/1.2
module load trimmomatic/0.33
module load python/3.5.0
module load hotspot/4.1
module load samblaster/0.1.22
