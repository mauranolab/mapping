#!/bin/bash
set -eu -o pipefail

#Boilerplate header
cat /vol/mauranolab/flowcells/src/SampleSheet.template.txt > SampleSheet.csv

#Now parse the sequencing sheet info from STDIN
awk -F "\t" 'BEGIN {OFS="\t"; parse=0} {print} $0!~/^#/ && parse==0 {parse=1; print "#Sample Name", "Sample #", "Sub-library", "Lab", "Made By", "Sample Type", "Barcode 1 (i7)", "Barcode 2 (i5)", "R1 Trim (P5)", "R2 Trim (P7)", "Sequencing primer R1", "Indexing primer BC1 (i7)", "Indexing primer BC2 (i5)", "Sequencing primer R2", "Library concentration (pM)", "", "Request Type", "Requested reads (M)", "Read format", "Scale factor", "Relative representation", "Amount put on FC (uL)", "Sequenced reads", "Actual representation"}' |
#Also creates info.txt
tee info.txt |
awk -F "\t" 'BEGIN {OFS=","; split("8,8", bclens, ",")} $1=="#Indices" && $2!="" {split($2, bclens, ",")} $0!~/^#/ && $1!="" {split($7, bc1, "_"); split($8, bc2, "_"); print "Sample_" $2 $3, $2 $3, "", "",  bc1[1], toupper(substr(bc1[2], 0, bclens[1])),  bc2[1], toupper(substr(bc2[2], 0, bclens[2])), "Project_" $4, "";}' >> SampleSheet.csv
