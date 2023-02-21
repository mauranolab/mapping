# mapping

This contains Maurano Lab sequencing processing pipelines:
* dnase - bwa mapping piepline for DNase-seq, ChIP-seq, or resequencing (including bcftools variant calling)
* dnase/bamintersect - analyzes multiple parallel mappings of the same reads to identify integration junctions
* dnase/trackhub - track hub creation for loading in UCSC Genome Browser
* flowcells - LIMS integration and demuxing support
* transposon - processing pipeline for amplicon sequencing of barcoded reporter libraries

These pipelines are not intended to be immediately portable but rather are intended to serve as reference for how the lab processes data.
