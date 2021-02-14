# mapping

This contains Maurano Lab sequencing processing pipelines:
* flowcells - LIMS integration and demuxing support
* dnase - bwa mapping piepline for DNase-seq, ChIP-seq, or resequencing (including bcftools variant calling)
* dnase/bamintersect - analyzes multiple parallel mappings of the same reads to identify integration junctions
* dnase/trackhub - track hub creation for loading in UCSC Genome Browser

These pipelines are not intended to be immediately portable but rather are intended to serve as reference for how the lab processes data.
