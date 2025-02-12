#!/bin/bash

#Wrapper around Hotspot 1

#BUGBUG can't handle % in file names


if [ "$#" -ne 6 ]; then
    echo "Wrong number of arguments"
    exit 1
fi


bam=${1}
dens=${2}
outdir=${3}
mappedgenome=${4}
mappableFile=${5}
chromsizes=${6}

date
echo "Running hotspot"
echo "bam=$bam"
echo "dens=$dens"
echo "outdir=$outdir"
echo "mappedgenome=$mappedgenome"
echo "mappableFile=$mappableFile"
echo "chromsizes=${chromsizes}"

name=`basename ${bam} | perl -pe 's/\.bam$//g;'`


chromFile=$TMPDIR/hotspots.${mappedgenome}.chromInfo.bed
awk 'BEGIN {OFS="\t"} {print $1, 0, $2, $1}' ${chromsizes} | grep -v 'alt\|random\|Un\|hap\|scaffold' | sort-bed - > ${chromFile} 

if [ ! -e "${mappableFile}" ]; then
    echo "WARNING can not find ${mappableFile}"
    mappableFile="${chromFile}"
fi

omit_satellite="/gpfs/data/isg/annotation/bed/${mappedgenome}/repeat_masker/Satellite.bed"
if [ ! -e "${omit_satellite}" ]; then
    echo "WARNING can not find ${omit_satellite}"
    omit_satellite=""
fi


if grep -q chrX ${chromsizes}; then
    check="T"
else
    check="F"
fi

HOTSPOT_DISTR=/cm/shared/apps/hotspot/4.1/hotspot-master/hotspot-distr/


cat <<EOF > runall.tokens.${mappedgenome}.txt
[script-tokenizer]

#######################################
## Notes:  If duplicate token definitions exist, the last definition in
##      the file will be used. Tokens can use .ini variables and declarations.
##      See http://docs.python.org/library/configparser.html
#######################################

#######################################
# Global tokens (used by most scripts)
#######################################

## Tags file in bam format (file extension .bam), or starched bed file
## (file extension .bed.starch).  If the tags file is in bed.starch
## format, and if it is in directory specified by _OUTDIR_, then it
## needs to be in the format
##
##     chr  5'start  5'start+1
##
## That is, the file should be three column bed (no ID field, etc.)
## containing the 1bp coordinates of the 5' ends of the tags. If
## _TAGS_ is a bam file, or a bed.starch file not in _OUTDIR_, then
## the bed.starch file in the above format will be generated (using
## the strand column if present), and put in _OUTDIR_.  NOTE: if you
## use run_badspot, then you must use bam files, or a bed.starch file
## not in _OUTDIR_.

_HOTSPOT_DISTR_ = $HOTSPOT_DISTR

_TAGS_ = $bam

## For ChIP data with an Input bam file, set _USE_INPUT_ to T, and set
## _INPUT_TAGS_ to the name of that bam file.

_USE_INPUT_ = F 
_INPUT_TAGS_ = 

## Genome 
_GENOME_ = $mappedgenome
## Tag length
_K_ = $Kreads
## Chromosome coordinates, bed format.
#$HOTSPOT_DISTR/hotspot-deploy/bin/writeChromInfoBed.pl /gpfs/data/isg/annotation/fasta/${mappedgenome}/${mappedgenome}all.fa && cat $TMPDIR/chromInfo.bed | grep -v hap | grep -v random | grep -v chrUn | grep -v alt | grep -v scaffold > $TMPDIR/hotspots.${mappedgenome}.chromInfo.bed && rm -f $TMPDIR/chromInfo.bed
_CHROM_FILE_ = $chromFile
## Location of uniquely mappable positions in the genome for this tag length.
_MAPPABLE_FILE_ = $mappableFile
#_MAPPABLE_FILE_ = /gpfs/data/isg/annotation/bed/%(_GENOME_)s/mappability/%(_GENOME_)s.K%(_K_)s.mappable_only.bed

## Set DUPOK to T for DNaseI data, F for ChIP-seq data (DUPOK = T means allow duplicate reads)
_DUPOK_ = T

## FDR levels, separated by spaces if more than one. Set to N if you
## do not want FDR thresholding (for example, if you just want SPOT
## score computed.)  
## _FDRS_ = "N"
_FDRS_ = "0.005 0.01 0.05 0"

## Tag density, 150bp window, sliding every 20bp, used for
## peak-finding.  Will be generated, based on the _TAGS_ file, if it
## does not exist. Assumed to be starched bed file, extension
## bed.starch.  Can be blank, in which case the density file will be
## generated in _OUTDIR_, and the the name will be the name of the
## tags file, minus the bam or bed.starch extension, with the added
## extension tagdensity.bed.starch.  
_DENS_ = $dens

## Output directories (can all be the same location).  Use full path names.
## _OUTDIR_ contains tags files in converted bed.starch and lib.txt formats (for hotspot
## program), and hotspot and peak results.
## _RANDIR_ contains generated random tags (for FDR thresholding) and hotspots called on random tags.
_OUTDIR_ = $outdir
_RANDIR_ = $outdir

## If there are any regions from which tags should be automatically
## omitted, include those here (only if you use run_badspot). May be
## left blank.
_OMIT_REGIONS_: "${omit_satellite} /gpfs/data/mauranolab/mapped/src/dnase/hotspots.omit.chrM.bed"

## Set to T if you want scripts to skip steps that have already been done.
_CHECK_ = ${check}

## If _CHECK_ = T, outputs are checked for completeness by searching
## for results for the following chromsome.
_CHKCHR_ = chrX

## Hotspot program binary
_HOTSPOT_ = %(_HOTSPOT_DISTR_)s/hotspot-deploy/bin/hotspot

## Clean up. Remove all intermediate files and directories if set to T.  See
## pipeline script run_final.
_CLEAN_ = T

## Peak-finding program.
_PKFIND_BIN_ = %(_HOTSPOT_DISTR_)s/hotspot-deploy/bin/wavePeaks
## Peak-finding smoothing level. If the resolution of the input file
## is x, then the results are smoothed out to a scale of (2^level)*x.
_PKFIND_SMTH_LVL_ = 3

## Random number seed, used for generating random tags for FDR thresholding.
_SEED_=101

## Hotspot program parameters
_THRESH_ = 2
_WIN_MIN_ = 200
_WIN_MAX_ = 300
_WIN_INCR_ = 50
_BACKGRD_WIN_ = 50000
_MERGE_DIST_ = 150
_MINSIZE_ = 10
EOF


#BUGBUG not customized for nonref genomes


scriptTokBin=$HOTSPOT_DISTR/ScriptTokenizer/src/script-tokenizer.py
pipeDir=$HOTSPOT_DISTR/pipeline-scripts
tokenFile=runall.tokens.${mappedgenome}.txt

## Do everything, including badspots and final cleanup
scripts="$pipeDir/run_badspot
    $pipeDir/run_make_lib
    $pipeDir/run_wavelet_peak_finding
    $pipeDir/run_10kb_counts
    $pipeDir/run_generate_random_lib
    $pipeDir/run_pass1_hotspot
    $pipeDir/run_pass1_merge_and_thresh_hotspots
    $pipeDir/run_pass2_hotspot
    $pipeDir/run_rescore_hotspot_passes
    $pipeDir/run_spot
    $pipeDir/run_thresh_hot.R
    $pipeDir/run_both-passes_merge_and_thresh_hotspots
    $pipeDir/run_add_peaks_per_hotspot
    $pipeDir/run_final"

$scriptTokBin \
    --clobber \
    --output-dir=`pwd`/src.${mappedgenome} \
    $tokenFile \
    $scripts

for script in $scripts
do
    src.${mappedgenome}/$(basename $script).tok
done


#BUGBUG bob's run_final cleans up pretty well but doesn't leave any peak file with scores
rm -f ${outdir}/${name}.bed.starch

#Output files are .bed instead of starch
for f in ${outdir}/${name}-final/*.bed; do
    base=`basename ${f} .bed`
    starch ${f} > ${outdir}/${name}-final/${base}.starch
    rm -f ${f}
done

#Output files are .txt instead of txt.gz
for f in ${outdir}/${name}-final/*.txt; do 
    #Getting intermittent pigz: abort: write error on >file> (Inappropriate ioctl for device). Tried switching to gzip to see if that fixes it.
    gzip -9 ${NSLOTS} ${f}
done


echo
date
echo "Done!"
