#!/bin/bash

#set -e -o pipefail

module load bedops

echo "Starting analysis" 
date 
#first input should be DNA (sample.barcode.counts.txt file)
DNA=$1
#Second input should be RNA (sample.barcode.counts.txt file)
RNA=$2
#Third input should be iPCR (sample.barcodes.coords.bed file)
iPCR=$3
#Fourth input should be output folder
OUTDIR=$4
#Tmp output
SampleName=`basename $OUTDIR`
#TMPDIR=/tmp
TMPOUT=$TMPDIR/${SampleName}.tmpout.txt


echo "Analyzing $SampleName"
echo "Input DNA $(basename $DNA)"
echo "Input RNA $(basename $RNA)"
echo "Input iPCR $(basename $iPCR)"


mkdir -p $OUTDIR
#Join all barcodes together and include NAs for each sample 
#Print all Barcodes to the AllBC.txt file
join -eNA -a1 -a2 $DNA -o 0 1.2 2.2 $RNA |
join -eNA -a1 -a2 -1 4 <(sort -k4,4 ${iPCR}) -o 1.1 1.2 1.3 0 1.6 2.2 2.3 1.5 -2 1 - |awk -F' ' 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5, $6, $7, $8, $9}' > $TMPOUT
header="chrom\tchromStart\tchromEnd\tBC\tstrand\tDNA\tRNA\tiPCR"
echo -e $header | cat - $TMPOUT |awk -F "\t" 'BEGIN {OFS="\t"} {print $4, $6, $7, $8}'> $OUTDIR/AllBCs.txt


#Remove iPCR insertions with more than one location
cat $iPCR| awk -F "\t" 'BEGIN {OFS="\t"} {print $4}'|sort |uniq -c|sort -nk1|awk '$1==1 {print $2}' > $TMPDIR/${SampleName}.singleIns.txt

#Create bed file of all inserted sites, including NAs
echo "Creating bed-file and removing barcodes with multiple insertions sites"
header="chrom\tchromStart\tchromEnd\tBC\tvalue\tstrand\tDNA\tRNA\tiPCR"
awk  -F "\t" 'BEGIN {OFS="\t"} $1!="NA" && $1!="chrY" && $1!="chrM" {print $1, $2, $3, $4, 0, $5, $6, $7, $8}' $TMPOUT|
#Retains barcodes from $TMPOUT that are present in $TMPDIR/${SampleName}.singleIns.txt
awk -F "\t" 'BEGIN {OFS="\t"} NR==FNR{a[$1];next} ($4) in a' $TMPDIR/${SampleName}.singleIns.txt - |
sort-bed - > $TMPOUT.new 
mv $TMPOUT.new $TMPOUT

#Sample name
echo "Adding sample name"
header="$header\tSampleName"
awk -F "\t" -v sampleName="$SampleName" 'BEGIN {OFS="\t"} {print $0, sampleName}' $TMPOUT > $TMPOUT.new
mv $TMPOUT.new $TMPOUT

###Add annotation
#Distance to TSS
echo 'doing distance to nearest TSS'
header="$header\tDistToTSS\tNearestRefseqGeneTSS\tNearestRefseqGeneStrand"
closest-features --delim '\t' --closest --dist --no-ref $TMPOUT /vol/isg/annotation/bed/hg38/refseq_gene/refGene.CombinedTxStarts.bed | 
awk -F "\t" 'BEGIN {OFS="\t"} {print $NF, $(NF-3), $(NF-1)}'|paste $TMPOUT - > $TMPOUT.new
mv $TMPOUT.new $TMPOUT


#echo 'doing nearest refseq TSS upstream'
#header="$header\tDistToTSS-up\tNearestRefseqGeneTSS-up\tNearestRefseqGeneStrand-up"
#closest-features  --delim '\t' --dist $TMPOUT /vol/isg/annotation/bed/hg38/refseq_gene/refGene.CombinedTxStarts.bed | awk -F "\t" 'BEGIN {OFS="\t"} {if ($6=="-") print $NF, $(NF-3), $(NF-1); else if  ($6=="+") print $(NF-7), $(NF-10), $(NF-8)}' |paste $TMPOUT - > $TMPOUT.new 
#mv $TMPOUT.new $TMPOUT
#
#echo 'doing nearest refseq TSS downstream'
#header="$header\tDistToTSS-dn\tNearestRefseqGeneTSS-dn\tNearestRefseqGeneStrand-dn"
#closest-features  --delim '\t' --dist $TMPOUT /vol/isg/annotation/bed/hg38/refseq_gene/refGene.CombinedTxStarts.bed | awk -F "\t" 'BEGIN {OFS="\t"} {if ($6=="-")  print $(NF-7), $(NF-10), $(NF-8); else if  ($6=="+") print $NF, $(NF-3), $(NF-1)}' |paste $TMPOUT - > $TMPOUT.new 
#mv $TMPOUT.new $TMPOUT
#

#Distance to DNase
echo 'doing distance to nearest DHS '
header="$header\tDistToNearestDHS"
sort-bed $TMPOUT |closest-features  --closest --dist --no-ref - /vol/isg/encode/dnase/mapped/K562-DS9764/hotspots/K562-DS9764.hg38_noalt-final/K562-DS9764.hg38_noalt.fdr0.01.pks.starch |
awk -F'|' 'BEGIN {OFS="\t"} {print $2}' |paste $TMPOUT - > $TMPOUT.new
mv $TMPOUT.new $TMPOUT

echo 'MCV of nearest DHS'
bedops -u /home/mauram01/scratch/hybridmice/dnase/lineagemcv/all.lineage.dhs.unmerged.starch <(awk 'BEGIN {OFS="\t"} {print $1, $2, $3, "K562-DS9764"}' /vol/isg/encode/dnase/mapped/K562-DS9764/hotspots/K562-DS9764.hg38_noalt-final/K562-DS9764.hg38_noalt.fdr0.01.pks.starch)|
awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4}'| sort-bed - | 
bedmap  --delim '\t' --skip-unmapped --echo --fraction-ref 1.0 --count  /vol/isg/encode/dnase/mapped/K562-DS9764/hotspots/K562-DS9764.hg38_noalt-final/K562-DS9764.hg38_noalt.fdr0.01.pks.starch - | 
awk 'BEGIN {OFS="\t"} {print $1, $2, $3, NR, $4}' > $TMPDIR/DHS_MCV.bed
header="$header\tMCVNearestDHS"
closest-features --delim '\t' --dist --closest $TMPOUT $TMPDIR/DHS_MCV.bed | awk '{print $(NF-1)}'| paste $TMPOUT - > $TMPOUT.new
mv $TMPOUT.new $TMPOUT


bedtools closest -d -k 20 -a $TMPOUT -b $TMPDIR/DHS_MCV.bed| awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $(NF-1), $NF}' > $OUTDIR/${SampleName}_k562_mcv_k20.bed

echo 'K562/Heme specific DHS'
unstarch /home/mauram01/scratch/hybridmice/dnase/lineagemcv/all.lineage.dhs.unmerged.starch |grep -v 'CD3\|CD34\|hESDCT\|hTH\|iPS' >$TMPDIR/DHS_nonHemeESC.bed
unstarch /home/mauram01/scratch/hybridmice/dnase/lineagemcv/all.lineage.dhs.unmerged.starch |grep  'CD3\|CD34\|hESDC\|hTH\|iPS' >$TMPDIR/DHS_HemeESC.bed

#bedops -m /vol/isg/encode/dnase/mapped/K562-DS9764/hotspots/K562-DS9764.hg38_noalt-final/K562-DS9764.hg38_noalt.fdr0.01.pks.starch $TMPDIR/DHS_HemeESC.bed |awk 'BEGIN {OFS="\t"} {print $1, $2, $3, "K562_Heme"}' >$TMPDIR/DHS_K562_heme.bed

bedmap --delim '\t' --echo --echo-map --count /vol/isg/encode/dnase/mapped/K562-DS9764/hotspots/K562-DS9764.hg38_noalt-final/K562-DS9764.hg38_noalt.fdr0.01.pks.starch  $TMPDIR/DHS_nonHemeESC.bed|awk '$NF==0' > $TMPDIR/K562HemeSpecificDHS.bed

bedmap --delim '\t' --echo --echo-map --count /vol/isg/encode/dnase/mapped/K562-DS9764/hotspots/K562-DS9764.hg38_noalt-final/K562-DS9764.hg38_noalt.fdr0.01.pks.starch $TMPDIR/DHS_nonHemeESC.bed|awk '$NF>=1' > $TMPDIR/K562NonSpecificDHS.bed
echo 'DistToNearest K562 Specific DHS'
header="$header\tDistToNearestK562SpecDHS"
closest-features --delim '\t' --dist --closest $TMPOUT  $TMPDIR/K562HemeSpecificDHS.bed| awk '{print $NF}'| paste $TMPOUT - > $TMPOUT.new
mv $TMPOUT.new $TMPOUT

echo 'n100kb K562 Specific DHS'
header="$header\tnK562SpecDHS"
bedops -u --range 100000 $TMPOUT| sort-bed -|bedmap --delim '\t' --echo --echo-map --count - $TMPDIR/K562HemeSpecificDHS.bed| awk '{print $NF}'| paste $TMPOUT - > $TMPOUT.new
mv $TMPOUT.new $TMPOUT

echo 'DistToNearest K562 Nonspecific DHS'
header="$header\tDistToNearestK562NonSpecDHS"
closest-features --delim '\t' --dist --closest $TMPOUT  $TMPDIR/K562NonSpecificDHS.bed| awk '{print $NF}'| paste $TMPOUT - > $TMPOUT.new
mv $TMPOUT.new $TMPOUT

echo 'n100kb K562 Nonspecific DHS'
header="$header\tnK562NonSpecDHS"
bedops -u --range 100000 $TMPOUT| sort-bed -|bedmap --delim '\t' --echo --echo-map --count - $TMPDIR/K562NonSpecificDHS.bed| awk '{print $NF}'| paste $TMPOUT - > $TMPOUT.new
mv $TMPOUT.new $TMPOUT


echo 'MCV of K562 Nonspecific DHS'
header="$header\tMCVK562NonSpecDHS"
bedmap  --delim '\t' --echo --fraction-ref 1.0 --count  /vol/isg/encode/dnase/mapped/K562-DS9764/hotspots/K562-DS9764.hg38_noalt-final/K562-DS9764.hg38_noalt.fdr0.01.pks.starch $TMPDIR/DHS_nonHemeESC.bed| awk '{print $1, $2, $3, $NF}' > $TMPDIR/K562NonSpecificDHS_MCV.bed
closest-features --delim '\t' --dist --closest $TMPOUT  $TMPDIR/K562NonSpecificDHS_MCV.bed| awk '{print $(NF-1)}'| paste $TMPOUT - > $TMPOUT.new
mv $TMPOUT.new $TMPOUT

echo 'doing distance to nearest hotspot2'
header="$header\tDistToNearestDHShot2"
sort-bed $TMPOUT |closest-features  --closest --dist --no-ref - /vol/isg/encode/dnase/mapped/K562-DS9764/hotspot2/K562-DS9764.hg38.hotspots.fdr0.05.starch |
awk -F'|' 'BEGIN {OFS="\t"} {print $2}' |paste $TMPOUT - > $TMPOUT.new
mv $TMPOUT.new $TMPOUT

#Distance to HeLa specific DNase
echo 'doing distance to nearest hela DHS '
header="$header\tHeLaDistToNearestDHS"
bedmap --delim '\t' --echo --echo-map --count /vol/isg/encode/dnase/mapped/HeLa_S3-DS10011/hotspots/HeLa_S3-DS10011.hg38-final/HeLa_S3-DS10011.hg38.fdr0.01.pks.bed \
/vol/isg/encode/dnase/mapped/K562-DS9767/hotspots/K562-DS9767.hg38-final/K562-DS9767.hg38.fdr0.01.pks.bed |awk 'BEGIN {OFS="\t"} $NF==0 {print $1, $2, $3}' > $TMPDIR/HeLa_DNase.bed

sort-bed $TMPOUT |closest-features  --closest --dist --no-ref - $TMPDIR/HeLa_DNase.bed |
awk -F'|' 'BEGIN {OFS="\t"} {print $2}' |paste $TMPOUT - > $TMPOUT.new
mv $TMPOUT.new $TMPOUT


#Nearest DNase peak density
echo 'doing neareast DNase peak density'
header="$header\tDNasePksDens"
paste /vol/isg/encode/dnase/mapped/K562-DS9764/hotspots/K562-DS9764.hg38_noalt-final/K562-DS9764.hg38_noalt.fdr0.01.pks.starch /vol/isg/encode/dnase/mapped/K562-DS9764/hotspots/K562-DS9764.hg38-final/K562-DS9764.hg38.fdr0.01.pks.dens.txt > $TMPDIR/${SampleName}.denspeaks.bed
sort-bed $TMPOUT |closest-features  --closest --dist --no-ref - $TMPDIR/${SampleName}.denspeaks.bed |
awk -F'|' 'BEGIN {OFS="\t"} {print $1}' |awk -F "\t" 'BEGIN {OFS="\t"} {print $NF}'|paste $TMPOUT - > $TMPOUT.new
mv $TMPOUT.new $TMPOUT


#DHS type
echo 'doing neareast DHS type'
header="$header\tNearestDHStype"
closest-features --delim '\t'  --closest --dist --no-ref $TMPOUT /home/maagj01/scratch/transposon/Analysis/K562_DNaseMasterList/DHStype_masterlist.starch |awk -F "\t" 'BEGIN {OFS="\t"} {print $(NF-3)}' | paste $TMPOUT - > $TMPOUT.new
mv $TMPOUT.new $TMPOUT


#DHS type
echo 'doing neareast DHS type strand'
header="$header\tNearestDHSstrand"
closest-features --delim '\t'  --closest --dist --no-ref $TMPOUT /home/maagj01/scratch/transposon/Analysis/K562_DNaseMasterList/DHStype_masterlist.starch |awk -F "\t" 'BEGIN {OFS="\t"} {print $(NF-1)}' | paste $TMPOUT - > $TMPOUT.new
mv $TMPOUT.new $TMPOUT


#Distance to DHS peak upstream an type
echo 'doing distance to nearest DHS - upstream '
header="$header\tDistToNearestDHS-up\tDHStype-up"
closest-features --delim '\t'  --dist  $TMPOUT /home/maagj01/scratch/transposon/Analysis/K562_DNaseMasterList/DHStype_masterlist.starch| awk -F "\t" 'BEGIN {OFS="\t"} {if ($6=="-") print $NF, $(NF-3); else if  ($6=="+") print $(NF-7), $(NF-10)}'|paste $TMPOUT - > $TMPOUT.new
mv $TMPOUT.new $TMPOUT


echo 'doing distance to nearest DHS - upstream Strand'
header="$header\tDHStype-upDirection"
closest-features --delim '\t'  --dist  $TMPOUT /home/maagj01/scratch/transposon/Analysis/K562_DNaseMasterList/DHStype_masterlist.starch| awk 'BEGIN {OFS="\t"} {if ($6=="-") print $(NF-1); else if  ($6=="+") print $(NF-8)}'|paste $TMPOUT - > $TMPOUT.new
mv $TMPOUT.new $TMPOUT


echo 'doing distance to nearest DHS - downstream '
header="$header\tDistToNearestDHS-dn\tDHStype-dn"
closest-features --delim '\t'  --dist  $TMPOUT /home/maagj01/scratch/transposon/Analysis/K562_DNaseMasterList/DHStype_masterlist.starch| awk -F "\t" 'BEGIN {OFS="\t"} {if ($6=="-") print $(NF-7), $(NF-10); else if  ($6=="+") print $NF, $(NF-3)}'|paste $TMPOUT - > $TMPOUT.new
mv $TMPOUT.new $TMPOUT


echo 'doing distance to nearest DHS - downstream Strand'
header="$header\tDHStype-dnDirection"
closest-features --delim '\t'  --dist  $TMPOUT /home/maagj01/scratch/transposon/Analysis/K562_DNaseMasterList/DHStype_masterlist.starch| awk 'BEGIN {OFS="\t"} {if ($6=="-") print $(NF-8); else if  ($6=="+") print $(NF-1)}'|paste $TMPOUT - > $TMPOUT.new
mv $TMPOUT.new $TMPOUT


#Distance to CTCF
echo 'doing distance to nearest CTCF peaks'
header="$header\tDistToNearestCTCF"
sort-bed $TMPOUT |closest-features  --closest --dist --no-ref - /home/maagj01/scratch/transposon/Analysis/K562_DNaseMasterList/Mauram2016cellReport.hg38.bed |
awk -F'|' 'BEGIN {OFS="\t"} {print $2}' |paste $TMPOUT - > $TMPOUT.new
mv $TMPOUT.new $TMPOUT




#Nearest CTCF peak density
echo 'doing neareast CTCF peak density'
header="$header\tCTCFPksDens"
gcat /vol/isg/encode/chipseq/mapped/K562-CTCF-DS11247/hotspots/K562-CTCF-DS11247.hg38_noalt-final/K562-CTCF-DS11247.hg38_noalt.fdr0.01.pks.starch | paste - /vol/isg/encode/chipseq/mapped/K562-CTCF-DS11247/hotspots/K562-CTCF-DS11247.hg38_noalt-final/K562-CTCF-DS11247.hg38_noalt.fdr0.01.pks.dens.txt  > $TMPDIR/${SampleName}.CTCF.denspeaks.bed
sort-bed $TMPOUT |closest-features  --closest --dist --no-ref - $TMPDIR/${SampleName}.CTCF.denspeaks.bed |
awk -F'|' 'BEGIN {OFS="\t"} {print $1}' |awk -F "\t" 'BEGIN {OFS="\t"} {print $NF}'|paste $TMPOUT - > $TMPOUT.new
mv $TMPOUT.new $TMPOUT

#echo 'doing neareast CTCF peak MCV score'
#bedmap  --delim '\t' --skip-unmapped --echo --fraction-ref 1.0 --count  /home/maagj01/scratch/transposon/Analysis/K562_DNaseMasterList/Mauram2016cellReport.hg38.bed /home/maagj01/scratch/transposon/Analysis/CTCF_MCV/Tissue_CTCF_unmerged.bed|
#awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $NF}' > $TMPDIR/CTCF_MCV_tissues.bed
#header="$header\tCTCFPksMCV"
#sort-bed $TMPOUT |closest-features --delim '\t'  --closest --dist --no-ref - $TMPDIR/CTCF_MCV_tissues.bed| awk '{print $(NF-1)}' |paste $TMPOUT - > $TMPOUT.new
#mv $TMPOUT.new $TMPOUT
#


#Nearest CTCF peak direction
echo 'doing neareast CTCF peak direction'
header="$header\tCTCFdirection"
sort-bed $TMPOUT |closest-features  --delim '\t' --closest --dist --no-ref - /home/maagj01/scratch/transposon/Analysis/K562_DNaseMasterList/CTCF_strand.bed |awk -F "\t" -v OFS="\t" '{print $(NF-1)}' |paste $TMPOUT - > $TMPOUT.new
mv $TMPOUT.new $TMPOUT


echo 'doing distance to nearest CTCF - upstream '
header="$header\tDistToNearestCTCF-up\tCTCFdirection-up"
closest-features --delim '\t'  --dist  $TMPOUT /home/maagj01/scratch/transposon/Analysis/K562_DNaseMasterList/CTCF_strand.bed | awk -F "\t" 'BEGIN {OFS="\t"} {if ($6=="-") print $NF, $(NF-1); else if  ($6=="+") print $(NF-7), $(NF-8)}'|paste $TMPOUT - > $TMPOUT.new
mv $TMPOUT.new $TMPOUT


echo 'doing distance to nearest CTCF - downstream '
header="$header\tDistToNearestCTCF-dn\tCTCFdirection-dn"
closest-features --delim '\t'  --dist   $TMPOUT /home/maagj01/scratch/transposon/Analysis/K562_DNaseMasterList/CTCF_strand.bed | awk -F "\t" 'BEGIN {OFS="\t"} {if ($6=="-") print $(NF-7), $(NF-8); else if  ($6=="+") print $NF, $(NF-1)}'|paste $TMPOUT - > $TMPOUT.new
mv $TMPOUT.new $TMPOUT

echo 'doing distance to loop anchor (Cohesin + CTCF)'
bedmap --delim '\t' --multidelim '\t'  --skip-unmapped --echo --echo-map /vol/isg/encode/chipseq/mapped/K562-RAD21-ENCLB559JFW/hotspots/K562-RAD21-ENCLB559JFW.hg38_noalt-final/K562-RAD21-ENCLB559JFW.hg38_noalt.fdr0.05.hot.starch /vol/isg/encode/chipseq/mapped/K562-SMC3-ENCLB209AOH/hotspots/K562-SMC3-ENCLB209AOH.hg38_noalt-final/K562-SMC3-ENCLB209AOH.hg38_noalt.fdr0.05.hot.starch |
bedmap --delim '\t' --multidelim '\t'  --skip-unmapped --echo --echo-map - /home/maagj01/scratch/transposon/Analysis/K562_DNaseMasterList/Mauram2016cellReport.hg38.bed |awk 'BEGIN {OFS="\t"} {print $1, $2, $3, "CohesinCTCF_"NR}' > $TMPDIR/Cohesin_CTCF.bed


bedtools closest -d -k 5 -a $TMPDIR/Cohesin_CTCF.bed -b $TMPDIR/Cohesin_CTCF.bed| sed 's/CohesinCTCF_//g'| awk '$4-$8==1 || $4-$8==-1'|
awk '$NF!=0' |awk 'BEGIN {OFS="\t"} {if($3<$(NF-3)) print $1, $2, $(NF-3), $4, $(NF-1), $NF; else if ($3>$(NF-3)) print $1, $(NF-3), $2, $(NF-1), $4, $NF }'| sort-bed -|
uniq |awk 'BEGIN {OFS="\t"} {print $1, $2, $3, "CohesinLoop_"NR, $NF}' > $TMPDIR/Cohesin_Loops.bed

header="$header\tcohesinLoop"
bedmap --delim '|'  --echo --echo-map $TMPOUT $TMPDIR/Cohesin_Loops.bed| awk -F '|' '{print $2}'| awk '{if ($NF=="") print "NA"; else print $(NF-1)}'|  paste $TMPOUT - > $TMPOUT.new
mv $TMPOUT.new $TMPOUT



#Distance to Dpn
echo 'doing distance to nearest Dpn sites'
header="$header\tDistToDpn"
sort-bed $TMPOUT |
closest-features --delim '\t'  --dist   - /home/maagj01/scratch/transposon/Analysis/Dpn_REsites/DpnPlusMinus.Sorted.bed |
awk -F "\t" 'BEGIN {OFS="\t"} function abs(value) {return (value<0?-value:value);} {if ($6=="-") print abs($(NF-7)); else if  ($6=="+") print  $NF}' |paste $TMPOUT -> $TMPOUT.new
mv $TMPOUT.new $TMPOUT


#Type of DNase peak
#bedops -u --range -2500:2500 /vol/isg/annotation/bed/hg38/refseq_gene/refGene.CombinedTxStarts.bed >$TMPDIR/${SampleName}.promoter.bed
#bedmap --echo --count  /vol/isg/encode/dnase/mapped/K562-DS9764/hotspots/K562-DS9764.hg38_noalt-final/K562-DS9764.hg38_noalt.fdr0.01.pks.starch /vol/isg/encode/chipseq/mapped/K562-CTCF-DS11508hotspots/K562-CTCF-DS11508.hg38_noalt-final/K562-CTCF-DS11508.hg38_noalt.fdr0.01.pks.starch |\
#awk -F'|' 'BEGIN {OFS="\t"} {print $1, $2}' |
#bedmap --echo --count - $TMPDIR/${SampleName}.promoter.bed|awk -F'|' 'BEGIN {OFS="\t"} {print $1, $2}'| 
#awk -v cutoff=1 '{if($4>=cutoff) {print $1, $2, $3, "CTCF", $5} else {print $1, $2, $3, "Distal", $5}}'|
#awk -v cutoff=1 '{if($5>=cutoff) {print $1, $2, $3, $4, "Promoter"} \
#else {print $1, $2, $3, $4, "Distal"}}'|
#awk '{if ($4=="Distal" && $5=="Distal") {print $1, $2, $3, "Distal"} \
#else {if ($4=="CTCF" && $5=="Distal") {print $1, $2, $3, "CTCF"} else {if ($5=="Promoter" && $4=="Distal") {print $1, $2, $3, "Promoter"} \
#else {if ($4=="CTCF" && $5=="Promoter") {print $1, $2, $3, "Promoter_CTCF"}}}}}'|awk -F' ' 'BEGIN {OFS="\t"} {print $1, $2, $3, $4}' > $TMPDIR/${SampleName}.DHStype.bed
#
#
#echo 'doing neareast DHS type'
#header="$header\tDHStype"
#closest-features  --closest --dist --no-ref $TMPOUT $TMPDIR/${SampleName}.DHStype.bed|
#awk -F'|' 'BEGIN {OFS="\t"} {print $1}' |awk -F "\t" 'BEGIN {OFS="\t"} {print $NF}'|paste $TMPOUT - > $TMPOUT.new
#mv $TMPOUT.new $TMPOUT

zcat /vol/isg/annotation/bed/hg38/gencodev24/src/gencode.v24.primary_assembly.annotation.gtf.gz|awk 'BEGIN {OFS="\t"} $3=="gene" && ($18=="1;"|| $18=="2;") {print $1, $4, $5, $6, $7, $10, $12}'|
sed 's/[;"]//g'| awk 'BEGIN {OFS="\t"} {split($6, id, ".")} {print $1, $2, $3, id[1], $4, $5, $7}'| 
bedtools flank -s -l 1 -r 0 -i - -g /vol/isg/annotation/fasta/hg38/hg38.chrom.sizes > $TMPDIR/GencodeV24TSS_all.bed

join -eNA -1 4 <(sort -k4,4 /tmp/GencodeV24TSS.bed) -2 4 <(sort -k4,4 /home/maagj01/scratch/transposon/Analysis/K562_GeneQuant/K562_GeneQuant_bedmapFormat.bed)| 
awk 'BEGIN {OFS="\t"} {print $2, $3, $4, $1, $5, $6, $7, $(NF-1)}'|sort-bed - > $TMPDIR/GencodeV24TSS_allTPM.bed



zcat /vol/isg/annotation/bed/hg38/gencodev24/src/gencode.v24.primary_assembly.annotation.gtf.gz|awk 'BEGIN {OFS="\t"} $3=="gene" && ($18=="1;"|| $18=="2;") {print $1, $4, $5, $6, $7, $10, $12}'|
sed 's/[;"]//g'| awk 'BEGIN {OFS="\t"} {split($6, id, ".")} {print $1, $2, $3, id[1], $4, $5, $7}'| 
bedtools flank -s -l 1 -r 0 -i - -g /vol/isg/annotation/fasta/hg38/hg38.chrom.sizes |grep 'protein_coding\|antisense\|lincRNA'> $TMPDIR/GencodeV24TSS_protLinc.bed

join -eNA -1 4 <(sort -k4,4 $TMPDIR/GencodeV24TSS_protLinc.bed) -2 4 <(sort -k4,4 /home/maagj01/scratch/transposon/Analysis/K562_GeneQuant/K562_GeneQuant_bedmapFormat.bed)| 
awk 'BEGIN {OFS="\t"} {print $2, $3, $4, $1, $5, $6, $7, $(NF-1)}'|sort-bed - > $TMPDIR/GencodeV24TSS_protLincTPM.bed
#Nearest Gencode TSS 
echo 'doing nearest Gencode.v24 gene'
header="$header\tDistToGencodeGene"
closest-features --closest --dist --no-ref $TMPOUT <( bedtools flank -s -l 1 -r 0 -i /vol/isg/annotation/bed/hg38/gencodev24/Gencodev24.gene.bed -g /vol/isg/annotation/fasta/hg38/hg38.chrom.sizes ) | awk -F'|' 'BEGIN {OFS="\t"} {print $2}' |paste $TMPOUT - > $TMPOUT.new 
mv $TMPOUT.new $TMPOUT


echo 'doing nearest Gencode.v24 TSS'
header="$header\tDistToGencodeTSS"
closest-features --closest --dist --no-ref $TMPOUT $TMPDIR/GencodeV24TSS_allTPM.bed | awk -F'|' 'BEGIN {OFS="\t"} {print $2}' |paste $TMPOUT - > $TMPOUT.new 
mv $TMPOUT.new $TMPOUT

echo 'doing nearest Gencode.v24 TSS ID'
header="$header\tToGencodeTSSID"
closest-features --closest --dist --no-ref $TMPOUT  $TMPDIR/GencodeV24TSS_allTPM.bed | awk  '{print $4}' |paste $TMPOUT - > $TMPOUT.new 
mv $TMPOUT.new $TMPOUT

echo 'doing nearest Gencode.v24 TSS TPM'
header="$header\tToGencodeTSSTPM"
closest-features --delim '\t' --closest --dist --no-ref $TMPOUT  $TMPDIR/GencodeV24TSS_allTPM.bed | awk  '{print $(NF-1)}'|paste $TMPOUT - > $TMPOUT.new 
mv $TMPOUT.new $TMPOUT


echo 'doing nearest Gencode.v24 TSS Protein and LncRNA'
header="$header\tDistToGencodeTSS_protLinc"
closest-features --closest --dist --no-ref $TMPOUT $TMPDIR/GencodeV24TSS_protLincTPM.bed | awk -F'|' 'BEGIN {OFS="\t"} {print $2}' |paste $TMPOUT - > $TMPOUT.new 
mv $TMPOUT.new $TMPOUT


echo 'doing nearest Gencode.v24 TSS ID Protein and LncRNA '
header="$header\tToGencodeTSSID_protLinc"
closest-features --closest --dist --no-ref $TMPOUT $TMPDIR/GencodeV24TSS_protLincTPM.bed| awk  '{print $4}' |paste $TMPOUT - > $TMPOUT.new 
mv $TMPOUT.new $TMPOUT

echo 'doing nearest Gencode.v24 TSS TPM expressed'
header="$header\tToGencodeTSSTPM_protLinc"
closest-features --delim '\t' --closest --dist --no-ref $TMPOUT $TMPDIR/GencodeV24TSS_protLincTPM.bed | awk  '{print $(NF-1)}'|paste $TMPOUT - > $TMPOUT.new 
mv $TMPOUT.new $TMPOUT
#
#echo 'doing nearest Gencode.v24 gene Name'
#header="$header\tdistGencodeGeneName"
#closest-features --closest --dist --no-ref $TMPOUT /home/maagj01/scratch/transposon/Analysis/K562_GeneQuant/GencodeV24.gene.bed | awk -F "\t" 'BEGIN {OFS="\t"} {print $4}' |paste $TMPOUT - > $TMPOUT.new 
#mv $TMPOUT.new $TMPOUT


#K562 expressed genes https://www.encodeproject.org/files/ENCFF104VTJ/@@download/ENCFF104VTJ.tsv,  Thomas Gingeras, CSHL, paired-ended, stranded, total RNA, biological replicate 1
#awk -F "\t" -v OFS='\t''{print $1, $6}' ENCFF104VTJ.tsv |grep ENSG >K562_GeneCount_TMP.tsv
#zcat /home/maagj01/scratch/transposon/Analysis/K562_GeneQuant/gencode.v24.primary_assembly.annotation.gtf.gz| 
#awk 'BEGIN{OFS="\t"}{if($3 == "gene"){split($0, a, " "); i=0; while(a[i] != "gene_id") i++; gsub(/[";]/, "", a[i+1]); print $1, $4-1, $5, a[i+1], ".", $7}}' - | 
#awk -F "\t" 'BEGIN {OFS="\t"} {split($4, name, "."); print $1, $2, $3, name[1], $5, $6}'|  
#sort -k4| 
#join -eNA -1 4 -2 1 -  -o 1.1 1.2 1.3 0 1.5 1.6 2.2 <(sort -k1 /home/maagj01/scratch/transposon/Analysis/K562_GeneQuant/K562_GeneCount_TMP.tsv| awk -F "\t" 'BEGIN {OFS="\t"} {split($1, name, "."); print name[1], $2}')| sort-bed - >  /home/maagj01/scratch/transposon/Analysis/K562_GeneQuant/K562_GeneQuant.bed

echo 'doing nearest Gencode.v24 gene expression'
header="$header\tDistToGencodeGeneTPM"
closest-features --closest --dist --no-ref $TMPOUT <(awk '$NF>=1' /home/maagj01/scratch/transposon/Analysis/K562_GeneQuant/K562_GeneQuant.bed| awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5, $6, $7}'|bedtools flank -s -l 10000 -r 0 -i - -g /vol/isg/annotation/fasta/hg38/hg38.chrom.sizes )|
awk -F'|' 'BEGIN {OFS="\t"} {print $2}' |awk  '{print $NF}' |paste $TMPOUT - > $TMPOUT.new 
mv $TMPOUT.new $TMPOUT



#Overlap with CGI
echo "doing Overlap with CGI"
header="$header\tCGIovr"
awk -F "\t" 'BEGIN {OFS="\t"} {$4="TRUE"; print}' /vol/isg/annotation/bed/hg38/cpg_islands/cpgIslands.bed |
bedmap --faster --delim '\t' --echo --echo-map-id $TMPOUT - | awk -F "\t" 'BEGIN {OFS="\t"} $NF~/TRUE/ {$NF="TRUE"} $NF=="" {$NF="FALSE"} {print}' > $TMPOUT.new
mv $TMPOUT.new $TMPOUT

#Nearest CGI
echo "doing distance nearest CGI"
header="$header\tDistToCGI"
sort-bed /vol/isg/annotation/bed/hg38/cpg_islands/cpgIslands.bed | closest-features --closest --dist --no-ref  $TMPOUT -| awk -F'|' 'BEGIN {OFS="\t"} {print $2}' |paste $TMPOUT  - > $TMPOUT.new
mv $TMPOUT.new $TMPOUT

bedmap --delim '\t'  --echo --echo-map  --mean /vol/isg/annotation/bed/hg38/cpg_islands/cpgIslands.bed /home/maagj01/scratch/transposon/Analysis/K562_reference_epigenome/WGBS/K562-CpGMethylation_Rep1_ENCFF827YGC.bed |
awk -F "\t" 'BEGIN {OFS="\t"} {if ($NF=="NAN")  print $1, $2, $3, $4, 0; else if  ($6!="NAN") print $1, $2, $3, $4, $NF}' > $TMPDIR/K562_CGI_methylation.bed
echo "Methylation of nearest CGI"
header="$header\tCGImethylation"
closest-features --delim '\t' --closest --dist --no-ref  $TMPOUT $TMPDIR/K562_CGI_methylation.bed| awk '{print $(NF-1)}' | paste $TMPOUT  - > $TMPOUT.new
mv $TMPOUT.new $TMPOUT


echo "Mean methylation +-2kb of insertion"
header="$header\tInsMethylation2kb"
bedops -u --range  2000 $TMPOUT|sort-bed -| bedmap --delim '\t' --echo --echo-map --mean  - /home/maagj01/scratch/transposon/Analysis/K562_reference_epigenome/WGBS/K562-CpGMethylation_Rep1_ENCFF827YGC.bed|awk '{print $NF}' | paste $TMPOUT  - > $TMPOUT.new
mv $TMPOUT.new $TMPOUT


echo "Number of CpG +-2kb of insertion"
header="$header\tInsCpG2kb"
bedops -u --range  2000 $TMPOUT|sort-bed -| bedmap --delim '\t' --echo --echo-map --count  - /home/maagj01/scratch/transposon/Analysis/K562_reference_epigenome/WGBS/K562-CpGMethylation_Rep1_ENCFF827YGC.bed|awk '{print $NF}' | paste $TMPOUT  - > $TMPOUT.new
mv $TMPOUT.new $TMPOUT


#Percentage G+C content 75bp around insertion
echo "doing G+C content +/-75bp"
header="$header\tPercentGC"
awk -F "\t" -v widen=75 'BEGIN {OFS="\t"} {$2=$2 > widen ? $2-widen : 0; $3+=widen; print;}' $TMPOUT | ~/bin/bed2fasta.pl - /vol/isg/annotation/fasta/hg38 2>/dev/null |
grep -v -e "^>" | tr '[a-z]' '[A-Z]' | perl -ne 'chomp; print length($_) != 0 ? ($_ =~ tr/[gcGC]//) / length($_) : "NA"; print "\n"' | paste $TMPOUT - > $TMPOUT.new
mv $TMPOUT.new $TMPOUT


echo "doing chromosome band"
header="$header\tchrArm"
gcat /vol/isg/annotation/bed/hg38/ucsc.chromosome.bands/cytoBand.txt.gz | awk -F "\t" 'BEGIN {OFS="\t"} {print $1, $2, $3, substr($1,4) $4}' | bedmap --faster --echo --echo-map-id $TMPOUT - |
awk -F'|' 'BEGIN {OFS="\t"} {print $2}'| paste $TMPOUT - > $TMPOUT.new 
mv $TMPOUT.new $TMPOUT


echo "doing cyto band type"
header="$header\tcytoBand"
gcat /vol/isg/annotation/bed/hg38/ucsc.chromosome.bands/cytoBand.txt.gz | cut -f1-3,5 | bedmap --faster --delim '\t' --echo --echo-map-id-uniq $TMPOUT - > $TMPOUT.new 
mv $TMPOUT.new $TMPOUT


echo "doing GenicLocation"
header="$header\tGenicLocation"
gcat /vol/isg/annotation/bed/hg38/refseq_paint/refGene_hg38.bed6.starch | grep -v promoter > $TMPDIR/refGene_hg38.bed6
bedmap --delim '\t' --echo-map-id $TMPOUT $TMPDIR/refGene_hg38.bed6 | awk -F "\t" 'BEGIN {OFS="\t"} {col=NF} $col~/coding/ {$col="coding"} $col~/TxS-1000/ {$col="TxS-1000"} $col~/TxStart-10kb/ {$col="TxStart-10kb"} $col~/3.UTR/ {$col="3UTR"} $col~/5.UTR/ {$col="5UTR"} $col~/intron/ {$col="intron"} $col=="" || $col~/3.proximal/ {$col="intergenic"} {print;}' | paste $TMPOUT - > $TMPOUT.new
mv $TMPOUT.new $TMPOUT


echo "doing GenicLocation.simple"
header="$header\tGenicLocation.simple"
gcat /vol/isg/annotation/bed/hg38/refseq_paint/refGene_hg38.bed6.starch | grep -v TxS-1000 | grep -v TxStart-10kb > $TMPDIR/refGene_hg38.simple.bed6
bedmap --delim '\t' --echo-map-id $TMPOUT $TMPDIR/refGene_hg38.simple.bed6 | awk -F "\t" 'BEGIN {OFS="\t"} {col=NF} $col~/coding/ {$col="coding"} $col~/promoter/ {$col="promoter"} $col~/3.UTR/ {$col="3UTR"} $col~/5.UTR/ {$col="5UTR"} $col~/intron/ {$col="intron"} $col=="" || $col~/3.proximal/ {$col="intergenic"} {print;}' | paste $TMPOUT - > $TMPOUT.new
mv $TMPOUT.new $TMPOUT


echo "doing number of DHS with 10kb of insertion"
header="$header\tnDHS10kb"
bedops -u --range 5000  $TMPOUT| bedmap --delim '\t' --echo --echo-map --count - /vol/isg/encode/dnase/mapped/K562-DS9764/hotspots/K562-DS9764.hg38_noalt-final/K562-DS9764.hg38_noalt.fdr0.01.pks.starch|awk '{print $NF}' |paste $TMPOUT -  > $TMPOUT.new 
mv $TMPOUT.new $TMPOUT


echo "doing number of DHS with 100kb of insertion"
header="$header\tnDHS100kb"
bedops -u --range 100000  $TMPOUT| bedmap --delim '\t' --echo --echo-map --count - /vol/isg/encode/dnase/mapped/K562-DS9764/hotspots/K562-DS9764.hg38_noalt-final/K562-DS9764.hg38_noalt.fdr0.01.pks.starch|awk '{print $NF}' |paste $TMPOUT -  > $TMPOUT.new 
mv $TMPOUT.new $TMPOUT

echo "doing number of DHS with 300kb +- of insertion"
header="$header\tnDHS300kb"
bedops -u --range 300000  $TMPOUT| bedmap --delim '\t' --echo --echo-map --count - /vol/isg/encode/dnase/mapped/K562-DS9764/hotspots/K562-DS9764.hg38_noalt-final/K562-DS9764.hg38_noalt.fdr0.01.pks.starch|awk '{print $NF}' |paste $TMPOUT -  > $TMPOUT.new 
mv $TMPOUT.new $TMPOUT


echo "doing number of DHS with 500kb +- of insertion"
header="$header\tnDHS500kb"
bedops -u --range 500000  $TMPOUT| bedmap --delim '\t' --echo --echo-map --count - /vol/isg/encode/dnase/mapped/K562-DS9764/hotspots/K562-DS9764.hg38_noalt-final/K562-DS9764.hg38_noalt.fdr0.01.pks.starch|awk '{print $NF}' |paste $TMPOUT -  > $TMPOUT.new 
mv $TMPOUT.new $TMPOUT


echo "doing number of DHS with 1Mb +- of insertion"
header="$header\tnDHS1Mb"
bedops -u --range 1000000  $TMPOUT| bedmap --delim '\t' --echo --echo-map --count - /vol/isg/encode/dnase/mapped/K562-DS9764/hotspots/K562-DS9764.hg38_noalt-final/K562-DS9764.hg38_noalt.fdr0.01.pks.starch|awk '{print $NF}' |paste $TMPOUT -  > $TMPOUT.new 
mv $TMPOUT.new $TMPOUT


paste /vol/isg/encode/dnase/mapped/K562-DS9764/hotspots/K562-DS9764.hg38_noalt-final/K562-DS9764.hg38_noalt.fdr0.01.pks.starch \
/vol/isg/encode/dnase/mapped/K562-DS9764/hotspots/K562-DS9764.hg38-final/K562-DS9764.hg38.fdr0.01.pks.dens.txt \
/vol/isg/encode/dnase/mapped/K562-DS9764/hotspots/K562-DS9764.hg38-final/K562-DS9764.hg38.fdr0.01.pks.pval.txt |
awk 'BEGIN {OFS="\t"} {print $1, $2, $3, "DHS_"NR, $4, $5}'> $TMPDIR/K562_DHS_dens_p.bed

echo "doing density of DHS with 100kb of insertion"
header="$header\tnDHS100kbDensity"
bedops -u --range 100000  $TMPOUT| bedmap --delim '\t' --echo --echo-map --sum - <(awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5}' $TMPDIR/K562_DHS_dens_p.bed)|awk '{print $NF}' |sed 's/NAN/0/g'|paste $TMPOUT -  > $TMPOUT.new 
mv $TMPOUT.new $TMPOUT


echo "doing density of DHS with 500kb +- of insertion"
header="$header\tnDHS500kbDensity"
bedops -u --range 500000  $TMPOUT| bedmap --delim '\t' --echo --echo-map --sum - <(awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5}' $TMPDIR/K562_DHS_dens_p.bed)|awk '{print $NF}' |sed 's/NAN/0/g'|paste $TMPOUT -  > $TMPOUT.new 
mv $TMPOUT.new $TMPOUT


echo "doing density of DHS with 1Mb +- of insertion"
header="$header\tnDHS1MbDensity"
bedops -u --range 1000000  $TMPOUT| bedmap --delim '\t' --echo --echo-map --sum - <(awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5}' $TMPDIR/K562_DHS_dens_p.bed)|awk '{print $NF}' |sed 's/NAN/0/g'|paste $TMPOUT -  > $TMPOUT.new 
mv $TMPOUT.new $TMPOUT


echo "doing number of HeLa DHS with 300kb +- of insertion"
header="$header\tnHeLaDHS300kb"
bedops -u --range 300000  $TMPOUT| bedmap --delim '\t' --echo --echo-map --count - /vol/isg/encode/dnase/mapped/HeLa_S3-DS10011/hotspots/HeLa_S3-DS10011.hg38-final/HeLa_S3-DS10011.hg38.fdr0.01.pks.bed|awk '{print $NF}' |paste $TMPOUT -  > $TMPOUT.new 
mv $TMPOUT.new $TMPOUT


echo "doing number of Gencode TSS within 100kb of insertion"
header="$header\tnTSS100kb"
awk '$4 ~/\-001/' /vol/isg/annotation/bed/hg38/gencodev24/Gencodev24.TxStarts.bed |sort-bed - > $TMPDIR/Gencodev24.TxStarts-001.bed
bedops -u --range 50000 $TMPOUT | bedmap --delim '\t' --echo --echo-map --count - $TMPDIR/Gencodev24.TxStarts-001.bed|awk '{print $NF}' |paste $TMPOUT -  > $TMPOUT.new 
mv $TMPOUT.new $TMPOUT



echo 'doing number of k562 Gencode v24 expressed genes (TPM>=1)  within 100kb of insertion'
header="$header\tnExpressed100kb"
#Had to filter for same gene names 
zcat /home/maagj01/scratch/transposon/Analysis/K562_GeneQuant/gencode.v24.primary_assembly.annotation.gtf.gz|awk '$3=="gene"'|
awk '{print $10,$16}'|sed 's/[";]//g'|awk -F' ' 'BEGIN {OFS="\t"} {split($1, name, "."); print name[1], $2}' > $TMPDIR/gencode.v24.primary_assembly.GeneNames.tsv


cat /home/maagj01/scratch/transposon/Analysis/K562_GeneQuant/K562_GeneQuant.bed| awk '$NF>=1'| awk '{print $4}'|
 fgrep -f - $TMPDIR/gencode.v24.primary_assembly.GeneNames.tsv|awk '{print $2}'|fgrep -wf - $TMPDIR/Gencodev24.TxStarts-001.bed >$TMPDIR/Gencodev24.TxStarts-001_Expressed.bed

#sed 's/-001//g' $TMPDIR/Gencodev24.TxStarts-001.bed|awk '{print $4}' |grep -v 'MT-\|-AS1\|-AS2'| fgrep -wf - $TMPDIR/gencode.v24.primary_assembly.GeneNames.tsv|
#awk '{print $1}'|fgrep -f - /home/maagj01/scratch/transposon/Analysis/K562_GeneQuant/K562_GeneQuant.bed| awk '$NF>=1'  >$TMPDIR/${SampleName}.k562expressedgenes.bed

bedops -u --range 50000 $TMPOUT |bedmap --delim '\t' --echo --echo-map --count  - $TMPDIR/Gencodev24.TxStarts-001_Expressed.bed|awk  '{print $NF}' |paste $TMPOUT -  > $TMPOUT.new 
mv $TMPOUT.new $TMPOUT


echo 'doing number of CGI within 100kb of insertion'
header="$header\tnCGI100kb"
bedops -u --range 50000 $TMPOUT |bedmap --delim '\t' --echo --echo-map --count - /vol/isg/annotation/bed/hg38/cpg_islands/cpgIslands.bed  |awk -F "\t" 'BEGIN {OFS="\t"} {print $NF}'|paste $TMPOUT -  > $TMPOUT.new 
mv $TMPOUT.new $TMPOUT


echo 'doing number of CTCF peaks in region'
header="$header\tnCTCFpeaks100kb"
bedops -u --range 50000 $TMPOUT |bedmap --delim '\t' --echo --echo-map --count  - /home/maagj01/scratch/transposon/Analysis/K562_DNaseMasterList/CTCF_strand.bed |awk -F "\t" 'BEGIN {OFS="\t"} {print $NF}'|paste $TMPOUT -  > $TMPOUT.new 
mv $TMPOUT.new $TMPOUT

#echo "doing ContactMap analysis"
#header="$header\tContactMap"
#bedmap --delim '\t' --echo --echo-map $TMPOUT /home/maagj01/scratch/transposon/Analysis/K562_reference_epigenome/Hi-C/DHSoverlappingOther_connections.starch|awk -F"\t" 'BEGIN {OFS="\t"}  $NF=="" {$NF="NA"} {print ;}' > $TMPOUT.new 
#mv $TMPOUT.new $TMPOUT
#


echo 'doing distance to nearest active enhancer -- defined by overlapped H3K4me1 and H3K27ac, surrounding a DHS within +- 10kb' 
header="$header\tDistToActiveEnhancer"
bedmap --delim '\t' --skip-unmapped --echo  /vol/isg/encode/chipseq/mapped/K562-H3K27ac-ENCLB695AJD/hotspots/K562-H3K27ac-ENCLB695AJD.hg38_noalt-final/K562-H3K27ac-ENCLB695AJD.hg38_noalt.fdr0.01.pks.starch /vol/isg/encode/chipseq/mapped/K562-H3K4me1-ENCLB695AOH/hotspots/K562-H3K4me1-ENCLB695AOH.hg38_noalt-final/K562-H3K4me1-ENCLB695AOH.hg38_noalt.fdr0.01.pks.starch |
closest-features --delim '\t' --dist /vol/isg/encode/dnase/mapped/K562-DS9764/hotspots/K562-DS9764.hg38_noalt-final/K562-DS9764.hg38_noalt.fdr0.01.pks.starch - | 
awk '$7>-10000 && $NF<10000 {print $1, $2, $3}' |
closest-features  --delim '\t' --closest --dist --no-ref $TMPOUT -|
awk '{print $NF}' |paste $TMPOUT - > $TMPOUT.new 

mv $TMPOUT.new $TMPOUT


echo 'doing distance to active promoter --- defined as overlapping H3K4me3 and H3K27ac, surrounding a DHS withing +-10kb (irrespective of gene annotation)'
header="$header\tDistToActivePromoter"
bedmap --delim '\t' --skip-unmapped --echo  /vol/isg/encode/chipseq/mapped/K562-H3K27ac-ENCLB695AJD/hotspots/K562-H3K27ac-ENCLB695AJD.hg38_noalt-final/K562-H3K27ac-ENCLB695AJD.hg38_noalt.fdr0.01.pks.starch /vol/isg/encode/chipseq/mapped/K562-H3K4me3-DS11507/hotspots/K562-H3K4me3-DS11507.hg38_noalt-final/K562-H3K4me3-DS11507.hg38_noalt.fdr0.01.pks.starch | 
closest-features --delim '\t' --dist /vol/isg/encode/dnase/mapped/K562-DS9764/hotspots/K562-DS9764.hg38_noalt-final/K562-DS9764.hg38_noalt.fdr0.01.pks.starch - |
awk '$7>-10000 && $NF<10000 {print $1, $2, $3}' | closest-features  --delim '\t' --closest --dist --no-ref $TMPOUT -|
awk '{print $NF}' |paste $TMPOUT - > $TMPOUT.new 

mv $TMPOUT.new $TMPOUT
####
#TADs
####
echo "doing Arrowhead analysis"
header="$header\tArrowhead"
bedmap --delim '\t' --echo --echo-map-id $TMPOUT /home/maagj01/scratch/transposon/Analysis/K562_reference_epigenome/Hi-C/Arrhead_K562_ID_collapsed_hg38.bed|awk -F"\t" 'BEGIN {OFS="\t"}  $NF=="" {$NF="NA"} {print ;}' > $TMPOUT.new 
mv $TMPOUT.new $TMPOUT


echo "doing TADs Armatus gamma 0.2 analysis"
header="$header\tTAD_armatus0.2"
bedmap --delim '\t' --echo --echo-map-id $TMPOUT /home/maagj01/scratch/transposon/Analysis/K562_reference_epigenome/Hi-C/K562/raodomains_5kb_KR/hg38/TADs/Id/Armatus_5kb_gamma.0.2.0.bed|awk -F"\t" 'BEGIN {OFS="\t"}  $NF=="" {$NF="NA"} {print ;}' > $TMPOUT.new 
mv $TMPOUT.new $TMPOUT

echo "doing TADs Armatus gamma 0.2 analysis --- DHS within TAD"
header="$header\tTAD_armatus0.2_DHS\tTAD_armatus0.2_length"
bedmap --delim '\t' --echo --echo-map --count /home/maagj01/scratch/transposon/Analysis/K562_reference_epigenome/Hi-C/K562/raodomains_5kb_KR/hg38/TADs/Id/Armatus_5kb_gamma.0.2.0.bed \
/vol/isg/encode/dnase/mapped/K562-DS9764/hotspots/K562-DS9764.hg38_noalt-final/K562-DS9764.hg38_noalt.fdr0.01.pks.starch|
awk '{print $NF}'| paste /home/maagj01/scratch/transposon/Analysis/K562_reference_epigenome/Hi-C/K562/raodomains_5kb_KR/hg38/TADs/Id/Armatus_5kb_gamma.0.2.0.bed - > $TMPDIR/Armatus_0.2.0_DHS.bed 
bedmap --delim '\t' --echo --echo-map --count $TMPOUT $TMPDIR/Armatus_0.2.0_DHS.bed  |awk -F"\t" 'BEGIN {OFS="\t"}  $NF==0 {$(NF-1)="NA"} {print $(NF-1), ($(NF-4)-$(NF-5))}' | paste $TMPOUT - > $TMPOUT.new 

mv $TMPOUT.new $TMPOUT

echo "doing TADs Armatus gamma 0.5 analysis"
header="$header\tTAD_armatus0.5"
bedmap --delim '\t' --echo --echo-map-id $TMPOUT /home/maagj01/scratch/transposon/Analysis/K562_reference_epigenome/Hi-C/K562/raodomains_5kb_KR/hg38/TADs/Id/Armatus_5kb_gamma.0.5.0.bed|awk -F"\t" 'BEGIN {OFS="\t"}  $NF=="" {$NF="NA"} {print ;}' > $TMPOUT.new 
mv $TMPOUT.new $TMPOUT

echo "doing TADs Armatus gamma 0.5 analysis --- DHS within TAD"
header="$header\tTAD_armatus0.5_DHS\tTAD_armatus0.5_length"
bedmap --delim '\t' --echo --echo-map --count /home/maagj01/scratch/transposon/Analysis/K562_reference_epigenome/Hi-C/K562/raodomains_5kb_KR/hg38/TADs/Id/Armatus_5kb_gamma.0.5.0.bed \
/vol/isg/encode/dnase/mapped/K562-DS9764/hotspots/K562-DS9764.hg38_noalt-final/K562-DS9764.hg38_noalt.fdr0.01.pks.starch|
awk '{print $NF}'| paste /home/maagj01/scratch/transposon/Analysis/K562_reference_epigenome/Hi-C/K562/raodomains_5kb_KR/hg38/TADs/Id/Armatus_5kb_gamma.0.5.0.bed - > $TMPDIR/Armatus_0.5.0_DHS.bed 
bedmap --delim '\t' --echo --echo-map --count $TMPOUT $TMPDIR/Armatus_0.5.0_DHS.bed  |awk -F"\t" 'BEGIN {OFS="\t"}  $NF==0 {$(NF-1)="NA"} {print $(NF-1), ($(NF-4)-$(NF-5))}' | paste $TMPOUT - > $TMPOUT.new 

mv $TMPOUT.new $TMPOUT


echo "doing TADs Armatus gamma 0.5 analysis --- distance to DHS in TAD"
header="$header\tTAD_nearestDHS"
bedmap --delim '\t' --echo --echo-map --skip-unmapped --count /vol/isg/encode/dnase/mapped/K562-DS9764/hotspots/K562-DS9764.hg38_noalt-final/K562-DS9764.hg38_noalt.fdr0.01.pks.starch \
/home/maagj01/scratch/transposon/Analysis/K562_reference_epigenome/Hi-C/K562/raodomains_5kb_KR/hg38/TADs/Id/Armatus_5kb_gamma.0.5.0.bed |
awk '$NF==1 {print $1, $2, $3, $7}' > $TMPDIR/DHS_in_TAD_0.05.bed
closest-features --delim '\t' --dist --closest $TMPOUT $TMPDIR/DHS_in_TAD_0.05.bed| awk '{print $(NF-1)}'| paste $TMPOUT - > $TMPOUT.new 

mv $TMPOUT.new $TMPOUT



echo "doing TADs Armatus gamma 0.8 analysis"
header="$header\tTAD_armatus0.8"
bedmap --delim '\t' --echo --echo-map-id $TMPOUT /home/maagj01/scratch/transposon/Analysis/K562_reference_epigenome/Hi-C/K562/raodomains_5kb_KR/hg38/TADs/Id/Armatus_5kb_gamma.0.8.0.bed|awk -F"\t" 'BEGIN {OFS="\t"}  $NF=="" {$NF="NA"} {print ;}' > $TMPOUT.new 
mv $TMPOUT.new $TMPOUT

echo "doing TADs Armatus gamma 0.8 analysis --- DHS within TAD"
header="$header\tTAD_armatus0.8_DHS\tTAD_armatus0.8_length"
bedmap --delim '\t' --echo --echo-map --count /home/maagj01/scratch/transposon/Analysis/K562_reference_epigenome/Hi-C/K562/raodomains_5kb_KR/hg38/TADs/Id/Armatus_5kb_gamma.0.8.0.bed \
/vol/isg/encode/dnase/mapped/K562-DS9764/hotspots/K562-DS9764.hg38_noalt-final/K562-DS9764.hg38_noalt.fdr0.01.pks.starch|
awk '{print $NF}'| paste /home/maagj01/scratch/transposon/Analysis/K562_reference_epigenome/Hi-C/K562/raodomains_5kb_KR/hg38/TADs/Id/Armatus_5kb_gamma.0.8.0.bed - > $TMPDIR/Armatus_0.8.0_DHS.bed 
bedmap --delim '\t' --echo --echo-map --count $TMPOUT $TMPDIR/Armatus_0.8.0_DHS.bed  |awk -F"\t" 'BEGIN {OFS="\t"}  $NF==0 {$(NF-1)="NA"} {print $(NF-1), ($(NF-4)-$(NF-5))}' | paste $TMPOUT - > $TMPOUT.new 

mv $TMPOUT.new $TMPOUT


echo "doing A/B compartment analysis --- value"
header="$header\tAB_compartment_value"
bedmap --echo-map $TMPOUT /home/maagj01/public_html/blog/2017Jun26/LAD/reverseValues/K562_ABcomparments_Reverse.bedGraph|awk -F"\t" 'BEGIN {OFS="\t"}  $NF=="" {$NF="NA"} {print $NF}'| paste $TMPOUT - > $TMPOUT.new
mv $TMPOUT.new $TMPOUT

echo "doing A/B compartment analysis --- Compartment"
header="$header\tABcomp"
bedmap --echo-map $TMPOUT /home/maagj01/public_html/blog/2017Jun26/LAD/reverseValues/K562_ABcomparments_Reverse.bedGraph|awk -F"\t" 'BEGIN {OFS="\t"}  $NF=="" {$NF="NA"} {print $NF}'| 
awk '{if ($1=="NA") print "NA"; else if ($1 <= 0 ) print "B"; else if ($1 > 0) print "A"; else print "NA"}'|paste $TMPOUT - > $TMPOUT.new
mv $TMPOUT.new $TMPOUT

echo "doing A/B compartment analysis --- Compartment number"
header="$header\tABcompNumber"
bedmap  --echo-map $TMPOUT /home/maagj01/public_html/blog/OverAllAnalysis/K562_geneExp_ABcomp/K562_ABcomparments_Reverse_Numbered.bedGraph|awk -F"\t" 'BEGIN {OFS="\t"}  $NF=="" {$NF="NA"} {print $NF}' |paste $TMPOUT - > $TMPOUT.new
mv $TMPOUT.new $TMPOUT



echo 'doing distance to opposite compartment'
cat /home/maagj01/public_html/blog/2017Jun26/LAD/reverseValues/K562_ABcomparments_Reverse.bedGraph| 
awk 'BEGIN {OFS="\t"} {if ($NF=="NA") print $1, $2, $3, "NA"; else if ($NF <= 0 ) print $1, $2, $3, "B"; else if ($NF > 0) print $1, $2, $3, "A"; else print $1, $2, $3, "NA"}' |awk '$NF=="A"' >$TMPDIR/Acomp.bed

cat /home/maagj01/public_html/blog/2017Jun26/LAD/reverseValues/K562_ABcomparments_Reverse.bedGraph| 
awk 'BEGIN {OFS="\t"} {if ($NF=="NA") print $1, $2, $3, "NA"; else if ($NF <= 0 ) print $1, $2, $3, "B"; else if ($NF > 0) print $1, $2, $3, "A"; else print $1, $2, $3, "NA"}' |awk '$NF=="B"' >$TMPDIR/Bcomp.bed

header="$header\tABcomp_Adist"
closest-features --ec --delim '\t' --dist --closest $TMPOUT $TMPDIR/Acomp.bed |awk '{print $NF}' |paste $TMPOUT - > $TMPOUT.new
mv $TMPOUT.new $TMPOUT

header="$header\tABcomp_Bdist"
closest-features --ec --delim '\t' --dist --closest $TMPOUT $TMPDIR/Bcomp.bed |awk '{print $NF}'| paste $TMPOUT - > $TMPOUT.new
mv $TMPOUT.new $TMPOUT

echo "doing ChromHMM"
header="$header\tchromHMM"
bedmap --echo-map $TMPOUT /home/maagj01/scratch/transposon/Analysis/K562_reference_epigenome/data/ChromHMM_18state_hg38_sorted.bed| awk -F"\t" 'BEGIN {OFS="\t"}  $4=="" {$4="NA"} {print $4}'| paste $TMPOUT - > $TMPOUT.new
mv $TMPOUT.new $TMPOUT


echo "doing genomicDivision"
header="$header\tgenomicDivision"
bedmap --echo-map $TMPOUT  /home/maagj01/public_html/blog/OverAllAnalysis/K562_geneExp_ABcomp/K562_genomicDivision.bed|
awk -F"\t" 'BEGIN {OFS="\t"}  $NF=="" {$NF="NA"} {print ;}'|awk '{print $NF}' |paste $TMPOUT - > $TMPOUT.new
mv $TMPOUT.new $TMPOUT



echo "doing genomicDivision No Expressed gene (TMP >0.1) (also removed heterochromatin)"
header="$header\tgenomicDivisionNoExp"
bedmap --echo-map $TMPOUT  /home/maagj01/public_html/blog/OverAllAnalysis/K562_geneExp_ABcomp/K562_genomicDivision_geneDesertGeneRichover0.1TPM.bed |
awk -F"\t" 'BEGIN {OFS="\t"}  $NF=="" {$NF="NA"} {print ;}'|awk '{print $NF}' |paste $TMPOUT - > $TMPOUT.new
mv $TMPOUT.new $TMPOUT



echo 'doing distance to opposite genomic regions'
cat /home/maagj01/public_html/blog/OverAllAnalysis/K562_geneExp_ABcomp/K562_genomicDivision_geneDesertGeneRichover0.1TPM.bed | 
awk '$4=="geneDesert"' > $TMPDIR/geneDesert.bed

cat /home/maagj01/public_html/blog/OverAllAnalysis/K562_geneExp_ABcomp/K562_genomicDivision_geneDesertGeneRichover0.1TPM.bed | 
awk '$4=="geneRich"' > $TMPDIR/generich.bed

header="$header\tgenomicDivision_desertDist"
closest-features --ec --delim '\t' --dist --closest $TMPOUT $TMPDIR/geneDesert.bed |awk '{print $NF}' |paste $TMPOUT - > $TMPOUT.new
mv $TMPOUT.new $TMPOUT

header="$header\tgenomicDivision_richDist"
closest-features --ec --delim '\t' --dist --closest $TMPOUT $TMPDIR/generich.bed|awk '{print $NF}'| paste $TMPOUT - > $TMPOUT.new
mv $TMPOUT.new $TMPOUT



#generating decile of random insertions
bedtools random -l 1 -n 1295 -g <( grep -v 'Un\|alt\|random' /vol/isg/annotation/fasta/hg38/hg38.chrom.sizes)| sort-bed - |
bedmap --echo-map -  /home/maagj01/public_html/blog/OverAllAnalysis/K562_geneExp_ABcomp/K562_genomicDivision.bed| awk -F"\t" 'BEGIN {OFS="\t"}  $NF=="" {$NF="NA"} {print ;}' |grep -v NA > $TMPDIR/randomGeneticDivision.bed


echo "doing distance to geneDesert"
header="$header\tDistToGeneDesert"
cat <(awk 'BEGIN {OFS="\t"} {print $1, $2, $2+1}' /home/maagj01/public_html/blog/OverAllAnalysis/K562_geneExp_ABcomp/geneDeserts.tsv) <(awk 'BEGIN {OFS="\t"} {print $1, $3, $3+1}' /home/maagj01/public_html/blog/OverAllAnalysis/K562_geneExp_ABcomp/geneDeserts.tsv)|sort-bed -|
closest-features --ec --delim '\t' --dist --closest $TMPOUT -| awk '{print $NF}' | paste $TMPOUT - > $TMPOUT.new
mv $TMPOUT.new $TMPOUT

echo "doing K562 specific gene desert"
header="$header\tK562geneDesert"
bedmap --echo-map $TMPOUT  /home/maagj01/public_html/blog/OverAllAnalysis/K562_geneExp_ABcomp/K562_noExpressed_geneDesert.tsv| 
awk -F"\t" 'BEGIN {OFS="\t"}  {if ($NF=="") print "NA"; else print "K562geneDesert"}'| paste $TMPOUT - > $TMPOUT.new
mv $TMPOUT.new $TMPOUT


echo "doing distance to K562geneDesert"
header="$header\tDistTok562GeneDesert"
cat <(awk 'BEGIN {OFS="\t"} {print $1, $2, $2+1}' /home/maagj01/public_html/blog/OverAllAnalysis/K562_geneExp_ABcomp/K562_noExpressed_geneDesert.tsv) <(awk 'BEGIN {OFS="\t"} {print $1, $3, $3+1}' /home/maagj01/public_html/blog/OverAllAnalysis/K562_geneExp_ABcomp/K562_noExpressed_geneDesert.tsv)|sort-bed -|
closest-features --ec --delim '\t' --dist --closest $TMPOUT -| awk '{print $NF}' | paste $TMPOUT - > $TMPOUT.new
mv $TMPOUT.new $TMPOUT


echo "doing Replication phase"
header="$header\trepliSeq"
bedmap --delim '\t' --echo --echo-map --count $TMPOUT /home/maagj01/public_html/blog/OverAllAnalysis/K562_geneExp_ABcomp/repliSeq/repliPhases.bed|
awk -F"\t" 'BEGIN {OFS="\t"}  {if ($NF==0 || $NF==2) print "NA"; else if  ($NF==1) print $(NF-6)}'| paste $TMPOUT - > $TMPOUT.new
mv $TMPOUT.new $TMPOUT


echo "doing distance to Ultra Conserved"
header="$header\tultraConserved"
closest-features --ec --delim '\t' --dist --closest $TMPOUT /home/maagj01/public_html/blog/OverAllAnalysis/K562_geneExp_ABcomp/UltraConservedElements.hg38.bed| awk '{print $NF}' | paste $TMPOUT - > $TMPOUT.new
mv $TMPOUT.new $TMPOUT


echo "doing distance to H3K27me3 NarroPeaks Encode"
header="$header\tH3K27me3NarrowPeaks"
closest-features --ec --delim '\t' --dist --closest $TMPOUT /home/maagj01/public_html/blog/2017Nov06/H3K27me3NarrowPeaks.bed| awk '{print $NF}' | paste $TMPOUT - > $TMPOUT.new
mv $TMPOUT.new $TMPOUT



echo "doing mappability 10kb +- around insertion"
header="$header\tmappability20kb"
bedops -u --range 10000 $TMPOUT | bedtools coverage -a - -b /home/maagj01/scratch/transposon/Analysis/GeneDeserts/hg38.K36.mappable_only.bed |awk '{print $NF}'| paste $TMPOUT - > $TMPOUT.new
mv $TMPOUT.new $TMPOUT



#echo "doing sum of POL2RA +- 50bp from insertion"
#header="$header\tsumPol2RA50bp"
#bedops -u --range 50 $TMPOUT| bedmap --delim '\t' --echo --echo-map --sum - /vol/isg/encode/chipseq/mapped/K562-POLR2A-ENCLB040QQO.hg38/K562-POLR2A-ENCLB040QQO.hg38.perBase.starch|
#awk ' {if ($NF=="NAN") print 0; else print $NF}' | paste $TMPOUT - > $TMPOUT.new
#mv $TMPOUT.new $TMPOUT


echo "doing distance to DHS hotspot clusters"
bedtools makewindows  -w 10000 -s 1000 -b <(awk 'BEGIN {OFS="\t"} {print $1, 0, $2}' /vol/isg/annotation/fasta/hg38/hg38.chrom.sizes|grep -v  'Un\|alt\|random')|
sort-bed -|bedtools intersect -wa -c -a - -b /vol/isg/encode/dnase/mapped/K562-DS9764/hotspots/K562-DS9764.hg38-final/K562-DS9764.hg38.fdr0.01.hot.bed | awk '$4>=4' |bedops -m - >/tmp/DHShotspots_10kb_1kb.bed

bedtools makewindows  -w 50000 -s 5000 -b <(awk 'BEGIN {OFS="\t"} {print $1, 0, $2}' /vol/isg/annotation/fasta/hg38/hg38.chrom.sizes|grep -v  'Un\|alt\|random')|
sort-bed -|bedtools intersect -wa -c -a - -b /vol/isg/encode/dnase/mapped/K562-DS9764/hotspots/K562-DS9764.hg38-final/K562-DS9764.hg38.fdr0.01.hot.bed  | awk '$4>=10'  |bedops -m ->/tmp/DHShotspots_50kb_1kb.bed

bedtools makewindows  -w 100000 -s 10000 -b <(awk 'BEGIN {OFS="\t"} {print $1, 0, $2}' /vol/isg/annotation/fasta/hg38/hg38.chrom.sizes|grep -v  'Un\|alt\|random')|
sort-bed -|bedtools intersect -wa -c -a - -b /vol/isg/encode/dnase/mapped/K562-DS9764/hotspots/K562-DS9764.hg38-final/K562-DS9764.hg38.fdr0.01.hot.bed  | awk '$4>=16'  |bedops -m ->/tmp/DHShotspots_100kb_1kb.bed


header="$header\tdistToDHScluster"
closest-features --delim '\t' --closest --dist --no-ref $TMPOUT /tmp/DHShotspots_50kb_1kb.bed|  awk '{print $NF}'| paste $TMPOUT - > $TMPOUT.new
mv $TMPOUT.new $TMPOUT

header="$header\tdistToDHScluster10kb"
closest-features --delim '\t' --closest --dist --no-ref $TMPOUT /tmp/DHShotspots_10kb_1kb.bed|  awk '{print $NF}'| paste $TMPOUT - > $TMPOUT.new
mv $TMPOUT.new $TMPOUT

header="$header\tdistToDHScluster100kb"
closest-features --delim '\t' --closest --dist --no-ref $TMPOUT /tmp/DHShotspots_100kb_1kb.bed|  awk '{print $NF}'| paste $TMPOUT - > $TMPOUT.new
mv $TMPOUT.new $TMPOUT



bedtools closest -d -k 5 -a $TMPOUT -b /tmp/DHShotspots_10kb_1kb.bed|awk 'BEGIN {OFS="\t"} {print $4, $5, $NF}' > $OUTDIR/${SampleName}_DHShotCluster10kb.bed
bedtools closest -d -k 5 -a $TMPOUT -b /tmp/DHShotspots_50kb_1kb.bed|awk 'BEGIN {OFS="\t"} {print $4, $5, $NF}' > $OUTDIR/${SampleName}_DHShotCluster50kb.bed
bedtools closest -d -k 5 -a $TMPOUT -b /tmp/DHShotspots_100kb_1kb.bed|awk 'BEGIN {OFS="\t"} {print $4, $5, $NF}' > $OUTDIR/${SampleName}_DHShotCluster100kb.bed



bedops -u /vol/isg/encode/dnase/mapped/K562-DS9764/hotspots/K562-DS9764.hg38-final/K562-DS9764.hg38.fdr0.01.hot.bed |
       bedmap --delim '\t' --echo --echo-map --count  - /home/maagj01/scratch/transposon/Analysis/K562_DNaseMasterList/Mauram2016cellReport.hg38.bed| awk '$NF==0'|
       bedmap --delim '\t' --echo --echo-map --count - /vol/isg/encode/chipseq/mapped/K562-H3K4me1-ENCLB060JFJ/hotspots/K562-H3K4me1-ENCLB060JFJ.hg38_noalt-final/K562-H3K4me1-ENCLB060JFJ.hg38_noalt.fdr0.01.hot.starch | awk '$NF==0'| 
       bedmap --delim '\t' --echo --echo-map --count - /vol/isg/encode/chipseq/mapped/K562-H3K4me3-DS11507/hotspots/K562-H3K4me3-DS11507.hg38_noalt-final/K562-H3K4me3-DS11507.hg38_noalt.fdr0.01.hot.starch |awk '$NF==0'| 
       bedmap --delim '\t' --echo --echo-map --count - /vol/isg/annotation/bed/hg38/refseq_gene/refGene.CombinedTxStarts.bed |awk '$NF==0' > $TMPDIR/${SampleName}.silencer.bed

header="$header\tdistToDHSsilencer"
closest-features --delim '\t' --closest --dist --no-ref $TMPOUT $TMPDIR/${SampleName}.silencer.bed|  awk '{print $NF}'| paste $TMPOUT - > $TMPOUT.new
mv $TMPOUT.new $TMPOUT
####
#Takes ages
####
#echo "doing promoter PhastCons +/- 10bp"
#header="$header\tphastCons"
#bedtools flank -g /vol/isg/annotation/bed/hg38/chromInfo.txt -s -i $TMPOUT  -l 1 -r 1 >$TMPDIR/${SampleName}.75bp.bed
#f="/vol/isg/annotation/bed/hg38/phastCons100way/phastCons100way.starch"
##f=/dev/null
#gcat $f | bedmap --faster --echo --mean --prec 3 $TMPDIR/${SampleName}.75bp.bed - | perl -pe 's/\tNAN$/\t0/g;'|
#awk -F'|' 'BEGIN {OFS="\t"} {print $2}'| paste $TMPOUT - > $TMPOUT.new 
#mv $TMPOUT.new $TMPOUT
#
#
#echo "doing promoter phyloP +/- 10bp"
#header="$header\tphyloP"
#f="/vol/isg/annotation/bed/hg38/phyloP100way/phyloP100way.starch"
##f=/dev/null
#gcat $f | bedmap --faster --bp-ovr 1 --echo --mean --prec 3 $TMPDIR/${SampleName}.75bp.bed - | perl -pe 's/\tNAN$/\tNA/g;' |
#awk -F'|' 'BEGIN {OFS="\t"} {print $2}'| paste $TMPOUT - > $TMPOUT.new 
#mv $TMPOUT.new $TMPOUT
#
#Print to file
echo -e $header | cat - <(sort-bed $TMPOUT)  > $OUTDIR/AllInsertions.annotated.bed
######
#4C data/CaptureC - Using normalised matrix from Rao et al. which I realigned with HiC-Pro and normalised with HiTC
#######
echo 'Creating normalised contactmaps in R. Only needs to run once for new data, comment out afterwards'
/home/maagj01/scratch/transposon/Analysis/InsertionRegion_analysis/src/CaptureCinsertionLoci.R  --SampleName ${SampleName} --OUTDIR ${OUTDIR}
echo 'Created contact map with Observed vs Expected contacts +-300kb around insertion'

####
#distance to nearest DHS and it's density
####
paste /vol/isg/encode/dnase/mapped/K562-DS9764/hotspots/K562-DS9764.hg38_noalt-final/K562-DS9764.hg38_noalt.fdr0.01.pks.starch \
/vol/isg/encode/dnase/mapped/K562-DS9764/hotspots/K562-DS9764.hg38-final/K562-DS9764.hg38.fdr0.01.pks.dens.txt \
/vol/isg/encode/dnase/mapped/K562-DS9764/hotspots/K562-DS9764.hg38-final/K562-DS9764.hg38.fdr0.01.pks.pval.txt |
awk 'BEGIN {OFS="\t"} {print $1, $2, $3, "DHS_"NR, $4, $5}'> $TMPDIR/K562_DHS_dens_p.bed



echo "doing summarised DHS density 300kb +- of insertion"
header="$header\tnDHS300kbsum"
bedops -u --range 300000  $TMPOUT| bedmap --delim '\t' --echo --echo-map --sum - $TMPDIR/K562_DHS_dens_p.bed|awk '{print $NF}' |sed 's/NAN/0/g' |paste $TMPOUT -  > $TMPOUT.new 
mv $TMPOUT.new $TMPOUT



awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5}' $OUTDIR/OBSvsExpcontactMap_raw.tsv > $OUTDIR/OBSvsExpcontactMap.tsv 
cat  $OUTDIR/OBSvsExpcontactMap.tsv|sort-bed -|
bedmap --echo --echo-map --delim '\t' --skip-unmapped --max - $TMPDIR/K562_DHS_dens_p.bed| awk 'BEGIN {OFS="\t"} $4>2 && $4!="NA" {print $1, $2, $3, $4, $5, $NF}' > $TMPDIR/${SampleName}_K562_DHS_dens_p_Obs2.bed


#Sum of DHS peaks in contact
cat  $OUTDIR/OBSvsExpcontactMap.tsv|sort-bed -|
bedmap --echo --echo-map --delim '\t' --skip-unmapped --sum - $TMPDIR/K562_DHS_dens_p.bed| awk 'BEGIN {OFS="\t"} $4>2 && $4!="NA" {print $1, $2, $3, $4, $5, $NF}' > $TMPDIR/${SampleName}_K562_DHS_SUM_Obs2.bed

#echo 'MVC of DHS in strongest contact with insert'
##First add K562 to the merged DHS list
bedops -u /home/mauram01/scratch/hybridmice/dnase/lineagemcv/all.lineage.dhs.unmerged.starch <(awk 'BEGIN {OFS="\t"} {print $1, $2, $3, "K562-DS9764"}' /vol/isg/encode/dnase/mapped/K562-DS9764/hotspots/K562-DS9764.hg38_noalt-final/K562-DS9764.hg38_noalt.fdr0.01.pks.starch)|
awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4}'| sort-bed - | 
bedmap  --delim '\t' --skip-unmapped --echo --fraction-ref 1.0 --count  /vol/isg/encode/dnase/mapped/K562-DS9764/hotspots/K562-DS9764.hg38_noalt-final/K562-DS9764.hg38_noalt.fdr0.01.pks.starch - | 
awk 'BEGIN {OFS="\t"} {print $1, $2, $3, NR, $4}' > $TMPDIR/DHS_MCV.bed

awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5}' $OUTDIR/OBSvsExpcontactMap_raw.tsv > $OUTDIR/OBSvsExpcontactMap.tsv 
cat  $OUTDIR/OBSvsExpcontactMap.tsv|sort-bed -|
bedmap  --delim '\t' --skip-unmapped --echo  --max  - $TMPDIR/DHS_MCV.bed| awk 'BEGIN {OFS="\t"} $4>2 && $4!="NA" {print $1, $2, $3, $4, $5, $NF}' > $TMPDIR/${SampleName}_DHS_MCV.bed



#####
#Create files with the distance to the nearest 10 DHS types for analysis
#####
echo 'creating distance to 10 nearest DHS type files'

paste /vol/isg/encode/dnase/mapped/K562-DS9764/hotspots/K562-DS9764.hg38_noalt-final/K562-DS9764.hg38_noalt.fdr0.01.pks.starch /vol/isg/encode/dnase/mapped/K562-DS9764/hotspots/K562-DS9764.hg38-final/K562-DS9764.hg38.fdr0.01.pks.dens.txt |
awk 'BEGIN {OFS="\t"} {print $1, $2, $3, "DHS_"NR, $4}' > $TMPDIR/K562DHSaccess.bed 
unstarch /home/maagj01/scratch/transposon/Analysis/K562_DNaseMasterList/DHStype_masterlist.starch > $TMPDIR/DHStype_masterlist.bed

bedtools closest -d -k 5 -a $TMPOUT -b <(grep Promoter $TMPDIR/DHStype_masterlist.bed) |awk 'BEGIN {OFS="\t"} {print $4, $5, $NF, "Promoter"}' > $OUTDIR/${SampleName}_DHSpromoter.bed
bedtools closest -d -k 5 -a $TMPOUT -b <(grep CTCF $TMPDIR/DHStype_masterlist.bed) |awk 'BEGIN {OFS="\t"} {print $4, $5, $NF, "CTCF"}' > $OUTDIR/${SampleName}_DHSctcf.bed
bedtools closest -d -k 5 -a $TMPOUT -b <(grep Distal $TMPDIR/DHStype_masterlist.bed) |awk 'BEGIN {OFS="\t"} {print $4, $5, $NF, "Distal"}' > $OUTDIR/${SampleName}_DHSdistal.bed
bedtools closest -d -k 5 -a $TMPOUT -b /vol/isg/encode/dnase/mapped/K562-DS9764/hotspots/K562-DS9764.hg38_noalt-final/K562-DS9764.hg38_noalt.fdr0.01.pks.starch |awk 'BEGIN {OFS="\t"} {print $4, $5, $NF, "all"}' > $OUTDIR/${SampleName}_DHSall.bed
bedtools closest -d -k 5 -a $TMPOUT -b <(grep -v CTCF $TMPDIR/DHStype_masterlist.bed) |awk 'BEGIN {OFS="\t"} {print $4, $5, $NF, "nonCTCF"}' > $OUTDIR/${SampleName}_nonCTCF.bed

bedtools closest -d -k 5 -a $TMPOUT -b  $TMPDIR/DHS_MCV.bed | awk 'BEGIN {OFS="\t"} {print  $4, $5, $(NF-1), $NF}' > $OUTDIR/${SampleName}_MCV.bed
bedtools closest -d -k 5 -a $TMPOUT -b  /vol/isg/encode/dnase/mapped/K562-DS9764/hotspots/K562-DS9764.hg38-final/K562-DS9764.hg38.fdr0.01.hot.bed | awk 'BEGIN {OFS="\t"} {print  $4, $5, $(NF-1), $NF}' > $OUTDIR/${SampleName}_Hotspots.bed
bedtools closest -d -k 21 -a $OUTDIR/AllInsertions.coords.bed -b  <(awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4 }' $OUTDIR/AllInsertions.coords.bed) | awk 'BEGIN {OFS="\t"} {print  $4, $(NF-1), $NF}' > $OUTDIR/${SampleName}_neighbours.bed




echo 'Nearest 20 neighbours and CTCF sites in between'
bedtools closest -d -k 21 -a $OUTDIR/AllInsertions.coords.bed -b  <(awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4 }' $OUTDIR/AllInsertions.coords.bed) |awk '$4!=$(NF-1) && $NF!=0'| 
awk 'BEGIN {OFS="\t"} {if($3<$(NF-3)) print $1, $2, $(NF-3), $4, $(NF-1), $NF; else if ($3>$(NF-3)) print $1, $(NF-3), $2, $(NF-1), $4, $NF }'|
sort-bed - |uniq| bedmap --delim '\t' --echo --echo-map --count - /home/maagj01/scratch/transposon/Analysis/K562_DNaseMasterList/Mauram2016cellReport.hg38.bed|
awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5, $6, $NF}' > $OUTDIR/AllNeighbour_CTCF.bed


#echo 'Nearest 20 neighbours and CTCF sites in between and MCV of CTCF'
#bedmap  --delim '\t' --skip-unmapped --echo --fraction-ref 1.0 --count  /home/maagj01/scratch/transposon/Analysis/K562_DNaseMasterList/Mauram2016cellReport.hg38.bed /home/maagj01/scratch/transposon/Analysis/CTCF_MCV/Tissue_CTCF_unmerged.bed|
#awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $NF}' > $TMPDIR/CTCF_MCV_tissues.bed
#
bedtools closest -d -k 21 -a $OUTDIR/AllInsertions.coords.bed -b  <(awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4 }' $OUTDIR/AllInsertions.coords.bed) |awk '$4!=$(NF-1) && $NF!=0'| 
awk 'BEGIN {OFS="\t"} {if($3<$(NF-3)) print $1, $2, $(NF-3), $4, $(NF-1), $NF; else if ($3>$(NF-3)) print $1, $(NF-3), $2, $(NF-1), $4, $NF }'|
sort-bed - |uniq| bedmap --delim '\t' --echo --echo-map --count - $TMPDIR/CTCF_MCV_tissues.bed| 
awk 'BEGIN {OFS="\t"} $NF<=1 {if ($NF==0) print $1, $2, $3, $4, $5, $6, 0, 0; else print $1, $2, $3, $4, $5, $6, $(NF-1), $NF}' > $OUTDIR/AllNeighbour_CTCF_MCV.bed


echo 'Nearest 20 neighbours and CTCF sites in between and around, including Distance to CTCF sites'
closest-features --delim '\t' --dist $OUTDIR/AllInsertions.coords.bed /home/maagj01/scratch/transposon/Analysis/K562_DNaseMasterList/Mauram2016cellReport.hg38.bed|
awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $(NF-11), $NF}'| sort-bed - > $TMPDIR/AllNeighbour_surroundingCTCF.bed

bedtools closest -d -k 21 -a $TMPDIR/AllNeighbour_surroundingCTCF.bed -b  $TMPDIR/AllNeighbour_surroundingCTCF.bed |
awk 'BEGIN {OFS="\t"} $4!=$10 && $NF!=0' |
awk 'BEGIN {OFS="\t"} {if($3<$(NF-5)) print $1, $2, $(NF-5), $4, $5, $6, $10, $11, $12, $NF; else if ($3>$(NF-4)) print $1, $(NF-5), $2,  $(NF-3), $11, $12,  $4, $5, $6, $NF}'|
sort-bed -| uniq| 
bedmap --delim '\t' --echo --echo-map --count - /home/maagj01/scratch/transposon/Analysis/K562_DNaseMasterList/Mauram2016cellReport.hg38.bed |
awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $NF}' > $TMPDIR/AllNeighbour_surroundingCTCF.annotated.bed

bedmap --delim '\t' --echo --echo-map --count $TMPDIR/AllNeighbour_surroundingCTCF.annotated.bed  /vol/isg/encode/dnase/mapped/K562-DS9764/hotspot2/K562-DS9764.hg38.hotspots.fdr0.05.starch| awk '{print $NF}'|
paste $TMPDIR/AllNeighbour_surroundingCTCF.annotated.bed - >  $TMPDIR/AllNeighbour_surroundingCTCF.bed 
echo -e 'chrom\tIns\tneighIns\tInsBC\tInsUpCTCF\tInsDownCTCF\tneighBC\tneighUpCTCF\tneighDownCTCF\tDistance\tnCTCF\tnDHS'|cat -  $TMPDIR/AllNeighbour_surroundingCTCF.bed >$OUTDIR/AllNeighbour_surroundingCTCF.bed 



neighheader="chrom\tchromStart\tchromEnd\tBC\tneighBC\tDistance\tnCTCF"
bedtools closest -d -k 2 -a $TMPOUT -b  <(awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4 }' $TMPOUT) |awk '$4!=$(NF-1) && $NF!=0'| 
awk 'BEGIN {OFS="\t"} {if($3<$(NF-3)) print $1, $2, $(NF-3), $4, $(NF-1), $NF; else if ($3>$(NF-3)) print $1, $(NF-3), $2, $(NF-1), $4, $NF }'|
uniq| sort-bed - |bedmap --delim '\t' --echo --echo-map --count - /home/maagj01/scratch/transposon/Analysis/K562_DNaseMasterList/Mauram2016cellReport.hg38.bed |
awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5, $6, $NF}' > $TMPDIR/Neighbour_CTCF.bed.new 

mv $TMPDIR/Neighbour_CTCF.bed.new $TMPDIR/Neighbour_CTCF.new2.bed

neighheader="$neighheader\tnDHS"
bedmap --delim '\t' --echo --echo-map --count  $TMPDIR/Neighbour_CTCF.new2.bed /vol/isg/encode/dnase/mapped/K562-DS9764/hotspots/K562-DS9764.hg38_noalt-final/K562-DS9764.hg38_noalt.fdr0.01.pks.starch|
awk '{print $NF}'|paste $TMPDIR/Neighbour_CTCF.new2.bed - >  $TMPDIR/Neighbour_CTCF.bed.new 
mv $TMPDIR/Neighbour_CTCF.bed.new $TMPDIR/Neighbour_CTCF.new2.bed

neighheader="$neighheader\tnCGI"
bedmap --delim '\t' --echo --echo-map --count  $TMPDIR/Neighbour_CTCF.new2.bed /vol/isg/annotation/bed/hg38/cpg_islands/cpgIslands.bed |
awk '{print $NF}'|paste $TMPDIR/Neighbour_CTCF.new2.bed - >  $TMPDIR/Neighbour_CTCF.bed.new 
mv $TMPDIR/Neighbour_CTCF.bed.new $TMPDIR/Neighbour_CTCF.new2.bed


echo -e $neighheader | cat - $TMPDIR/Neighbour_CTCF.new2.bed > $OUTDIR/Neighbour_CTCF.bed




neighheader="chrom\tchromStart\tchromEnd\tBC\tneighBC\tDistance\tnCTCF"
bedtools closest -d -k 21 -a $TMPOUT -b  <(awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4 }' $TMPOUT) |awk '$4!=$(NF-1) && $NF!=0'| 
awk 'BEGIN {OFS="\t"} {if($3<$(NF-3)) print $1, $2, $(NF-3), $4, $(NF-1), $NF; else if ($3>$(NF-3)) print $1, $(NF-3), $2, $(NF-1), $4, $NF }'|
 sort-bed - |uniq| bedmap --delim '\t' --echo --echo-map --count - /home/maagj01/scratch/transposon/Analysis/K562_DNaseMasterList/Mauram2016cellReport.hg38.bed |
awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5, $6, $NF}' > $TMPDIR/Neighbour_CTCF.bed.new 

mv $TMPDIR/Neighbour_CTCF.bed.new $TMPDIR/Neighbour_CTCF.new2.bed

neighheader="$neighheader\tnDHS"
bedmap --delim '\t' --echo --echo-map --count  $TMPDIR/Neighbour_CTCF.new2.bed /vol/isg/encode/dnase/mapped/K562-DS9764/hotspots/K562-DS9764.hg38_noalt-final/K562-DS9764.hg38_noalt.fdr0.01.pks.starch|
awk '{print $NF}'|paste $TMPDIR/Neighbour_CTCF.new2.bed - >  $TMPDIR/Neighbour_CTCF.bed.new 
mv $TMPDIR/Neighbour_CTCF.bed.new $TMPDIR/Neighbour_CTCF.new2.bed

neighheader="$neighheader\tnCGI"
bedmap --delim '\t' --echo --echo-map --count  $TMPDIR/Neighbour_CTCF.new2.bed /vol/isg/annotation/bed/hg38/cpg_islands/cpgIslands.bed |
awk '{print $NF}'|paste $TMPDIR/Neighbour_CTCF.new2.bed - >  $TMPDIR/Neighbour_CTCF.bed.new 
mv $TMPDIR/Neighbour_CTCF.bed.new $TMPDIR/Neighbour_CTCF.new2.bed


echo -e $neighheader | cat - $TMPDIR/Neighbour_CTCF.new2.bed > $OUTDIR/Neighbour_20neigh_CTCF.bed




echo 'creating 10kb windows over 100kb and count nDHSbins, tagsDHSbins, maxTagsDHSbins, totalTagsbins'
echo '       ---nDHS'
bedtools intersect -wa -wb -a <(bedtools flank -i $TMPOUT -l 50000 -r 50000 -g /vol/isg/annotation/fasta/hg38/hg38.chrom.sizes|sort-bed -) -b <(bedtools flank -i $TMPOUT -l 50000 -r 50000 -g /vol/isg/annotation/fasta/hg38/hg38.chrom.sizes|sort-bed -| bedtools makewindows -b - -w 10000 -i src)|
awk '$4==$NF' |
awk 'BEGIN {OFS="\t"} {print $(NF-3), $(NF-2), $(NF-1), $NF, $5}'|sort-bed -| 
bedmap --delim '\t'  --echo --echo-map --count - /vol/isg/encode/dnase/mapped/K562-DS9764/hotspots/K562-DS9764.hg38_noalt-final/K562-DS9764.hg38_noalt.fdr0.01.pks.starch |
awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5, $NF}'  >${OUTDIR}/${SampleName}_linearDistanceDHS_bedtools.bed

echo '       ---TagsDHS'
bedtools intersect -wa -wb -a <(bedtools flank -i $TMPOUT -l 50000 -r 50000 -g /vol/isg/annotation/fasta/hg38/hg38.chrom.sizes|sort-bed -) -b <(bedtools flank -i $TMPOUT -l 50000 -r 50000 -g /vol/isg/annotation/fasta/hg38/hg38.chrom.sizes|sort-bed -| bedtools makewindows -b - -w 10000 -i src)|
awk '$4==$NF' |
awk 'BEGIN {OFS="\t"} {print $(NF-3), $(NF-2), $(NF-1), $NF, $5}'|sort-bed -| 
bedmap --delim '\t'  --echo --echo-map --sum - $TMPDIR/K562DHSaccess.bed  |
awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5, $NF}'  > ${OUTDIR}/${SampleName}_DistDensity_bedtools.bed


echo '       ---maxTagsDHS'
bedtools intersect -wa -wb -a <(bedtools flank -i $TMPOUT -l 50000 -r 50000 -g /vol/isg/annotation/fasta/hg38/hg38.chrom.sizes|sort-bed -) -b <(bedtools flank -i $TMPOUT -l 50000 -r 50000 -g /vol/isg/annotation/fasta/hg38/hg38.chrom.sizes|sort-bed -| bedtools makewindows -b - -w 10000 -i src)|
awk '$4==$NF' |
awk 'BEGIN {OFS="\t"} {print $(NF-3), $(NF-2), $(NF-1), $NF, $5}'|sort-bed -| 
bedmap --delim '\t'  --echo --echo-map --max - $TMPDIR/K562DHSaccess.bed  |
awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5, $NF}'  > ${OUTDIR}/${SampleName}_DistDensityMax_bedtools.bed
       
#####
#Compare window activity by DHS from hostposts or directly from bam file
#####
echo '       ---totalTags'
echo '       ---Total DNase tags'
bedtools intersect -wa -wb -a <(bedtools flank -i $TMPOUT -l 50000 -r 50000 -g /vol/isg/annotation/fasta/hg38/hg38.chrom.sizes|sort-bed -) -b <(bedtools flank -i $TMPOUT -l 50000 -r 50000 -g /vol/isg/annotation/fasta/hg38/hg38.chrom.sizes|sort-bed -| bedtools makewindows -b - -w 10000 -i src)|
awk '$4==$NF' |
awk 'BEGIN {OFS="\t"} {print $(NF-3), $(NF-2), $(NF-1), $NF, $5}'|sort-bed - > ${OUTDIR}/${SampleName}_10kbWindows.bed

echo 'bam to bed'
samtools view /vol/isg/encode/dnase/mapped/K562-DS9764/K562-DS9764.hg38.bam  | convert2bed --input=sam -  > $TMPDIR/K562-DS9764.hg38.bed 
awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4}'  $TMPDIR/K562-DS9764.hg38.bed > $TMPDIR/K562-DS9764.hg38.clean.bed  
bedtools intersect -sorted -c -wa -a ${OUTDIR}/${SampleName}_10kbWindows.bed -b $TMPDIR/K562-DS9764.hg38.clean.bed  >${OUTDIR}/${SampleName}_10kbWindows_tagsFromBam.bed


echo '       ---MCV'
bedtools intersect -wa -wb -a <(bedtools flank -i $TMPOUT -l 50000 -r 50000 -g /vol/isg/annotation/fasta/hg38/hg38.chrom.sizes|sort-bed -) -b <(bedtools flank -i $TMPOUT -l 50000 -r 50000 -g /vol/isg/annotation/fasta/hg38/hg38.chrom.sizes|sort-bed -| bedtools makewindows -b - -w 10000 -i src)|
awk '$4==$NF' |
awk 'BEGIN {OFS="\t"} {print $(NF-3), $(NF-2), $(NF-1),  $5, $NF}'|sort-bed -| 
bedmap  --delim '\t' --echo  --echo-map  --max -  <( awk 'BEGIN {OFS="\t"} {print $1, $2, $3, "DHS_"NR, $NF}' $TMPDIR/DHS_MCV.bed)  |
awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5, $NF}'  > ${OUTDIR}/${SampleName}_Dist_MCV.bed



#echo '       ---maxTagsDHSnoCTCF'
#bedmap --delim '\t' --echo --echo-map --skip-unmapped /tmp/K562_DHS_dens_p.bed <(grep -v CTCF $TMPDIR/DHStype_masterlist.bed)| awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5}' > $TMPDIR/DHSdens_nonCTCF.bed
#
#
#bedtools intersect -wa -wb -a <(bedtools flank -i $TMPOUT -l 50000 -r 50000 -g /vol/isg/annotation/fasta/hg38/hg38.chrom.sizes|sort-bed -) -b <(bedtools flank -i $TMPOUT -l 50000 -r 50000 -g /vol/isg/annotation/fasta/hg38/hg38.chrom.sizes|sort-bed -| bedtools makewindows -b - -w 10000 -i src)|
#awk '$4==$NF' |
#awk 'BEGIN {OFS="\t"} {print $(NF-3), $(NF-2), $(NF-1), $NF, $5}'|sort-bed -| 
#bedmap --delim '\t'  --echo --echo-map --max - $TMPDIR/DHSdens_nonCTCF.bed |
#awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5, $NF}'  > ${OUTDIR}/${SampleName}_Dist_maxDHSnonCTCF.bed
#
#
#bedtools intersect -wa -wb -a <(bedtools flank -i $TMPOUT -l 50000 -r 50000 -g /vol/isg/annotation/fasta/hg38/hg38.chrom.sizes|sort-bed -) -b <(bedtools flank -i $TMPOUT -l 50000 -r 50000 -g /vol/isg/annotation/fasta/hg38/hg38.chrom.sizes|sort-bed -| bedtools makewindows -b - -w 10000 -i src)|
#awk '$4==$NF' |
#awk 'BEGIN {OFS="\t"} {print $(NF-3), $(NF-2), $(NF-1), $NF, $5}'|sort-bed -| 
#bedmap --delim '\t'  --echo --echo-map --count - $TMPDIR/DHSdens_nonCTCF.bed |
#awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5, $NF}'  > ${OUTDIR}/${SampleName}_Dist_nDHSnonCTCF.bed
#
#
#bedtools intersect -wa -wb -a <(bedtools flank -i $TMPOUT -l 50000 -r 50000 -g /vol/isg/annotation/fasta/hg38/hg38.chrom.sizes|sort-bed -) -b <(bedtools flank -i $TMPOUT -l 50000 -r 50000 -g /vol/isg/annotation/fasta/hg38/hg38.chrom.sizes|sort-bed -| bedtools makewindows -b - -w 10000 -i src)|
#awk '$4==$NF' |
#awk 'BEGIN {OFS="\t"} {print $(NF-3), $(NF-2), $(NF-1), $NF, $5}'|sort-bed -| 
#bedmap --delim '\t'  --echo --echo-map --sum - $TMPDIR/DHSdens_nonCTCF.bed |
#awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5, $NF}'  > ${OUTDIR}/${SampleName}_Dist_sumDHSnonCTCF.bed


echo 'creating 10kb windows over 100kb and count nDHSbins, tagsDHSbins, maxTagsDHSbins,totalTagsbins for hotspots instead peaks'
echo '       ---nDHS'
bedtools intersect -wa -wb -a <(bedtools flank -i $TMPOUT -l 50000 -r 50000 -g /vol/isg/annotation/fasta/hg38/hg38.chrom.sizes|sort-bed -) -b <(bedtools flank -i $TMPOUT -l 50000 -r 50000 -g /vol/isg/annotation/fasta/hg38/hg38.chrom.sizes|sort-bed -| bedtools makewindows -b - -w 10000 -i src)|
awk '$4==$NF' |
awk 'BEGIN {OFS="\t"} {print $(NF-3), $(NF-2), $(NF-1), $NF, $5}'|sort-bed -| 
bedmap --delim '\t'  --echo --echo-map --count - /vol/isg/encode/dnase/mapped/K562-DS9764/hotspots/K562-DS9764.hg38-final/K562-DS9764.hg38.fdr0.01.hot.bed |
awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5, $NF}'  >${OUTDIR}/${SampleName}_nHotspots.bed



#echo 'creating 2kb +- window from insertion and counting Tags from histone marks and from DHS'
##bedtools intersect -wa -wb -a <(bedtools flank -i $TMPOUT -l 2000 -r 2000 -g /vol/isg/annotation/fasta/hg38/hg38.chrom.sizes|sort-bed -) -b <(bedtools flank -i $TMPOUT -l 2000 -r 2000 -g /vol/isg/annotation/fasta/hg38/hg38.chrom.sizes|sort-bed -| 
##bedtools makewindows -b - -w 50 -i src)| awk 'BEGIN {OFS="\t"} $4==$NF {print $(NF-3), $(NF-2), $(NF-1), $NF}' |sort-bed - > $TMPDIR/${SampleName}_promoter.bed
#
#echo 'bam to bed'
#samtools view /vol/isg/encode/dnase/mapped/K562-DS9764/K562-DS9764.hg38.bam  | convert2bed --input=sam -  > $TMPDIR/K562-DS9764.hg38.bed 
#awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4}'  $TMPDIR/K562-DS9764.hg38.bed > $TMPDIR/K562-DS9764.hg38.clean.bed  
#
#samtools view /vol/isg/encode/chipseq/mapped/K562-H3K9ac-ENCLB695ASG/K562-H3K9ac-ENCLB695ASG.hg38.bam | convert2bed --input=sam -  > $TMPDIR/K562-H3K9ac-ENCLB695ASG.hg38.bed 
#awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4}'  $TMPDIR/K562-H3K9ac-ENCLB695ASG.hg38.bed > $TMPDIR/K562-H3K9ac-ENCLB695ASG.hg38.clean.bed  
#
#samtools view /vol/isg/encode/chipseq/mapped/K562-H3K4me3-DS11507/K562-H3K4me3-DS11507.hg38_noalt.bam | convert2bed --input=sam -  > $TMPDIR/K562-H3K4me3-DS11507.hg38.bed 
#awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4}'  $TMPDIR/K562-H3K4me3-DS11507.hg38.bed > $TMPDIR/K562-H3K4me3-DS11507.hg38.clean.bed  
#
#samtools view /vol/isg/encode/chipseq/mapped/K562-H3K4me2-ENCLB695AGB/K562-H3K4me2-ENCLB695AGB.hg38_noalt.bam | convert2bed --input=sam -  > $TMPDIR/K562-H3K4me2-ENCLB695AGB.hg38.bed 
#awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4}'  $TMPDIR/K562-H3K4me2-ENCLB695AGB.hg38.bed > $TMPDIR/K562-H3K4me2-ENCLB695AGBG.hg38.clean.bed  
#
#samtools view /vol/isg/encode/chipseq/mapped/K562-H3K79me2-ENCLB695AHH/K562-H3K79me2-ENCLB695AHH.hg38_noalt.bam | convert2bed --input=sam -  > $TMPDIR/K562-H3K79me2-ENCLB695AHH.hg38.bed 
#awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4}'  $TMPDIR/K562-H3K79me2-ENCLB695AHH.hg38.bed > $TMPDIR/K562-H3K79me2-ENCLB695AHH.hg38.clean.bed  
#
#samtools view /vol/isg/encode/chipseq/mapped/K562-H3K27ac-ENCLB695AJD/K562-H3K27ac-ENCLB695AJD.hg38_noalt.bam| convert2bed --input=sam -  > $TMPDIR/K562-H3K27ac-ENCLB695AJD.hg38.bed 
#awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4}'  $TMPDIR/K562-H3K27ac-ENCLB695AJD.hg38.bed > $TMPDIR/K562-H3K27ac-ENCLB695AJD.hg38.clean.bed  
#
#echo 'Count Tags'
#tagHeader="chrom\tstart\tend\tBC"
#bedops -u --range 2000 $TMPOUT|awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4}' >$TMPDIR/${SampleName}_promoter.bed
#
#tagHeader="$tagHeader\tDNaseTags"
#bedmap --delim '\t'  --echo --echo-map --count  $TMPDIR/${SampleName}_promoter.bed $TMPDIR/K562-DS9764.hg38.clean.bed |
#awk '{print $NF}'  |paste $TMPDIR/${SampleName}_promoter.bed - > $TMPDIR/${SampleName}_promoter.bed.new 
#mv $TMPDIR/${SampleName}_promoter.bed.new $TMPDIR/${SampleName}_promoter.bed
#
#tagHeader="$tagHeader\tH3K9acTags"
#bedmap --delim '\t'  --echo --echo-map --count  $TMPDIR/${SampleName}_promoter.bed $TMPDIR/K562-H3K9ac-ENCLB695ASG.hg38.clean.bed |
#awk '{print $NF}'  |paste $TMPDIR/${SampleName}_promoter.bed - > $TMPDIR/${SampleName}_promoter.bed.new 
#mv $TMPDIR/${SampleName}_promoter.bed.new $TMPDIR/${SampleName}_promoter.bed
#
#
#tagHeader="$tagHeader\tH3K4me3Tags"
#bedmap --delim '\t'  --echo --echo-map --count  $TMPDIR/${SampleName}_promoter.bed $TMPDIR/K562-H3K4me3-DS11507.hg38.clean.bed |
#awk '{print $NF}'  |paste $TMPDIR/${SampleName}_promoter.bed - > $TMPDIR/${SampleName}_promoter.bed.new 
#mv $TMPDIR/${SampleName}_promoter.bed.new $TMPDIR/${SampleName}_promoter.bed
#
#tagHeader="$tagHeader\tH3K4me2Tags"
#bedmap --delim '\t'  --echo --echo-map --count  $TMPDIR/${SampleName}_promoter.bed $TMPDIR/K562-H3K4me2-ENCLB695AGBG.hg38.clean.bed  |
#awk '{print $NF}'  |paste $TMPDIR/${SampleName}_promoter.bed - > $TMPDIR/${SampleName}_promoter.bed.new 
#mv $TMPDIR/${SampleName}_promoter.bed.new $TMPDIR/${SampleName}_promoter.bed
#
#tagHeader="$tagHeader\tH3K79me2Tags"
#bedmap --delim '\t'  --echo --echo-map --count  $TMPDIR/${SampleName}_promoter.bed $TMPDIR/K562-H3K79me2-ENCLB695AHH.hg38.clean.bed |
#awk '{print $NF}'  |paste $TMPDIR/${SampleName}_promoter.bed - > $TMPDIR/${SampleName}_promoter.bed.new 
#mv $TMPDIR/${SampleName}_promoter.bed.new $TMPDIR/${SampleName}_promoter.bed
#
#tagHeader="$tagHeader\tH3K27acTags"
#bedmap --delim '\t'  --echo --echo-map --count  $TMPDIR/${SampleName}_promoter.bed $TMPDIR/K562-H3K27ac-ENCLB695AJD.hg38.clean.bed |
#awk '{print $NF}'  |paste $TMPDIR/${SampleName}_promoter.bed - > $TMPDIR/${SampleName}_promoter.bed.new 
#
#
#echo -e $tagHeader | cat - <(sort-bed $TMPDIR/${SampleName}_promoter.bed.new )  >$OUTDIR/histoneDNaseTags.bed



samtools view  /vol/isg/encode/chipseq/mapped/K562-H3K27me3-DS12066/K562-H3K27me3-DS12066.hg38_noalt.bam| convert2bed --input=sam -  > $TMPDIR/K562-H3K27me3-DS12066.hg38.bed 
awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4}'  $TMPDIR/K562-H3K27me3-DS12066.hg38.bed > $TMPDIR/K562-H3K27me3-DS12066.hg38.clean.bed  

echo 'Count Tags'
bedops -u --range 10000 $TMPOUT|awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4}' >$TMPDIR/${SampleName}_promoter.bed

header="$header\tH3K27me3Tags10kb"
bedops -u --range 50000 $TMPOUT| bedmap --delim '\t'  --echo --echo-map --count  - $TMPDIR/K562-H3K27me3-DS12066.hg38.clean.bed  |
awk '{print $NF}'  |paste $TMPOUT - > $TMPOUT.new 
mv $TMPOUT.new $TMPOUT


####
#Find max DHS in contact with insertion
####
echo "doing DNase peak in 5kb windows in contact based on highest DHS density"

R --quiet --no-save << EOF
library(dplyr)
library(data.table)

DHSdens <- fread("$TMPDIR/${SampleName}_K562_DHS_dens_p_Obs2.bed")
#Highest DHS density
DHSdensMax <- DHSdens %>%
       group_by(V5) %>%
       filter(V6 == max(V6)) %>% #Get window with max DHS density 
       filter(V4 == max(V4)) # if max DHS is the same for multiple windows, get the one with strongest contact
write.table(DHSdensMax, file="$TMPDIR/${SampleName}_K562_DHS_dens_p_Obs2_maxDHS.bed", sep='\t', col.names=F, row.names=F, quote=F)

#Strongest contact DHS density
DHScontactMax <- DHSdens %>%
       group_by(V5) %>%
       filter(V4 == max(V4)) %>% #Get window with strongest HiC contact
       filter(V6 == max(V6)) # if max DHS is the same for multiple windows, get the one with strongest contact
write.table(DHScontactMax, file="$TMPDIR/${SampleName}_K562_DHS_dens_p_Obs2_maxContact.bed", sep='\t', col.names=F, row.names=F, quote=F)

##MCV
DHSmcv <- fread("$TMPDIR/${SampleName}_DHS_MCV.bed")
DHScontactMCV <- DHSmcv %>%
       group_by(V5) %>%
       filter(V4 == max(V4)) %>% #Get window with strongest HiC contact
       filter(V6 == max(V6)) # if max DHS is the same for multiple windows, get the one with strongest contact
write.table(DHScontactMCV, file="$TMPDIR/${SampleName}_K562_MCV.bed", sep='\t', col.names=F, row.names=F, quote=F)

#Sum of contact DHS density
DHSsum <- fread("$TMPDIR/${SampleName}_K562_DHS_SUM_Obs2.bed")
DHScontactSum <- DHSsum %>%
       group_by(V5) %>%
       summarise (V6 = sum(V6))
write.table(DHScontactSum, file="$TMPDIR/${SampleName}_K562_DHS_SUM.bed", sep='\t', col.names=F, row.names=F, quote=F)

EOF

echo 'finding distance to Max DHS peak'
join -eNA -a1 -a2 -1 4 <(sort -k4,4 $TMPOUT) -2 5 <(sort -k5,5 $TMPDIR/${SampleName}_K562_DHS_dens_p_Obs2_maxDHS.bed) -o 1.1,1.2,1.3,1.4,2.1,2.2,2.3,2.4,2.5,2.6|grep -v NA |
awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}'| 
#awk 'BEGIN {OFS="\t"} $(NF-4)> $3 && $(NF-3) >$2 {print $1, $2, $3, $4, $(NF-4)-$3, $(NF-2), $(NF-1), $NF, $(NF-4), $(NF-3)}'
awk 'BEGIN {OFS="\t"} { if ($(NF-4)> $3 && $(NF-3) >$2) print $1, $2, $(NF-3), $4, $(NF-4)-$3, $(NF-2), $(NF-1), $NF; else if ($(NF-4)< $3 && $(NF-3) <$2) print $1, $(NF-4), $3, $4, $2-$(NF-3), $(NF-2), $(NF-1), $NF}'|
sort-bed - >$OUTDIR/${SampleName}_DHS_dens_max_distance.bed


header="$header\tmaxDHScontact\tmaxDHSdensity\tmaxDHSdistance"
join -eNA -a1 -a2 -1 4 <(sort -k4,4 $TMPOUT) -2 4 <(sort -k4,4 $OUTDIR/${SampleName}_DHS_dens_max_distance.bed) -o 2.6,2.8,2.5|
awk 'BEGIN {OFS="\t"} {print $1, $2, $3}' |paste <(sort -k4,4 $TMPOUT) -| sort-bed - >  $TMPOUT.new
mv $TMPOUT.new $TMPOUT

echo "doing DNase peak in 5kb windows in contact based on strongest HiC contact"

echo 'finding distance to strongest HiC contact'
join -eNA -a1 -a2 -1 4 <(sort -k4,4 $TMPOUT) -2 5 <(sort -k5,5 $TMPDIR/${SampleName}_K562_DHS_dens_p_Obs2_maxContact.bed) -o 1.1,1.2,1.3,1.4,2.1,2.2,2.3,2.4,2.5,2.6|grep -v NA |
awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}'| 
#awk 'BEGIN {OFS="\t"} $(NF-4)> $3 && $(NF-3) >$2 {print $1, $2, $3, $4, $(NF-4)-$3, $(NF-2), $(NF-1), $NF, $(NF-4), $(NF-3)}'
awk 'BEGIN {OFS="\t"} { if ($(NF-4)> $3 && $(NF-3) >$2) print $1, $2, $(NF-3), $4, $(NF-4)-$3, $(NF-2), $(NF-1), $NF; else if ($(NF-4)< $3 && $(NF-3) <$2) print $1, $(NF-4), $3, $4, $2-$(NF-3), $(NF-2), $(NF-1), $NF}'|
sort-bed - > $OUTDIR/${SampleName}_DHS_contact_max_distance.bed


header="$header\tmaxContactContact\tmaxContactDHSdensity\tmaxConcatdistance"
join -eNA -a1 -a2 -1 4 <(sort -k4,4 $TMPOUT) -2 4 <(sort -k4,4 $OUTDIR/${SampleName}_DHS_contact_max_distance.bed) -o 2.6,2.8,2.5|
awk 'BEGIN {OFS="\t"} {print $1, $2, $3}' |paste <(sort -k4,4 $TMPOUT) - |sort-bed - >  $TMPOUT.new
mv $TMPOUT.new $TMPOUT



header="$header\tmaxContactMCVcontact\tmaxContactMCV"
join -eNA -a1 -a2 -1 4 <(sort -k4,4 $TMPOUT) -2 5 <(sort -k5,5 $TMPDIR/${SampleName}_K562_MCV.bed) -o 2.4,2.6|
awk 'BEGIN {OFS="\t"} {print $1, $2}' |paste <(sort -k4,4 $TMPOUT) - |sort-bed - >  $TMPOUT.new
mv $TMPOUT.new $TMPOUT


header="$header\tmaxContactDHSsum"
join -eNA -a1 -a2 -1 4 <(sort -k4,4 $TMPOUT) -2 1 <(sort -k1,1 $TMPDIR/${SampleName}_K562_DHS_SUM.bed) -o 2.2|
awk 'BEGIN {OFS="\t"} {print $1}' |paste <(sort -k4,4 $TMPOUT) - |sort-bed - >  $TMPOUT.new
mv $TMPOUT.new $TMPOUT

#R --quiet --no-save << EOF
#library(dplyr)
#library(data.table)
#DHSdens <- fread("$TMPDIR/${SampleName}_DHS_MCV.bed")
#DHSdensMax <- DHSdens %>%
#       group_by(V5) %>%
#       filter(V6 == max(V6)) %>% #Get window with max DHS density 
#       filter(V4 == max(V4)) # if max DHS is the same for multiple windows, get the one with strongest contact
#write.table(DHSdensMax, file="$TMPDIR/${SampleName}_DHS_MCV_max.bed", sep='\t', col.names=F, row.names=F, quote=F)
#EOF
######
#DNase  contact maps TODO fix this code
######
echo "doing DNase peaks in contact"
header="$header\tcontactDNase"
awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5}' $OUTDIR/OBSvsExpcontactMap_raw.tsv > $OUTDIR/OBSvsExpcontactMap.tsv 
cat  $OUTDIR/OBSvsExpcontactMap.tsv|sort-bed -|
bedmap --echo --echo-map --delim '\t' --skip-unmapped --count - /vol/isg/encode/dnase/mapped/K562-DS9764/hotspots/K562-DS9764.hg38_noalt-final/K562-DS9764.hg38_noalt.fdr0.01.pks.starch |
awk '$4>2 {print $5, $NF}'|sort | awk '{arr[$1]+=$2} END {for (i in arr) {print i,arr[i]}}' |sort -k1,1 |awk 'BEGIN {OFS="\t"} {print $1, $2}' > $TMPDIR/${SampleName}_DNAcontact.tsv

join -eNA -a1 -a2 -1 4 <(sort -k4,4 $TMPOUT) -2 1 <(sort -k1,1 $TMPDIR/${SampleName}_DNAcontact.tsv) -o 1.1,1.2,1.3,1.4,2.1,2.2|
awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5, $6}'|awk '$1!="NA"'| sort-bed -|awk '{print $NF}' | sed 's/NA/0/g'| paste $TMPOUT - > $TMPOUT.new

mv $TMPOUT.new $TMPOUT


echo "doing DNase peak in 5kb windows in contact"
header="$header\tcontactDNase_window"
awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5}' $OUTDIR/OBSvsExpcontactMap_raw.tsv > $OUTDIR/OBSvsExpcontactMap.tsv 
cat  $OUTDIR/OBSvsExpcontactMap.tsv|sort-bed -|
bedmap --echo --echo-map --delim '\t' --skip-unmapped --count - /vol/isg/encode/dnase/mapped/K562-DS9764/hotspots/K562-DS9764.hg38_noalt-final/K562-DS9764.hg38_noalt.fdr0.01.pks.starch |
awk '$4>2 {print $5}'|sort | uniq -c |awk 'BEGIN {OFS="\t"} {print $2, $1}' > $TMPDIR/${SampleName}_DNAcontact_window.tsv

join -eNA -a1 -a2 -1 4 <(sort -k4,4 $TMPOUT) -2 1 <(sort -k1,1 $TMPDIR/${SampleName}_DNAcontact_window.tsv) -o 1.1,1.2,1.3,1.4,2.1,2.2|
awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5, $6}'|awk '$1!="NA"'| sort-bed -|awk '{print $NF}' | sed 's/NA/0/g'| paste $TMPOUT - > $TMPOUT.new

mv $TMPOUT.new $TMPOUT

#Hela DNase peaks
echo "doing HeLa DNase peaks in contact"
header="$header\tHeLacontactDNase"
awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5}' $OUTDIR/OBSvsExpcontactMap_raw.tsv > $OUTDIR/OBSvsExpcontactMap.tsv 
cat  $OUTDIR/OBSvsExpcontactMap.tsv|sort-bed -|
bedmap --echo --echo-map --delim '\t' --skip-unmapped --count - $TMPDIR/HeLa_DNase.bed|
awk '$4>2 {print $5, $NF}'|sort | awk '{arr[$1]+=$2} END {for (i in arr) {print i,arr[i]}}' |awk 'BEGIN {OFS="\t"} {print $1, $2}' > $TMPDIR/${SampleName}_HeLa_DNAcontact.tsv

join -eNA -a1 -a2 -1 4 <(sort -k4,4 $TMPOUT) -2 1 <(sort -k1,1 $TMPDIR/${SampleName}_HeLa_DNAcontact.tsv) -o 1.1,1.2,1.3,1.4,2.1,2.2|
awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5, $6}'|awk '$1!="NA"'| sort-bed -|awk '{print $NF}' | sed 's/NA/0/g'| paste $TMPOUT - > $TMPOUT.new

mv $TMPOUT.new $TMPOUT

echo "doing HeLa DNase peak in 5kb windows in contact"
header="$header\tHeLacontactDNase_window"
awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5}' $OUTDIR/OBSvsExpcontactMap_raw.tsv > $OUTDIR/OBSvsExpcontactMap.tsv 
cat  $OUTDIR/OBSvsExpcontactMap.tsv|sort-bed -|
bedmap --echo --echo-map --delim '\t' --skip-unmapped --count - $TMPDIR/HeLa_DNase.bed|
awk '$4>2 {print $5}'|sort | uniq -c |awk 'BEGIN {OFS="\t"} {print $2, $1}' >   $TMPDIR/${SampleName}_HeLa_DNAcontact_window.tsv

join -eNA -a1 -a2 -1 4 <(sort -k4,4 $TMPOUT) -2 1 <(sort -k1,1 $TMPDIR/${SampleName}_HeLa_DNAcontact_window.tsv) -o 1.1,1.2,1.3,1.4,2.1,2.2|
awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5, $6}'|awk '$1!="NA"'| sort-bed -|awk '{print $NF}' | sed 's/NA/0/g'| paste $TMPOUT - > $TMPOUT.new

mv $TMPOUT.new $TMPOUT

echo "doing HeLa DNase peak in 5kb windows in contact - Remove windows with K562 peaks"
header="$header\tHeLacontactDNase_window_noK562"
awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5}' $OUTDIR/OBSvsExpcontactMap_raw.tsv > $OUTDIR/OBSvsExpcontactMap.tsv 
cat  $OUTDIR/OBSvsExpcontactMap.tsv|sort-bed -|
bedmap --echo --echo-map --delim '\t' --skip-unmapped --count - $TMPDIR/HeLa_DNase.bed| 
bedmap --echo --echo-map --delim '\t' --count - /vol/isg/encode/dnase/mapped/K562-DS9764/hotspots/K562-DS9764.hg38_noalt-final/K562-DS9764.hg38_noalt.fdr0.01.pks.starch|
 awk '$4>2 && $NF==0 {print $5}'|sort | uniq -c |awk 'BEGIN {OFS="\t"} {print $2, $1}' >   $TMPDIR/${SampleName}_HeLa_DNAcontact_window_noK562.tsv

join -eNA -a1 -a2 -1 4 <(sort -k4,4 $TMPOUT) -2 1 <(sort -k1,1 $TMPDIR/${SampleName}_HeLa_DNAcontact_window_noK562.tsv) -o 1.1,1.2,1.3,1.4,2.1,2.2|
awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5, $6}'|awk '$1!="NA"'| sort-bed -|awk '{print $NF}' | sed 's/NA/0/g'| paste $TMPOUT - > $TMPOUT.new

mv $TMPOUT.new $TMPOUT


#echo "doing number of DNase peaks in contact with insertion"
#for BC in `cat $TMPOUT|awk '{print $4}'`; do
#       NumDNase=$(cat $TMPOUT|awk -v Barcode="$BC" '$4==Barcode' |bedops -u --range 300000 -| #Get a range of +-300kb from the insertion
#       bedops --chop 5000 -|paste - $OUTDIR/contactMaps/${BC}_contactMap.tsv 2> /dev/null|  #Chop it into 5kb windows i.e. the resolution used to calculate contacts
#        #Only use 5kb windows with twice the expected contact
#       bedmap --echo --echo-map - /vol/isg/encode/dnase/mapped/K562-DS9764/hotspots/K562-DS9764.hg38_noalt-final/K562-DS9764.hg38_noalt.fdr0.01.pks.starch | awk -F'|' '$2!=""'|awk -F "\t" 'BEGIN {OFS="\t"} $5>2'| wc -l);
#       echo $BC $NumDNase  #Overlap DNase peaks
#done > $OUTDIR/DnaseContact2.txt


echo "doing Contact maps for TFs/Histone peaks in contact"
#Marks picked from Mercer et al. 2013 Nature Gen. DNase I-hypersensitive exons colocalize with promoters and distal regulatory elements.
for chip in `cat /vol/isg/encode/chipseq/SamplesForTrackhub.tsv |grep K562| grep rep1|grep -w 'E2F6\|H3K27me3\|H3K4me3\|H3K9ac\|CCNT2\|HMGN3\|E2F4\|POLR2A\|TBP\|MYC\|MXI1\|BHLH40\|EP300\|CEPBP\|TAL1\|GATA1\|GATA2\|MAFK\|H3K4me2\|H3K4me1\|SMC3\|RAD21\|CTCF\|ZNF143'|awk '{print $1"-"$5"-"$2}'`; do 
       mark=$(echo $chip|awk -F '-' '{print $2}')
       header="$header\tcontact"${mark}
       cat  $OUTDIR/OBSvsExpcontactMap.tsv|sort-bed -|
       bedmap --echo --echo-map --delim '\t' --skip-unmapped --count -  /vol/isg/encode/chipseq/mapped/${chip}/hotspots/${chip}.hg38_noalt-final/${chip}.hg38_noalt.fdr0.01.pks.starch |
       awk '$4>2 {print $5, $NF}'|sort  | awk '{arr[$1]+=$2} END {for (i in arr) {print i,arr[i]}}'|awk 'BEGIN {OFS="\t"} {print $1, $2}'> $TMPDIR/${SampleName}_${mark}contact.tsv
       
       join -eNA -a1 -a2 -1 4 <(sort -k4,4 $TMPOUT) -2 1 <(sort -k1,1  $TMPDIR/${SampleName}_${mark}contact.tsv) -o 1.1,1.2,1.3,1.4,2.1,2.2|awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5, $6}'|awk '$1!="NA"'| sort-bed -|awk '{print $NF}'  |sed 's/NA/0/g'| paste $TMPOUT - > $TMPOUT.new
       mv $TMPOUT.new $TMPOUT
done





######
##Looking at 300kb surrounding the insertion  
######
#echo "doing number of DNase peaks in +- 300kb of insertion "
#for BC in `cat $TMPOUT|awk '{print $4}'`; do
#       NumDNase=$(cat $TMPOUT|awk -v Barcode="$BC" '$4==Barcode' |bedops -u --range 300000 -| #Get a range of +-300kb from the insertion
#       bedops --chop 5000 -|paste - $OUTDIR/contactMaps/${BC}_contactMap.tsv 2> /dev/null|  #Chop it into 5kb windows i.e. the resolution used to calculate contacts
#        #Only use 5kb windows with twice the expected contact
#       bedmap --echo --echo-map - /vol/isg/encode/dnase/mapped/K562-DS9764/hotspots/K562-DS9764.hg38_noalt-final/K562-DS9764.hg38_noalt.fdr0.01.pks.starch | awk -F'|' '$2!=""'| wc -l);
#       echo $BC $NumDNase  #Overlap DNase peaks
#done > $OUTDIR/DnaseContact_all.txt
  


######
#4C data - Using normalised matrix from Rao et al. which I realigned with HiC-Pro and normalised with HiTC
#######
#mkdir -p $OUTDIR/Contacts/
#echo "doing number of CTCF peaks in contact with insertion"
#for BC in `cat $TMPOUT|awk '{print $4}'`; do
#       NumDNase=$(cat $TMPOUT|awk -v Barcode="$BC" '$4==Barcode' |bedops -u --range 300000 -| #Get a range of +-300kb from the insertion
#       bedops --chop 5000 -|paste - $OUTDIR/contactMaps/${BC}_contactMap.tsv 2> /dev/null|  #Chop it into 5kb windows i.e. the resolution used to calculate contacts
#        #Only use 5kb windows with twice the expected contact
#       bedmap --echo --echo-map - /vol/isg/encode/chipseq/mapped/K562-CTCF-DS11247/hotspots/K562-CTCF-DS11247.hg38_noalt-final/K562-CTCF-DS11247.hg38_noalt.fdr0.01.pks.starch  | awk -F'|' '$2!=""'|awk -F "\t" 'BEGIN {OFS="\t"} $5>2'| wc -l);
#       echo $BC $NumDNase  #Overlap DNase peaks
#done > $OUTDIR/Contacts/CTCF_Contact.txt
#  
#
#
#
######
##Looking at 300kb surrounding the insertion  
######
#echo "doing number of CTCF peaks in +- 300kb of insertion "
#for BC in `cat $TMPOUT|awk '{print $4}'`; do
#       NumDNase=$(cat $TMPOUT|awk -v Barcode="$BC" '$4==Barcode' |bedops -u --range 300000 -| #Get a range of +-300kb from the insertion
#       bedops --chop 5000 -|paste - $OUTDIR/contactMaps/${BC}_contactMap.tsv 2> /dev/null|  #Chop it into 5kb windows i.e. the resolution used to calculate contacts
#        #Only use 5kb windows with twice the expected contact
#       bedmap --echo --echo-map - /vol/isg/encode/chipseq/mapped/K562-CTCF-DS11247/hotspots/K562-CTCF-DS11247.hg38_noalt-final/K562-CTCF-DS11247.hg38_noalt.fdr0.01.pks.starch | awk -F'|' '$2!=""'| wc -l);
#       echo $BC $NumDNase  #Overlap DNase peaks
#done > $OUTDIR/Contacts/CTCF_Contact_all.txt
  
  


#######
#TF analysis
#######
#TODO fix working from normalised insertion data
echo "doing TF biding analysis"
#cat $OUTDIR/AllInsertions.coords.bed > $TMPDIR/${SampleName}_TF_analysisTest.tsv

cat /vol/isg/encode/chipseq/SamplesForTrackhub.tsv |grep K562| grep rep1|grep 'Cell_lines_TFs\|Cell_lines_Histones'|awk '{print $5}'| grep -v 'H3K\|H4K' > $OUTDIR/TFmarks
cat /vol/isg/encode/chipseq/SamplesForTrackhub.tsv |grep K562| grep rep1|grep 'Cell_lines_TFs\|Cell_lines_Histones'|awk '{print $5}'| grep  'H3K\|H4K' > $OUTDIR/HistoneMarks

#TFOUT=$TMPDIR/${SampleName}_TF_analysisTest.tsv


#TFheader="chrom\tchromStart\tchromEnd\tBC\texpression\tstrand"
for chip in `cat /vol/isg/encode/chipseq/SamplesForTrackhub.tsv |grep K562| grep rep1|grep 'Cell_lines_TFs\|Cell_lines_Histones'|awk '{print $1"-"$5"-"$2".hg38"}'`; do 
       mark=$(echo $chip|awk -F '-' '{print $2}')
       header="$header\t"${mark}
       #Closest upstream
       #closest-features --delim '\t'  --dist   $TMPOUT /vol/isg/encode/chipseq/mapped/${chip}/hotspots/${chip}-final/${chip}.fdr0.01.pks.bed |  awk 'BEGIN {OFS="\t"} {if ($6=="-") print $(NF); else if  ($6=="+") print $(NF-4)}'|paste $TMPOUT - > $TMPOUT.new
       #Closest 
       closest-features --delim '\t'  --dist --closest  $TMPOUT /vol/isg/encode/chipseq/mapped/${chip}/hotspots/${chip}-final/${chip}.fdr0.01.pks.bed |  awk -F '\t' 'BEGIN {OFS="\t"} {print $NF}'|paste $TMPOUT - > $TMPOUT.new
       mv $TMPOUT.new $TMPOUT
done
#

if [ -f "$OUTDIR/AllInsertions.coords.bed" ]; then
       echo "doing ChIP enrichment over expressed and nonexpressed insertions"
       numInsertions=$(cat $OUTDIR/AllInsertions.coords.bed| wc -l)
       Insquartile=$(echo "$(($numInsertions/4))")
       sort -rnk5 $OUTDIR/AllInsertions.coords.bed| head -$Insquartile| sort-bed - > $TMPDIR/${SampleName}_highExpressed.bed
       sort -rnk5 $OUTDIR/AllInsertions.coords.bed| tail -$Insquartile| sort-bed - > $TMPDIR/${SampleName}_lowExpressed.bed
       
       for j in {2000,5000,20000,50000,100000,200000}; do 
              for chip in `cat /vol/isg/encode/chipseq/SamplesForTrackhub.tsv |grep K562| grep rep1|grep 'Cell_lines_TFs'|awk '{print $1"-"$5"-"$2".hg38"}'`; do 
                     bash /home/maagj01/scratch/transposon/Analysis/TAD_CTCFenrichment/Enrichment.sh $TMPDIR/${SampleName}_highExpressed.bed\
                     /vol/isg/encode/chipseq/mapped/${chip}/hotspots/${chip}-final/${chip}.fdr0.01.hot.bed ${j}; 
              done; 
       done| awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5, "High"}' > $OUTDIR/${SampleName}_ChIPhighExpressed.tsv

       for j in {2000,5000,20000,50000,100000,200000}; do 
              for chip in `cat /vol/isg/encode/chipseq/SamplesForTrackhub.tsv |grep K562| grep rep1|grep 'Cell_lines_TFs'|awk '{print $1"-"$5"-"$2".hg38"}'`; do 
                     bash /home/maagj01/scratch/transposon/Analysis/TAD_CTCFenrichment/Enrichment.sh $TMPDIR/${SampleName}_lowExpressed.bed\
                     /vol/isg/encode/chipseq/mapped/${chip}/hotspots/${chip}-final/${chip}.fdr0.01.hot.bed ${j}; 
              done; 
       done| awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5, "Low"}' > $OUTDIR/${SampleName}_ChIPlowExpressed.tsv
       
       
       sort -rnk5 $OUTDIR/AllInsertions.coords.bed| head -$Insquartile| 
              sort-bed - |bedops -u --range 10000 -|
              bedmap --delim '\t' --echo --echo-map --count - /vol/isg/annotation/bed/hg38/refseq_gene/refGene.CombinedTxStarts.bed |
              awk 'BEGIN {OFS="\t"} $NF>=1 {print $1, $2, $3, $4}'|bedops -u --range -10000 - > $TMPDIR/${SampleName}_highExpressed_Promoter.bed
              
              for j in {2000,5000,20000,50000,100000,200000}; do 
                     for chip in `cat /vol/isg/encode/chipseq/SamplesForTrackhub.tsv |grep K562| grep rep1|grep 'Cell_lines_TFs'|awk '{print $1"-"$5"-"$2".hg38"}'`; do 
                            bash /home/maagj01/scratch/transposon/Analysis/TAD_CTCFenrichment/Enrichment.sh $TMPDIR/${SampleName}_highExpressed_Promoter.bed\
                            /vol/isg/encode/chipseq/mapped/${chip}/hotspots/${chip}-final/${chip}.fdr0.01.hot.bed ${j}; 
                     done; 
              done| awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5, "High", "Promoter"}' > $OUTDIR/${SampleName}_ChIPhighExpressed_Promoter.tsv
              
              
              
       sort -rnk5 $OUTDIR/AllInsertions.coords.bed| head -$Insquartile| 
              sort-bed - |bedops -u --range 10000 -|
              bedmap --delim '\t' --echo --echo-map --count - /vol/isg/annotation/bed/hg38/refseq_gene/refGene.CombinedTxStarts.bed |
              awk 'BEGIN {OFS="\t"} $NF==0 {print $1, $2, $3, $4}'|bedops -u --range -10000 - > $TMPDIR/${SampleName}_highExpressed_Distal.bed
              
              for j in {2000,5000,20000,50000,100000,200000}; do 
                     for chip in `cat /vol/isg/encode/chipseq/SamplesForTrackhub.tsv |grep K562| grep rep1|grep 'Cell_lines_TFs'|awk '{print $1"-"$5"-"$2".hg38"}'`; do 
                            bash /home/maagj01/scratch/transposon/Analysis/TAD_CTCFenrichment/Enrichment.sh $TMPDIR/${SampleName}_highExpressed_Distal.bed\
                            /vol/isg/encode/chipseq/mapped/${chip}/hotspots/${chip}-final/${chip}.fdr0.01.hot.bed ${j}; 
                     done; 
              done| awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5, "High", "Distal"}' > $OUTDIR/${SampleName}_ChIPhighExpressed_Distal.tsv
fi




############
#Bedfile with all HiC contacts +- 10Mb
############
#First find all insertions in contact with something 
#TODO multiple insertions overlapping same 5kb window...
#First step is time consuming for all chromosomes
##bedmap --delim '\t' --echo --echo-map-id  --skip-unmapped  /home/maagj01/scratch/transposon/Analysis/CaptureHiC/K562/Results_K562_HiC069-074/hic_results/matrix/K562/K562_5kb_intraChromContacts.starch  $TMPOUT > $TMPDIR/${SampleName}_allcontactMappedID.bed
#echo 'Creating +- 1Mb contact map and overlapping with DHS peaks in contact'
#echo 'No filtering here just report sum of DHS and CTCF peaks'
#bedtools intersect -wa -wb -a /home/maagj01/scratch/transposon/Analysis/CaptureHiC/K562/Results_K562_HiC069-074/hic_results/matrix/K562/K562_5kb_intraChromContacts.bed  -b $TMPOUT |
#awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5, $6, $7, $11}'> $TMPDIR/${SampleName}_allcontactMappedIDtools.bed
#
##Create 5kb windows around each insertions +- 1Mb away with the insertion BC
#awk 'BEGIN {OFS="\t"} {print $1, 0, $2}' /vol/isg/annotation/fasta/hg38/hg38.chrom.sizes|grep -v 'Un\|alt\|random'| sort-bed - |
#bedops -w 5000 -| 
#bedtools intersect -wa -wb -a - -b <(bedops -u --range 1000000 $TMPOUT)| 
#awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $7}' > $TMPDIR/${SampleName}_5kbwindows.bed
# 
# 
#bedtools intersect -wao -b <(awk 'BEGIN {OFS="\t"} {print $4, $5, $6, $7, $8}' $TMPDIR/${SampleName}_allcontactMappedIDtools.bed|sort-bed -) -a $TMPDIR/${SampleName}_5kbwindows.bed|
#awk 'BEGIN {OFS="\t"} {if ($4==$(NF-1)) print $1, $2, $3, $4, $(NF-2); else if ($4!=$(NF-1)) print $1, $2, $3, $4, 0}' | #####Adds extra rows with zero whenever a window overlaps another insertion.
#uniq > $TMPDIR/${SampleName}_1MbContacts.bed
#
#R --quiet --no-save << EOF
#library(dplyr)
#library(data.table)
#library(parallel)
#SampleContacts <-  fread("$TMPDIR/${SampleName}_1MbContacts.bed", sep='\t')
#samplesClean <- SampleContacts %>% 
#       group_by(V1, V2, V3, V4) %>%
#       filter(V5 == max(V5))
#
######
##Loess normalise all connections with full length
######
##Will loose 173 insertions close to the edge of the chromosomes
#options(mc.cores = 16)
#samplesClean <- as.data.frame(samplesClean)
#contactList <- split(samplesClean, samplesClean[,"V4"])
#
#output4 <- mclapply(contactList, function(loci) {
#       if ( nrow(loci) ==401) {
#              loci[,"distance"] <-  seq(from=-1e6, to=1e6, by=5000)
#              loci[,"percentageOfMax"] <- loci[,"V5"]/max(loci[,"V5"])*100
#       } else {loci[,"distance"] <- NA
#              loci[,"percentageOfMax"] <- NA}
#       loci
#})
#
#contactListDistance <- as.data.frame(rbindlist(output4))
#
#contactListDistance2 <- contactListDistance[!is.na(contactListDistance[,"distance"]),]
#delta <- 0.005*diff(range(contactListDistance2[,"distance"]))
#
#lowess.fit <-lowess(x=contactListDistance2[,"distance"], y=contactListDistance2[,"V5"], f=0.01, delta=delta)[["y"]]
#
#
#contactListDistance2[,"lowess"] <- rep(unique(lowess.fit), length(unique(contactListDistance2[,"V4"])))
#contactListDistance2[,"ObsExp"] <- contactListDistance2[,"V5"]/contactListDistance2[,"lowess"]
#
#lowess.fit.Per <-lowess(x=contactListDistance2[,"distance"], y=contactListDistance2[,"percentageOfMax"], f=0.01, delta=delta)[["y"]]
#contactListDistance2[,"lowessPercentage"] <- rep(unique(lowess.fit.Per),length(unique(contactListDistance2[,"V4"])))
#contactListDistance2[,"ObsExpPercentage"] <- contactListDistance2[,"percentageOfMax"]/contactListDistance2[,"lowessPercentage"]
#
#
#write.table(contactListDistance2, file="$OUTDIR/${SampleName}_1Mb_contact_loess.bed", sep='\t', col.names=F, row.names=F, quote=F)
#
##Colnames "chr"       "start"       "end"       "BC"       "ICE-normalised contact"       distance       percentageOfMax    lowess   ObsExp lowessPercentage     ObsExpPercentage
#EOF
#
#
#
#oneMBheader="chrom\tstart\tend\tBC\tnormContact\tdistance\tpercentageOfMax\tlowess\tObseExp\tlowessPercentage\tObsExpPercentage"
#oneMBheader="$oneMBheader\tDHSdensity"
#sort-bed $OUTDIR/${SampleName}_1Mb_contact_loess.bed| bedmap  --delim '\t' --echo --echo-map --sum - <( awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5}' $TMPDIR/K562_DHS_dens_p.bed)|
#awk '{print $NF}' |sed 's/NAN/0/g'| paste <(sort-bed $OUTDIR/${SampleName}_1Mb_contact_loess.bed) -  > $TMPDIR/${SampleName}_1Mb_contact_loess.bed.new
#
#
#mv $TMPDIR/${SampleName}_1Mb_contact_loess.bed.new $OUTDIR/${SampleName}_1Mb_contact_loess_annotated.bed
#
#
#oneMBheader="$oneMBheader\tCTCFdensity"
#sort-bed $OUTDIR/${SampleName}_1Mb_contact_loess_annotated.bed| bedmap  --delim '\t' --echo --echo-map --sum - <( awk 'BEGIN {OFS="\t"} {print $1, $2, $3, NR, $4}' $TMPDIR/${SampleName}.CTCF.denspeaks.bed)|
#awk '{print $NF}' |sed 's/NAN/0/g'| paste <(sort-bed $OUTDIR/${SampleName}_1Mb_contact_loess_annotated.bed) -  > $TMPDIR/${SampleName}_1Mb_contact_loess.bed.new
#
#
#echo -e $oneMBheader | cat - $TMPDIR/${SampleName}_1Mb_contact_loess.bed.new > $OUTDIR/${SampleName}_1Mb_contact_loess_annotated.bed
#

######
#Big bed with all contacts
#######
echo 'creating bigBed with all contacts'
bedtools intersect -sorted -wa -wb -a <(awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4}' $TMPOUT| sort-bed -) -b /home/maagj01/scratch/transposon/Analysis/CaptureHiC/K562/Results_K562_HiC069-074/hic_results/matrix/K562/K562_5kb_intraChromContacts_Loess_2Mb.bed |
awk 'BEGIN {OFS="\t"} $11<=1000000 && $11>=-1000000 {print $5, $8, $9, $4, $10, $11, $12, $13, $14, $15, $16}' >$OUTDIR/${SampleName}_1Mb_contact_loess.bed


oneMBheader="chrom\tstart\tend\tBC\tnormContact\tdistance\tpercentageOfMax\tlowess\tObseExp\tlowessPercentage\tObsExpPercentage"
oneMBheader="$oneMBheader\tDHSdensity"
sort-bed $OUTDIR/${SampleName}_1Mb_contact_loess.bed| bedmap  --delim '\t' --echo --echo-map --sum - <( awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5}' $TMPDIR/K562_DHS_dens_p.bed)|
awk '{print $NF}' |sed 's/NAN/0/g'| paste <(sort-bed $OUTDIR/${SampleName}_1Mb_contact_loess.bed) -  > $TMPDIR/${SampleName}_1Mb_contact_loess.bed.new


mv $TMPDIR/${SampleName}_1Mb_contact_loess.bed.new $OUTDIR/${SampleName}_1Mb_contact_loess_annotated.bed


oneMBheader="$oneMBheader\tCTCFdensity"
sort-bed $OUTDIR/${SampleName}_1Mb_contact_loess_annotated.bed| bedmap  --delim '\t' --echo --echo-map --sum - <( awk 'BEGIN {OFS="\t"} {print $1, $2, $3, NR, $4}' $TMPDIR/${SampleName}.CTCF.denspeaks.bed)|
awk '{print $NF}' |sed 's/NAN/0/g'| paste <(sort-bed $OUTDIR/${SampleName}_1Mb_contact_loess_annotated.bed) -  > $TMPDIR/${SampleName}_1Mb_contact_loess.bed.new


echo -e $oneMBheader | cat - $TMPDIR/${SampleName}_1Mb_contact_loess.bed.new > $OUTDIR/${SampleName}_1Mb_contact_loess_annotated.bed

#bedtools intersect -wa -wb -a  $TMPOUT -b <(awk '$9>5' $OUTDIR/${SampleName}_1Mb_contact_loess_annotated.bed|tail -n +2| sort-bed - )|  awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $(NF-9), $(NF-8)}' > $TMPDIR/${SampleName}_InsContacts.tsv


####
#TODO Add GTex tissue specificity for different expression levels  
####
#echo "doing expression of GTEx genes - all expressed"
#header="$header\tGTExAll"
#
#echo "doing neighbor expression analysis --- All expressed GTEx "
#awk '$5=="-"'  /vol/mauranolab/maagj01/ExAc/GTEx/GTExAllexpressed.sorted.bed |sort-bed - >$TMPDIR/${SampleName}.GTExAllexpressed.NegStrand.bed
#awk '$5=="+"'  /vol/mauranolab/maagj01/ExAc/GTEx/GTExAllexpressed.sorted.bed |sort-bed - >$TMPDIR/${SampleName}.GTExAllexpressed.PosStrand.bed
#closest-features --closest --dist --no-ref $TMPOUT  /vol/isg/annotation/bed/hg38/gencodev25/Gencodev25.gene.bed|awk -F'|' 'BEGIN {OFS="\t"} {print $1}'|paste - $TMPOUT|
#awk '$6=="+"'| sort-bed -|  bedmap  --delim '\t' --fraction-ref 1.0 --echo --count - $TMPDIR/${SampleName}.GTExAllexpressed.PosStrand.bed >$TMPDIR/${SampleName}.PosStrand 
# 
#closest-features --closest --dist --no-ref $TMPOUT /vol/isg/annotation/bed/hg38/gencodev25/Gencodev25.gene.bed|awk -F'|' 'BEGIN {OFS="\t"} {print $1}'|paste - $TMPOUT|
#awk '$6=="-"'| sort-bed -|  bedmap  --delim '\t' --fraction-ref 1.0 --echo --count - $TMPDIR/${SampleName}.GTExAllexpressed.NegStrand.bed >$TMPDIR/${SampleName}.NegStrand 
#
#header="$header\tnearGencodeGTEx"
#cat $TMPDIR/${SampleName}.PosStrand  $TMPDIR/${SampleName}.NegStrand |cut -f7- > $TMPOUT.new
#mv $TMPOUT.new $TMPOUT




#Print to file
echo -e $header | cat - <(sort-bed $TMPOUT)  > $OUTDIR/AllInsertions.annotated.bed

echo -e $header|tr ' ' '\n'
mkdir -p $OUTDIR/graphs
date


echo 'Plotting data in R'
/vol/mauranolab/mapped/src/transposon/AnalyseOverlappedBCs_plotting.R --SampleName ${SampleName} --OUTDIR ${OUTDIR}
echo 'Finished plotting in R'


echo 'Converting pdf to png'
for i in `ls $OUTDIR/graphs/*pdf|sed 's/.pdf//g'`; do  convert -density 75 $i.pdf -quality 50 $i.png ;done


echo 'sodaPlots of geneDeserts with expressed genes'
bedmap --echo --echo-map ${OUTDIR}/Sodaplot/geneDeserts.bed /home/maagj01/public_html/blog/2017Jul10/genomeWalking/geneDeserts_10percente.bed |awk -F '|' '{print $2}'|sort|uniq >${OUTDIR}/Sodaplot/geneDeserts_whole_chunk.bed

#module load python/2.7.10 
#if [ -d "${OUTDIR}/Sodaplot/geneDeserts" ]; then
#  echo 'geneDesert folder exists. Removing and recreating'
#  rm -r ${OUTDIR}/Sodaplot/geneDeserts
#fi
#~/src/soda/soda.py -v -r ${OUTDIR}/Sodaplot/geneDeserts_whole_chunk.bed -g http://genome.isg.med.nyu.edu/ -j 300 -b hg38 -s 873_nhPQZL28pa3NVqsIRS0C5opzrGL9 -o ${OUTDIR}/Sodaplot/geneDeserts
#
#echo 'Creating UCSC track of insertions'
#awk -v sample=$sample -F "\t" 'BEGIN {OFS="\t"; print "track name=" sample " description=" sample "-integrations"} {print}' $OUTDIR/AllInsertions.annotated.bed > $OUTDIR/$sample.barcodes.coords.ucsc.bed
#####
#Analysis of regions in between insertion sites
#####

echo "doing analysis of regions between insertions points"
RegionHeader="chrom\tchromStart\tchromEnd\tleftSignal\trightSignal"
REGOUT=$TMPDIR/${SampleName}.regions.bed
tail -n +2 $OUTDIR/AllRegions.tsv|sort-bed - > $REGOUT


echo "length between insertions sites"
RegionHeader="$RegionHeader\tregionWidth"
awk "{print \$3-\$2}"  $REGOUT|paste $REGOUT - >$REGOUT.new
mv $REGOUT.new $REGOUT


echo 'doing number of CTCF peaks in region'
RegionHeader="$RegionHeader\tnCTCFpeaks"
bedmap --delim '\t' --echo --echo-map --count  $REGOUT /home/maagj01/scratch/transposon/Analysis/K562_DNaseMasterList/CTCF_strand.bed |awk -F "\t" 'BEGIN {OFS="\t"} {print $NF}'|paste $REGOUT - > $REGOUT.new
mv $REGOUT.new $REGOUT


echo 'doing number of Genes in region'
RegionHeader="$RegionHeader\tnGencodeGenes"
bedmap --delim '\t' --echo --echo-map --count  $REGOUT /vol/isg/annotation/bed/hg38/gencodev24/Gencodev24.gene.bed |awk -F "\t" 'BEGIN {OFS="\t"} {print $NF}'|paste $REGOUT - > $REGOUT.new
mv $REGOUT.new $REGOUT


echo 'doing number of DHS in region'
RegionHeader="$RegionHeader\tnDHS"
bedmap --delim '\t' --echo --echo-map --count  $REGOUT /vol/isg/encode/dnase/mapped/K562-DS9764/hotspots/K562-DS9764.hg38_noalt-final/K562-DS9764.hg38_noalt.fdr0.01.pks.starch |awk -F "\t" 'BEGIN {OFS="\t"} {print $NF}'|paste $REGOUT - > $REGOUT.new
mv $REGOUT.new $REGOUT


echo 'doing number of CGI in region'
RegionHeader="$RegionHeader\tnCGI"
bedmap --delim '\t' --echo --echo-map --count  $REGOUT  /vol/isg/annotation/bed/hg38/cpg_islands/cpgIslands.bed  |awk -F "\t" 'BEGIN {OFS="\t"} {print $NF}'|paste $REGOUT - > $REGOUT.new
mv $REGOUT.new $REGOUT


echo 'doing number of Dpn sites in region'
RegionHeader="$RegionHeader\tnDpn"
bedmap --delim '\t' --echo --echo-map --count  $REGOUT  /home/maagj01/scratch/transposon/Analysis/Dpn_REsites/DpnPlusMinus.Sorted.bed |awk -F "\t" 'BEGIN {OFS="\t"} {print $NF}'|paste $REGOUT - > $REGOUT.new
mv $REGOUT.new $REGOUT


echo 'doing number of k562 expressed genes (TPM>=1) in Region'
RegionHeader="$RegionHeader\tnExpressedGenes"
awk '$NF>=1' /home/maagj01/scratch/transposon/Analysis/K562_GeneQuant/K562_GeneQuant.bed >$TMPDIR/${SampleName}.k562expressedgenes.bed

bedmap --delim '\t' --echo --echo-map --count  $REGOUT $TMPDIR/${SampleName}.k562expressedgenes.bed| awk -F'|' 'BEGIN {OFS="\t"} {print $1}' |awk  '{print $NF}' |paste $REGOUT - > $REGOUT.new
mv $REGOUT.new $REGOUT


echo -e $RegionHeader | cat - <(sort-bed $REGOUT) > $OUTDIR/AllRegions.annotated.bed

date
mkdir -p $OUTDIR/Regiongraphs/
echo 'Plotting data in R'
/vol/mauranolab/mapped/src/transposon/AnalyseOverlappedBCs_Region_plotting.R --SampleName ${SampleName} --OUTDIR ${OUTDIR}
echo 'Finished plotting in R'


echo 'Converting pdf to png'
for i in `ls $OUTDIR/Regiongraphs/*pdf|sed 's/.pdf//g'`; do  convert -density 75 $i.pdf -quality 50 $i.png ;done

echo -e $header|tr ' ' '\n'
echo 'Done!!!'
date


