#!/bin/bash
set -e -o pipefail

Dir=$1
Flow=$(basename `pwd`|sed 's/aligned.//g')
echo $Flow
OUTDIR=$Dir/$Flow
mkdir -p $OUTDIR


#####
#Weblogos
######
mkdir -p $OUTDIR/Weblogos
#For raw sequence
mkdir -p $OUTDIR/Weblogos/raw
find -name *.R2.raw.png |grep -v bak| xargs cp -t $OUTDIR/Weblogos/raw
find -name *.R1.raw.png |grep -v bak| xargs cp -t $OUTDIR/Weblogos/raw
for i in `ls $OUTDIR/Weblogos/raw/*R1.raw.png|sort|awk -F'/' '{print $NF}'`;do echo '<a href="'${i}'"><img src="'${i}'" height="120"></a>' ;done  >$OUTDIR/Weblogos/raw/R1index.html
for i in `ls $OUTDIR/Weblogos/raw/*R2.raw.png|sort|awk -F'/' '{print $NF}'`;do echo '<a href="'${i}'"><img src="'${i}'"  style= "position:absolute; LEFT:900px"; height="120"></a>' ;done  >$OUTDIR/Weblogos/raw/R2index.html
cat $OUTDIR/Weblogos/raw/R1index.html $OUTDIR/Weblogos/raw/R2index.html |sort| awk ' {print;} NR % 2 == 0 { print "<br>"; }'> $OUTDIR/Weblogos/raw/index.html



#For processed sequence
mkdir -p $OUTDIR/Weblogos/processed
find -name *.BC.processed.png |grep -v bak| xargs cp -t $OUTDIR/Weblogos/processed
find -name *.plasmid.processed.png |grep -v bak| xargs cp -t $OUTDIR/Weblogos/processed
for i in `ls $OUTDIR/Weblogos/processed/*BC.processed.png|sort|awk -F'/' '{print $NF}'`;do echo '<a href="'${i}'"><img src="'${i}'" height="120"></a>' ;done  >$OUTDIR/Weblogos/processed/BCindex.html
for i in `ls $OUTDIR/Weblogos/processed/*plasmid.processed.png|sort|awk -F'/' '{print $NF}'`;do echo '<a href="'${i}'"><img src="'${i}'"  style= "position:absolute; LEFT:900px"; height="120"></a>' ;done  >$OUTDIR/Weblogos/processed/plasmidindex.html
cat $OUTDIR/Weblogos/processed/BCindex.html $OUTDIR/Weblogos/processed/plasmidindex.html |sort| awk ' {print;} NR % 2 == 0 { print "<br>"; }'> $OUTDIR/Weblogos/processed/index.html


#For Barcode sequence
mkdir -p $OUTDIR/Weblogos/Barcode
find -name *barcodes.png |grep -v bak| xargs cp -t $OUTDIR/Weblogos/Barcode
for i in `ls $OUTDIR/Weblogos/Barcode/*barcodes.png|sort|awk -F'/' '{print $NF}'`;do echo '<a href="'${i}'"><img src="'${i}'" height="120"></a>' ;done  |awk ' {print;} NR % 1 == 0 { print "<br>"; }'> $OUTDIR/Weblogos/Barcode/index.html

#For UMI sequence
mkdir -p $OUTDIR/Weblogos/UMI
find -name *UMIs.png |grep -v bak| xargs cp -t $OUTDIR/Weblogos/UMI
for i in `ls $OUTDIR/Weblogos/UMI/*UMIs.png|sort|awk -F'/' '{print $NF}'`;do echo '<a href="'${i}'"><img src="'${i}'" height="120"></a>' ;done  |awk ' {print;} NR % 1 == 0 { print "<br>"; }'> $OUTDIR/Weblogos/UMI/index.html




######
#Levenstein distance per sample
#######
mkdir -p $OUTDIR/EditDist

#
for i in `find -type d| grep BS|grep -v 'bak'`;do  grep 'BC levensthein distance' ${i}/extract*|awk -v OFS='\t' -F'[' '{print $1, $3}'|sed 's/]//g';done| perl -pe 's/, /\t/g' > $OUTDIR/EditDist/BC_EditDist
for i in `find -type d| grep BS|grep -v 'bak'`;do  grep 'Plasmid levensthein distance' ${i}/extract*|awk -v OFS='\t' -F'[' '{print $1, $3}'|sed 's/]//g';done| perl -pe 's/, /\t/g' > $OUTDIR/EditDist/Plasmid_EditDist

R  --quiet --no-save << EOF

library(stringr)
library(reshape2)
BCleven<-read("$OUTDIR/EditDist/BC_EditDist",stringsAsFactors=F)

BCleven<-summaryBy(.~V1, data=BCleven,FUN=sum)
colnames(BCleven)<-c('Sample',0:18)
BCleven<-melt(BCleven)
BCleven[,4]<-"Barcode"

Plasmidleven<-read("$OUTDIR/EditDist/Plasmid_EditDist",stringsAsFactors=F)
Plasmidleven<-summaryBy(.~V1, data=Plasmidleven,FUN=sum)
colnames(Plasmidleven)<-c('Sample',0:(ncol(Plasmidleven)-2))
Plasmidleven<-melt(Plasmidleven)
Plasmidleven[,4]<-'Plasmid'

BCleven<-rbind(BCleven,Plasmidleven)
colnames(BCleven)[1:4]<-c('Sample','levenDist','Reads','Read')
BCleven[,2]<-as.numeric(as.character(BCleven[,2]))
BCleven[,3]<-as.numeric(as.character(BCleven[,3]))

colnames(BCleven)[1:4]<-c('Sample','levenDist','Reads','Read')
#plot
BCleven[,1]<-gsub('\\\.o.*','',gsub('extract.','',(gsub('.*/','',BCleven[,1]))))
pdf("$OUTDIR/EditDist/LevenDistance_BC.pdf",height=10,width=20)
ggplot(BCleven[grep('Barcode',BCleven[,4]),], aes(levenDist,Reads)) + 
geom_bar(stat = "identity") +
#ggtitle('GGlo')+
facet_wrap(~Sample+Read,scales = "free_y")+ 
theme_classic()+
#scale_fill_brewer('Set1')+
theme(legend.position = "bottom",legend.text=element_text(size=8))+
#scale_x_continuous(label= fancy_scientific)+
scale_y_continuous(label= fancy_scientific)
dev.off()


#plot
pdf("$OUTDIR/EditDist/LevenDistance_Plasmid.pdf",height=10,width=20)
ggplot(BCleven[grep('Plasmid',BCleven[,4]),], aes(levenDist,Reads)) + 
geom_bar(stat = "identity") +
#ggtitle('GGlo')+
facet_wrap(~Sample+Read,scales = "free_y")+ 
theme_classic()+
#scale_fill_brewer('Set1')+
theme(legend.position = "bottom",legend.text=element_text(size=8))+
#scale_x_continuous(label= fancy_scientific)+
scale_y_continuous(label= fancy_scientific)
dev.off()

EOF

convert -density 300 $OUTDIR/EditDist/LevenDistance_BC.pdf -quality 100 $OUTDIR/EditDist/LevenDistance_BC.png
convert -density 300 $OUTDIR/EditDist/LevenDistance_Plasmid.pdf -quality 100 $OUTDIR/EditDist/LevenDistance_Plasmid.png

for i in `ls $OUTDIR/EditDist/*.png|sort|sed 's/.png//g' |awk -F'/' '{print $NF}'`;do echo '<a href="'${i}.pdf'"><img src="'${i}.png'" height="700"></a>' ;done  >$OUTDIR/EditDist/index.html



#####
#Barcode frequency
#####
mkdir -p $OUTDIR/BarcodeFreq
for i in `find -type d|grep BS00|grep -v bak|sed 's/^..//g'`; do sort -nk2 ${i}/${i}.barcode.counts.txt| awk -v OFS='\t' -F'\t' -v sample="$i" '{print $1, $2, sample}';done > $OUTDIR/BarcodeFreq/BarcodeFreq.txt

R  --quiet --no-save << EOF
BC <- read("$OUTDIR/BarcodeFreq/BarcodeFreq.txt")
colnames(BC) <- c("Barcodes","barcodeFreq", "Sample")
BC\$BC.bin <- cut(BC\$barcodeFreq, breaks=c( 0, 1,5, 10, 50, 100, 1000, 10000, 100000, 1000000), right=F, include.lowest=T, labels=c("+1", "+5", "+10", "+50", "+100", "+1000", "+10000", "+100000","+1000000"))


t<-ggplot(BC, aes(x=BC.bin)) +
geom_bar(color="black", size=0.25) +
theme_classic()+
facet_wrap(~Sample,scales="free_y")+ 
theme_classic()+
theme(axis.title.x = element_text(size=14),axis.text.x  = element_text(size=12,angle=60,hjust=1),axis.title.y=element_text(size=14),axis.text.y=element_text(size=12))
ggsave(t,file="$OUTDIR/BarcodeFreq/BarcodeFreq.pdf", width=20, height=10)

EOF

convert -density 300 $OUTDIR/BarcodeFreq/BarcodeFreq.pdf -quality 100 $OUTDIR/BarcodeFreq/BarcodeFreq.png

for i in `ls $OUTDIR/BarcodeFreq/*.png|sort|sed 's/.png//g' |awk -F'/' '{print $NF}'`;do echo '<a href="'${i}.pdf'"><img src="'${i}.png'" height="700"></a>' ;done  >$OUTDIR/BarcodeFreq/index.html

######
#iPCR specific 
######

if [[ `find -type d |grep BS|grep -v bak|grep iPCR| wc -l` -gt 1 ]]
then
       echo 'iPCR samples found'
       mkdir -p $OUTDIR/iPCR/
       for iPCR in `find -type d|grep iPCR|grep -v bak|sed 's/^..//g'`; do awk -v OFS='\t' -F'\t' -v sample="$iPCR" '{print $1, $2, $3, $4, $5, sample}' ${iPCR}/DistDpn.bed;done > $OUTDIR/iPCR/DistDpn.bed
       
       for iPCR in `find -type d|grep iPCR|grep -v bak|sed 's/^..//g'`; do awk -v OFS='\t' -F'\t' -v sample="$iPCR" '{print $1, sample}' ${iPCR}/DistToTSS.txt;done >$OUTDIR/iPCR/DistToTSS.bed
       
       for iPCR in `find -type d|grep iPCR|grep -v bak|sed 's/^..//g'`; do awk -v OFS='\t' -F'\t' -v sample="$iPCR"  '{print $1, $2, $3, $4, sample}' ${iPCR}/DistToDNase.bed;done >$OUTDIR/iPCR/DistToDNase.bed
       
       
       R --quiet --no-save << EOF
       library(stringr)
       library(reshape2)
       Dpn<-read("$OUTDIR/iPCR/DistDpn.bed",stringsAsFactors=F)
       colnames(Dpn)[1:6]<-c('chr','start','end','dpnsite','Distance','Sample')
       Dpn[,5]<-as.numeric(as.character(Dpn[,5]))
       
       p<-ggplot(Dpn,aes(x=Sample,y=Distance))+
       geom_boxplot(outlier.shape = NA)+
       #scale_fill_brewer(palette ='Set1')+
       theme_classic()+
       theme(axis.title.x = element_text(size=14),axis.text.x  = element_text(size=12,angle=90,hjust=1),axis.title.y=element_text(size=14),axis.text.y=element_text(size=12))+
       ylab('Distance to Dpn sites')+
       xlab('Samples')+
       scale_y_continuous(limits = c(0, 750))+
       coord_flip()
       ggsave(p, file="$OUTDIR/iPCR/DistDpn.pdf", width=12, height=10) 


       DNase<-read("$OUTDIR/iPCR/DistToDNase.bed",stringsAsFactors=F)
       colnames(DNase)[1:5]<-c('chr','start','end','Distance','Sample')
       DNase[,4]<-as.numeric(as.character(DNase[,4]))

       d<-ggplot(DNase, aes(x=abs(Distance+1))) +
       geom_histogram(color="black", size=0.25)+ 
       theme_classic()+
       facet_wrap(~Sample,scales = "free_y")+ 
       theme(axis.title.x = element_text(size=14),axis.text.x  = element_text(size=12),axis.title.y=element_text(size=14),axis.text.y=element_text(size=12))+
       #xlim(c(0,1500))+
       xlab('Distance to DNase site')+
       scale_x_log10(  breaks = scales::trans_breaks("log10", function(x) 10^x),  labels = scales::trans_format("log10", scales::math_format(10^.x)))+
       annotation_logticks(sides = "b",long = unit(0.2, "cm"))    
       ggsave(d, file="$OUTDIR/iPCR/DistDNase.pdf", width=20, height=10) 

       
       TSS <- read("$OUTDIR/iPCR/DistToTSS.bed")
       colnames(TSS) <- c("DistToTSS", "Sample")
       TSS\$DistToTSS.bin <- cut(TSS\$DistToTSS, breaks=c(-10^7, -10^6, -10^5, -10^4, -2500, 0, 5000, 10^4, 10^5, 10^6, 10^7), right=F, include.lowest=T, labels=c("-10 Mb", "-1 Mb", "-100 kb", "-10 kb", "-2.5 kb", "+5 kb", "+10 kb", "+100 kb", "+1 Mb", "+10 Mb"))
       
       
       t<-ggplot(TSS, aes(x=DistToTSS.bin)) +
       geom_bar(color="black", size=0.25) +
       theme_classic()+
       facet_wrap(~Sample,scales="free_y")+ 
       theme_classic()+
       theme(axis.title.x = element_text(size=14),axis.text.x  = element_text(size=12,angle=90,hjust=1),axis.title.y=element_text(size=14),axis.text.y=element_text(size=12))
       ggsave(t,file="$OUTDIR/iPCR/DistToTSS.pdf", width=20, height=10)

       
EOF

       convert -density 300 $OUTDIR/iPCR/DistDpn.pdf -quality 100 $OUTDIR/iPCR/DistDpn.png
       convert -density 300 $OUTDIR/iPCR/DistDNase.pdf -quality 100 $OUTDIR/iPCR/DistDNase.png
       convert -density 300 $OUTDIR/iPCR/DistToTSS.pdf -quality 100 $OUTDIR/iPCR/DistToTSS.png
       
       for i in `ls $OUTDIR/iPCR/*.png|sort|sed 's/.png//g' |awk -F'/' '{print $NF}'`;do echo '<a href="'${i}.pdf'"><img src="'${i}.png'" height="700"></a>' ;done |awk ' {print;} NR % 1 == 0 { print "<br>"; }' >$OUTDIR/iPCR/index.html
else 
       echo 'No iPCR samples found'     
       echo 'Exit'
fi

echo 'Done!!!'
date