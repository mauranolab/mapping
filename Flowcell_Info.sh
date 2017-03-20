#!/bin/bash

Dir=$1
Flow=$(basename `pwd`|sed 's/aligned.//g')
echo $Flow
OUTDIR=$Dir/$Flow
mkdir -p $OUTDIR


#####
#Weblogos
######
#For raw sequence



if [ "$Flow" != "Merged" ]
then
       mkdir -p $OUTDIR/Weblogos_raw
       find -name *.R2.raw.png |grep -v bak| xargs cp -t $OUTDIR/Weblogos_raw
       find -name *.R1.raw.png |grep -v bak| xargs cp -t $OUTDIR/Weblogos_raw
       for i in `ls $OUTDIR/Weblogos_raw/*R1.raw.png|sort|awk -F'/' '{print $NF}'`;do echo '<a href="'${i}'"><img src="'${i}'" height="120"></a>' ;done  >$OUTDIR/Weblogos_raw/R1index.html
       for i in `ls $OUTDIR/Weblogos_raw/*R2.raw.png|sort|awk -F'/' '{print $NF}'`;do echo '<a href="'${i}'"><img src="'${i}'"  style= "position:absolute; LEFT:900px"; height="120"></a>' ;done  >$OUTDIR/Weblogos_raw/R2index.html
       cat $OUTDIR/Weblogos_raw/R1index.html $OUTDIR/Weblogos_raw/R2index.html |sort| awk ' {print;} NR % 2 == 0 { print "<br>"; }'> $OUTDIR/Weblogos_raw/index.html
       
       
       
       #For processed sequence
       mkdir -p $OUTDIR/Weblogos_processed
       find -name *.BC.processed.png |grep -v bak| xargs cp -t $OUTDIR/Weblogos_processed
       find -name *.plasmid.processed.png |grep -v bak| xargs cp -t $OUTDIR/Weblogos_processed
       for i in `ls $OUTDIR/Weblogos_processed/*BC.processed.png|sort|awk -F'/' '{print $NF}'`;do echo '<a href="'${i}'"><img src="'${i}'" height="120"></a>' ;done  >$OUTDIR/Weblogos_processed/BCindex.html
       for i in `ls $OUTDIR/Weblogos_processed/*plasmid.processed.png|sort|awk -F'/' '{print $NF}'`;do echo '<a href="'${i}'"><img src="'${i}'"  style= "position:absolute; LEFT:900px"; height="120"></a>' ;done  >$OUTDIR/Weblogos_processed/plasmidindex.html
       cat $OUTDIR/Weblogos_processed/BCindex.html $OUTDIR/Weblogos_processed/plasmidindex.html |sort| awk ' {print;} NR % 2 == 0 { print "<br>"; }'> $OUTDIR/Weblogos_processed/index.html
       
fi

#For Barcode sequence
mkdir -p $OUTDIR/Weblogos_Barcode
find -name *barcodes.png |grep -v bak| xargs cp -t $OUTDIR/Weblogos_Barcode
for i in `ls $OUTDIR/Weblogos_Barcode/*barcodes.png|sort|awk -F'/' '{print $NF}'`;do echo '<a href="'${i}'"><img src="'${i}'" height="120"></a>' ;done  |awk ' {print;} NR % 1 == 0 { print "<br>"; }'> $OUTDIR/Weblogos_Barcode/index.html

#For UMI sequence
mkdir -p $OUTDIR/Weblogos_UMI
find -name *UMIs.png |grep -v bak| xargs cp -t $OUTDIR/Weblogos_UMI
for i in `ls $OUTDIR/Weblogos_UMI/*UMIs.png|sort|awk -F'/' '{print $NF}'`;do echo '<a href="'${i}'"><img src="'${i}'" height="120"></a>' ;done  |awk ' {print;} NR % 1 == 0 { print "<br>"; }'> $OUTDIR/Weblogos_UMI/index.html




#####
#Summarize flowcell info
#####
mkdir -p $OUTDIR/FlowcellSummary
R  --quiet --no-save << EOF

Barcodes<-list.files('./',pattern='^BS00',include.dirs = FALSE)
Barcodes<-Barcodes[grep('bak|minReads|BS00118A_19_20_121A|iPCRvsDNA',Barcodes,invert=T)]#REMOVE WHEN FILES ARE FINISHED
Barcodes<-Barcodes[grep('\\\.o[0-9]', Barcodes)]
Barcodes<-Barcodes[file.info(Barcodes)\$size>200]#Removes files with errors
sumBarcode<-as.data.frame(matrix(ncol=19,nrow=length(Barcodes))) 
colnames(sumBarcode)<-c('Flowcell','BSnumber','Name','sampleType','Total reads','Total barcodes','BC+UMI','UMI length','Unique BC','BC length','Complexity','PF reads','QC mapped reads','Prop PF','Unique sites','Unique sites +-5bp','#BC 1 site','#BC 2+ sites','Prop 2+')
for (files in 1:length(Barcodes)){cat(Barcodes[files],'\n')
    Sample<-readLines(Barcodes[files])
    Sample<-sub("^\\\s+","",Sample)
   sumBarcode[files,]\$BSnumber<-gsub('-.*','',Barcodes[files])
   sumBarcode[files,]\$Name<-gsub('_RNA\$|_DNA\$|_iPCR\$','',gsub('\\\.o.*|.*-','',Barcodes[files]))
   if(length(grep('Merged',Barcodes[files]))==1){sumBarcode[files,]\$sampleType<-gsub('.*_|\\\.o.*','',gsub('_Merged','',Barcodes[files]))}
   else{sumBarcode[files,]\$sampleType<-gsub('.*_|\\\.o.*','',Barcodes[files])}
   sumBarcode[files,]\$Flowcell<-gsub('.*/','',getwd())
   sumBarcode[files,][,5]<-splitLines('Number of total reads','\t')\$X3
   sumBarcode[files,][,6]<-splitLines('Number of total read barcodes','\t')\$X3
   if (length(grep('No UMIs found',Sample))==0){
       sumBarcode[files,][,7]<-splitLines('Number of unique barcodes+UMI','\t')\$X3
       sumBarcode[files,][,9]<-splitLines('Number of unique barcodes','\t')[2,]\$X3
       if(!is.na(getLength('UMI lengths')\$X2)) {sumBarcode[files,][,8]<-getLength('UMI lengths')\$X2}
          else if (getLength('UMI lengths')\$X1!=sumBarcode[files,][,7]){
                 sumBarcode[files,][,6] <- strsplit(Sample[grep(sumBarcode[files,][,5],Sample,fixed=T)][2], "\\\s+")[[1]][2]
                 if (is.na(strsplit(Sample[grep(sumBarcode[files,][,7],Sample,fixed=T)][2], "\\\s+")[[1]][2]) && length(Sample[grep(as.numeric(sumBarcode[files,][,7])-1,Sample,fixed=T)])==0){
                 sumBarcode[files,][,8]<-NA
                 }else sumBarcode[files,][,8] <- strsplit(Sample[grep(as.numeric(sumBarcode[files,][,7])-1,Sample,fixed=T)], "\\\s+")[[1]][2]
       }else sumBarcode[files,][,8]<-getLength('UMI lengths')\$X2
   }else if (length(splitLines('Number of unique barcodes','\t')\$X3)==2){sumBarcode[files,][,8]<-NA
   }else sumBarcode[files,][,9]<-splitLines('Number of unique barcodes','\t')\$X3
   if (getLength('Barcode lengths')\$X1!=sumBarcode[files,][,9]){
       sumBarcode[files,][,10] <- strsplit(Sample[grep(as.character(as.numeric(sumBarcode[files,][,9])-1),Sample,fixed=T)], "\\\s+")[[1]][2]
   }else sumBarcode[files,][,10] <- getLength('Barcode lengths')\$X2
   
    if (length(grep('Total PF tags',Sample))==1){ 
        sumBarcode[files,][,12]<-splitLines('Total PF tags','\t')\$X3
        sumBarcode[files,][,13]<-splitLines('Number of tags passing all filters and having barcodes assigned','\t')\$X2
        sumBarcode[files,][,15]<-splitLines('Total uniq sites','\t')\$X3[1]
        sumBarcode[files,][,16]<-splitLines('Total uniq sites (within 5 bp, ignoring strand)','\t')\$X3[1]
        sumBarcode[files,][17]<-getInsert(Sample)[1]
        sumBarcode[files,][18]<-getInsert(Sample)[2]
    
    }
        
}   

SummarizeFlowcell<-sumBarcode
SummarizeFlowcell\$Complexity<-0
for (i in 1:nrow(SummarizeFlowcell)){
    #Complexity
    if (is.na(SummarizeFlowcell[,7][i])) {SummarizeFlowcell\$Complexity[i]<-1
    } else SummarizeFlowcell\$Complexity[i]<-as.numeric(SummarizeFlowcell[,9][i])/as.numeric(SummarizeFlowcell[,7][i])
    cat(as.numeric(SummarizeFlowcell[,9][i])/as.numeric(SummarizeFlowcell[,7][i]),'\n')
    
    #Prop PF
    if(!is.na(SummarizeFlowcell[,12][i])) {SummarizeFlowcell[,14][i]<-as.numeric(SummarizeFlowcell[,13][i])/(as.numeric(SummarizeFlowcell[,12][i])+as.numeric(SummarizeFlowcell[,13][i]))}
    if(!is.na(SummarizeFlowcell[,12][i])) {SummarizeFlowcell[,19][i]<-as.numeric(SummarizeFlowcell[,18][i])/(as.numeric(SummarizeFlowcell[,17][i])+as.numeric(SummarizeFlowcell[,18][i]))}
}
  
  
  
for (i in 5:ncol(SummarizeFlowcell)){
SummarizeFlowcell[,i]<-format(as.numeric(SummarizeFlowcell[,i]),big.mark=",", trim=TRUE)
}  
library(tableHTML)
SumHTMLtable<-tableHTML(SummarizeFlowcell) %>%  add_css_row(css = list('background-color', 'lightblue'),rows = odd(1:nrow(SummarizeFlowcell)))
write_tableHTML(SumHTMLtable, file = "$OUTDIR/FlowcellSummary/index.html")

EOF


######
#Levenstein distance per sample
#######
if [ "$Flow" != "Merged" ]
then
       mkdir -p $OUTDIR/EditDist
       
       for i in `find \( -name "extract*" -o -name "map*" \)|grep -v bak`; do Dist=$(grep 'BC Hamming distance' $i); echo $i $Dist;done |awk -v OFS='\t' -F'[' '{print $1, $3}'|sed 's/]//g'| perl -pe 's/, /\t/g'|grep 'Hamming distance' > $OUTDIR/EditDist/BC_EditDist
       for i in `find \( -name "extract*" -o -name "map*" \)|grep -v bak`; do Dist=$(grep 'Plasmid Hamming distance' $i); echo $i $Dist;done |awk -v OFS='\t' -F'[' '{print $1, $3}'|sed 's/]//g'| perl -pe 's/, /\t/g'|grep 'Hamming distance' > $OUTDIR/EditDist/Plasmid_EditDist
       
       R  --quiet --no-save << EOF
       
       library(stringr)
       library(reshape2)
       BCleven<-read("$OUTDIR/EditDist/BC_EditDist",stringsAsFactors=F)
       BCleven[,1]<-gsub('\\\.o.*','',gsub('extract.|map.','',(gsub('.*/','',BCleven[,1]))))
       
       BCleven<-summaryBy(.~V1, data=BCleven,FUN=sum)
       colnames(BCleven)<-c('Sample',0:18)
       BCleven<-melt(BCleven)
       BCleven[,4]<-"Barcode"
       
       Plasmidleven<-read("$OUTDIR/EditDist/Plasmid_EditDist",stringsAsFactors=F)
       Plasmidleven[,1]<-gsub('\\\.o.*','',gsub('extract.|map.','',(gsub('.*/','',Plasmidleven[,1]))))
       
       Plasmidleven<-summaryBy(.~V1, data=Plasmidleven,FUN=sum)
       colnames(Plasmidleven)<-c('Sample',0:(ncol(Plasmidleven)-2))
       Plasmidleven<-melt(Plasmidleven)
       Plasmidleven[,4]<-'Plasmid'
       
       BCleven<-rbind(BCleven,Plasmidleven)
       colnames(BCleven)[1:4]<-c('Sample','Hamming','Reads','Read')
       BCleven[,2]<-as.numeric(as.character(BCleven[,2]))
       BCleven[,3]<-as.numeric(as.character(BCleven[,3]))
       
       colnames(BCleven)[1:4]<-c('Sample','Hamming','Reads','Read')
       
       
       myColors <- brewer.pal(4,"Set1")
       names(myColors)<-factor(c('RNA','iPCR','Plasmid','DNA'),levels=c('RNA','iPCR','Plasmid','DNA'))
       colScale <- scale_fill_manual(name = "Type",values = myColors)
       
       
       #plot
       BCleven\$Type<-gsub(".*_",'',BCleven\$Sample)
       BCleven\$Sample<- substr(BCleven\$Sample, 0, min(nchar(BCleven\$Sample)))
       #pdf("$OUTDIR/EditDist/LevenDistance_BC.pdf",height=10,width=20)
       bcl<-ggplot(BCleven[grep('Barcode',BCleven[,4]),], aes(Hamming,Reads,fill=Type)) + 
       geom_bar(stat="identity",color='black',size=0.25,alpha=0.4) +
       #ggtitle('GGlo')+
       facet_wrap(~Sample+Read,scales = "free_y")+ 
       theme_classic()+
       colScale+
       theme(legend.position = "bottom",legend.text=element_text(size=8))+
       #scale_x_continuous(label= fancy_scientific)+
       scale_y_continuous(label= fancy_scientific)
       #dev.off()
       gp = ggplotGrob(bcl)
       gp = editGrob(grid.force(gp), gPath("GRID.stripGrob", "GRID.text"), grep = TRUE, global = TRUE, just = "left", x = unit(0.1,"npc"))
       ggsave(gp ,file="$OUTDIR/EditDist/LevenDistance_BC.pdf", width=20, height=10)
       
       
       #plot
       #pdf("$OUTDIR/EditDist/LevenDistance_Plasmid.pdf",height=10,width=20)
       pll<-ggplot(BCleven[grep('Plasmid',BCleven[,4]),], aes(Hamming,Reads,fill=Type)) + 
       geom_bar(stat="identity",color="black", size=0.1,alpha=0.4) +
       #ggtitle('GGlo')+
       facet_wrap(~Sample+Read,scales = "free_y")+ 
       theme_classic()+
       colScale+
       #scale_fill_brewer('Set1')+
       theme(legend.position = "bottom",legend.text=element_text(size=8))+
       #scale_x_continuous(label= fancy_scientific)+
       scale_y_continuous(label= fancy_scientific)
       #dev.off()
       gp = ggplotGrob(pll)
       gp = editGrob(grid.force(gp), gPath("GRID.stripGrob", "GRID.text"), grep = TRUE, global = TRUE, just = "left", x = unit(0.1,"npc"))
       ggsave(gp ,file="$OUTDIR/EditDist/LevenDistance_Plasmid.pdf", width=20, height=10)
       
EOF
       
       convert -density 300 $OUTDIR/EditDist/LevenDistance_BC.pdf -quality 100 $OUTDIR/EditDist/LevenDistance_BC.png
       convert -density 300 $OUTDIR/EditDist/LevenDistance_Plasmid.pdf -quality 100 $OUTDIR/EditDist/LevenDistance_Plasmid.png
       
       for i in `ls $OUTDIR/EditDist/*.png|sort|sed 's/.png//g' |awk -F'/' '{print $NF}'`;do echo '<a href="'${i}.pdf'"><img src="'${i}.png'" height="700"></a>' ;done  >$OUTDIR/EditDist/index.html
       

fi
#####
#Barcode frequency
#####
mkdir -p $OUTDIR/BarcodeFreq
for i in `find -name *barcode.counts.txt|sed 's/^..//g'|sed 's/\/.*//g'|grep -v bak`; do sort -nk2 ${i}/${i}.barcode.counts.txt| awk -v OFS='\t' -F'\t' -v sample="$i" '{print $1, $2, sample}';done > $OUTDIR/BarcodeFreq/BarcodeFreq.txt


R  --quiet --no-save << EOF
BC <- read("$OUTDIR/BarcodeFreq/BarcodeFreq.txt")
colnames(BC) <- c("Barcodes","barcodeFreq", "Sample")
BC\$BC.bin <- cut(BC\$barcodeFreq, breaks=c( 0, 1.1, 5.1, 10.1, 50.1, 100.1, 1000.1, 10000.1, 100000.1, 1000000.1), right=F, include.lowest=T, labels=c("0-1", "1-5", "5-10", "10-50", "50-100", "100-1000", "100-10000", "10000-100000","+100000"))
BC\$Type<-gsub(".*_",'',gsub('_Merged','',BC\$Sample))

myColors <- brewer.pal(4,"Set1")
names(myColors)<-factor(c('RNA','iPCR','Plasmid','DNA'),levels=c('RNA','iPCR','Plasmid','DNA'))
colScale <- scale_fill_manual(name = "Type",values = myColors)
BC\$Sample<- substr(BC\$Sample, 0, min(nchar(BC\$Sample)))

t<-ggplot(BC, aes(x=BC.bin,fill=Type)) +
geom_bar(color="black", size=0.25,alpha=0.4) +
theme_classic()+
colScale+
facet_wrap(~Sample,scales="free_y")+ 
theme_classic()+
theme(axis.title.x = element_text(size=14),axis.text.x  = element_text(size=12,angle=60,hjust=1),axis.title.y=element_text(size=14),axis.text.y=element_text(size=12))
#ggsave(t,file="$OUTDIR/BarcodeFreq/BarcodeFreq.pdf", width=20, height=10)

gp = ggplotGrob(t)
gp = editGrob(grid.force(gp), gPath("GRID.stripGrob", "GRID.text"), grep = TRUE, global = TRUE, just = "left", x = unit(0.1,"npc"))
ggsave(gp ,file="$OUTDIR/BarcodeFreq/BarcodeFreq.pdf", width=20, height=10)


EOF

rm $OUTDIR/BarcodeFreq/BarcodeFreq.txt
convert -density 300 $OUTDIR/BarcodeFreq/BarcodeFreq.pdf -quality 100 $OUTDIR/BarcodeFreq/BarcodeFreq.png

for i in `ls $OUTDIR/BarcodeFreq/*.png|sort|sed 's/.png//g' |awk -F'/' '{print $NF}'`;do echo '<a href="'${i}.pdf'"><img src="'${i}.png'" height="700"></a>' ;done  >$OUTDIR/BarcodeFreq/index.html



#####
#Saturation curve
#####
if [[ `find -name *Saturation*|sed 's/^..//g'|sed 's/\/.*//g'|grep -v bak| wc -l` -ge 1 ]]
       then
       mkdir -p $OUTDIR/SaturationCurve
       for i in `find -name *Saturation*|sed 's/^..//g'|sed 's/\/.*//g'|grep -v bak`; do cat ${i}/${i}.Saturation*| awk -v OFS='\t' -F'\t' -v sample="$i" '{print $1, $2, sample}';done|perl -pe 's/ /\t/g' > $OUTDIR/SaturationCurve/SaturationCurve.txt
       
       
       R  --quiet --no-save << EOF
       SaturationCurve <- read("$OUTDIR/SaturationCurve/SaturationCurve.txt")
       
       #SaturationCurve<-SaturationCurve[!is.na(SaturationCurve$\Unique_BC),]
       SaturationCurve[,2]<-as.numeric(SaturationCurve[,2])
       SaturationCurve[,1]<-as.numeric(SaturationCurve[,1])
       SaturationCurve[,3]<-gsub('minreads','',SaturationCurve[,3])
       SaturationCurve\$Type<-gsub(".*_",'',SaturationCurve[,5])
       colnames(SaturationCurve)[1:6] <- c("Reads","Unique_BC","minReads","NA","Sample",'Type')
       
       SaturationCurve\$Sample<- substr(SaturationCurve\$Sample, 0, min(nchar(SaturationCurve\$Sample)))
       myColors <- brewer.pal(4,"Set1")
       names(myColors)<-factor(c('RNA','iPCR','Plasmid','DNA'),levels=c('RNA','iPCR','Plasmid','DNA'))
       colScale <- scale_colour_manual(name = "Type",values = myColors)
       
       
       s1 <- ggplot(SaturationCurve, aes(Reads, Unique_BC,color=Type)) + 
       geom_line(aes(linetype=minReads)) +
       #geom_line()+
       #ggtitle('GGlo')+
       facet_wrap(~Sample,scales = "free")+ 
       theme_classic()+
       colScale+
       theme(legend.position = "bottom",legend.text=element_text(size=8))+
       scale_x_continuous(label= fancy_scientific)+
       scale_y_continuous(label= fancy_scientific)
       
       gp = ggplotGrob(s1)
       gp = editGrob(grid.force(gp), gPath("GRID.stripGrob", "GRID.text"), grep = TRUE, global = TRUE, just = "left", x = unit(0.1,"npc"))
       ggsave(gp ,file="$OUTDIR/SaturationCurve/SaturationCurve.pdf", width=20, height=10)
       
       
EOF
       
       convert -density 300 $OUTDIR/SaturationCurve/SaturationCurve.pdf -quality 100 $OUTDIR/SaturationCurve/SaturationCurve.png
       
       for i in `ls $OUTDIR/SaturationCurve/*.png|sort|sed 's/.png//g' |awk -F'/' '{print $NF}'`;do echo '<a href="'${i}.pdf'"><img src="'${i}.png'" height="700"></a>' ;done  >$OUTDIR/SaturationCurve/index.html
fi

######
#iPCR specific 
######

if [[ `find -type d |grep BS|grep -v bak|grep iPCR| wc -l` -ge 1 ]]
then
       echo 'iPCR samples found'
       mkdir -p $OUTDIR/iPCR/
       if [[ `find -name *DistDpn.bed|sed 's/^..//g'|sed 's/\/.*//g'|grep -v bak| wc -l` -ge 1 ]]
       then
              for iPCR in `find -name *DistDpn.bed|sed 's/^..//g'|sed 's/\/.*//g'|grep -v bak`; do awk -v OFS='\t' -F'\t' -v sample="$iPCR" '{print $1, $2, $3, $4, $5, sample}' ${iPCR}/DistDpn.bed;done > $OUTDIR/iPCR/DistDpn.bed
       fi
       
       if [[ `find -name *DistToTSS.txt|sed 's/^..//g'|sed 's/\/.*//g'|grep -v bak| wc -l` -ge 1 ]]
       then
              for iPCR in `find -name *DistToTSS.txt|sed 's/^..//g'|sed 's/\/.*//g'|grep -v bak`; do awk -v OFS='\t' -F'\t' -v sample="$iPCR" '{print $1, sample}' ${iPCR}/DistToTSS.txt;done >$OUTDIR/iPCR/DistToTSS.bed
       fi
       
       if [[ `find -name *ObsvsExp.txt|sed 's/^..//g'|sed 's/\/.*//g'|grep -v bak| wc -l` -ge 1 ]]
       then
              for iPCR in `find -name *ObsvsExp.txt|sed 's/^..//g'|sed 's/\/.*//g'|grep -v bak`; do awk -v OFS='\t' -F'\t' -v sample="$iPCR" '{print sample, $0}' ${iPCR}/${iPCR}.ObsvsExp.txt;done >$OUTDIR/iPCR/ObsvsExp.txt
       fi
       
       if [[ `find -name *DistToDNase.bed|sed 's/^..//g'|sed 's/\/.*//g'|grep -v bak| wc -l` -ge 1 ]]
       then
              for iPCR in `find -name *DistToDNase.bed|sed 's/^..//g'|sed 's/\/.*//g'|grep -v bak`; do awk -v OFS='\t' -F'\t' -v sample="$iPCR"  '{print $1, $2, $3, $4, sample}' ${iPCR}/DistToDNase.bed;done >$OUTDIR/iPCR/DistToDNase.bed
       
              
              R --quiet --no-save << EOF
              library(stringr)
              library(reshape2)
              myColors <- brewer.pal(4,"Set1")
              names(myColors)<-factor(c('RNA','iPCR','Plasmid','DNA'),levels=c('RNA','iPCR','Plasmid','DNA'))
              colScale <- scale_fill_manual(name = "Type",values = myColors)

              
              Dpn<-read("$OUTDIR/iPCR/DistDpn.bed",stringsAsFactors=F)
              colnames(Dpn)[1:6]<-c('chr','start','end','dpnsite','Distance','Sample')
              Dpn[,5]<-as.numeric(as.character(Dpn[,5]))
              Dpn\$Type<-gsub(".*_",'',gsub('_Merged','',Dpn\$Sample))

              p<-ggplot(Dpn,aes(x=Sample,y=Distance,fill=Type))+
              geom_boxplot(outlier.shape = NA,alpha=0.4)+
              #scale_fill_brewer(palette ='Set1')+
              theme_classic()+
              theme(axis.title.x = element_text(size=14),axis.text.x  = element_text(size=12,angle=90,hjust=1),axis.title.y=element_text(size=14),axis.text.y=element_text(size=12))+
              ylab('Distance to Dpn sites')+
              xlab('Samples')+
              colScale+
              scale_y_continuous(limits = c(0, 750))+
              coord_flip()
              ggsave(p, file="$OUTDIR/iPCR/DistDpn.pdf", width=12, height=10) 
       
       
              DNase<-read("$OUTDIR/iPCR/DistToDNase.bed",stringsAsFactors=F)
              colnames(DNase)[1:5]<-c('chr','start','end','Distance','Sample')
              DNase[,4]<-as.numeric(as.character(DNase[,4]))
              DNase\$Type<-gsub(".*_",'',gsub('_Merged','',DNase\$Sample))
              
              d<-ggplot(DNase, aes(x=abs(Distance+1),fill=Type)) +
              geom_histogram(color="black", size=0.25,alpha=0.4)+ 
              theme_classic()+
              facet_wrap(~Sample,scales = "free_y")+ 
              theme(axis.title.x = element_text(size=14),axis.text.x  = element_text(size=12),axis.title.y=element_text(size=14),axis.text.y=element_text(size=12))+
              #xlim(c(0,1500))+
              xlab('Distance to DNase site')+
              colScale+
              scale_x_log10(  breaks = scales::trans_breaks("log10", function(x) 10^x),  labels = scales::trans_format("log10", scales::math_format(10^.x)))+
              annotation_logticks(sides = "b",long = unit(0.2, "cm"))    
              ggsave(d, file="$OUTDIR/iPCR/DistDNase.pdf", width=20, height=10) 
       
              
              TSS <- read("$OUTDIR/iPCR/DistToTSS.bed")
              colnames(TSS) <- c("DistToTSS", "Sample")
              TSS\$DistToTSS.bin <- cut(TSS\$DistToTSS, breaks=c(-10^7, -10^6, -10^5, -10^4, -2500, 0, 5000, 10^4, 10^5, 10^6, 10^7), right=F, include.lowest=T, labels=c("-10 Mb", "-1 Mb", "-100 kb", "-10 kb", "-2.5 kb", "+5 kb", "+10 kb", "+100 kb", "+1 Mb", "+10 Mb"))
              
              TSS\$Type<-gsub(".*_",'',gsub('_Merged','',TSS\$Sample))
              t<-ggplot(TSS, aes(x=DistToTSS.bin,fill=Type)) +
              geom_bar(color="black", size=0.25,alpha=0.4) +
              theme_classic()+
              colScale+
              facet_wrap(~Sample,scales="free_y")+ 
              theme_classic()+
              theme(axis.title.x = element_text(size=14),axis.text.x  = element_text(size=12,angle=90,hjust=1),axis.title.y=element_text(size=14),axis.text.y=element_text(size=12))
              ggsave(t,file="$OUTDIR/iPCR/DistToTSS.pdf", width=20, height=10)
              
              
              
              ObsvsExp<-read('$OUTDIR/iPCR/ObsvsExp.txt',header=T)
              colnames(ObsvsExp)[1]<-c('Sample')
              ObsvsExp<-ObsvsExp[grep('l2fold',ObsvsExp\$l2fold,invert=T),]
              ObsvsExp\$l2fold<-as.numeric(ObsvsExp\$l2fold)
              ObsvsExp\$significant<-as.numeric(ObsvsExp\$qvalue)<0.05
              
              
              myColors <- c('#a6cee3','#e41a1c')
              names(myColors)<-factor(c('FALSE','TRUE'),levels=c('FALSE','TRUE'))
              colScale <- scale_fill_manual(name = "significant",values = myColors)

              
              o<-ggplot(ObsvsExp, aes(x=annotation,y=l2fold,fill=significant)) +
              geom_bar(stat='identity',color="black", size=0.25,alpha=0.4)+ 
              theme_classic()+
              colScale+
              #scale_fill_manual(values=c('#a6cee3','#e41a1c'))+
              theme(axis.title.x = element_text(size=14),axis.text.x  = element_text(size=12,angle=90,hjust=1),axis.title.y=element_text(size=14),axis.text.y=element_text(size=12))+
              facet_wrap(~Sample)+ 
              ggtitle("Obs/Exp")+
              ylab('log2FC(Observed/Expected)')
              ggsave(o,file="$OUTDIR/iPCR/ObsvsExp.pdf", width=20, height=10)

              
              
EOF

              convert -density 300 $OUTDIR/iPCR/DistDpn.pdf -quality 100 $OUTDIR/iPCR/DistDpn.png
              convert -density 300 $OUTDIR/iPCR/DistDNase.pdf -quality 100 $OUTDIR/iPCR/DistDNase.png
              convert -density 300 $OUTDIR/iPCR/DistToTSS.pdf -quality 100 $OUTDIR/iPCR/DistToTSS.png
              convert -density 300 $OUTDIR/iPCR/ObsvsExp.pdf -quality 100 $OUTDIR/iPCR/ObsvsExp.png
              for i in `ls $OUTDIR/iPCR/*.png|sort|sed 's/.png//g' |awk -F'/' '{print $NF}'`;do echo '<a href="'${i}.pdf'"><img src="'${i}.png'" height="700"></a>' ;done |awk ' {print;} NR % 1 == 0 { print "<br>"; }' >$OUTDIR/iPCR/index.html
       fi
       



else 
       echo 'No iPCR samples found'     
       echo 'Exit'
fi

echo 'Done!!!'
date