#!/bin/bash
set -eu  -o pipefail

src="/vol/mauranolab/mapped/src/transposon"

module load samtools/1.9


fc=$(basename `pwd`)

OUTDIR=$1/${fc}
mkdir -p $OUTDIR

echo "Output to $OUTDIR"


dirs=`find -maxdepth 1 -not -path "*/bak*" -not -path "*/trash*" -type d`


#####
#TODO fix if only HiC samples are running
#####
#######
#HiC data
######
if [[ `find ${dirs} -maxdepth 1 -type d | grep 'HiC\|CapC\|dsDNA\|3C'| wc -l` -ge 1 ]]; then
    echo 'HiC samples found'
    mkdir -p $OUTDIR/HiC/
    
    
    find ${dirs} -maxdepth 1 -name *.R2.raw.png | grep 'HiC\|CapC\|dsDNA\|3C'| xargs --no-run-if-empty cp -t $OUTDIR/HiC/
    find ${dirs} -maxdepth 1 -name *.R1.raw.png | grep 'HiC\|CapC\|dsDNA\|3C'| xargs --no-run-if-empty cp -t $OUTDIR/HiC/
    for i in `ls $OUTDIR/HiC/*R1.raw.png | sort | awk -F "/" '{print $NF}'`; do echo '<a href="'${i}'"><img src="'${i}'" height="120"></a>' ; done  >$OUTDIR/HiC/R1index.html
    for i in `ls $OUTDIR/HiC/*R2.raw.png | sort | awk -F "/" '{print $NF}'`; do echo '<a href="'${i}'"><img src="'${i}'"  style= "position:absolute; LEFT:1500px"; height="120"></a>' ; done  >$OUTDIR/HiC/R2index.html
    cat $OUTDIR//HiC/R1index.html $OUTDIR/HiC//R2index.html |sort| awk ' {print;} NR % 2 == 0 { print "<br>"; }'> $OUTDIR/HiC/Weblogoindex.html
    
#    if [[ `find ${dirs} -name *R1.mmapstat |sed 's/^..//g' | sed 's/\/.*//g' | wc -l` -ge 1 ]]
#    then
#        for HiC in `find ${dirs} -name *R1.mmapstat | sed 's/^..//g' | sed 's/\/.*//g' `; do awk -v OFS='\t' -F'\t' -v HiC="$HiC" '{print $1, $2, "R1", HiC}' ${HiC}/bowtie_results/bwt2/${HiC}/*R1.mmapstat |
#        grep -v "#" ; done > $OUTDIR/HiC/mmapstatR1.tsv
#        for HiC in `find ${dirs} -name *R2.mmapstat | sed 's/^..//g' | sed 's/\/.*//g' `; do awk -v OFS='\t' -F'\t' -v HiC="$HiC" '{print $1, $2, "R2", HiC}' ${HiC}/bowtie_results/bwt2/${HiC}/*R2.mmapstat |
#        grep -v "#" ; done >$OUTDIR/HiC/mmapstatR2.tsv
#    fi
#    
#    if [[ `find ${dirs} -name *mpairstat |sed 's/^..//g' | sed 's/\/.*//g' | wc -l` -ge 1 ]]
#    then
#        for HiC in `find ${dirs} -name *mmapstat | sed 's/^..//g' | sed 's/\/.*//g' `; do awk -v OFS='\t' -F'\t' -v HiC="$HiC" '{print $1, $2, $3, HiC}' ${HiC}/bowtie_results/bwt2/${HiC}/*mpairstat |
#        grep -v "#" ; done > $OUTDIR/HiC/mpairstat.tsv
#    fi
#    
#    if [[ `find ${dirs} -maxdepth 2 -name *_R1.bam |grep -v bak| wc -l` -ge 1 ]]
#    then
#        for HiC in `find ${dirs} -maxdepth 1 -type d | grep 'HiC\|CapC\|dsDNA\|3C'|sed 's/^..//g'`; do samtools view -F4 -q30 ${HiC}/${HiC}_R1.bam |
#        cut -f3 | sort | uniq -c | awk -v OFS='\t' -v name="$HiC" '{print $2, $1, "R1", name}' ; done > $OUTDIR/HiC/SE_mappingR1.tsv
#        for HiC in `find ${dirs} -maxdepth 1 -type d | grep 'HiC\|CapC\|dsDNA\|3C'|sed 's/^..//g'`; do samtools view -F4 -q30 ${HiC}/${HiC}_R2.bam |
#        cut -f3 | sort | uniq -c | awk -v OFS='\t' -v name="$HiC" '{print $2, $1, "R2", name}' ; done > $OUTDIR/HiC/SE_mappingR2.tsv
#    fi
#    
#    if [[ `find ${dirs} -name *.bwt2glob.unmap_trimmed.fastq | wc -l` -ge 1 ]]
#    then
#        for HiC in `find ${dirs} -name *.bwt2glob.unmap_trimmed.fastq | awk -F '/' '{print $2}'  | sort | uniq`; do cat ${HiC}/bowtie_results/bwt2_global/${HiC}/*R1*.bwt2glob.unmap_trimmed.fastq |
#        awk 'NR%4==2 {print length($1)}' | sort | uniq -c | awk -v OFS='\t' -v sample="$HiC" '{print $2, $1, "R1", sample}' ; done > $OUTDIR/HiC/trimmedReadsR1.tsv
#        for HiC in `find ${dirs} -name *.bwt2glob.unmap_trimmed.fastq | awk -F '/' '{print $2}'  | sort | uniq`; do cat ${HiC}/bowtie_results/bwt2_global/${HiC}/*R2*.bwt2glob.unmap_trimmed.fastq |
#        awk 'NR%4==2 {print length($1)}' | sort | uniq -c | awk -v OFS='\t' -v sample="$HiC" '{print $2, $1, "R2", sample}'; done > $OUTDIR/HiC/trimmedReadsR2.tsv
#    fi
#    
#    if [[ `find ${dirs} -name *mergestat | sed 's/^..//g' | sed 's/\/.*//g' | wc -l` -ge 1 ]]
#    then
#        for HiC in `find ${dirs} -name *mergestat | sed 's/^..//g' | sed 's/\/.*//g' `; do awk -v OFS='\t' -F'\t' -v HiC="$HiC" '{print $1, $2, HiC}' ${HiC}/hic_results/data/${HiC}/*mergestat |
#        grep -v "#" ; done > $OUTDIR/HiC/mergestat.tsv
#    fi
#    
#    if [[ `find ${dirs} -name *mRSstat | sed 's/^..//g' | sed 's/\/.*//g' | wc -l` -ge 1 ]]
#    then
#        for HiC in `find ${dirs} -name *mRSstat | sed 's/^..//g' | sed 's/\/.*//g' `; do awk -v OFS='\t' -F'\t' -v HiC="$HiC" '{print $1, $2, HiC}' ${HiC}/hic_results/data/${HiC}/*mRSstat|
#        grep -v "#" ; done > $OUTDIR/HiC/mRSstat.tsv
#    fi
#    
#     for HiC in `find ${dirs} -maxdepth 1 -type d | grep 'HiC\|CapC\|dsDNA\|3C' | sed 's/^..//g'`; do reads=$(grep Surviving submit.${HiC}* | awk '{print $4}'); \
#        adapter=$(grep Surviving submit.${HiC}* | awk '{print $7}'); \
#        validPairs=$(cat ${HiC}/hic_results/data/${HiC}/${HiC}_allValidPairs | wc -l); \
#        validPairsSameChrom=$(cat ${HiC}/hic_results/data/${HiC}/${HiC}_allValidPairs | awk '$2==$5' | wc -l); \
#        gatcR1=$(zcat ${HiC}/${HiC}_R1.fastq.gz| grep 'GATC' | wc -l); \
#        gatcR2=$(zcat ${HiC}/${HiC}_R2.fastq.gz| grep 'GATC' | wc -l); \
#        numLinkR1=$(zcat ${HiC}/${HiC}_R1.fastq.gz| grep 'CGCGATATCTTATCTGA\|TCAGATAAGATATCGC' | wc -l); \
#        numLinkR2=$(zcat ${HiC}/${HiC}_R2.fastq.gz | grep 'CGCGATATCTTATCTGA\|TCAGATAAGATATCGC' | wc -l); \
#        mappedR1=$(samtools view -F4 -q30 ${HiC}/${HiC}_R1.sorted.bam | wc -l); \
#        mappedR2=$(samtools view -F4 -q30 ${HiC}/${HiC}_R2.sorted.bam | wc -l); \
#        dangling=$(grep Dangling ${HiC}/hic_results/data/${HiC}/${HiC}.mRSstat | awk '{print $2}'); \
#        echo ${fc} $HiC $reads $adapter $mappedR1 $mappedR2 $gatcR1 $gatcR2 $numLinkR1 $numLinkR2 $dangling $validPairs $validPairsSameChrom  ; done | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13}' > $OUTDIR/HiC/HiC_summary.tsv
#          
#
#    
#    R --quiet --no-save << EOF
#    library(stringr)
#    library(reshape2)
#    library(data.table)
#    library(dplyr)
#    library(gridExtra)
#    library(gtools)
#    library(tableHTML)
#    
#    
#    #####
#    #Create HiC summary table
#    #####
#    HiCsummary <- read.table("$OUTDIR/HiC/HiC_summary.tsv", header=F,stringsAsFactors=F)
#    
#    for (i in 3:ncol(HiCsummary)){
#        HiCsummary[,i]<-format(as.numeric(HiCsummary[,i]),big.mark=",", trim=TRUE)
#    
#    }
#    HiCsummary <- HiCsummary[order(HiCsummary[,"V2"]),]
#    rownames(HiCsummary) <-1:nrow(HiCsummary)
#    colnames(HiCsummary) <- c('Flowcell', 'Sample', 'Reads', 'Survived trimming', 'SE map R1 (MAPQ>=30)', 'SE map R2(MAPQ>=30)', 'gatcR1', 'gatcR2', 'linkerR1', 'linkerR2', 'dangling ends', 'validPairs', 'cis-validPairs')
#    SumHTMLtable<-tableHTML(HiCsummary) %>%  add_css_row(css = list('background-color', 'lightblue'),rows = odd(1:nrow(HiCsummary)))
#    write_tableHTML(SumHTMLtable, file = "$OUTDIR/HiC/HiCsummary.html")
#    
#    
#    #####
#    #Plot SE mapping
#    #####
#    R1 <- fread("$OUTDIR/HiC/SE_mappingR1.tsv")
#    R2 <- fread("$OUTDIR/HiC/SE_mappingR2.tsv")
#    SEmapped <- rbind(R1,R2)
#    SEmapped <- as.data.frame(SEmapped)
#    SEmapped <- SEmapped[grep('consensus|chrM|Un|alt|random|chrM|chrY',SEmapped[,1], invert=T),]
#    hg38 <- read.table('/vol/isg/annotation/fasta/hg38/hg38.chrom.sizes')
#    hg38 <- hg38[grep('Un|alt|random|chrM|chrY',hg38[,"V1"],invert=T),]
#    hg38[,"Mb"] <- hg38[,"V2"]/1e6
#    colnames(hg38) <-c('V1', 'size','Mb')
#
#    SEmapped <- SEmapped %>% 
#        arrange(V1) %>%
#        left_join(hg38, by='V1') %>%
#        mutate(normReads=V2/Mb)
#    head(SEmapped)
#
#    
#    height <- ceiling(length(unique(SEmapped[,"V4"]))/3)*2
#    
#    ggSEmapped <- SEmapped %>%
#        mutate(V1 = factor(V1, levels=unique(mixedsort(V1)))) %>%
#        ggplot(aes(x=V1, y=normReads,fill=V3))+
#        geom_bar(stat='identity') +
#        facet_wrap(~V4, ncol=3)+
#        theme_classic() +
#        xlab('Chromosomes') +
#        ylab('Mapped reads/Mb') +
#        scale_fill_manual(values=rev(brewer.pal(n=3,"Set1")[1:2]))+
#        ggtitle('Normalised reads per chromosome (MAPQ >=30)')+
#        theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12,angle=60,hjust=1), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))
#    ggsave(ggSEmapped,file="$OUTDIR/HiC/1.SE_perChrom_mapping.pdf", width=13, height=height)
#
#    
#    
#    #####
#    #Plot mapping stat
#    #####
#    R1 <- read.table("$OUTDIR/HiC/mmapstatR1.tsv", sep='\t', header=F, stringsAsFactors=F)
#    R2 <- read.table("$OUTDIR/HiC/mmapstatR2.tsv", sep='\t', header=F, stringsAsFactors=F)
#    
#    mapStat <- rbind(R1,R2)
#    mapped <- mapStat[grep('mapped',mapStat\$V1, invert=T),]
#    mappedUnmelt <- dcast(mapped,formula=V1~V4+V3,value.var='V2')
#    mappedUnmelt[4,] <- c('Unmapped', as.numeric(mappedUnmelt[,-1][3,])- (as.numeric(mappedUnmelt[,-1][2,])+as.numeric(mappedUnmelt[,-1][1,])))
#    mapped<- melt(mappedUnmelt,id='V1')
#    mapped <- mapped[grep('total',mapped\$V1, invert=T),]
#    mapped\$Read <- substr(as.character(mapped\$variable),nchar(as.character(mapped\$variable))-1, nchar(as.character(mapped\$variable)))
#    mapped\$variable <- substr(as.character(mapped\$variable),1, nchar(as.character(mapped\$variable))-3)
#    
#    mapped[grep('global',mapped\$V1),]\$V1 <-'Full_read_mapped'
#    mapped[grep('local',mapped\$V1),]\$V1 <-'Trimmed_read_mapped'
#    mapped\$value <-as.numeric(as.character(mapped\$value))
#
#    colnames(mapped) <- c('Mapping', 'Sample','numberReads', 'Read')
#      
#    mapped\$Mapping <- factor(mapped\$Mapping,levels=c( 'Unmapped', 'Trimmed_read_mapped','Full_read_mapped' ))
#    #mapped\$Sample <- gsub('MSH_061517_GM12865~DpnII~','',mapped\$Sample)
#    mapped\$Sample <- substr(mapped\$Sample,1, 8)
#    t<-ggplot(mapped, aes(x=Read,y=numberReads, fill=Mapping)) +
#    geom_bar(stat='identity',color="black", size=0.25) +
#    theme_classic()+
#    scale_fill_manual(values=c('#bdbdbd','#9ecae1','#4292c6'))+
#    facet_wrap(~Sample, ncol=3)+ 
#    theme_classic()+
#    ggtitle('Mapped reads (MAPQ>=0)')+
#    scale_y_continuous(label= fancy_scientific) +
#    theme(axis.title.x = element_text(size=14),axis.text.x  = element_text(size=12,angle=90,hjust=1),axis.title.y=element_text(size=14),axis.text.y=element_text(size=12))
#    
#    pdf(NULL)
#    gp = ggplotGrob(t)
#    gp = editGrob(grid.force(gp), gPath("GRID.stripGrob", "GRID.text"), grep = TRUE, global = TRUE, just = "left", x = unit(0.1,"npc"))
#    
#    ggsave(gp,file="$OUTDIR/HiC/2.mapping.pdf", width=13, height=height)
#    
#    
#    
#    #####
#    #Plot merged interaction stat
#    #####
#    mergeStat <- read.table("$OUTDIR/HiC/mergestat.tsv", header=F, sep='\t', stringsAsFactors=F)
#    
#    mergeStat <- mergeStat[grep('valid_interaction|Range', mergeStat\$V1, invert=T),]
#    colnames(mergeStat) <- c('Interaction', 'numberReads', 'Sample')
#    mergeStat\$Sample <- substr(mergeStat\$Sample,1, 8)
#    mergeStat\$Interaction <- gsub("_interaction",'',mergeStat\$Interaction)
#    
#    t<-ggplot(mergeStat, aes(x=Interaction,y=numberReads, fill=Interaction)) +
#    geom_bar(stat='identity',color="black", size=0.25) +
#    theme_classic()+
#    scale_fill_manual(values=c('#9ecae1','#4292c6'))+
#    facet_wrap(~Sample, ncol=3)+ 
#    theme_classic()+
#    scale_y_continuous(label= fancy_scientific) +
#    theme(axis.title.x = element_text(size=14),axis.text.x  = element_text(size=12,angle=90,hjust=1),axis.title.y=element_text(size=14),axis.text.y=element_text(size=12))
#    
#    ggsave(t,file="$OUTDIR/HiC/3.Interactions.pdf", width=13, height=height)
#    
#    
#    #####
#    #Plot valid pairs
#    #####
#    mRSstat <- read.table("$OUTDIR/HiC/mRSstat.tsv", header=F, stringsAsFactors=F, sep='\t')
#    
#    mRSstat\$Passed <- ''
#    mRSstat[grep('Valid', mRSstat\$V1),]\$Passed <- 'Valid'
#    mRSstat[grep('Valid', mRSstat\$V1, invert=T),]\$Passed <- 'Dumped'
#    
#    mRSstat<- mRSstat[!mRSstat\$V1 %in% 'Valid_interaction_pairs',]
#    colnames(mRSstat) <- c('Interaction', 'Pairs', 'Sample','Passed')
#    mRSstat\$Sample <- substr(mRSstat\$Sample,1, 8)
#    
#    t<-ggplot(mRSstat, aes(x=Sample,y=Pairs, fill=Interaction)) +
#    geom_bar(stat='identity',color="black", size=0.25) +
#    theme_classic()+
#    #scale_fill_manual(values=c('#9ecae1','#4292c6'))+
#    facet_wrap(~Passed, ncol=3)+ 
#    theme_classic()+
#    scale_y_continuous(label= fancy_scientific) +
#    theme(axis.title.x = element_text(size=14),axis.text.x  = element_text(size=12,angle=90,hjust=1),axis.title.y=element_text(size=14),axis.text.y=element_text(size=12))
#    
#    ggsave(t,file="$OUTDIR/HiC/4.ValidPairs.pdf", width=13, height=height)
#    
#EOF
#        convert -density 300 $OUTDIR/HiC/1.SE_perChrom_mapping.pdf -quality 100 $OUTDIR/HiC/1.SE_perChrom_mapping.png
#        convert -density 300 $OUTDIR/HiC/2.mapping.pdf -quality 100 $OUTDIR/HiC/2.mapping.png
#        convert -density 300 $OUTDIR/HiC/3.Interactions.pdf -quality 100 $OUTDIR/HiC/3.Interactions.png
#        convert -density 300 $OUTDIR/HiC/4.ValidPairs.pdf -quality 100 $OUTDIR/HiC/4.ValidPairs.png
#        for i in `ls $OUTDIR/HiC/*.png | grep 'SE_perChrom_mapping\|mapping\|Interactions\|ValidPairs' |sort | sed 's/.png//g' |awk -F "/" '{print $NF}'`; do echo '<a href="'${i}.pdf'"><img src="'${i}.png'" width="1200"></a>' ; done | awk ' {print;} NR % 1 == 0 { print "<br>"; }' >$OUTDIR/HiC/Images.html
#        cat $OUTDIR/HiC/HiCsummary.html <(echo "") <(echo "<br>") $OUTDIR/HiC/Weblogoindex.html <(echo "") <(echo "<br>") $OUTDIR/HiC/Images.html >$OUTDIR/HiC/index.html
else
    echo 'No HiC samples found'
fi


######
echo
echo "Weblogos"
#######
##For raw sequence
if [[ ! `pwd` =~ ^\/vol\/mauranolab\/transposon\/aggregations\/ ]]; then
    echo "Copying raw weblogos"
    
    mkdir -p $OUTDIR/Weblogos_raw
    find ${dirs} -maxdepth 1 -name "*.R2.raw.png" | xargs --no-run-if-empty cp -t $OUTDIR/Weblogos_raw
    find ${dirs} -maxdepth 1 -name "*.R1.raw.png" | xargs --no-run-if-empty cp -t $OUTDIR/Weblogos_raw
    for i in `ls $OUTDIR/Weblogos_raw/*R1.raw.png | sort | awk -F "/" '{print $NF}'`; do echo '<a href="'${i}'"><img src="'${i}'" height="120"></a>' ; done  >$OUTDIR/Weblogos_raw/R1index.html
    for i in `ls $OUTDIR/Weblogos_raw/*R2.raw.png | sort | awk -F "/" '{print $NF}'`; do echo '<a href="'${i}'"><img src="'${i}'"  style= "position:absolute; LEFT:1200px"; height="120"></a>' ; done  >$OUTDIR/Weblogos_raw/R2index.html
    cat $OUTDIR/Weblogos_raw/R1index.html $OUTDIR/Weblogos_raw/R2index.html | sort| awk ' {print;} NR % 2 == 0 { print "<br>"; }'> $OUTDIR/Weblogos_raw/index.html
    
    
    #For processed sequence
    mkdir -p $OUTDIR/Weblogos_processed
    find ${dirs} -maxdepth 1 -name "*.BC.processed.png" | xargs --no-run-if-empty cp -t $OUTDIR/Weblogos_processed
    find ${dirs} -maxdepth 1 -name "*.plasmid.processed.png" | xargs --no-run-if-empty cp -t $OUTDIR/Weblogos_processed
    for i in `ls $OUTDIR/Weblogos_processed/*BC.processed.png | sort | awk -F "/" '{print $NF}'`; do echo '<a href="'${i}'"><img src="'${i}'" height="120"></a>' ; done  >$OUTDIR/Weblogos_processed/BCindex.html
    for i in `ls $OUTDIR/Weblogos_processed/*plasmid.processed.png | sort | awk -F "/" '{print $NF}'`; do echo '<a href="'${i}'"><img src="'${i}'" style= "position:absolute; LEFT:1200px"; height="120"></a>' ; done  >$OUTDIR/Weblogos_processed/plasmidindex.html
    cat $OUTDIR/Weblogos_processed/BCindex.html $OUTDIR/Weblogos_processed/plasmidindex.html | sort| awk ' {print;} NR % 2 == 0 { print "<br>"; }'> $OUTDIR/Weblogos_processed/index.html
    
    echo "weblogo for UMI sequence"
    mkdir -p $OUTDIR/Weblogos_UMI
    find ${dirs} -maxdepth 1 -name "*.UMIs.png" | xargs --no-run-if-empty cp -t $OUTDIR/Weblogos_UMI
    if find ${dirs} -maxdepth 1 -name "*.UMIs.png" | grep -q "." ; then
        for i in `ls $OUTDIR/Weblogos_UMI/*.UMIs.png | sort | awk -F "/" '{print $NF}'`; do echo '<a href="'${i}'"><img src="'${i}'" height="120"></a>' ; done | awk ' {print;} NR % 1 == 0 { print "<br>"; }'> $OUTDIR/Weblogos_UMI/index.html
    fi
fi

echo "weblogo for genomic reads"
mkdir -p $OUTDIR/Weblogos_genomic
find ${dirs} -maxdepth 1 -name "*.genomic.png" | xargs --no-run-if-empty cp -t $OUTDIR/Weblogos_genomic
for i in `ls $OUTDIR/Weblogos_genomic/*.genomic.png | sort | awk -F "/" '{print $NF}'`; do echo '<a href="'${i}'"><img src="'${i}'" height="120"></a>' ; done | awk ' {print;} NR % 1 == 0 { print "<br>"; }'> $OUTDIR/Weblogos_genomic/index.html

echo "weblogo for Barcode sequence"
mkdir -p $OUTDIR/Weblogos_Barcode
find ${dirs} -maxdepth 1 -name "*barcodes.png" | xargs --no-run-if-empty cp -t $OUTDIR/Weblogos_Barcode
for i in `ls $OUTDIR/Weblogos_Barcode/*barcodes.png | sort | awk -F "/" '{print $NF}'`; do echo '<a href="'${i}'"><img src="'${i}'" height="120"></a>' ; done | awk ' {print;} NR % 1 == 0 { print "<br>"; }'> $OUTDIR/Weblogos_Barcode/index.html


echo "weblogo for Cell Barcode sequence"
mkdir -p $OUTDIR/Weblogos_CellBC
find ${dirs} -maxdepth 1 -name "*cellBCs.png" | xargs --no-run-if-empty cp -t $OUTDIR/Weblogos_CellBC
for i in `ls $OUTDIR/Weblogos_CellBC/*cellBCs.png | sort | awk -F "/" '{print $NF}'`; do echo '<a href="'${i}'"><img src="'${i}'" height="120"></a>' ; done | awk ' {print;} NR % 1 == 0 { print "<br>"; }'> $OUTDIR/Weblogos_CellBC/index.html


#####
echo
echo "Summarize flowcell info"
#####
mkdir -p $OUTDIR/FlowcellSummary
${src}/FlowcellSummary.R $OUTDIR

#Convert to Excel
PYTHON_ARG="${fc}" INPUT_ARG="$OUTDIR/FlowcellSummary/summaryflowcell.tsv" OUTPUT_ARG="$OUTDIR/FlowcellSummary/${fc}.xlsx" python - <<END
import csv
import os
from xlsxwriter.workbook import Workbook
tsv_file = os.environ['INPUT_ARG']
xlsx_file = os.environ['OUTPUT_ARG']

# Create an XlsxWriter workbook object and add a worksheet.
workbook = Workbook(xlsx_file)
worksheet = workbook.add_worksheet()

# Create a TSV file reader.
tsv_reader = csv.reader(open(tsv_file, 'rt'), delimiter='\t')

# Read the row data from the TSV file and write it to the XLSX file.
for row, data in enumerate(tsv_reader):
    worksheet.write_row(row, 0, data)

# Close the XLSX file.
workbook.close()

END


echo '<a href="'${fc}.xlsx'"><font size="6">FlowcellSummary</font></a>'| cat $OUTDIR/FlowcellSummary/index.html - > $OUTDIR/FlowcellSummary/index1
mv $OUTDIR/FlowcellSummary/index1 $OUTDIR/FlowcellSummary/index.html


######
#Edit distance per sample
#######
#BUGBUG hardcoded number of columns, can't handle mixed runs with different number of columns
if [[ ! `pwd` =~ ^\/vol\/mauranolab\/transposon\/aggregations\/ ]]; then
    mkdir -p $OUTDIR/EditDist
    
    for i in `find ${dirs} -maxdepth 1 \( -name "extract*" -o -name "map*" \)`; do Dist=$(grep 'BC read edit distances' $i); echo $i $Dist; done | awk -v OFS='\t' -F'[' '{print $1, $3}' | sed 's/]//g' | perl -pe 's/, /\t/g' > $OUTDIR/EditDist/BC_EditDist
    for i in `find ${dirs} -maxdepth 1 \( -name "extract*" -o -name "map*" \)`; do Dist=$(grep 'Plasmid read edit distances' $i); echo $i $Dist; done | awk -v OFS='\t' -F'[' '{print $1, $3}' | sed 's/]//g' | perl -pe 's/, /\t/g' > $OUTDIR/EditDist/Plasmid_EditDist
    
    
    R --quiet --no-save << EOF
    source("${src}/Flowcell_Info.R")
    
    library(stringr)
    library(reshape2)
    BCEditDist<-read("$OUTDIR/EditDist/BC_EditDist",stringsAsFactors=F)
    BCEditDist[,1]<-gsub('\\\.o.*','',gsub('extract.|map.','',(gsub('.*/','',BCEditDist[,1]))))
    
    BCEditDist<-summaryBy(.~V1, data=BCEditDist,FUN=sum)
    colnames(BCEditDist)<-c('Sample',0:18)
    BCEditDist<-melt(BCEditDist)
    BCEditDist[,4]<-"Barcode"
    
    PlasmidEditDist<-read("$OUTDIR/EditDist/Plasmid_EditDist",stringsAsFactors=F)
    PlasmidEditDist[,1]<-gsub('\\\.o.*','',gsub('extract.|map.','',(gsub('.*/','',PlasmidEditDist[,1]))))
    
    PlasmidEditDist<-summaryBy(.~V1, data=PlasmidEditDist,FUN=sum)
    colnames(PlasmidEditDist)<-c('Sample',0:(ncol(PlasmidEditDist)-2))
    PlasmidEditDist<-melt(PlasmidEditDist)
    PlasmidEditDist[,4]<-'Plasmid'
    
    if (length(grep('iPCR\$',PlasmidEditDist[,1]))>0) {
        PlasmidEditDist[grep('iPCR\$',PlasmidEditDist[,1]),][,4] <-'Primer'
    }
    
    BCEditDist<-rbind(BCEditDist,PlasmidEditDist)
    colnames(BCEditDist)[1:4]<-c('Sample','Hamming','Reads','Read')
    BCEditDist[,2]<-as.numeric(as.character(BCEditDist[,2]))
    BCEditDist[,3]<-as.numeric(as.character(BCEditDist[,3]))
    
    colnames(BCEditDist)[1:4]<-c('Sample','Hamming','Reads','Read')
    
    
    myColors <- brewer.pal(5,"Set1")
    names(myColors)<-factor(c('RNA','iPCR','Plasmid','DNA','10xRNA'),levels=c('RNA','iPCR','Plasmid','DNA','10xRNA'))
    colScale <- scale_fill_manual(name = "Type",values = myColors)
    height <- ceiling(length(unique(BCEditDist[,"Sample"]))/3)*2
    
    #plot
    BCEditDist\$Type<-gsub(".*_",'',BCEditDist\$Sample)
    BCEditDist\$Sample<- substr(BCEditDist\$Sample, 0, min(nchar(BCEditDist\$Sample)))
    #pdf("$OUTDIR/EditDist/EditDistDistance_BC.pdf",height=10,width=20)
    bcl<-ggplot(BCEditDist[grep('Barcode',BCEditDist[,4]),], aes(Hamming,Reads,fill=Type)) + 
    geom_bar(stat="identity",color='black',size=0.25,alpha=0.4) +
    #ggtitle('GGlo')+
    facet_wrap(~Sample+Read,scales = "free_y", ncol=3)+ 
    theme_classic()+
    colScale+
    xlab('Hamming distance')+
    theme(legend.position = "bottom",legend.text=element_text(size=8))+
    theme(axis.title.x = element_text(size=14),axis.text.x  = element_text(size=12),axis.title.y=element_text(size=14),axis.text.y=element_text(size=12))+
    #scale_x_continuous(label= fancy_scientific)+
    scale_y_continuous(label= fancy_scientific)
    #dev.off()
    pdf(NULL)
    gp = ggplotGrob(bcl)
    gp = editGrob(grid.force(gp), gPath("GRID.stripGrob", "GRID.text"), grep = TRUE, global = TRUE, just = "left", x = unit(0.1,"npc"))
    ggsave(gp ,file="$OUTDIR/EditDist/EditDistDistance_BC.pdf", width=13, height=height)
    
    
    #plot
    #pdf("$OUTDIR/EditDist/EditDistDistance_Plasmid.pdf",height=10,width=20)
    pll<-ggplot(BCEditDist[grep('Plasmid|Primer',BCEditDist[,4]),], aes(Hamming,Reads,fill=Type)) + 
    geom_bar(stat="identity",color="black", size=0.1,alpha=0.4) +
    #ggtitle('GGlo')+
    facet_wrap(~Sample+Read,scales = "free_y", ncol=3)+ 
    theme_classic()+
    colScale+
    xlab('Hamming distance')+
    #scale_fill_brewer('Set1')+
    theme(legend.position = "bottom",legend.text=element_text(size=8))+
    theme(axis.title.x = element_text(size=14),axis.text.x  = element_text(size=12),axis.title.y=element_text(size=14),axis.text.y=element_text(size=12))+
    #scale_x_continuous(label= fancy_scientific)+
    scale_y_continuous(label= fancy_scientific)
    #dev.off()
    pdf(NULL)
    gp = ggplotGrob(pll)
    gp = editGrob(grid.force(gp), gPath("GRID.stripGrob", "GRID.text"), grep = TRUE, global = TRUE, just = "left", x = unit(0.1,"npc"))
    ggsave(gp ,file="$OUTDIR/EditDist/EditDistDistance_Plasmid.pdf", width=13, height=height)
    
EOF
    
    convert -density 300 $OUTDIR/EditDist/EditDistDistance_BC.pdf -quality 100 $OUTDIR/EditDist/EditDistDistance_BC.png
    convert -density 300 $OUTDIR/EditDist/EditDistDistance_Plasmid.pdf -quality 100 $OUTDIR/EditDist/EditDistDistance_Plasmid.png
    
    for i in `ls $OUTDIR/EditDist/*.png | sort | sed 's/.png//g' | awk -F "/" '{print $NF}'`; do echo '<a href="'${i}.pdf'"><img src="'${i}.png'" width="1200"></a>' ; done  >$OUTDIR/EditDist/index.html
fi


#####
#Barcode frequency
#####
mkdir -p $OUTDIR/BarcodeFreq
for i in `find ${dirs} -maxdepth 1 -name *barcode.counts.txt | sed 's/^..//g' | sed 's/\/.*//g'`; do reads=$(cut -f2 ${i}/${i}.barcode.counts.txt | awk '{sum+=$1} END {print sum}'); sort -nk2 ${i}/${i}.barcode.counts.txt| awk -v OFS='\t' -F'\t' -v sample="$i" -v reads="$reads" '{print $1, $2, sample, reads}'; done > $OUTDIR/BarcodeFreq/BarcodeFreq.txt


#R --quiet --no-save << EOF
#library(dplyr)
#BC <- read("$OUTDIR/BarcodeFreq/BarcodeFreq.txt")
#colnames(BC) <- c("Barcodes","barcodeFreq", "Sample", 'TotalBCcount')
#BC[,"Type"]<-gsub(".*_",'',gsub('_Merged','',BC[,"Sample"]))
#BC[,"cpm"] <- BC[,"barcodeFreq"]/(BC[,"TotalBCcount"]/1e6)
#BC[,"BC.bin"] <- cut(BC[,"cpm"], breaks=c( 0, 1.1, 5.1, 10.1, 50.1, 100.1, 1000.1, 10000.1, 100000.1, 1000000.1), right=F, include.lowest=T, labels=c("0-1", "1-5", "5-10", "10-50", "50-100", "100-1000", "100-10000", "10000-100000","+100000"))
##BC <- BC %>%
##    group_by(Sample, cpm) %>%
##    mutate(quantile = ntile(cpm, 10))
##quantBCs <- lsit()
##for (i in 1:length(unique(BC\$Sample))) {
##    quantile(BC[BC\$Sample %in% 'BS493A-MSH_K562~A1~GGlo~HS2~A1~1Linear~ExoI_Rep2_DNA',]\$cpm, prob = seq(0, 1, length = 11), type = 5
##
#myColors <- brewer.pal(5,"Set1")
#names(myColors) <- factor(c('RNA','iPCR','Plasmid','DNA','10xRNA'),levels=c('RNA','iPCR','Plasmid','DNA','10xRNA'))
#colScale <- scale_fill_manual(name = "Type",values = myColors)
#BC[,"Sample"]<- substr(BC[,"Sample"], 0, min(nchar(BC[,"Sample"])))
#
#height <- ceiling(length(unique(BC[,"Sample"]))/3)*2
#
#ggHist <- BC %>%
#    ggplot(aes(x=log10(cpm+1),fill=Type)) +
#    geom_histogram(binwidth=0.1, colour='black', alpha=0.4) +
#    theme_classic()+
#    colScale+
#    facet_wrap(~Sample,scales="free_y", ncol=3)+ 
#    theme_classic()+
#    coord_cartesian(xlim = c(0, 4))+
#    theme(axis.title.x = element_text(size=14),axis.text.x  = element_text(size=12),axis.title.y=element_text(size=14),axis.text.y=element_text(size=12))
#
#pdf(NULL)
#gp = ggplotGrob(ggHist)
#gp = editGrob(grid.force(gp), gPath("GRID.stripGrob", "GRID.text"), grep = TRUE, global = TRUE, just = "left", x = unit(0.1,"npc"))
#ggsave(gp ,file="$OUTDIR/BarcodeFreq/BarcodeHist.pdf", width=13, height=height)
#
#
###ggsave(t,file="$OUTDIR/BarcodeFreq/BarcodeFreq.pdf", width=20, height=10)
##t<-ggplot(BC, aes(x=BC.bin,fill=Type)) +
##geom_bar(color="black", size=0.25) +
##theme_classic()+
##colScale+
##xlab('BC bin (cpm)') +
##facet_wrap(~Sample,scales="free_y")+ 
##theme_classic()+
##theme(axis.title.x = element_text(size=14),axis.text.x  = element_text(size=12,angle=60,hjust=1),axis.title.y=element_text(size=14),axis.text.y=element_text(size=12))
###ggsave(t,file="$OUTDIR/BarcodeFreq/BarcodeFreq.pdf", width=20, height=10)
##
##pdf(NULL)
##gp = ggplotGrob(t)
##gp = editGrob(grid.force(gp), gPath("GRID.stripGrob", "GRID.text"), grep = TRUE, global = TRUE, just = "left", x = unit(0.1,"npc"))
##ggsave(gp ,file="$OUTDIR/BarcodeFreq/BarcodeFreq.pdf", width=20, height=10)
##
#EOF
#
#rm $OUTDIR/BarcodeFreq/BarcodeFreq.txt
#for i in `ls  $OUTDIR/BarcodeFreq | grep pdf | sed 's/.pdf//g'`; do convert -density 300 $OUTDIR/BarcodeFreq//${i}.pdf -quality 100 $OUTDIR/BarcodeFreq/${i}.png; done
#
#for i in `ls $OUTDIR/BarcodeFreq/*.png | sort | sed 's/.png//g' | awk -F "/" '{print $NF}'`; do echo '<a href="'${i}.pdf'"><img src="'${i}.png'" width="1200"></a> <br>' ; done  >$OUTDIR/BarcodeFreq/index.html


#####
#Barcode frequency with UMI
#####
if [[ `find ${dirs} -maxdepth 1 -name "*.barcode.counts.withUMI.txt.gz" | sed 's/^..//g' | sed 's/\/.*//g' | wc -l` -ge 1 ]]; then
    mkdir -p $OUTDIR/BarcodeFreq_UMI
    for i in `find ${dirs} -maxdepth 1 -name "*.barcode.counts.withUMI.txt.gz" | sed 's/^..//g' | sed 's/\/.*//g'`; do uniqueMol=$(zcat -f ${i}/${i}.barcode.counts.withUMI.txt.gz | wc -l); zcat -f ${i}/${i}.barcode.counts.withUMI.txt.gz | cut -f1 | sort | uniq -c | sort -nk1| awk 'BEGIN {OFS="\t"} {print $2, $1}' | awk -v OFS='\t' -F'\t' -v sample="$i" -v uniqMol="$uniqueMol" '{print $1, $2, sample, uniqMol}'; done > $OUTDIR/BarcodeFreq_UMI/BarcodeFreq_UMI.txt
    
    R --quiet --no-save << EOF
    BC <- read("$OUTDIR/BarcodeFreq_UMI/BarcodeFreq_UMI.txt")
    colnames(BC) <- c("Barcodes","BarcodeFreq_UMI", "Sample", "TotalUMIcount")
    BC[,"Type"]<-gsub(".*_",'',gsub('_Merged','',BC[,"Sample"]))
    BC[,"cpm"] <- BC[,"BarcodeFreq_UMI"]/(BC[,"TotalUMIcount"]/1e6)
    BC[,"BC.bin"] <- cut(BC[,"cpm"], breaks=c( 0, 1.1, 5.1, 10.1, 50.1, 100.1, 1000.1, 10000.1, 100000.1, 1000000.1), right=F, include.lowest=T, labels=c("0-1", "1-5", "5-10", "10-50", "50-100", "100-1000", "100-10000", "10000-100000","+100000"))
    #BC <- BC %>%
      # 
      # BC\$BC.bin <- cut(BC\$BarcodeFreq_UMI, breaks=c( 0, 1.1, 5.1, 10.1, 50.1, 100.1, 1000.1, 10000.1, 100000.1, 1000000.1), right=F, include.lowest=T, labels=c("0-1", "1-5", "5-10", "10-50", "50-100", "100-1000", "100-10000", "10000-100000","+100000"))
      # BC\$Type<-gsub(".*_",'',gsub('_Merged','',BC\$Sample))
    
    myColors <- brewer.pal(5,"Set1")
    names(myColors)<-factor(c('RNA','iPCR','Plasmid','DNA','10xRNA'),levels=c('RNA','iPCR','Plasmid','DNA','10xRNA'))
    colScale <- scale_fill_manual(name = "Type",values = myColors)
    BC\$Sample<- substr(BC\$Sample, 0, min(nchar(BC\$Sample)))
    
    height <- ceiling(length(unique(BC[,"Sample"]))/3)*2
    
    ggHist <- ggplot(BC, aes(x=log10(cpm+1),fill=Type)) +
    geom_histogram(binwidth=0.1, colour='black', alpha=0.4) +
    theme_classic()+
    colScale+
    facet_wrap(~Sample,scales="free_y", ncol=3)+ 
    theme_classic()+
    coord_cartesian(xlim = c(0, 4))+
    theme(axis.title.x = element_text(size=14),axis.text.x  = element_text(size=12),axis.title.y=element_text(size=14),axis.text.y=element_text(size=12))
    
    pdf(NULL)
    gp = ggplotGrob(ggHist)
    gp = editGrob(grid.force(gp), gPath("GRID.stripGrob", "GRID.text"), grep = TRUE, global = TRUE, just = "left", x = unit(0.1,"npc"))
    ggsave(gp ,file="$OUTDIR/BarcodeFreq_UMI/BarcodeFreq_UMI_hist.pdf", width=13, height=height)
    
    #t<-ggplot(BC, aes(x=BC.bin,fill=Type)) +
    #geom_bar(color="black", size=0.25,alpha=0.4) +
    #theme_classic()+
    #colScale+
    #facet_wrap(~Sample,scales="free_y")+ 
    #xlab('Number of barcodes') +
    #theme_classic()+
    #theme(axis.title.x = element_text(size=14),axis.text.x  = element_text(size=12,angle=60,hjust=1),axis.title.y=element_text(size=14),axis.text.y=element_text(size=12))
    ##ggsave(t,file="$OUTDIR/BarcodeFreq_UMI/BarcodeFreq_UMI.pdf", width=20, height=10)
    #
    #pdf(NULL)
    #gp = ggplotGrob(t)
    #gp = editGrob(grid.force(gp), gPath("GRID.stripGrob", "GRID.text"), grep = TRUE, global = TRUE, just = "left", x = unit(0.1,"npc"))
    #ggsave(gp ,file="$OUTDIR/BarcodeFreq_UMI/BarcodeFreq_UMI.pdf", width=20, height=10)
    
    
EOF
    
    rm $OUTDIR/BarcodeFreq_UMI/BarcodeFreq_UMI.txt
      # convert -density 300 $OUTDIR/BarcodeFreq_UMI/BarcodeFreq_UMI.pdf -quality 100 $OUTDIR/BarcodeFreq_UMI/BarcodeFreq_UMI.png
    for i in `ls  $OUTDIR/BarcodeFreq_UMI | grep pdf | sed 's/.pdf//g'`; do convert -density 300 $OUTDIR/BarcodeFreq_UMI//${i}.pdf -quality 100 $OUTDIR/BarcodeFreq_UMI/${i}.png; done

    for i in `ls $OUTDIR/BarcodeFreq_UMI/*.png | sort | sed 's/.png//g' | awk -F "/" '{print $NF}'`; do echo '<a href="'${i}.pdf'"><img src="'${i}.png'" width="1200"></a>' ; done  >$OUTDIR/BarcodeFreq_UMI/index.html
    
    
    mkdir -p $OUTDIR/UMI_distribution
    for i in `find ${dirs} -maxdepth 1 -name "*.barcode.counts.withUMI.txt.gz" | sed 's/^..//g' | sed 's/\/.*//g'`; do  zcat -f $i/$i.barcode.counts.withUMI.txt.gz | awk -v OFS='\t' -F'\t' -v sample="$i" '{print $1, $2, $3, sample}' ; done > $OUTDIR/UMI_distribution/UMI_distribution.txt


    R --quiet --no-save << EOF
    
    library(data.table)
    
    maxUMI <- fread("$OUTDIR/UMI_distribution/UMI_distribution.txt",header=F,stringsAsFactors=F)
    maxUMI <- as.data.frame(maxUMI)
    #maxUMI <- fread("/home/maagj01/public_html/blog/flowcellInfo/FCHVTCCBGX2/UMI_distribution/UMI_distribution.txt",header=F,stringsAsFactors=F)
    height <- ceiling(length(unique(maxUMI[,"V4"]))/3)*2
    cat(height, '\n')
    maxUMI\$Type<-gsub(".*_",'',gsub('_Merged','',maxUMI\$V4))
    
    myColors <- brewer.pal(5,"Set1")
    names(myColors)<-factor(c('RNA','iPCR','Plasmid','DNA','10xRNA'),levels=c('RNA','iPCR','Plasmid','DNA','10xRNA'))
    colScale <- scale_fill_manual(name = "Type",values = myColors)
    maxUMI\$V4<- substr(maxUMI\$V4, 0, min(nchar(maxUMI\$V4)))
    library(dplyr)
    
    ggmaxUMI <- maxUMI %>%
    select(V1,V3, V4, Type) %>%
    group_by(V1,V4, Type) %>%
    summarise_each(funs(max)) %>%
    ggplot(aes(x=log10(V3+1), fill=Type)) +
    geom_histogram(color="black", size=0.25,alpha=0.4)+
    facet_wrap(~V4, scale='free', ncol=3)+
    theme_classic()+
    colScale+
    xlab('log10(Max UMI copies per barcode)')+
    ylab('Counts')+
    theme(axis.title.x = element_text(size=14),axis.text.x  = element_text(size=12,angle=60,hjust=1),axis.title.y=element_text(size=14),axis.text.y=element_text(size=12))
    
    pdf(NULL)
    gp = ggplotGrob(ggmaxUMI)
    gp = editGrob(grid.force(gp), gPath("GRID.stripGrob", "GRID.text"), grep = TRUE, global = TRUE, just = "left", x = unit(0.1,"npc"))
    ggsave(gp ,file="$OUTDIR/UMI_distribution/UMI_distribution.pdf", width=13, height=height)

EOF

    rm $OUTDIR/UMI_distribution/UMI_distribution.txt
    convert -density 300 $OUTDIR/UMI_distribution/UMI_distribution.pdf -quality 100 $OUTDIR/UMI_distribution/UMI_distribution.png
    
    for i in `ls $OUTDIR/UMI_distribution/*.png | sort | sed 's/.png//g' | awk -F "/" '{print $NF}'`; do echo '<a href="'${i}.pdf'"><img src="'${i}.png'" width="1200"></a>' ; done  >$OUTDIR/UMI_distribution/index.html
fi


#####
#Saturation curve
#####
if [[ `find ${dirs} -maxdepth 1 -name "*Saturation*.txt" | sed 's/^..//g' | sed 's/\/.*//g' | wc -l` -ge 1 ]]; then
    mkdir -p $OUTDIR/SaturationCurve
    for i in `find ${dirs} -maxdepth 1 -name "*Saturation*.txt" | sed 's/^..//g' | sed 's/\/.*//g'`; do cat ${i}/${i}.Saturation*| awk -v OFS='\t' -F'\t' -v sample="$i" '{print $1, $2, sample}'; done | perl -pe 's/ /\t/g' > $OUTDIR/SaturationCurve/SaturationCurve.txt
    
    rm -f $OUTDIR/SaturationCurve/*pdf
    rm -f $OUTDIR/SaturationCurve/*png
    R --quiet --no-save << EOF
    source("${src}/Flowcell_Info.R")
    
    SaturationCurve <- read("$OUTDIR/SaturationCurve/SaturationCurve.txt")
    
    #SaturationCurve<-SaturationCurve[!is.na(SaturationCurve$\Unique_BC),]
    SaturationCurve[,2]<-as.numeric(SaturationCurve[,2])
    SaturationCurve[,1]<-as.numeric(SaturationCurve[,1])
    SaturationCurve[,3]<-gsub('minreads','',SaturationCurve[,3])
    #SaturationCurve\$Type<-gsub(".*_",'',SaturationCurve[,5])
    SaturationCurve\$Type<-gsub(".*_",'',gsub('_Merged','',SaturationCurve[,5]))
    
    colnames(SaturationCurve)[1:6] <- c("Reads","Unique_BC","minReads","NA","Sample",'Type')
    
    SaturationCurve\$Sample<- substr(SaturationCurve\$Sample, 0, min(nchar(SaturationCurve\$Sample)))
    myColors <- brewer.pal(5,"Set1")
    names(myColors) <- factor(c('RNA','iPCR','Plasmid','DNA','10xRNA'),levels=c('RNA','iPCR','Plasmid','DNA','10xRNA'))
    colScale <- scale_colour_manual(name = "Type", values = myColors)
    height <- ceiling(length(unique(SaturationCurve[,"Sample"]))/3)*2
    
    s1 <- ggplot(SaturationCurve, aes(Reads, Unique_BC,color=Type)) + 
    geom_line(aes(linetype=minReads)) +
    #geom_line()+
    #ggtitle('GGlo')+
    facet_wrap(~Sample,scales = "free", ncol=3)+ 
    theme_classic()+
    colScale+
    theme(legend.position = "bottom",legend.text=element_text(size=8))+
    theme(axis.title.x = element_text(size=14),axis.text.x  = element_text(size=12),axis.title.y=element_text(size=14),axis.text.y=element_text(size=12))+
    scale_x_continuous(label= fancy_scientific)+
    scale_y_continuous(label= fancy_scientific)
    
    pdf(NULL)
    gp = ggplotGrob(s1)
    gp = editGrob(grid.force(gp), gPath("GRID.stripGrob", "GRID.text"), grep = TRUE, global = TRUE, just = "left", x = unit(0.1,"npc"))
    ggsave(gp ,file="$OUTDIR/SaturationCurve/SaturationCurve.pdf", width=13, height=height)
EOF
    
    convert -density 300 $OUTDIR/SaturationCurve/SaturationCurve.pdf -quality 100 $OUTDIR/SaturationCurve/SaturationCurve.png
    
    for i in `ls $OUTDIR/SaturationCurve/*.png | sort | sed 's/.png//g' | awk -F "/" '{print $NF}'`; do echo '<a href="'${i}.pdf'"><img src="'${i}.png'" width="1200"></a>' ; done  >$OUTDIR/SaturationCurve/index.html
fi


######
#iPCR specific 
######
if [[ `find ${dirs} -maxdepth 1 -type d | grep BS | grep iPCR| wc -l` -ge 1 ]]; then
    echo 'iPCR samples found'
    mkdir -p $OUTDIR/iPCR/
    if [[ `find ${dirs} -maxdepth 1 -name *DistDpn.bed | sed 's/^..//g' | sed 's/\/.*//g'| wc -l` -ge 1 ]]; then
        for iPCR in `find ${dirs} -maxdepth 1 -name *DistDpn.bed | sed 's/^..//g' | sed 's/\/.*//g' `; do awk -v OFS='\t' -F'\t' -v sample="$iPCR" '{print $1, $2, $3, $4, $5, sample}' ${iPCR}/DistDpn.bed; done > $OUTDIR/iPCR/DistDpn.bed
    fi
    
    if [[ `find ${dirs} -maxdepth 1 -name *DistToTSS.txt | sed 's/^..//g' | sed 's/\/.*//g' | wc -l` -ge 1 ]]; then
        for iPCR in `find ${dirs} -maxdepth 1 -name *DistToTSS.txt | sed 's/^..//g' | sed 's/\/.*//g' `; do awk -v OFS='\t' -F'\t' -v sample="$iPCR" '{print $1, sample}' ${iPCR}/DistToTSS.txt; done >$OUTDIR/iPCR/DistToTSS.bed
    fi
    
    if [[ `find ${dirs} -maxdepth 1 -name *ObsvsExp.txt | sed 's/^..//g' | sed 's/\/.*//g' | wc -l` -ge 1 ]]; then
        for iPCR in `find ${dirs} -maxdepth 1 -name *ObsvsExp.txt | sed 's/^..//g' | sed 's/\/.*//g' `; do awk -v OFS='\t' -F'\t' -v sample="$iPCR" '{print sample, $0}' ${iPCR}/${iPCR}.ObsvsExp.txt; done >$OUTDIR/iPCR/ObsvsExp.txt
    fi
    
    if [[ `find ${dirs} -maxdepth 1 -name *DistToDNase.bed | sed 's/^..//g' | sed 's/\/.*//g' | wc -l` -ge 1 ]]; then
        for iPCR in `find ${dirs} -maxdepth 1 -name *DistToDNase.bed | sed 's/^..//g' | sed 's/\/.*//g' `; do awk -v OFS='\t' -F'\t' -v sample="$iPCR"  '{print $1, $2, $3, $4, sample}' ${iPCR}/DistToDNase.bed; done >$OUTDIR/iPCR/DistToDNase.bed
    fi
    
    if [[ `find ${dirs} -maxdepth 1 -name *DistMsp1.bed | sed 's/^..//g' | sed 's/\/.*//g' | wc -l` -ge 1 ]]; then
        for iPCR in `find ${dirs} -maxdepth 1 -name *DistMsp1.bed | sed 's/^..//g' | sed 's/\/.*//g' `; do awk -v OFS='\t' -F'\t' -v sample="$iPCR" '{print $1, $2, $3, $4, $5, sample}' ${iPCR}/DistMsp1.bed; done > $OUTDIR/iPCR/DistMsp1.bed
    fi
    
    
    
    R --quiet --no-save << EOF
    library(stringr)
    library(reshape2)
    myColors <- brewer.pal(5,"Set1")
    names(myColors)<-factor(c('RNA','iPCR','Plasmid','DNA','10xRNA'),levels=c('RNA','iPCR','Plasmid','DNA','10xRNA'))
    colScale <- scale_fill_manual(name = "Type",values = myColors)
    
    
    Dpn<-read("$OUTDIR/iPCR/DistDpn.bed",stringsAsFactors=F)
    colnames(Dpn)[1:6]<-c('chr','start','end','dpnsite','Distance','Sample')
    Dpn[,5]<-as.numeric(as.character(Dpn[,5]))
    Dpn\$Type<-gsub(".*_",'',gsub('_Merged','',Dpn\$Sample))
    height <- ceiling(length(unique(Dpn[,"Sample"]))/3)*2
    
    p<-ggplot(Dpn,aes(x=Sample,y=Distance,fill=Type))+
    geom_boxplot(outlier.shape = NA,alpha=0.4)+
    #scale_fill_brewer(palette ='Set1')+
    theme_classic()+
    theme(axis.title.x = element_text(size=14),axis.text.x  = element_text(size=12,angle=90,hjust=1),axis.title.y=element_text(size=14),axis.text.y=element_text(size=12))+
    ylab('Distance to Dpn sites')+
    xlab('Samples')+
    colScale+
    coord_flip(ylim = c(0, 750))
    ggsave(p, file="$OUTDIR/iPCR/DistDpn.pdf", width=13, height=10) 
    
    
    DNase<-read("$OUTDIR/iPCR/DistToDNase.bed",stringsAsFactors=F)
    colnames(DNase)[1:5]<-c('chr','start','end','Distance','Sample')
    DNase[,4]<-as.numeric(as.character(DNase[,4]))
    DNase\$Type<-gsub(".*_",'',gsub('_Merged','',DNase\$Sample))
    
    d<-ggplot(DNase, aes(x=abs(Distance+1),fill=Type)) +
    geom_histogram(color="black", size=0.25,alpha=0.4)+ 
    theme_classic()+
    facet_wrap(~Sample,scales = "free_y", ncol=3)+ 
    theme(axis.title.x = element_text(size=14),axis.text.x  = element_text(size=12),axis.title.y=element_text(size=14),axis.text.y=element_text(size=12))+
    #xlim(c(0,1500))+
    xlab('Distance to DNase site')+
    colScale+
    scale_x_log10(  breaks = scales::trans_breaks("log10", function(x) 10^x),  labels = scales::trans_format("log10", scales::math_format(10^.x)))+
    annotation_logticks(sides = "b",long = unit(0.2, "cm"))    
    ggsave(d, file="$OUTDIR/iPCR/DistDNase.pdf", width=13, height=height) 
    
    
    TSS <- read("$OUTDIR/iPCR/DistToTSS.bed")
    colnames(TSS) <- c("DistToTSS", "Sample")
    #TSS\$DistToTSS.bin <- cut(TSS\$DistToTSS, breaks=c(-10^7, -10^6, -10^5, -10^4, -2500, 0, 5000, 10^4, 10^5, 10^6, 10^7), right=F, include.lowest=T, labels=c("-10 Mb", "-1 Mb", "-100 kb", "-10 kb", "-2.5 kb", "+5 kb", "+10 kb", "+100 kb", "+1 Mb", "+10 Mb"))
    
    TSS\$Type<-gsub(".*_",'',gsub('_Merged','',TSS\$Sample))
    t<-ggplot(TSS, aes(x=abs(DistToTSS)+1,fill=Type)) +
    #geom_bar(color="black", size=0.25,alpha=0.4) +
    geom_histogram(binwidth=0.1, colour="black", size=0.25, alpha=0.4) + 
    theme_classic()+
    colScale+
    facet_wrap(~Sample,scales="free_y", ncol=3)+ 
    theme_classic()+
    theme(axis.title.x = element_text(size=14),axis.text.x  = element_text(size=12),axis.title.y=element_text(size=14),axis.text.y=element_text(size=12)) +
    scale_x_log10(  breaks = scales::trans_breaks("log10", function(x) 10^x),  labels = scales::trans_format("log10", scales::math_format(10^.x)))+
    annotation_logticks(sides = "b",long = unit(0.2, "cm"))    
    ggsave(t,file="$OUTDIR/iPCR/DistToTSS.pdf", width=13, height=height)
    
    
    # 
    # Msp1<-read("$OUTDIR/iPCR/DistMsp1.bed",stringsAsFactors=F)
    # colnames(Msp1)[1:6]<-c('chr','start','end','Msp1site','Distance','Sample')
    # Msp1[,5]<-as.numeric(as.character(Msp1[,5]))
    # Msp1[,"Type"]<-gsub(".*_",'',gsub('_Merged','',Msp1[,"Sample"]))
    #
    # p<-ggplot(Msp1,aes(x=Sample,y=Distance,fill=Type))+
    # geom_boxplot(outlier.shape = NA,alpha=0.4)+
    # #scale_fill_brewer(palette ='Set1')+
    # theme_classic()+
    # theme(axis.title.x = element_text(size=14),axis.text.x  = element_text(size=12,angle=90,hjust=1),axis.title.y=element_text(size=14),axis.text.y=element_text(size=12))+
    # ylab('Distance to Msp1 sites')+
    # xlab('Samples')+
    # colScale+
    # coord_flip(ylim = c(0, 1500))
    # ggsave(p, file="$OUTDIR/iPCR/DistMsp1.pdf", width=13, height=10) 
    # 
    # ObsvsExp<-read('$OUTDIR/iPCR/ObsvsExp.txt',header=T)
    # colnames(ObsvsExp)[1]<-c('Sample')
    # ObsvsExp<-ObsvsExp[grep('l2fold',ObsvsExp\$l2fold,invert=T),]
    # ObsvsExp\$l2fold<-as.numeric(ObsvsExp\$l2fold)
    # ObsvsExp\$significant<-as.numeric(ObsvsExp\$qvalue)<0.05
    # 
    # 
    # myColors <- c('#a6cee3','#e41a1c')
    # names(myColors)<-factor(c('FALSE','TRUE'),levels=c('FALSE','TRUE'))
    # colScale <- scale_fill_manual(name = "significant",values = myColors)
    #
    # 
    # o<-ggplot(ObsvsExp, aes(x=annotation,y=l2fold,fill=significant)) +
    # geom_bar(stat='identity',color="black", size=0.25,alpha=0.4)+ 
    # theme_classic()+
    # colScale+
    # #scale_fill_manual(values=c('#a6cee3','#e41a1c'))+
    # theme(axis.title.x = element_text(size=14),axis.text.x  = element_text(size=12,angle=90,hjust=1),axis.title.y=element_text(size=14),axis.text.y=element_text(size=12))+
    # facet_wrap(~Sample, ncol=3)+ 
    # ggtitle("Obs/Exp")+
    # ylab('log2FC(Observed/Expected)')
    # ggsave(o,file="$OUTDIR/iPCR/ObsvsExp.pdf", width=13, height=height)

    
    
EOF
    for i in `ls  $OUTDIR/iPCR | grep pdf | sed 's/.pdf//g'`; do convert -density 300 $OUTDIR/iPCR//${i}.pdf -quality 100 $OUTDIR/iPCR/${i}.png; done

    #convert -density 300 $OUTDIR/iPCR/DistDpn.pdf -quality 100 $OUTDIR/iPCR/DistDpn.png
    #convert -density 300 $OUTDIR/iPCR/DistDNase.pdf -quality 100 $OUTDIR/iPCR/DistDNase.png
    #convert -density 300 $OUTDIR/iPCR/DistToTSS.pdf -quality 100 $OUTDIR/iPCR/DistToTSS.png
    #convert -density 300 $OUTDIR/iPCR/ObsvsExp.pdf -quality 100 $OUTDIR/iPCR/ObsvsExp.png
    #convert -density 300 $OUTDIR/iPCR/DistMsp1.pdf -quality 100 $OUTDIR/iPCR/DistMsp1.png
    for i in `ls $OUTDIR/iPCR/*.png | sort | sed 's/.png//g' | awk -F "/" '{print $NF}'`; do echo '<a href="'${i}.pdf'"><img src="'${i}.png'" width="1200"></a>' ; done | awk ' {print;} NR % 1 == 0 { print "<br>"; }' >$OUTDIR/iPCR/index.html
else 
    echo 'No iPCR samples found'     
fi



echo
echo 'Done!!!'
date
