library(reshape)
library(ggplot2)
library(stringr)
library(RColorBrewer)
#for i in {10,50,100,500,1000,5000,10000,50000,100000}; do echo head -${i} BS00067A-5xBGlo_K562d4_2hDpn_iPCR.barcodes.txt| py -m cProfile DeDup_adjacencyOutput.py /dev/fd/0 -o testTime >Profiling/first_${i}_time.txt;done

#Read all cProfile files
timeProfile<-list.files('./',recursive=T)
timeProfile<-timeProfile[grep('complex',timeProfile)]
timeProfile<-timeProfile[file.info(list.files(path="./", full.names=TRUE))$size>1000]
timeProfile<-timeProfile[grep('Edit2_|EditLoop_|EditSwap_|EditLoopList2_|ListCountEditEqual_|EditLoopAppend_|ListCompReverse_|ListComp_',timeProfile)]
#Skip the file not yet completed

#Functions of interest
DeDup<-c('.py:2(','dedup_dir_adj','edit_distance',
'get_adj_list_directional_adjacency','levenshtein','listcomp','dictcomp',"method 'encode" )

DeDup<-c('.py:2(')
#Create list with all variables
timeProfile<-expand.grid(DeDup,timeProfile)

#Create list for results
cumTime<-as.data.frame(matrix(0, ncol = 5, nrow = nrow(timeProfile)))
colnames(cumTime)<-c('CumulativeTime','nrows','Function','NumberCalls','Complexity')


for (file in 1:nrow(timeProfile)){cat(as.character(timeProfile$Var1[file]),' ',as.character(timeProfile$Var2[file]),'\n')
	Time <-readLines(paste0('./',timeProfile$Var2[file]),skip=5)[-c(1:5)]
	#sub("^\\s+","",readLines(paste0('./',timeProfile$Var2[file]),skip=5)[1])
	#Time<-Time[grep(timeProfile$Var1[file],Time,,fixed=TRUE)]
	#Time<-head(Time,1)
	Time<-sub("^\\s+","",Time)
	Time<-gsub("\\s+", " ", str_trim(Time))
	Time<-Time[c(1:(length(Time)-7))]
	Time<-Time[grep('Shadow.py:2\\(',Time,invert=T)]
	Time<-Time[grep('.py:2\\(',Time)]
	#cat(Time)
	df1<-data.frame(do.call('rbind', strsplit(as.character(Time),' ',fixed=TRUE)),stringsAsFactors=F)
	##df1<-data.frame(do.call('rbind',as.data.frame(t(unlist(str_split(Time," ",n=6))))),stringsAsFactors=F)
	colnames(df1)<-c('ncalls','tottime','percall','cumtime','percall','function1')
	cumTime[,1][file]<-as.numeric(df1$cumtime)
	cumTime[,2][file]<-as.numeric(gsub('High_complex_|Low_complex_|EditSwap_|ListCountEditEqual_|EditLoopList2_|EditLoopList1_|EditLoopAppend_|EditLoopList|EditLoop_|ListCompReverse_|ListComp_|Edit2_|Cython_|\\..*','',timeProfile$Var2[file]))
	cumTime[,3][file]<-as.character(df1$function1)
	cumTime[,4][file]<-as.numeric(as.character(df1$ncalls))
	cumTime[,5][file]<-as.character(gsub('EditSwap_|EditLoopList2_|EditLoopList1_|ListCountEditEqual_|EditLoopList|EditLoop_|EditLoopAppend_|ListComp_|ListCompReverse_|Edit2_|Cython_|_.*','',timeProfile$Var2[file]))
	cumTime[,5][file]
	#cumTime[,6][file]<-as.character(gsub('EditSwap_|EditLoopList2_|EditLoopList1_|EditLoopList|EditLoop_|EditLoopEncode_|Edit2_|Cython_|_.*','',timeProfile$Var2[file]))
}

#cumTime<-cumTime[grep("DeDup_adjacencyOutput_edit_distance_Swap.py:2(<module>)",cumTime$Function,fixed=TRUE,invert=T),]
ReNameFunction<-as.data.frame(c('EdistCount_encodeBoth_OLDloop','EqualCountEdist_PreEncode_append_innerloop','EqualCountEdist_encodeUmi2_append_inneerloop',
'CountsEdist_encodeBoth_OLDloop','EqualEdistCount_PreEncoded_listComp','CountEdistEqual_PreEncoded_listComp_countsAremadeTolist'),
c("DeDup_adjacencyOutput_edit_distance.py:2(<module>)", "DeDup_adjacencyOutput_edit_distance_LoopAppend.py:2(<module>)", 
"DeDup_adjacencyOutput_edit_distance_Loop_Listcomp.py:2(<module>)", 
"DeDup_adjacencyOutput_edit_distance_Swap.py:2(<module>)", "DeDup_adjacencyOutput_Loop_ListComp.py:2(<module>)", 
"DeDup_adjacencyOutput_Loop_ListComp_ReverseCount.py:2(<module>)"))


for (i in 1:6){
cumTime$Function<-gsub(rownames(ReNameFunction)[i],as.character(ReNameFunction[i,][1]),cumTime$Function,fixed=TRUE)
}


pdf('../TimeProfiling_cumTime.pdf',height=4,width=12)
ggplot(cumTime,aes(x=nrows,y=CumulativeTime,colour=Function))+
geom_line()+
geom_point(size=0.5)+
facet_wrap(~Complexity)+
scale_colour_brewer(palette ='Dark2')+
theme_classic()+
ylab('Cumulative time (s)')+
xlab('Number of input lines')+
scale_x_log10(  breaks = scales::trans_breaks("log10", function(x) 10^x),  labels = scales::trans_format("log10", scales::math_format(10^.x)))+
annotation_logticks(sides = "b",long = unit(0.2, "cm"))    
dev.off()



#####
#Check order of loop
#####
orderFunction<-list.files('../',recursive=F)
allOrders<-c('cdupListCountEditEqual','cdupListEdist','cdupListEditCount','cdupListEditCountEqual','cdupListEqualCountEdit','cdupListEqualEditCount')
orderFunction<-timeProfile[timeProfile%in%allOrders]

functionOrder<-as.data.frame(matrix(0, ncol = 2, nrow = length(orderFunction)))
colnames(functionOrder)<-c('CumulativeTime','Order')

library(stringr)
for (file in 1:length(orderFunction)){cat(as.character(orderFunction[file]),'\n')
	Time <-readLines(paste0('../',orderFunction[file]))[1]
	Time<-sub("^\\s+","",Time)
	TimeAll<-as.data.frame(t(unlist(str_split(Time," ",n=9))))$V8
	functionOrder[,1][file]<-as.numeric(as.character(TimeAll))
	functionOrder[,2][file]<-gsub('cdupList','',as.character(orderFunction[file]))
}

library(ggplot2)
pdf('../OrderOfLoop.pdf',height=4,width=4)
ggplot(functionOrder,aes(x=Order,y=CumulativeTime))+
geom_bar(stat = "identity")+
theme_classic()+
ylab('Cumulative time (s)')+
xlab('Order')+
ggtitle(paste0('Argument order effect on time for', '\n',' 1M low-complex reads'))+
theme(axis.title.x = element_text(size=14),axis.text.x  = element_text(size=12,angle=45,hjust=1),axis.title.y=element_text(size=14),axis.text.y=element_text(size=12))
#scale_x_log10(  breaks = scales::trans_breaks("log10", function(x) 10^x),  labels = scales::trans_format("log10", scales::math_format(10^.x)))+
##scale_y_log10(  breaks = scales::trans_breaks("log10", function(x) 10^x),  labels = scales::trans_format("log10", scales::math_format(10^.x)))+
#annotation_logticks(sides = "b",long = unit(0.2, "cm"))    
dev.off()






pdf('../TimeProfiling_nCalls.pdf',height=4,width=8)
ggplot(cumTime,aes(x=nrows,y=NumberCalls,colour=Complexity))+
geom_line()+
geom_point(size=0.5)+
facet_wrap(~Function,scale='free')+
scale_colour_brewer(palette ='Set1')+
theme_classic()+
ylab('Number of calls per function')+
xlab('Number of input lines')+
scale_x_log10(  breaks = scales::trans_breaks("log10", function(x) 10^x),  labels = scales::trans_format("log10", scales::math_format(10^.x)))+
#scale_y_log10(  breaks = scales::trans_breaks("log10", function(x) 10^x),  labels = scales::trans_format("log10", scales::math_format(10^.x)))+
annotation_logticks(sides = "b",long = unit(0.2, "cm"))    
dev.off()
