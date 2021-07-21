rm(list=ls())

#######
#### prepared inputs of ABSOLUTE

    #### CNVs


setwd("~/MPLC/cnvkit/")

file<-list.files(pattern=".seg$")

mclapply(1:47,function(x){
  
  na<-sub(".seg","",file[x])
  
  a<-data.table::fread(file[x])
  
  a<-a[,-1]
  
  colnames(a)<-c("Chromosome","Start","End","Num_Probes","Segment_Mean")

  a$Chromosome<-sub("chrX","chr23",a$Chromosome)
  
  chrom<-paste("chr",1:23,sep='')

  a<-a[a$Chromosome%in%chrom,]
  
  write.table(a,file = paste("~/MPLC/ABSOLUTE/input_cnv/",na,".cnv.txt",sep=""),sep="\t",row.names = F)
  
  },mc.cores = 47)






           #### SNVs

setwd("~/MPLC/maf/")

file<-list.files(pattern=".maf$")

mclapply(1:47,function(x){
  
  na<-sub(".maf","",file[x])
  
  a<-data.table::fread(file[x])
  
  a<-a[a$FILTER=="PASS",]
  
  #a<-a[,c(41,42,15,6,16,1,5)]
  
  colnames(a)[grep("Start_Position",colnames(a))]<-"Start_position"
  
  #colnames(a)<-c("t_ref_count","t_alt_count","dbSNP_Val_Status","Start_position","Tumor_Sample_Barcode","Hugo_Symbol","Chromosome")
  
  a$Chromosome<-sub("chr","",a$Chromosome)
  
  a$Chromosome<-sub("X","23",a$Chromosome)
  
  write.table(a,file = paste("~/MPLC/ABSOLUTE/input_snv/",na,".snv.maf",sep=""),sep="\t",row.names = F)
  
},mc.cores = 47)




####################
#### run




library(ABSOLUTE)
#library(foreach)
#genome<-"hg38"




setwd("~/MPLC/maf/")

file<-list.files(pattern=".maf$")
file<-sub(".maf","",file)

cnv_dir="~/MAPLC/ABSOLUTE/input_cnv/"
snv_dir="~/MAPLC/ABSOLUTE/input_snv/"
out_main_dir = "~/MAPLC/ABSOLUTE/output/"

log.dir<-paste(out_main_dir,"log/",sep="")
if(!file.exists(log.dir)){dir.create(log.dir,recursive = T)}

mclapply(1:47,function(x){
  
  
  seg.dat.fn<-paste(cnv_dir,file[x],".cnv.txt",sep="")
  maf.fn<-paste(snv_dir,file[x],".snv.maf",sep="")
  #a<-data.table::fread(seg.dat.fn)
  sample.name<-file[x]
  platform<-"Illumina_WES"
  primary.disease<-c("Lung Cancer")
  sigma.p<-0
  max.sigma.h<-0.015
  min.ploidy<-0.95
  max.ploidy<-10
  max.as.seg.count<-1500
  max.non.clonal<-0.05 
  max.neg.genome<- 0.005
  copy_num_type<-"total"
  min.mut.af<-0
  output.fn.base<-file[x]

    results.dir<-paste(out_main_dir,"run/",file[x],"/",sep="")
  
    if(!file.exists(results.dir)){dir.create(results.dir,recursive = T)}
    
    
    
    sink(paste("~/MAPLC/ABSOLUTE/output/sink/",file[x],".run",sep=""))
      RunAbsolute(seg.dat.fn,sigma.p,max.sigma.h,min.ploidy,max.ploidy,primary.disease,platform,sample.name,results.dir,max.as.seg.count,max.non.clonal,
              max.neg.genome,copy_num_type,maf.fn,min.mut.af,output.fn.base,verbose=T)
    sink()  
      
},mc.cores = 47)
  

#####
#####
#####

obj.name<-"MPLC_WES"

setwd("~/MAPLC/ABSOLUTE/output/run/")
absolute.files<-list.files(pattern=".RData$",recursive = T)
indv.results.dir<-"~/MAPLC/ABSOLUTE/output/summarize"

if(!file.exists(indv.results.dir)){dir.create(indv.results.dir,recursive = T)}

CreateReviewObject(obj.name, absolute.files, indv.results.dir, "total", verbose=TRUE)

#####
####
###

reviewed.pp.calls.fn<-"~/MAPLC/ABSOLUTE/output/summarize/MPLC_WES.PP-calls_tab.txt"
modes.fn<-"~/MAPLC/ABSOLUTE/output/summarize/MPLC_WES.PP-modes.data.RData"
out.dir.base<-"~/MAPLC/ABSOLUTE/output/review/"
if(!file.exists(out.dir.base)){dir.create(out.dir.base,recursive = T)}


analyst.id<-"Laveny"
ExtractReviewedResults(reviewed.pp.calls.fn, analyst.id, modes.fn, out.dir.base, obj.name, "total")
