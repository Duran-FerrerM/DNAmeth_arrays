
###########################################################
## Define dirs and necessary variables to start analysis
###########################################################

# call script on termial like this

"R CMD BATCH" 
"--args 
project='project_name'

samp.sheet.path.450k='path'
baseDir.450k='path'

samp.sheet.path.EPIC=NULL
baseDir.EPIC=NULL

cpus=2 
GetOnlyQualityControlReport=FALSE
DoQualityControl=TRUE
addRGsets='Endothelial_RGset.EC.LN.RData'
addRGsetsOutType='IlluminaHumanMethylation450k'
ArrayOutType='IlluminaHumanMethylation450k'

RemoveSNPs=TRUE
MAF=0
RemoveSexualChrs=TRUE
TypeNorm='SWAN'

"
"GetMethArrays.v.2.R"


options(stringsAsFactors = F,error=NULL)
if(!interactive()){
  args=(commandArgs(TRUE))
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}else{
  #User defined parameters interactively
  project <- "/Volumes/Marti_SSD/WORK/CLL_ERIC_Subsets/CLL_subsets_450k"
  
  #450k
  samp.sheet.path.450k <- "/Volumes/Marti_SSD/WORK/CLL_ERIC_Subsets/Extended_Cohort_CLL_Subsets_Illumina-450k.csv"
  baseDir.450k <- "/Volumes/Marti_SSD/WORK/IDATS/IDATS_450k/"
  
  #EPIC
  samp.sheet.path.EPIC <- NULL#"/media/mduran/5A129A25129A06631/Marti/EXTENDEND_ANALYSIS/REVISION/EPIC_DLBCLs/sampleSheet_IMS_31_m.csv"
  baseDir.EPIC <- NULL#"/media/mduran/5A129A25129A06631/Marti/IDATS/IDATS_850k/"
  
  # controls and calculations
  cpus <- 1
  DoQualityControl <- T
  GetOnlyQualityControlReport <- F
  addRGsets<- NULL#"/Volumes/Marti_SSD/WORK/pipeline/Metadata_each_Cancer/RGset_Creations/Endothelial_RGset.EC.LN.RData"
  if(!is.null(addRGsets)){
    addRGsetsOutType <- "IlluminaHumanMethylation450k"#"IlluminaHumanMethylation450k" #"IlluminaHumanMethylationEPIC"
    ArrayOutType <- "IlluminaHumanMethylation450k"#"IlluminaHumanMethylationEPIC"#"IlluminaHumanMethylation450k"
  }
  
  # NORMALIZATION
  RemoveSNPs <- TRUE
  MAF <- 0
  RemoveSexualChrs <- TRUE
  TypeNorm <- "SWAN" #

  
}



###########################################################
## Define dirs and necessary variables to start analysis
###########################################################

writeLines("Starting the analysis")
writeLines(paste0(rep("_",100),collapse = ""))

writeLines("Loading requires libraries...")
writeLines(paste0(rep("_",100),collapse = ""))
suppressMessages(library(minfi))
suppressMessages(library(IlluminaHumanMethylation450kmanifest))
suppressMessages(library(IlluminaHumanMethylation450kanno.ilmn12.hg19))
suppressMessages(library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19))
suppressMessages(library(IlluminaHumanMethylationEPICmanifest))
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))
suppressMessages(library(reshape2))
suppressMessages(library(data.table))


print(sessionInfo())
###########################################################
## Define dirs and necessary variables to start analysis
###########################################################

## mine arguments to extract things needed
projectName <- unlist(strsplit(x = project,split = "\\/"))[length(unlist(strsplit(x = project,split = "\\/")))]
dir.project <- unlist(strsplit(x = project,split = projectName))[[1]]
excludingCpGs <- as.character(read.table("/Volumes/Marti_SSD/WORK/pipeline/Pipelines/filters450k/ExcludingCpGs")[,1])

###########################################################
## Reading IDAT files and pheno data
###########################################################

writeLines("1. Reading IDAT files")
writeLines(paste0(rep("_",100),collapse = ""))
setwd(dir.project)
dir.create(project)
setwd(project)



options(scipen=999)
CreateRGsetFromIDATS <- function(IDATS.path = NULL,Sample.sheet.path = NULL){
  stopifnot(!is.null(IDATS.path) && !is.null(Sample.sheet.path))
  baseDir <- IDATS.path
  samp.sheet.path <- Sample.sheet.path
  samp.sheet <- unlist(strsplit(x = samp.sheet.path,split = "\\/"))[length(unlist(strsplit(x = samp.sheet.path,split = "\\/")))]
  
  # checking
  if(!dir.exists(dir.project)){
    stop(paste0("The following directory do not exist!:\n",dir.project))
  }
  if(!file.exists(samp.sheet.path)){
    stop(paste0("File do not exits!:\n",samp.sheet.path))
  }
  if(!dir.exists(baseDir)){
    stop(paste0("Directory of IDAT files do no exist:\n",baseDir))
  }
  if(file.exists(file.path(baseDir,samp.sheet))){
    invisible(file.remove(file.path(baseDir,samp.sheet)))
  }
  FileCopied <- file.copy(from = samp.sheet.path,to = baseDir)
  if(!FileCopied){
    stop(paste0("Could not start analysis. Cannot arrange necessary directories!"))
  }
  targets <- read.metharray.sheet(base = baseDir,pattern = samp.sheet,ignore.case = F,verbose = T)
  stopifnot(all(!is.na(targets$Basename)))
  # print(head(targets,20))
  # print(tail(targets,20))
  ##little check of Sample sheet
  if(!all(c("Sample_Name","Sample_Group","Sample_SubGroup")%in%colnames(targets))){
    stop("You shoul specify 'Sample_Name', 'Sample_Group','Sample_SubGroup' in you Sample Sheet")
  }
  write.table(targets,paste0("targets.",samp.sheet,".csv"),sep="\t",row.names = F,quote = F)
  targets <- targets[!is.na(targets$Slide),] # new line of code to allow mergind other data, eg endothelial cells
  RGset <- read.metharray.exp(targets = targets,verbose = T,force = T)
  invisible(gc())
  return(RGset)
}

## Read IDATS from 450k arrays
if(!is.null(baseDir.450k) && !is.null(samp.sheet.path.450k)){
  RGset450k <- CreateRGsetFromIDATS(IDATS.path = baseDir.450k,Sample.sheet.path = samp.sheet.path.450k)
  writeLines("450k data loaded!")
}
invisible(gc())

## Read IDATS from EPIC arrays
if(!is.null(baseDir.EPIC) && !is.null(samp.sheet.path.EPIC)){
  RGsetEPIC <- CreateRGsetFromIDATS(IDATS.path = baseDir.EPIC,Sample.sheet.path = samp.sheet.path.EPIC,Give.CorrectBetas.params = F)
  writeLines("EPIC data loaded!")
}
invisible(gc())

if(exists("RGset450k") && exists("RGsetEPIC")){
  RGset <- combineArrays(object1 = RGset450k,object2 = RGsetEPIC,outType = ArrayOutType,verbose=T)
  rm(RGset450k)
  rm(RGsetEPIC)
}else if(exists("RGset450k")){
  RGset <- RGset450k
  rm(RGset450k)
}else{
  RGset <- RGsetEPIC
  rm(RGsetEPIC)
}
invisible(gc())


## get phenoData of RGset for analysis
pd <- pData(RGset)
pd$Sample_SubGroup <- make.unique(pd$Sample_SubGroup)
rownames(pd) <- pd$Sample_SubGroup
rownames(pData(RGset)) <- pd$Sample_SubGroup

###########################################################
## Quality Control
###########################################################

if(DoQualityControl){
  
  writeLines("2. Performing quality control of samples...")
  writeLines(paste0(rep("_",100),collapse = ""))
  
  invisible(gc())
  dir.create("QC")
  setwd("QC")
  
  qcReport(rgSet = RGset)
  invisible(gc())
  
  pdf("Signal_Intensities.pdf",width = 10,height = 10)
  plotQC(getQC(preprocessRaw(rgSet = RGset)))
  grid()
  title("Signal intensities per sample")
  invisible(dev.off())
  invisible(gc())
}

if(GetOnlyQualityControlReport){
  writeLines("Quality control done!")
  quit()
}


## Add RGset for endothelial cells
##
## Merging possibly added new RGset where signal insted of IDATS are available.
##

if(!is.null(addRGsets)){
  writeLines("1. Merging RGset channels...")
  load(addRGsets,verbose = T)
  RGset <- combineArrays(object1 = RGset,object2 = RGset.endo,outType = addRGsetsOutType,verbose=T)
  invisible(gc())
}


###########################################################
## Start normalization
###########################################################

invisible(gc())
options(scipen=0)

writeLines("3. Starting normalization and filtering")
writeLines(paste0(rep("_",100),collapse = ""))

if(TypeNorm=="Functional"){
  set.seed(6)
  writeLines("Using 'Functional Normalization' to normalize data...")
  GRset <- preprocessFunnorm(rgSet = RGset,verbose = T)
}else if(TypeNorm=="Quantile"){
  writeLines("Using 'Quantile normalization' to normalize data...")
  GRset <- preprocessQuantile(object = RGset,verbose = T)
}else if(TypeNorm=="SWAN"){
  set.seed(6)
  writeLines("Using 'SWAN normalization' to normalize data...")
  GRset <- preprocessSWAN(rgSet = RGset,verbose = T)
  invisible(gc())
  GRset <- ratioConvert(mapToGenome(GRset))
}else if(TypeNorm=="ssNob"){
  set.seed(6)
  writeLines("Using 'ssNob normalization' to normalize data...")
  GRset <- preprocessNoob(rgSet = RGset,verbose = T)
  invisible(gc())
  GRset <- ratioConvert(mapToGenome(GRset))
}else{
  writeLines("Using 'Illumina normalization' to normalize data...")
  GRset <- preprocessIllumina(object = RGset)
  invisible(gc())
  GRset <- ratioConvert(mapToGenome(GRset))
}
invisible(gc())

##get full annotation for reporting info of removed CpGs
Full.CpGs.annot <- getAnnotation(GRset)

## Filering based on pvalues
writeLines("Finding probes with bad signal intensities...")
excludingCpGsPvals <- detectionP(rgSet = RGset)
invisible(gc())
excludingCpGsPvals <- excludingCpGsPvals>0.01
excludingCpGsPvals <- names(which(rowMeans(excludingCpGsPvals,na.rm = T)>0.1)) # Remove positions that do not have good pvalue in more than 10% of the samples 
invisible(gc())

if(length(excludingCpGsPvals)>0){
  writeLines(paste0(length(excludingCpGsPvals)," CpGs did not pass pvalue filter. Removing them."))
  GRset <- GRset[which(!rownames(GRset)%in%excludingCpGsPvals),]
}
invisible(gc())

## Obtain sex based on X and Y probes before removing further CpGs
Sex <- getSex(GRset)
GRset <- addSex(object = GRset,sex = Sex)
pSex <- ggplot(data.frame(pData(GRset)),
       aes(x = xMed,
           y= yMed,
           fill=predictedSex
           )
       )+
  geom_point(pch=21,size=3)+
  xlab("X chr, median total intensity (log2)")+ 
  ylab("Y chr, median total intensity (log2)")+
  theme_bw()
ggsave(filename = "Sex.predictions.pdf",plot = pSex,height = 5,width = 6,useDingbats = F)

##  get out SNPs
if(RemoveSNPs){
  writeLines(paste0("3.1 Getting out SNPs"))
  GRset <- dropMethylationLoci(object = GRset,dropRS = T,dropCH = T)
  GRset <- dropLociWithSnps(object = GRset,snps=c("SBE","CpG"), maf=MAF)
}
invisible(gc())

##  get out  and individual specific CpGs 
writeLines("3.2 Getting out individual-specific CpGs...")
GRset <- GRset[which(!rownames(GRset)%in%excludingCpGs),]

## Filter sexual chormoomes
if(RemoveSexualChrs){
  writeLines("3.3 Removing sexual chromosomes...")
  GRset <- GRset[which(as.character(seqnames(GRset))!="chrX" & as.character(seqnames(GRset))!="chrY"),]
  seqlevels(GRset) <- as.character(unique(seqnames(GRset)))
}
invisible(gc())

##report removed CpGs
Removed.CpGs <- Full.CpGs.annot[setdiff(rownames(Full.CpGs.annot),rownames(GRset)),c("chr","pos","strand","Name","Probe_rs","Probe_maf","CpG_rs","CpG_maf","SBE_rs","SBE_maf","Relation_to_Island","UCSC_RefGene_Group","UCSC_RefGene_Name")]
Removed.CpGs$Individual.specific.methylation.Bcells <- F
Removed.CpGs$Individual.specific.methylation.Bcells[which(Removed.CpGs$Name%in%excludingCpGs)] <- T

fwrite(as.data.frame(Removed.CpGs),"Removed.CpGs.tsv",sep="\t",na = "NA")

###########################################################
## 3. Get beta values
###########################################################

writeLines("4. Getting beta values...")
## Get betas
betas <- as.data.frame(getBeta(GRset))
colnames(betas) <- pd$Sample_SubGroup
invisible(gc())

## Gett matrix without converting NA's
betasNAs <- betas
whichBetasNAs <- is.na(getBeta(RGset))
whichBetasNAs <- whichBetasNAs[rownames(betas),]
betasNAs[whichBetasNAs] <- NA
rm("whichBetasNAs")
invisible(gc())

if(DoQualityControl){
  ## Quality control after normalization 
  ## Density of methylation
  dData <- mclapply(1:ncol(betasNAs),function(i){
    print(i)
    d <- density(betasNAs[,i],na.rm=T,from=0,to=1)
    return(data.frame(x=d$x,y=d$y,Samples=pd$Sample_SubGroup[i]))
  },mc.cores = cpus)
  dData <- do.call(rbind,dData)
  invisible(gc())
  
  ## Density of beta values
  dPlot <- ggplot(data = dData,aes(x = x,y = y,colour=Samples))
  dPlot <- dPlot + geom_line() +
    scale_x_continuous(breaks = seq(0,1,0.25)) +  xlab("Betas") + ylab("Density") + 
    ggtitle(paste0("Density of Beta values normalized with method: ",TypeNorm)) + 
    theme_bw()
  Legend <- cowplot::get_legend(dPlot)
  dPlot <- dPlot + theme(legend.position="none")
  ggsave(filename = "QualityControl1_Normalized.pdf",plot = dPlot,width = 8,height = 5)
  ggsave(filename = "QualityControl1_Normalized.Legend.png",plot = cowplot::plot_grid(Legend),width = 49,height = 5,units = "in")
  rm("dData")
  invisible(gc())
}



###########################################################
## 4. Export betas 
###########################################################

writeLines("4. Exporting betas without purity adjustments...")
writeLines(paste0(rep("_",100),collapse = ""))

dir.create("../Results")
setwd("../Results")
data.table::fwrite(x = cbind(CpGs=rownames(betas),as.data.frame(betas)),
                   file = paste0(projectName,"_Betas.txt"),sep = "\t")
data.table::fwrite(x = cbind(CpGs=rownames(betasNAs),as.data.frame(betasNAs)),
                   file = paste0(projectName,"_BetasNAs.txt"),sep = "\t")

if(!is.null(samp.sheet.path.450k))
  invisible(file.remove(paste0(baseDir.450k,"/",unlist(strsplit(x = samp.sheet.path.450k,split = "\\/"))[length(unlist(strsplit(x = samp.sheet.path.450k,split = "\\/")))])))
if(!is.null(samp.sheet.path.EPIC))
  invisible(file.remove(paste0(baseDir.EPIC,"/",unlist(strsplit(x = samp.sheet.path.EPIC,split = "\\/"))[length(unlist(strsplit(x = samp.sheet.path.EPIC,split = "\\/")))])))

save(RGset,GRset,file="RGset.GRset.RData")
invisible(gc())

## save phenoData of RGset for analysis
fwrite(data.frame(pData(GRset)),paste0("../PhenoData_",projectName,".tsv"),sep = "\t",na = "NA")

writeLines("ANALYSIS DONE! Congrats! :-)")





