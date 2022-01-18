
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
"GetMethArrays.v.1.1.R"


options(stringsAsFactors = F,error=NULL)

if(!interactive()){
  args=(commandArgs(TRUE))
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}else{
  #User defined parameters interactively
  project <- "/media/mduran/Marti_SSD1/WORK/ALCL/Control_Tcells_TALL"
  
  #450k
  samp.sheet.path.450k <- NULL#"/Volumes/Marti_SSD/WORK/CLL_ERIC_Subsets/Extended_Cohort_CLL_Subsets_Illumina-450k.csv"
  baseDir.450k <- NULL#"/Volumes/Marti_SSD/WORK/IDATS/IDATS_450k/"
  
  #EPIC
  samp.sheet.path.EPIC <- "/media/mduran/Marti_SSD1/WORK/ALCL/Data_Marti/GSE155333_Phenodata_Complete_for_analysis.csv"
  baseDir.EPIC <- "/media/mduran/Marti_SSD1/WORK/ALCL/Data_Marti/GSE155333/GSE155333_RAW/"
  
  # controls and calculations
  excludingCpGs <- NULL#as.character(read.table("/Volumes/Marti_SSD/WORK/pipeline/Pipelines/filters450k/ExcludingCpGs")[,1])
  cpus <- 1
  DoQualityControl <- T
  PlotDensityNormBetas <- F

  # NORMALIZATION
  RemoveSNPs <- TRUE
  MAF <- 0
  RemoveSexualChrs <- TRUE
  TypeNorm <- "ssNob"#SWAN" #

  
}



###########################################################
## Define dirs and necessary variables to start analysis
###########################################################

writeLines("Starting the analysis")
writeLines(paste0(rep("_",100),collapse = ""))

writeLines("Loading requires libraries")
writeLines(paste0(rep("_",100),collapse = ""))
suppressMessages(library(minfi))
suppressMessages(library(IlluminaHumanMethylation450kmanifest))
suppressMessages(library(IlluminaHumanMethylationEPICmanifest))
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

###########################################################
## Reading IDAT files and pheno data
###########################################################

writeLines("Reading IDAT files")
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
  if(!all(c("Sample_Name","Sample_Group","Sample_Name_Analysis")%in%colnames(targets))){
    stop("You should specify 'Sample_Name', 'Sample_Group','Sample_Name_Analysis' in you Sample Sheet")
  }
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
  RGsetEPIC <- CreateRGsetFromIDATS(IDATS.path = baseDir.EPIC,Sample.sheet.path = samp.sheet.path.EPIC)
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
pd$Sample_Name_Analysis <- make.unique(pd$Sample_Name_Analysis)
rownames(pd) <- pd$Sample_Name_Analysis
rownames(pData(RGset)) <- pd$Sample_Name_Analysis

###########################################################
## Quality Control
###########################################################

if(DoQualityControl){
  
  writeLines("Performing quality control of samples")
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


###########################################################
## Start normalization
###########################################################

invisible(gc())
options(scipen=0)

writeLines("Normalizing and filtering")
writeLines(paste0(rep("_",100),collapse = ""))

if(TypeNorm=="Functional"){
  set.seed(6)
  warning("Using 'Functional Normalization' to normalize data. Between array normalization, beta values will change upon batches! Returning an grSet.")
  RGset.norm <- preprocessFunnorm(rgSet = RGset,verbose = T)
}else if(TypeNorm=="Quantile"){
  warning("Using 'Quantile normalization' to normalize data. Between array normalization, beta values will change upon batches! Returning an grSet.")
  RGset.norm <- preprocessQuantile(object = RGset,verbose = T)
}else if(TypeNorm=="SWAN"){
  set.seed(6)
  writeLines("Using 'SWAN normalization' to normalize data")
  RGset.norm <- preprocessSWAN(rgSet = RGset,verbose = T)
  invisible(gc())
}else if(TypeNorm=="ssNob"){
  set.seed(6)
  writeLines("Using 'ssNob normalization' to normalize data")
  RGset.norm <- preprocessNoob(rgSet = RGset,verbose = T)
  invisible(gc())
}else{
  writeLines("Using 'Illumina normalization' to normalize data")
  RGset.norm <- preprocessIllumina(object = RGset)
  invisible(gc())
}
invisible(gc())


##get full annotation for reporting info of removed CpGs
Full.CpGs.annot <- getAnnotation(RGset.norm)

if(class(RGset.norm)=="MethylSet"){
  grSet <- mapToGenome(ratioConvert(RGset.norm))
}else{
  grSet <- RGset.norm
}

## Obtain sex based on X and Y probes before removing any CpGs
Sex <- getSex(grSet)
grSet <- addSex(object = grSet,sex = Sex)
pSex <- ggplot(data.frame(pData(grSet)),
               aes(x = xMed,
                   y= yMed,
                   fill=predictedSex
               )
)+
  geom_point(pch=21,size=3)+
  ggtitle(paste0("Intensity values normalized by ",TypeNorm))+
  xlab("X chr, median total intensity (log2)")+ 
  ylab("Y chr, median total intensity (log2)")+
  theme_bw()
ggsave(filename = "Sex.predictions.pdf",plot = pSex,height = 5,width = 6,useDingbats = F)

##EXPORT RGset AND normalized RGset data
pData(RGset.norm) <- pData(grSet)
save(RGset,RGset.norm,file=paste0("RGset.RGset.norm.",TypeNorm,".RData"))

## Filering based on pvalues
writeLines("Finding probes with bad signal intensities")
CpGs.pval.ok <- detectionP(rgSet = RGset)
invisible(gc())
CpGs.pval.ok <- CpGs.pval.ok<=1e-16 ##extreme pvalue increases data quality!
CpGs.pval.ok <- rownames(CpGs.pval.ok)[which(rowAlls(x = CpGs.pval.ok,value = TRUE))]
invisible(gc())

##release some RAM
rm(RGset.norm)
invisible(gc())

if(length(CpGs.pval.ok)<length(rownames(grSet))){
  writeLines("Removing CpGs with low signal intensity")
  grSet <- grSet[rownames(grSet)[which(rownames(grSet) %in% CpGs.pval.ok)],] ## in b4 annotation, some cg are not mapped to annotation package.
}
invisible(gc())


##  get out SNPs
if(RemoveSNPs){
  writeLines(paste0("Getting out SNPs"))
  grSet <- dropMethylationLoci(object = grSet,dropRS = T,dropCH = T)
  grSet <- dropLociWithSnps(object = grSet,snps=c("SBE","CpG"), maf=MAF)
}
invisible(gc())

##  get out  and individual specific CpGs
if(!is.null(excludingCpGs)){
  writeLines("Getting out individual-specific CpGs...")
  grSet <- grSet[which(!rownames(grSet)%in%excludingCpGs),]
}

## Filter sexual chormoomes
if(RemoveSexualChrs){
  writeLines("Removing sexual chromosomes")
  grSet <- grSet[which(as.character(seqnames(grSet))!="chrX" & as.character(seqnames(grSet))!="chrY"),]
  seqlevels(grSet) <- as.character(unique(seqnames(grSet)))
}
invisible(gc())

##report removed CpGs
Removed.CpGs <- Full.CpGs.annot[setdiff(rownames(Full.CpGs.annot),rownames(grSet)),c("chr","pos","strand","Name","Probe_rs","Probe_maf","CpG_rs","CpG_maf","SBE_rs","SBE_maf","Relation_to_Island","UCSC_RefGene_Group","UCSC_RefGene_Name")]
if(!is.null(excludingCpGs)){
  Removed.CpGs$Individual.specific.methylation.Bcells <- F
  Removed.CpGs$Individual.specific.methylation.Bcells[which(Removed.CpGs$Name%in%excludingCpGs)] <- T  
}

fwrite(as.data.frame(Removed.CpGs),"Removed.CpGs.tsv",sep="\t",na = "NA")

###########################################################
## 3. Get beta values
###########################################################

writeLines("Getting beta values")
## Get betas
betas <- as.data.frame(getBeta(grSet))
invisible(gc())

## Gett matrix without converting NA's
betasNAs <- betas
whichBetasNAs <- is.na(getBeta(RGset))
whichBetasNAs <- whichBetasNAs[rownames(betas),]
betasNAs[whichBetasNAs] <- NA
rm("whichBetasNAs")
invisible(gc())

##MDS plot
pdf("mdsPlot.pdf",height = 6,width = 6)
mdsPlot(dat = as.matrix(betas),sampGroups = pd$Sample_Group,sampNames = pd$Sample_Name_Analysis)
dev.off()

if(PlotDensityNormBetas){
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

## save phenoData of RGset for analysis
fwrite(data.frame(pData(grSet)),paste0("../PhenoData_",projectName,".tsv"),sep = "\t",na = "NA")

writeLines("ANALYSIS DONE! Congrats! :-)")





