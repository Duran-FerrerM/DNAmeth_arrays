
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
ArrayOutType='IlluminaHumanMethylation450k'

excludingCpGs=NULL
RemoveSNPs=TRUE
Detection.Pval.Mean.Sample.cutoff=0.01
Detection.Pval.CpGs.cutoff=1e-16
MAF=0
RemoveSexualChrs=TRUE
Normalization.type='SWAN'

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
  project <- "name_project"
  
  #450k
  samp.sheet.path.450k <- NULL#"path"
  baseDir.450k <- NULL#"path"
  
  #EPIC
  samp.sheet.path.EPIC <- "path"
  baseDir.EPIC <- "path"
  
  # controls and calculations
  excludingCpGs <- NULL#as.character(read.table("path")[,1]) ## B-cell related stuff!
  cpus <- 1

  # NORMALIZATION
  ArrayOutType<-NULL
  Detection.Pval.Mean.Sample.cutoff <- 0.01
  Detection.Pval.CpGs.cutoff <- 1e-16
  RemoveSNPs <- TRUE
  MAF <- 0
  RemoveSexualChrs <- TRUE
  Normalization.type <- "ssNob"#SWAN" #2 recomended, batch-independent normalizations.

  
}



###########################################################
## Define dirs and necessary variables to start analysis
###########################################################

writeLines("Starting the analysis")
writeLines(paste0(rep("_",100),collapse = ""))

writeLines("Loading requires libraries")
writeLines(paste0(rep("_",100),collapse = ""))
suppressMessages(library(minfi))
# suppressMessages(library(IlluminaHumanMethylation450kmanifest))
# suppressMessages(library(IlluminaHumanMethylationEPICmanifest))
# suppressMessages(library(IlluminaHumanMethylation450kanno.ilmn12.hg19))
# suppressMessages(library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19))
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))
suppressMessages(library(data.table))


print(sessionInfo())
###########################################################
## Define dirs and necessary variables to start analysis
###########################################################

## mine arguments to extract things needed
project.name <- unlist(strsplit(x = project,split = "\\/"))[length(unlist(strsplit(x = project,split = "\\/")))]
dir.project <- unlist(strsplit(x = project,split = project.name))[[1]]

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
rownames(pData(RGset)) <- make.unique(as.character(pData(RGset)$Sample_Name_Analysis))

##save RGSet
writeLines("Saving RGSet object")
dir.create("RData")
save(RGset,file=paste0("RData/RGset.RData"))

###########################################################
## Quality Control
###########################################################



writeLines("Performing quality control of samples")
writeLines(paste0(rep("_",100),collapse = ""))

invisible(gc())
dir.create("QC")

##general QC report
qcReport(rgSet = RGset,pdf = "QC/qcReport.pdf",
         sampNames = pData(RGset)$Sample_Name_Analysis,
         sampGroups = factor(pData(RGset)$Sample_Group,levels = unique(pData(RGset)$Sample_Group)),
         maxSamplesPerPage = ifelse(nrow(pData(RGset))>100,100,nrow(pData(RGset)))
)
invisible(gc())

##general mdsplot
##MDS plot
set.seed(6)
pdf("QC/mdsPlot.RGset.pdf",height = 6,width = 6)
mdsPlot(dat = RGset,
        sampGroups = factor(pData(RGset)$Sample_Group,levels = unique(pData(RGset)$Sample_Group)),
        sampNames = pData(RGset)$Sample_Name_Analysis
        )
dev.off()

##calculate intensity signal
RGset.preprocessedRaw <- preprocessRaw(rgSet = RGset)
invisible(gc())
qc <- getQC(object = RGset.preprocessedRaw)
qc$Signal.Intensity.ok <- ifelse(((qc$mMed + qc$uMed)/2)>=10.5,TRUE,FALSE) ##minfi cutoff

pSignals <- ggplot(data.frame(qc),
                   aes(x = mMed,
                       y= uMed,
                       fill=factor(Signal.Intensity.ok,levels = c("FALSE","TRUE"),labels = c("Bad","Good"))
                   )
)+
  geom_point(pch=21,size=2)+
  geom_text_repel(aes(label=ifelse(!qc$Signal.Intensity.ok,rownames(qc),NA)),
                  size=1,max.overlaps = Inf,min.segment.length = 0,segment.size=0.2
  )+
  guides(fill=guide_legend(title = "Signal intensity"))+
  xlab("Meth median intensity (log2)")+ 
  ylab("Unmeth median intensity (log2)")+
  ggtitle(label = "Samples signal intensities",
          subtitle = paste0("Good=",length(which(qc$Signal.Intensity.ok)),
                            ", Bad=",length(which(!qc$Signal.Intensity.ok))
          )
  )+
  coord_cartesian(xlim = c(8,14),ylim = c(8,14))+
  geom_abline(intercept = 10.5*2,slope = -1,lty=2)+
  theme_bw()
ggsave(filename = "QC/Signal_Intensities.pdf",plot = pSignals,height = 5,width = 6,useDingbats = F)

##MDS plot with bad and good signal samples highlighted
set.seed(6)
pdf("QC/mdsPlot.RGset.Signal.intensity.ok.pdf",height = 6,width = 6)
mdsPlot(dat = RGset,
        sampGroups = qc$Signal.Intensity.ok,
        sampNames = rownames(qc)
        )
dev.off()
invisible(gc())

## Obtain sex based on X and Y probes before removing any CpGs
Sex <- getSex(object = mapToGenome(RGset.preprocessedRaw))
pSex <- ggplot(data.frame(Sex),
               aes(x = xMed,
                   y= yMed,
                   fill=predictedSex
               )
)+
  geom_point(pch=21,size=3)+
  ggtitle(paste0("Raw intensity values"))+
  xlab("X chr, median total intensity (log2)")+ 
  ylab("Y chr, median total intensity (log2)")+
  theme_bw()
ggsave(filename = "QC/Sex.predictions.pdf",plot = pSex,height = 5,width = 6,useDingbats = F)

##get NAs from CpGs before normalization
whichBetasNAs <- getBeta(object = RGset.preprocessedRaw) ##get CpGs with NAs
invisible(gc())
whichBetasNAs <- is.na(whichBetasNAs)
rm(RGset.preprocessedRaw)
invisible(gc())




###########################################################
## Start normalization
###########################################################

invisible(gc())
options(scipen=0)

writeLines("Normalizing and filtering")
writeLines(paste0(rep("_",100),collapse = ""))

if(Normalization.type=="Functional"){
  set.seed(6)
  warning("Using 'Functional Normalization' to normalize data. Between array normalization, beta values will change upon batches! Returning an grSet.")
  RGset.norm <- preprocessFunnorm(rgSet = RGset,verbose = T)
}else if(Normalization.type=="Quantile"){
  warning("Using 'Quantile normalization' to normalize data. Between array normalization, beta values will change upon batches! Returning an grSet.")
  RGset.norm <- preprocessQuantile(object = RGset,verbose = T)
}else if(Normalization.type=="SWAN"){
  set.seed(6)
  writeLines("Using 'SWAN normalization' to normalize data")
  RGset.norm <- preprocessSWAN(rgSet = RGset,verbose = T)
  invisible(gc())
}else if(Normalization.type=="ssNob"){
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
  gc()
}else{
  grSet <- RGset.norm
}

## add sex predictions and intensity signal to phenodata of normalized RGSET
grSet <- addSex(object = grSet,sex = Sex)
grSet <- addQC(object = grSet,qc = qc)

##EXPORT RGset AND normalized RGset data
pData(RGset.norm) <- pData(grSet)
save(RGset.norm,file=paste0("RData/RGset.norm.",Normalization.type,".RData"))

##release some RAM
rm(RGset.norm)
invisible(gc())

## Filering based on pvalues
writeLines("Finding probes with bad signal intensities")
CpGs.detection.pvals <- detectionP(rgSet = RGset)
invisible(gc())
save(CpGs.detection.pvals,file = "RData/CpGs.detection.pvals.RData")

## visualize detection pvalues from all samples
pData(object = grSet)$Mean.detection.pval <- colMeans(CpGs.detection.pvals)

#plot
pMean.detection.pval <- ggplot(data.frame(pData(grSet)[order(pData(grSet)$Mean.detection.pval),]),
       aes(x = factor(Sample_Name_Analysis,levels = unique(Sample_Name_Analysis)),
           y = Mean.detection.pval,
           fill=Signal.Intensity.ok
           )
       )+
  geom_col()+
  geom_hline(yintercept = Detection.Pval.Mean.Sample.cutoff,color="red",linetype="dashed",lwd=0.2)+
  coord_cartesian(ylim = c(0,0.05))+
  facet_grid(.~Sample_Group,scales = "free_x",space = "free_x")+
  xlab("Samples")+
  ylab("Mean detection pvalue")+
  theme_bw()+
  theme(axis.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.text.x = element_text(size = 5,angle = 90)
        )
ggsave(filename = "QC/Mean.Detection.Pvals.Samples.pdf",plot = pMean.detection.pval,height = 3,width = 10,useDingbats=F)


##FILTER CPGS BASED ON SAMPLES WITH GOOD INTENSITY SIGNAL DEFINED ABOVE ACCORDING TO MINFI CUTOFF,extreme pvalue increases data quality! CpGs nearby deletions are removed
Samples.Signal.Pvals.ok <- pData(grSet)$Sample_Name_Analysis[which(pData(grSet)$Signal.Intensity.ok & pData(grSet)$Mean.detection.pval<=Detection.Pval.Mean.Sample.cutoff)]
pData(grSet)$Sample.good.quality <- ifelse(pData(grSet)$Sample_Name_Analysis %in% Samples.Signal.Pvals.ok,TRUE,FALSE)

CpGs.pval.ok <- rownames(CpGs.detection.pvals)[which(rowAlls(x = CpGs.detection.pvals[,which(colnames(CpGs.detection.pvals) %in% Samples.Signal.Pvals.ok)]<=Detection.Pval.CpGs.cutoff,value = TRUE))] ## all samples with lower pval
length(CpGs.pval.ok)

##Release RAM
rm(CpGs.detection.pvals)
invisible(gc())


if(length(CpGs.pval.ok)<length(rownames(grSet))){
  writeLines("Removing CpGs with bad signal/nose ratio intensity")
  grSet <- grSet[rownames(grSet)[which(rownames(grSet) %in% CpGs.pval.ok)],] ## in b4 annotation, some cg are not mapped to annotation package.
}
invisible(gc())


##  get out SNPs
if(RemoveSNPs){
  writeLines(paste0("Removing SNPs"))
  grSet <- dropMethylationLoci(object = grSet,dropRS = T,dropCH = T)
  grSet <- dropLociWithSnps(object = grSet,snps=c("SBE","CpG"), maf=MAF)
}
invisible(gc())

##  get out  and individual specific CpGs
if(!is.null(excludingCpGs)){
  writeLines("Getting out individual-specific CpGs...")
  grSet <- grSet[which(!rownames(grSet)%in%excludingCpGs),]
}
invisible(gc())


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
invisible(gc())

fwrite(as.data.frame(Removed.CpGs),"QC/Removed.CpGs.tsv",sep="\t",na = "NA")

###########################################################
## 3. Get beta values
###########################################################

writeLines("Getting beta values")
## Get betas
betas <- as.data.frame(getBeta(grSet))
invisible(gc())

# ## Gett matrix without converting NA's by normalizations
betasNAs <- betas
whichBetasNAs <- whichBetasNAs[rownames(betas),]
betasNAs[whichBetasNAs] <- NA
rm("whichBetasNAs")
invisible(gc())

##MDS plot with normalized data
set.seed(6)
pdf(paste0("QC/mdsPlot.normalized.",Normalization.type,".pdf"),height = 6,width = 6)
mdsPlot(dat = as.matrix(betas),
        sampGroups = factor(pData(grSet)$Sample_Group,levels = unique(pData(grSet)$Sample_Group)),
        sampNames = pData(grSet)$Sample_Name_Analysis
        )
dev.off()
invisible(gc())

set.seed(6)
pdf(paste0("QC/mdsPlot.normalized.",Normalization.type,".Samples.good.quality.pdf"),height = 6,width = 6)
mdsPlot(dat = as.matrix(betas),
        sampGroups = pData(grSet)$Sample.good.quality,
        sampNames = pData(grSet)$Sample_Name_Analysis
        )
dev.off()
invisible(gc())


set.seed(6)
pdf(paste0("QC/densityPlot.",Normalization.type,".Samples.good.quality.pdf"),height = 6,width = 6)
densityPlot(dat = as.matrix(betas[,pData(grSet)$Sample_Name_Analysis[which(pData(grSet)$Sample.good.quality)]]),
                sampGroups = factor(pData(grSet)$Sample_Group[which(pData(grSet)$Sample.good.quality)],levels = unique(pData(grSet)$Sample_Group[which(pData(grSet)$Sample.good.quality)])),
                sampNames = pData(grSet)$Sample_Name_analysis[which(pData(grSet)$Sample.good.quality)]
)
dev.off()
invisible(gc())


if(any(!pData(grSet)$Sample.good.quality)){
  set.seed(6)
  pdf(paste0("QC/densityPlot.",Normalization.type,".Samples.bad.quality.pdf"),height = 6,width = 6)
  densityPlot(dat = as.matrix(betas[,pData(grSet)$Sample_Name_Analysis[which(!pData(grSet)$Sample.good.quality)]]),
              sampGroups = factor(pData(grSet)$Sample_Group[which(pData(grSet)$Sample.good.quality)],levels = unique(pData(grSet)$Sample_Group[which(pData(grSet)$Sample.good.quality)])),
              sampNames = pData(grSet)$Sample_Name_analysis[which(!pData(grSet)$Sample.good.quality)],ylim = c(0,5)
  )
  dev.off()
  invisible(gc())
}




set.seed(6)
pdf(paste0("QC/densityBeanPlot.",Normalization.type,".Samples.good.quality.pdf"),height = 6,width = 6)
densityBeanPlot(dat = as.matrix(betas[,pData(grSet)$Sample_Name_Analysis[which(pData(grSet)$Sample.good.quality)]]),
                sampGroups = factor(pData(grSet)$Sample_Group[which(pData(grSet)$Sample.good.quality)],levels = unique(pData(grSet)$Sample_Group[which(pData(grSet)$Sample.good.quality)])),
                sampNames = pData(grSet)$Sample_Name_analysis[which(pData(grSet)$Sample.good.quality)]
            )
dev.off()
invisible(gc())

###########################################################
## 4. Export betas 
###########################################################

writeLines("Exporting normalized betas")
writeLines(paste0(rep("_",100),collapse = ""))

dir.create("Normalized.Betas")
data.table::fwrite(x = cbind(CpGs=rownames(betas),as.data.frame(betas)),
                   file = paste0("Normalized.Betas/",project.name,"_Betas_",Normalization.type,".tsv"),sep = "\t")
data.table::fwrite(x = cbind(CpGs=rownames(betasNAs),as.data.frame(betasNAs)),
                   file = paste0("Normalized.Betas/",project.name,"_Betas.w.NA_",Normalization.type,".tsv"),sep = "\t")

## export all phenodata associated
data.table::fwrite(x = as.data.frame(pData(grSet)),
                   file = paste0("Normalized.Betas/Phenodata_",project.name,".tsv"),sep = "\t",na = "NA")



if(!is.null(samp.sheet.path.450k))
  invisible(file.remove(paste0(baseDir.450k,"/",unlist(strsplit(x = samp.sheet.path.450k,split = "\\/"))[length(unlist(strsplit(x = samp.sheet.path.450k,split = "\\/")))])))
if(!is.null(samp.sheet.path.EPIC))
  invisible(file.remove(paste0(baseDir.EPIC,"/",unlist(strsplit(x = samp.sheet.path.EPIC,split = "\\/"))[length(unlist(strsplit(x = samp.sheet.path.EPIC,split = "\\/")))])))

writeLines("ANALYSIS DONE! Congrats! :-)")










