.libPaths(c(.libPaths(), "/home/eviun/R/x86_64-unknown-linux-gnu-library/2.15/"))

library(edgeR)
library(limma)

args=(commandArgs(TRUE))
if(length(args)==0){
  print("No arguments supplied.")
  ##supply default values
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}

#arguments include 1) path 2) priordf 3) spaired

setwd(path)

path
priordf
spaired

#data preparation
targets = read.delim(file="targets",stringsAsFactors=FALSE,row.names="label")
targets$tissue <- factor(targets$tissue)
targets$patient <- factor(targets$patient)

d <- readDGE(targets,comment.char="#")
colnames(d) = row.names(targets)
keep <- rowSums(cpm(d)) >= 1 #at least one count per million in at least one sample
d <- d[keep,]


#filter some uniq expression
#ff <- function (x, y) {
#  j = 0
#  for (i in 1:length(x)){
#    if (x[i] > y) j = j+1
#  }
# j
#}

#normalization
d$samples$lib.size <- colSums(d$counts) #recalculate sample size
d <- calcNormFactors(d)

#png(filename="MDS_bayer.png", width=600, height=600)
#plotMDS.dge(d2, main="MDS")
#dev.off()


#estimating dispersion (both common and tagwise)
if (spaired == 1) {
  design <- model.matrix(~patient+tissue, data = targets)
  d <- estimateGLMCommonDisp(d, design, verbose=TRUE)
  d <- estimateGLMTrendedDisp(d, design)
  d <- estimateGLMTagwiseDisp(d, design, prior.df = priordf)
} else {
  d$samples$patient <- NULL
  d$samples$files <- NULL
  colnames(d$samples) = c("group","lib.size","norm.factors")
  d <- estimateCommonDisp(d, verbose=TRUE)
  d <- estimateTagwiseDisp(d, prior.df = priordf)
}

#d$CR.common.dispersion
#sqrt(d$CR.common.dispersion)


#differential expression analysis

if (spaired == 1) {
  d.fit <- glmFit(d, design)
  d.lrt <- glmLRT(d.fit)
  d.top <- topTags(d.lrt, n=dim(d)[1])
  d.ifDE <- decideTestsDGE(d.lrt)
  d.ifDE.table <- cbind(d.lrt$table, d.ifDE)
} else {
  d.et <- exactTest(d, dispersion = "tagwise")
  d.top <- topTags(d.et, n=dim(d)[1])
  d.ifDE <- decideTestsDGE(d.et)
  d.ifDE.table <- cbind(d.et$table, d.ifDE)
}


write.table(d.top$table, file="topDE.txt", sep="\t", quote=F)
write.table(d.ifDE.table, file="ifDE.txt", sep="\t", quote=F)
