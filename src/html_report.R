library(R2HTML)
library(RColorBrewer)
library(fields)
library(KernSmooth)
library(lattice)

args=(commandArgs(TRUE))
if(length(args)==0){
  print("No arguments supplied.")
  ##supply default values
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}

#arguments include 1).lane 2).path
setwd(path)

stats.dir  = paste(path, "03_STATS", sep = "/")
reads.dir  = paste(path, "01_READS", sep = "/")
name.rep   = paste(lane, "report", sep = ".")
dir.html   = paste(stats.dir, name.rep, sep ="/")
dir.create(dir.html)

target <- HTMLInitFile(dir.html, name.rep, BackGroundColor="FFFFFF",Title = lane)
HTML(as.title(paste("Mapping report for", lane)),file=target)


##The table of mapping statistics
HTML("<br> General statistics on the mapping:", file = target)
mapping.stats.file = paste(lane, "mapping.stats", sep = ".")
mapping.stats.file.real = paste(stats.dir, mapping.stats.file, sep = "/")
mapping.stats = read.table(mapping.stats.file.real)
pN = function(x,...) {return(prettyNum(x, big.mark = ",", ...))} #prettyNum
pP = function(x,...) {return(paste(round(x*100,0),"%",sep=""))}  #percentage
row.infos =  c("N read-pairs sequenced ",
               "N read-pairs mapped ",
               "N read-pairs unmappable ",
               "N read-pairs mapped uniquely ",
               "N read-pairs mapped as singleton ",
               "N read-pairs mapped as proper pair ",
               "N read-pairs mapped as wrong pair ",
               "N read-pairs mapped as exon junction ",
               "N read-pairs mapped to multi-locations "
              )

row.values = c(pN(mapping.stats[13,2]),
               pN(mapping.stats[1,2]),
               pN(mapping.stats[13,2]-mapping.stats[1,2]),
               pN(mapping.stats[4,2]),
               pN(mapping.stats[8,2]),
               pN(mapping.stats[9,2]),
               pN(mapping.stats[10,2]),
               pN(mapping.stats[11,2]),
               pN(mapping.stats[12,2])    
              )

row.fracs  = c(pP(mapping.stats[13,2]/mapping.stats[13,2]),
               pP(mapping.stats[1,2]/mapping.stats[13,2]),
               pP((mapping.stats[13,2]-mapping.stats[1,2])/mapping.stats[13,2]),
               pP(mapping.stats[4,2]/mapping.stats[13,2]),
               pP(mapping.stats[8,2]/mapping.stats[13,2]),
               pP(mapping.stats[9,2]/mapping.stats[13,2]),
               pP(mapping.stats[10,2]/mapping.stats[13,2]),
               pP(mapping.stats[11,2]/mapping.stats[13,2]),
               pP(mapping.stats[12,2]/mapping.stats[13,2])    
              )
mapping.stats.summary = data.frame(Property=row.infos, Count=row.values, Percent=row.fracs)
HTML(mapping.stats.summary, innerBorder = 2, file= target)



#quality distribution of reads
HTML("<br> Read base quality distributions:", file = target)
qual.plot = "base_qualities.png"
mate1.qc = paste(lane,"_1.fq.qc",sep="")
mate2.qc = paste(lane,"_2.fq.qc",sep="")
reads.qual.mate1 = read.table(paste(reads.dir, mate1.qc, sep="/"), header=T)
reads.qual.mate2 = read.table(paste(reads.dir, mate2.qc, sep="/"), header=T)

cols.qual = brewer.pal(3, "Set1")[1:2]
quality.profile1 = cbind(reads.qual.mate1$mean, reads.qual.mate1$med)
quality.profile2 = cbind(reads.qual.mate2$mean, reads.qual.mate2$med)

png(file = paste(dir.html, qual.plot, sep ="/"), width = 1200, height = 400)
layout(matrix(1:2,1,2))
matplot(quality.profile1, type = "l" , col = cols.qual , lty = 1, lwd = 2,
        main = "Quality values (mate 1)", ylab = "Phred Quality Score", xlab = "position")
legend(x = "topright", col = cols.qual,legend = c("mean","median"), pch = 19)

matplot(quality.profile2, type = "l" , col = cols.qual , lty = 1, lwd = 2,
        main = "Quality values (mate 2)", ylab = "Phred Quality Score", xlab = "position")
legend(x = "topright", col = cols.qual,legend = c("mean","median"), pch = 19)
dev.off()

HTMLInsertGraph(Caption = "profile of average and median quality values",
                GraphFileName = qual.plot, Width = 1200,file = target)


HTML("<br> Read base composition distributions:", file = target)
compos.plot = "base_composition.png"
png(file = paste(dir.html, compos.plot, sep ="/"), width = 1200, height = 400)
layout(matrix(1:2,1,2))
plot(reads.qual.mate1$column, reads.qual.mate1$A_Count/reads.qual.mate1$count,type="l",ylim=c(0,0.5),col="#008B8B",lwd=2,main="base composition (mate 1)", xlab="base", ylab="fraction")
lines(reads.qual.mate1$column, reads.qual.mate1$C_Count/reads.qual.mate1$count,col="#CD8500",lwd=2)
lines(reads.qual.mate1$column, reads.qual.mate1$G_Count/reads.qual.mate1$count,col="#B03060",lwd=2)
lines(reads.qual.mate1$column, reads.qual.mate1$T_Count/reads.qual.mate1$count,col="#436EEE",lwd=2)
lines(reads.qual.mate1$column, reads.qual.mate1$N_Count/reads.qual.mate1$count,col="#030303",lwd=2)
legend("topright",legend=c("A","C","G","T","N"),col=c("#008B8B","#CD8500","#B03060","#436EEE","#030303"),bty="n",lwd=2)

plot(reads.qual.mate2$column, reads.qual.mate2$A_Count/reads.qual.mate2$count,type="l",ylim=c(0,0.5),col="#008B8B",lwd=2,main="base composition (mate 2)", xlab="base", ylab="fraction")
lines(reads.qual.mate2$column, reads.qual.mate2$C_Count/reads.qual.mate2$count,col="#CD8500",lwd=2)
lines(reads.qual.mate2$column, reads.qual.mate2$G_Count/reads.qual.mate2$count,col="#B03060",lwd=2)
lines(reads.qual.mate2$column, reads.qual.mate2$T_Count/reads.qual.mate2$count,col="#436EEE",lwd=2)
lines(reads.qual.mate2$column, reads.qual.mate2$N_Count/reads.qual.mate2$count,col="#030303",lwd=2)
legend("topright",legend=c("A","C","G","T","N"),col=c("#008B8B","#CD8500","#B03060","#436EEE","#030303"),bty="n",lwd=2)
dev.off()

HTMLInsertGraph(Caption = "base composition plot",
                GraphFileName = compos.plot, Width = 1200,file = target)


#insert size distributions
HTML("<br> Estimated fragment size based on the mapping result:", file = target)
insert.file = paste(lane,"ins",sep=".")
insert.size = read.table(paste(stats.dir,insert.file,sep="/"))
insert.plot = "insert_size.png"
png(file = paste(dir.html, insert.plot, sep ="/"), width = 600, height = 600)
plot(density(insert.size[,1]), xlab="size (bp)", ylab="Frenquency", main="Estimated Fragment Size")
dev.off()

HTMLInsertGraph(Caption = "predicted fragment size based on the mapping",
                GraphFileName = insert.plot, Width = 600,file = target)


#chromosome distribution of uniquely mapped reads
HTML("<br> Uniquely Mappable Reads (UMR) on Chromosomes:", file = target)
transcriptome.length.hg19 = paste(anno, "transcriptome_length_HG19", sep="/")
tmp1 <- pipe(paste("cut -f 2", transcriptome.length.hg19, sep=" "))
tmp2 <- pipe(paste("cut -f 1", transcriptome.length.hg19, sep=" "))
transcriptome.length <- as.list(scan(tmp1,comment.char = "#"))
names(transcriptome.length) = scan(tmp2,comment.char = "#", what="")
close(tmp1)
close(tmp2)

chrmap.file = paste(lane, "chrmap", sep=".")
chrmap<-read.table(paste(stats.dir,chrmap.file,sep="/"))
chrmap.sort1<-chrmap[order(-chrmap$V2),]
chrmap.sort2<-chrmap[order(-chrmap$V3),]
tr.cov = lapply(as.list(chrmap.sort1[,1]),function(key) chrmap.sort1$V2[chrmap.sort1[,1] == key]*readlen/transcriptome.length[[key]])
names(tr.cov) = chrmap.sort1[,1]
tr.cov = t(data.frame(tr.cov))

cols1 = brewer.pal(12, "Paired")[1:12]
cols2 = c("gray","darkred")

chrmap.plot = "chrmap.png"
png(file = paste(dir.html, chrmap.plot, sep ="/"), width = 800, height = 400)
barplot(chrmap.sort1[,2],names.arg=chrmap.sort1[,1], col=cols1,las=3, main="Reads on chromosomes",ylab="Count of uniquely mapped reads")
dev.off()
HTMLInsertGraph(Caption = "Reads distributed on chromosomes",
                GraphFileName = chrmap.plot, Width = 800,file = target)

chrcov.plot = "chrcov.png"
png(file = paste(dir.html, chrcov.plot, sep ="/"), width = 800, height = 400)
barplot(log2(tr.cov[,1]), names.arg=chrmap.sort1[,1], col=cols1, las=3,
        xlab="Chromosome", ylab="log2(Coverage per base in transcriptome)", main="base coverage in transcriptome (exon union)")
dev.off()
HTMLInsertGraph(Caption = "base coverage in transcriptome (exon union)",
                GraphFileName = chrcov.plot, Width = 800,file = target)


#start position distribution
HTML("<br> Unique read start positions on Chromosomes:", file = target)
chrpos.plot = "chrpos.png"
par(cex.lab=2,cex.axis=2)
t.chr.position<-t(as.matrix(chrmap.sort2[,3:4]))
png(file = paste(dir.html, chrpos.plot, sep ="/"), width = 800, height = 400)
barplot(t.chr.position, names.arg=chrmap.sort2[,1], col=cols2, beside=T, las=3, ylab="Count of start positions")
legend(x="topright",legend=c("total start positions","start positions with more than 2 reads"),col=cols2,pch=15, bty="n",cex=2)
dev.off()
HTMLInsertGraph(Caption = "unique start positions on chromosomes",
                GraphFileName = chrpos.plot, Width = 800,file = target)


#RPKM distribution (Reads Per Kilobase per Million of mapped reads)
HTML("<br> RPKM distribution (Reads Per Kilobase per Million of mapped reads)", file = target)
expr.file = paste(lane, "expr", sep=".")
expr <- read.table(paste(stats.dir, expr.file, sep="/"))
RPKM.plot = "RPKM.png"
png(file = paste(dir.html, RPKM.plot, sep ="/"), width = 600, height = 400)
hist(log2(expr$V5[which(expr$V5>0)]*10^9/(mapping.stats[4,2]*2-mapping.stats[8,2])), breaks=200, col="lightblue", xlab="log2(RPKM)", main="")
dev.off()
HTMLInsertGraph(Caption = "RPKM distribution",
                GraphFileName = RPKM.plot, Width = 600,file = target)

#position ~ RPKM
source(paste(src, "smkey.R", sep="/"))
postag.plot = "pos_tag.png"
png(file = paste(dir.html, postag.plot, sep ="/"), width = 600, height = 500)
smkey(log2(expr$V5[which(expr$V5>0)]*10^9/(mapping.stats[4,2]*2-mapping.stats[8,2])),expr$V6[which(expr$V5>0)],xlab="log2(RPKM)",ylab="fraction of covered positions", main="positions vs tags")
abline(v=log2(mean(expr$V5[which(expr$V5>0)]*10^9/(mapping.stats[4,2]*2-mapping.stats[8,2]))), lwd = 2, col="red")
abline(h=mean(expr$V6[which(expr$V5>0)]), lwd=2, col="black")
legend("topleft", legend=c("mean RPKM of expressed transcripts", "mean covered fraction of expressed transcripts"), col=c("red","black"), pch="-", lwd=2, text.col="white", bty="n")
dev.off()
HTMLInsertGraph(Caption = "positions vs tags",
                GraphFileName = postag.plot, Width = 600,file = target)


#read coverage distribution
readcov.plot = "readcov.png"
expr.file = paste(lane, "RefSeq.expr", sep=".")
expr <- read.table(paste(stats.dir, expr.file, sep="/"))
sel.cov.range <- quantile(log2(expr$V7[which(expr[,7] >= 1)]), c(0,1))
ntrans = length(expr$V7[which(expr[,7] >= 1)])
x.cov <- seq(0, sel.cov.range[2], by = 1)
cov.index = findInterval( x.cov, sort(log2(expr$V7[which(expr[,7] >= 1)])), rightmost.closed = TRUE)
cov.above = ntrans - cov.index
y.above = cov.above/ntrans
png(file = paste(dir.html, readcov.plot, sep ="/"), width = 600, height = 500)
plot(sel.cov.range, c(0,1), type = "n", xlab = "log(Number of Reads covered for each expressed RefSeq genes)", ylab = "Fraction of all expressed RefSeq genes", main=paste("Read coverage on expressed RefSeq Genes N=", ntrans, sep=""))
points(x.cov, y.above, type = "l", lty = 1, lwd = 3 ,col = rgb(1,0,0,alpha=0.5))
points(density(log2(expr$V7[which(expr[,7] >= 1)]), n=300), type = "l", lty = 1, lwd = 3, col= rgb(0,0,1,alpha=0.5))
legend("topright", legend = c("cumulative fraction", "density"), col = c(rgb(1,0,0,alpha=0.5),rgb(0,0,1,alpha=0.5)), pch = 19, bty="n")
dev.off()
HTMLInsertGraph(Caption = "read coverage of each expressed RefSeq genes (with at least one read)",
                GraphFileName = readcov.plot, Width = 600,file = target)


#% refgenes contributing
trans = expr$V7[which(expr[,7] >= 1)]
trans.sort = sort(trans, decreasing=T)
ntrans = length(expr$V7[which(expr[,7] >= 1)])
perc = (seq(0,100))*0.01
index.perc = round(ntrans*perc)
cum.perc = rep(0,101)
for (i in 1:101) {
  current = sum(trans.sort[1:index.perc[i]])/sum(trans.sort)
  cum.perc[i] = current
}
cum.perc = cum.perc*100
perc = perc*100
percov.plot = "percov.png"
png(file = paste(dir.html, percov.plot, sep ="/"), width = 600, height = 500)
plot(perc, cum.perc, xlab = "% RefSeq genes contributing", ylab = "% of total counts of all RefSeq genes" ,type = "l", lwd = 4, col=rgb(0,0,1,alpha=0.7))
dev.off()
HTMLInsertGraph(Caption = "Cumulative percentage of total read count, starting with the RefSeq gene with the highest read count",
                GraphFileName = percov.plot, Width = 600,file = target)


#covered fraction distribution
frac = expr$V6[which(expr[,6] > 0)]
fscale = seq(0,1,0.01)
fperc = rep(0,101)
for (i in 1:101){
  fperc[i] = length(which(frac >= fscale[i]))/length(which(frac > 0))
}
frac.hist = hist(expr$V6[which(expr[,6] > 0)], breaks=50, plot=F)
frachist.plot = "frachist.png"
png(file = paste(dir.html, frachist.plot, sep ="/"), width = 1000, height = 500)
layout(matrix(1:2,1,2))
plot(frac.hist$mids, frac.hist$counts, log="y", xlim=c(0,1), pch=20, col="blue", xlab= "fraction of positions with starting reads", ylab="count of refseq genes", main="fraction histogram 1")
plot(fscale, fperc, xlim=c(0,1), ylim=c(0,1), xlab=">= fraction of positions with starting reads", ylab="fraction of refseq genes", main="fraction histogram 2", pch=20, col="darkgreen")
dev.off()
HTMLInsertGraph(Caption = "Fraction of all the positions with starting reads of Refseq genes",
                GraphFileName = frachist.plot, Width = 1000,file = target)


#locus_bias
HTML("<br> Read distribution along transcripts (start site+1000; stop site-1000)", file = target)
lbias.file = paste(lane, "lbias", sep=".")
lbias <- read.table(paste(stats.dir, lbias.file, sep="/"))
lbias.plot = "lbias.png"
png(file = paste(dir.html, lbias.plot, sep ="/"), width = 800, height = 400)
par(mfrow=c(1,2))
plot(lbias[,1],log(lbias[,2]), pch=3, col="red", xlab="downstream the start of a transcript", ylab="log(total counts)", ylim=c(8,15))
plot(-lbias[,1], log(lbias[,4]),pch=4, col="blue", xlab="upstream the end of a transcript", ylab="", ylim=c(8,15))
dev.off()
HTMLInsertGraph(Caption = "locus bias plot",
                GraphFileName = lbias.plot, Width = 800,file = target)


#category distribution
#HTML("<br> Read coverages of transcripts belonging to different categories", file = target)
cate.file = paste(lane, "cate", sep=".")
cate = read.table(paste(stats.dir, cate.file, sep="/"))
#cate.plot = "cate.png"
#png(file = paste(dir.html, cate.plot, sep ="/"), width = 800, height = 800)
#histogram(~log(cate$V2[which(cate[,2] != 0)])|cate$V3[which(cate[,2] != 0)], breaks=50, scales=list(x=list(log = 2)), xlab="log2(transcript_coverage)")
#dev.off()
#HTMLInsertGraph(Caption = "coverage distribution of different kinds of transcripts",
#                GraphFileName = cate.plot, Width = 800,file = target)


HTML("<br> Boxplot of log2-scale coverage for transcripts in each category. Red dot: mean / Black dot : median.", file = target)
trc.plot = "trc.png"
png(file = paste(dir.html, trc.plot, sep ="/"), width = 600, height = 800)
trc<-data.frame(x=cate$V2, y=cate$V3)
bwplot(y~(x+1), scales=list(x=list(log = 2)), data=trc, xlab="transcript coverage", panel=function(x,y,...) {
  panel.bwplot(x, y, ...)
  panel.points(x=log2(unlist(as.list(by(trc$x,trc$y,function(z) mean(z))))+1), y=1:12, col="red",pch=19)
})
dev.off()
HTMLInsertGraph(Caption = "coverage sum of transcripts belonging to different categories",
                GraphFileName = trc.plot, Width = 600,file = target)


#end
HTMLEndFile(target)

