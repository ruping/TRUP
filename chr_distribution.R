#load local R library
.libPaths(c(.libPaths(), "/scratch/ngsvin/R_libraries/"))

require(R2HTML, quiet =  TRUE)
library(lattice)
library(geneplotter)
library(RColorBrewer)
cols1 = brewer.pal(9, "Pastel1")
cols2 = brewer.pal(3, "Paired")[1:2]
cols3 = brewer.pal(3, "Paired")
cols4 = brewer.pal(3, "Accent")
cols5 = brewer.pal(3, "Set1")

#accept command line arguments
args=(commandArgs(TRUE))
if(length(args)==0){
    print("No arguments supplied.")
    ##supply default values
}else{
    for(i in 1:length(args)){
      eval(parse(text=args[[i]]))
    }
}

setwd(path)
report_html = paste(lane, "report", sep=".")
target <- HTMLInitFile(path, report_html, BackGroundColor="FFFFFF", Title = lane)
HTML(as.title(paste("Mapping report for", lane)), file=target)

#read transcriptome length and put it into a list
tmp1 <- pipe("cut -f 2 /scratch/ngsvin/RNA-seq/MPI-NF/ANNOTATION/transcriptome_length_HG18")
tmp2 <- pipe("cut -f 1 /scratch/ngsvin/RNA-seq/MPI-NF/ANNOTATION/transcriptome_length_HG18")
transcriptome.length <- as.list(scan(tmp1,comment.char = "#"))
names(transcriptome.length) = scan(tmp2,comment.char = "#", what="")
close(tmp1)
close(tmp2)

chr_map_stat = paste(lane,"chr_map_stat",sep=".")
d<-read.table(chr_map_stat)
d$V1 = paste("chr", d$V1, sep="")
dsort1<-d[order(-d$V2),]
dsort2<-d[order(-d$V3),]


#generate the transcriptme coverage
tr.cov = lapply(as.list(dsort1$V1),function(key) dsort1$V2[dsort1$V1 == key]*80/transcriptome.length[[key]])
names(tr.cov) = dsort1$V1
tr.cov = t(data.frame(tr.cov))

chr_counts=paste(lane,"chr_counts.png",sep=".")
png(filename=chr_counts,width=1200,height=1200,res=110)
layout(matrix(1:2, 2, 1))
barplot(dsort1$V2,names.arg=dsort1$V1,col=cols1,las=3, xlab="",ylab="Count of UMR")
barplot(log2(tr.cov[,1]),names.arg=dsort1$V1,col=cols1,las=3, xlab="Chromosome",ylab="log2(Coverage per base in transcriptome)")
dev.off()


t.chr.position<-t(as.matrix(dsort2[,3:4]))
chr_position=paste(lane,"chr_position.png",sep=".")
png(filename=chr_position,width=1200,height=600,res=110)
barplot(t.chr.position, names.arg=dsort2$V1 , col=cols2, beside=T, las=3, xlab="Chromosome", ylab="Count of Positions")
legend(x="topright",legend=c("total positions","piling up positions"),col=cols2,pch=15, bty="n")
dev.off()


tot_stat=paste(lane,"tot_stat",sep=".")
tot<-read.table(tot_stat, header=T)
UMR_whole = tot$UMR-tot$UMR_partial
cols3 = brewer.pal(4, "Paired")
pie<-c(UMR_whole,tot[1,2],tot[1,3],tot[1,4])
car_labels <- round(pie/sum(pie) * 100, 1)
car_labels <- paste(c("Uniquely Mappable Reads (completely)", "Uniquely Mappable Reads (Partially)" ,"Multi-location Mappable Reads","Unmappable Reads"), car_labels, sep = "-")
map_pie= paste(lane,"map_pie.png",sep=".")
png(filename=map_pie,width=1200,height=600,res=110)
pie(pie,labels=car_labels,col=cols3,clockwise=T,border=F)
#legend(1, 0.5, c("Uniquely Mappable Reads (completely)", "Uniquely Mappable Reads (Partially)" ,"Multi-location Mappable Reads","Unmappable Reads"), cex=1, fill=cols3)
dev.off()

#data.frame of general statistics of the mapping
options(digits = 7)
pN = function(x,...) {return(prettyNum(x, big.mark = ",", ...))}
coverage_in_transcriptome=paste(lane,"coverage_in_transcriptome",sep=".")
d.cov.tr<-read.table(coverage_in_transcriptome)
mean.cov = d.cov.tr[1,4]
N_reads_sequenced = tot$UMR+tot$MMR+tot$unmapped
row.infos = c("N reads sequenced",
              "N reads mapped uniquely (UMR)",
              "N reads mapped uniquely (completely)",
              "N reads mapped uniquely (partially)",
              "N reads mapped to multiple locations",
              "Total number of positions of UMR",
              "Total number of piling up positions",
              "mean coverage in transcriptome"
              )
row.values = c(pN(N_reads_sequenced),
               pN(tot$UMR),
               pN(UMR_whole),
               pN(tot$UMR_partial),
               pN(tot$MMR),
               pN(tot$total_pos),
               pN(tot$pileup_pos),
               pN(mean.cov,digits=3)
              )
lane.summary = data.frame(Property = row.infos,value=row.values);


UMR_in_gene_ratio = paste(lane,"UMR_in_gene_ratio",sep=".")
d<-read.table(UMR_in_gene_ratio)
dsort<-d[order(-d$V2),]
t.ratio<-t(as.matrix(dsort[,2:4]))
UMR_in_gene_ratio_png = paste(lane,"UMR_in_gene_ratio.png",sep=".")
png(filename=UMR_in_gene_ratio_png,width=1200,height=600,res=110)
barplot(t.ratio, names.arg=dsort$V1 , col=cols4, beside=T, las=3, xlab="Chromosome", ylab="Loci with no piled up reads in a certion region", ylim = c(0,1))
legend(x="toprigh",legend=c("in mRNA","in exon","in CDS"),col=cols4,pch=15,bty="n")
dev.off()

#generate the coverage stuff
coverage_in_transcriptome_perbase = paste(lane, "coverage_in_transcriptome.perbase", sep=".")
cov <- scan(pipe(paste("cut -f 3", coverage_in_transcriptome_perbase, sep=" ")),comment.char = "#")
sel.cov.range <- quantile(cov, c(0,0.99))
nbases = length(cov)
x.cov <- seq(0, sel.cov.range[2], by = 1)
cov.index = findInterval( x.cov, sort(cov), rightmost.closed = TRUE)
cov.above = nbases - cov.index
y.above = cov.above/nbases
coverage_in_transcriptome_png = paste(lane, "coverage_in_transcriptome.png", sep=".")
png(filename=coverage_in_transcriptome_png,width=1200,height=600,res=110)
plot(sel.cov.range, c(0,1), type = "n", xlab = "Coverage", ylab = "Cumulative fraction or density")
points(x.cov, y.above, type = "l", lty = 1, lwd = 2 ,col = cols5[1])
points(density(cov, n=300000), type = "l", lty = 1, lwd = 2, col=cols5[2])
legend(x = "topright", legend = c("cumulative fraction", "density"), col = cols5, pch = 19, bty="n")
dev.off()


#the base composition stuff
base_distribution = paste(lane, "base_distribution", sep=".")
base.d<-read.table(base_distribution)
base_distribution_png = paste(lane, "base_distribution.png", sep=".")
png(filename=base_distribution_png,width=1200,height=600,res=110)
plot(base.d$V1,base.d$V2,ylim=c(0.1,0.4),type="l",lwd=3,xaxt="n",yaxt="n",xlab="position from 5'end of the UMR", ylab="base composition",col="red")
lines(base.d$V1,base.d$V3,type="l",lwd=3,col="green")
lines(base.d$V1,base.d$V4,type="l",lwd=3,col="orange")
lines(base.d$V1,base.d$V5,type="l",lwd=3,col="blue")
axis(side=2,c(0.1,0.25,0.4))
axis(side=1,c(1,20,100),labels=c("-20","1","80"))
legend(x="bottomright", legend=c("A","C","G","T"),col=c("red","green","orange","blue"),lwd=3,bty="n")
dev.off()


HTML("<br> General statistics on the mapping:", file = target)
HTML(lane.summary, innerBorder = 2)

qual.av.file = paste(lane, "quality.average", sep=".")
qual.med.file = paste(lane, "quality.median", sep=".")
qual.av = read.table(qual.av.file, col.names = "average.quality" )
qual.med = read.table(qual.med.file, col.names = "median.quality" )
quality.profile = cbind(qual.av, qual.med)
cols.qual = brewer.pal(3, "Set1")[1:2]

qual.plot = paste(lane,"basequality.png",sep=".")
png(filename = qual.plot, width = 1200, height = 600,res=110)
matplot(quality.profile, type = "l" , col = cols.qual , lty = 1, lwd = 2, main = "Quality values", ylab = "Phred Quality Score", xlab = "position" )
legend(x = "topright", col = cols.qual, legend = names(quality.profile), pch = 19)
dev.off()

#expression stuff
SRPKM = paste(lane, "SRPKM", sep=".")
SRP<-read.table(SRPKM)
log_RPKM = paste(lane, "log_RPKM.png", sep=".")
png(filename = log_RPKM, width = 1200, height = 600,res=110)
hist(log(SRP$V5),breaks=200, col=cols2[2],xlab="log(RPKM)",main="")
dev.off()
position_tags = paste(lane, "position_tags.png", sep=".")
png(filename = position_tags, width = 1200, height = 600, res=110)
smoothScatter(log(SRP$V5),SRP$V6,xlab="log(RPKM)",ylab="covered position fraction")
dev.off()

HTMLInsertGraph(Caption = "base quality of the reads", GraphFileName = qual.plot, Width = 900, file = target)
HTMLInsertGraph(Caption = "mapping pie", GraphFileName = map_pie, Width = 900, file = target)
HTMLInsertGraph(Caption = "base distribution of UMR", GraphFileName = base_distribution_png, Width = 900, file = target)
HTMLInsertGraph(Caption = "chromosome distribution of UMR and coverage in transcriptome", GraphFileName = chr_counts, Width = 900, file = target)
HTMLInsertGraph(Caption = "Coverage density and Cumulative coverage", GraphFileName = coverage_in_transcriptome_png, Width = 900, file = target)
HTMLInsertGraph(Caption = "position information", GraphFileName = chr_position, Width = 900, file = target)
HTMLInsertGraph(Caption = "UMR in gene ratio", GraphFileName = UMR_in_gene_ratio_png, Width = 900, file = target)
HTMLInsertGraph(Caption = "Reads Per Kilo base of gene model per Million mapped reads (log RPKM)", GraphFileName = log_RPKM, Width = 900, file = target)
HTMLInsertGraph(Caption = "Covered position fraction vs RPKM", GraphFileName = position_tags, Width = 900, file = target)


HTMLEndFile(target)
