#parse arguments

args=(commandArgs(TRUE))

if(length(args)==0){
    print("No arguments supplied.")
}else{
    for(i in 1:length(args)){
      eval(parse(text=args[[i]]))
    }
}
        

setwd(dir)
library(RColorBrewer)
qual.av = read.table(qual.av.file, col.names = "average.quality" )
qual.med = read.table(qual.med.file, col.names = "median.quality" )
quality.profile = cbind(qual.av, qual.med)
cols.qual = brewer.pal(3, "Set1")[1:2]
png(filename = qual.plot, width = 1200, height = 800, res = 150)
matplot(quality.profile, type = "l" , col = cols.qual , lty = 1, lwd = 2, main = "Quality values", ylab = "Phred Quality Score", xlab = "position" )
legend(x = "topright", col = cols.qual, legend = names(quality.profile), pch = 19)
dev.off()

