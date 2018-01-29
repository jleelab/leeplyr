# plot(log10(rev(seq_along(table(mapped$gene.symbol)))) , log10(sort(as.integer(table(mapped$gene.symbol)) )), ylab='Amplicons per gene', xlab='Gene ranked by expression', axes=F, type='s')
# points(log10(rev(seq_along(table(mapped$gene.symbol)))), log10(sort(as.integer(table(mapped$gene.symbol)) )), pch=16, cex=0.5)

# axis(2, at=log10(c(seq(1,10), seq(20, 100, by=10), seq(200, 1000, by=100), 2000)), labels=FALSE )
# axis(2, at=log10(c(1,10,100,1000)), labels=c(1,10,100,1000) , las=1)
# axis(1, at=log10(c(seq(1,10), seq(20, 100, by=10), seq(200, 1000, by=100), 2000)), labels=FALSE)
# axis(1, at=log10(c(1,10,100,1000)), labels=c(1,10,100,1000) )
# box()

# text(log10(c(1:5)), log10(rev(sort(as.integer(table(mapped$gene.symbol))))[1:5]), c("45S", "18S", "PTOV1-AS2", "LOC105372073", "EZR") , pos=4, offset=0.4, cex=0.5)
# text(log10(which(names(sort(table(mapped$gene.symbol), decreasing=TRUE))=='RBFOX3')), log10(sort(table(mapped$gene.symbol), decreasing=TRUE)[55]), names(table(mapped$gene.symbol)[2055]), pos=4, offset=0.4, cex=0.5 )

