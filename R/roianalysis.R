#' Get ROIs
#'
#' This function reads in a set of ImageJ ROIs and returns a polygon data frame.
#' @param folder character describing the system file path for the folder containing ROIs (can't be a zip-file) so unzip before.
#' @keywords ROI
#' @export
#' @examples
#' cytosol.folder<-system.file('data/roi/RoiSet_cytosol', package='leeplyr')
#' nuclei.folder<-system.file('data/roi/RoiSet_cytosol', package='leeplyr')
#' cytosol<-get.rois(cytosol.folder)
#' nuclei<-get.rois(nuclei.folder)
get.rois<-function(folder){
  roi.cytosol<-dir(folder, full.names=TRUE)
  
  polygon<-read.ijroi(roi.cytosol[1])
  poly<-data.frame(x = polygon$coords[,1], y = polygon$coords[,2], id = rep(1, nrow(polygon$coords)))
  for(i in 2:length(roi.cytosol)){
    polygon<-read.ijroi(roi.cytosol[i])
    poly.tmp<-data.frame(x = polygon$coords[,1], y = polygon$coords[,2], id = rep(i, nrow(polygon$coords)))
    poly<-rbind(poly, poly.tmp)
  }
  
  cat(paste('Number of ROIs loaded:', length(unique(poly$id)), '\n'))
  return(poly)
  
}

#' Plot polygon
#'
#' Draws polygon for a polygon data frame.
#' @param dataframe data frame, a polygon data frame with list items x, y, and id.
#' @param col character, fill color for polygon. Default, NA, is to leave polygons unfilled.
#' @param border character, the color to draw the border. The default is 'black'.
#' @keywords polygon
#' @export
#' @examples
#' cytosol<-get.rois(cytosol.folder)
#' plot(cytosol[,1:2], col = 0)
#' plot.polygon(cytosol, col='pink')
plot.polygon<-function(dataframe, col = NA, border = 'black'){
  lapply(unique(dataframe$id), function(x){polygon(dataframe$x[dataframe$id==x], dataframe$y[dataframe$id==x], col=col, border = border)})
}


#' Map transcripts to ROI 
#'
#' Maps each individual transcripts to one or more region of interests (ROIs).
#' @param transcripts FISSEQ readr output tsv file path.
#' @param roi.folder character vector, each element is the location of folder containing ImageJ roi files.
#' @param roi.labels character vector, .
#' @param OR.cutoff numeric, absolute OR value on log2 coordinates to use as cut-off. Default is 5.
#' @param p.value.cutoff numeric, p-value on negative log10 coordinates to use as cut-off. Default is 20.
#' @param plot boolean, plot output graph. Default is TRUE.
#' @param roi.col character vector, Color of regions of interest. Default is c('#C8EAEB', '#FCDAD1').
#' @param transcript.col character, color or transcripts in plot. Default is 'black'.
#' @param cex numeric, Size of transcripts on plot showing alll transcripts. Default is 0.2. 
#' @export
#' @examples
#' cytosol.folder<-system.file('roi/RoiSet_cytosol', package='leeplyr')
#' nuclei.folder<-system.file('roi/RoiSet_nuclei', package='leeplyr')
#' transcripts<-system.file('fisseq/res_001_001FISSEQ.out', package='leeplyr')
#' map.to.roi(transcripts, roi.folder = c(cytosol.folder, nuclei.folder), roi.labels =  c('cytoplasm', 'nucleus'), OR.cutoff = 5, p.value.cutoff = 20)
map.to.roi<-function(transcripts, roi.folder = c(cytosol.folder, nuclei.folder), roi.labels =  c('cytoplasm', 'nucleus'), OR.cutoff = 5, p.value.cutoff = 20, plot = TRUE, roi.col = c('#C8EAEB', '#FCDAD1'), transcript.col = 'black', cex = 0.2){
  #read ROIs
  rois<-list()
  ylim<-numeric()
  xlim<-numeric()
  for(i in seq_along(roi.folder) ){
    rois[[i]]<-get.rois(roi.folder[i])
    xlim<-append(xlim, range(rois[[i]]$x) )
    ylim<-append(ylim, range(rois[[i]]$y) )
  }
  #read transcripts
  if(class(transcripts)=='character'){
    amplicons <- read_delim(transcripts, "\t", escape_double = FALSE, trim_ws = TRUE)
  }else{
    amplicons <- transcripts
  }
  xlim<-append(xlim, range(amplicons$centroid_y) )
  ylim<-append(ylim, range(amplicons$centroid_x) )
  ylim<-range(ylim)
  xlim<-range(xlim)
  cat('Assigning transcripts to polygons... \n')
  #check in polygon
  check.inside<-list()
  for(i in seq_along(roi.folder) ){
    check.inside[[i]]<-rep(0, nrow(amplicons))
    for(j in unique(rois[[i]]$id) ){
      check.inside[[i]]<-check.inside[[i]]+point.in.polygon(amplicons$centroid_y, amplicons$centroid_x, rois[[i]]$x[which(rois[[i]]$id==j)], rois[[i]]$y[which(rois[[i]]$id==j)])
    }
  }
  
  
  roi.output<-data.frame(do.call("cbind", check.inside))
  names(roi.output)<-paste0('ROI.', roi.labels)  
  roi.output<-cbind( data.frame(gene.symbol = amplicons$string_gene_symbols, strand = amplicons$strand, x = amplicons$centroid_y, y = amplicons$centroid_x), roi.output)
  
  inside.roi<-list()
  for(i in seq_along(roi.folder)[-length(roi.folder)] )
    inside.roi[[i]]<-( (roi.output[,4+i] > 0)&(roi.output[,4+i+1] == 0) )
  inside.roi[[i+1]]<-(roi.output[,4+i+1] > 0)
  inside.roi<-do.call("cbind", inside.roi)
  
  gene.symbols<-unique(roi.output$gene.symbol)
  OR<-rep(NA, length(gene.symbols))
  p.value<-rep(NA, length(gene.symbols))
  for(k in unique(roi.output$gene.symbol)){
    N<-table(roi.output$gene.symbol == k, 1-inside.roi[,2])
    if(!is.null(dim(N))){
        OR[which(gene.symbols == k)] <- fisher.test(N)$estimate
        p.value[which(gene.symbols == k)] <- fisher.test(N)$p.value
    }
  }
  
  p.value[which(p.value==0)]<-as.numeric( noquote(unlist(format(.Machine)))[3] )
  
  interesting<-which(abs(log2(OR))> OR.cutoff & -log10(p.value)> p.value.cutoff)

  print.output<-data.frame(gene = gene.symbols, OR, p.value)[interesting, ]
  print.output$roi<- log2(print.output$OR)>0
  print.output$OR<-log2(print.output$OR)
  print.output$p.value<- -log10( print.output$p.value)
  print.output$roi<-rev(roi.labels)[print.output$roi+1]  
  print.output<-print.output[order(print.output$OR),]
  print.output<-print.output[is.finite(print.output$OR),]

  #plot results
  if(plot){
    par(mfrow=c(1,3))
    plot(rois[[1]]$x, rois[[1]]$y, col=0, asp=1, ylim=rev(ylim), xlim = xlim, ylab='', xlab='')
    for(i in seq_along(roi.folder)){
      plot.polygon(rois[[i]], col = roi.col[i])
    }
    points(amplicons$centroid_y, amplicons$centroid_x, pch = 16, cex = cex, col = transcript.col)
    
    
    plot(rois[[1]]$x, rois[[1]]$y, col=0, asp=1, ylim=rev(ylim), xlim = xlim, ylab='', xlab='')
    for(i in seq_along(roi.folder)){
      plot.polygon(rois[[i]], col = roi.col[i])
    }
    max.roi1<-print.output$gene[1]
    max.roi2<-print.output$gene[2]
    
    points(amplicons$centroid_y, amplicons$centroid_x, pch = 16, cex = 0.7, col = c(NA, 'purple', 'green3', 'red')[(amplicons$string_gene_symbols == max.roi1)*3 + (amplicons$string_gene_symbols == max.roi2)*2 + (amplicons$string_gene_symbols== 'RNA18S5,RNA45S5')*1 + 1] )
    legend('bottomleft', c('18S', as.character(max.roi2), as.character(max.roi1)), pch=16, col=c('purple', 'green3', 'red'), bg='white')
    
    par(yaxs='i', xaxs='r')
    ylim<-rep(max(abs(range(log2(OR), na.rm=TRUE, finite = TRUE))), 2)*c(-1,1)
    plot(-log10(p.value), log2(OR), col=0, ylim=ylim, ylab=expression(log[2]*(OR)), xlab=expression(-log[10]*(P-value)) )
    polygon(c(-20, max(-log10(p.value), na.rm=TRUE)+20,  max(-log10(p.value), na.rm=TRUE)+20, -20), c(0, 0, rep(ylim[1]*1.1, 2)), col='#FCDAD1' )
    polygon(c(-20, max(-log10(p.value), na.rm=TRUE)+20,  max(-log10(p.value), na.rm=TRUE)+20, -20), c(0, 0, -rep(ylim[1]*1.1, 2)), col='#C8EAEB' )
    
    points(-log10(p.value), log2(OR), pch=16, col=rgb(0,0,0,0.2))
    
    
    for(i in interesting){
      symb<-gene.symbols[i]
        if(symb == 'RNA18S5,RNA45S5')
          symb<-'18S'
        if(symb == 'RNA45S5,RNA28S5')
          symb<-'28S'
      text(-log10(p.value)[i], log2(OR)[i], symb, pos = 2, offset = 0.75)
    }
    
   
  }
  

  print(print.output)

  return(roi.output)
  
  
}
  
  
  
  
  
  

