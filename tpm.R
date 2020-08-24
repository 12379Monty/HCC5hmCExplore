
 #suppressMessages(reqduire(spatstat))

 BCV_mtx <- do.call('cbind', lapply(unique(sampDescA$group), 
 function(GRP) {
  GRP_dgel <- edgeR::DGEList(counts=featureCountsA[, sampDescA$group==GRP])
  GRP_dgel <- edgeR::estimateDisp(GRP_dgel)

  sqrt(GRP_dgel$tagwise.dispersion)
  }))

 colnames(BCV_mtx) <- unique(sampDescA$group)


 # Plot
 plot(spatstat::CDF(density(BCV_mtx[,1])), col=groupCol[colnames(BCV_mtx)[1]], lwd=2, ylab='Prob(BCV<x)',
   xlim=c(0, 0.5))
 for(JJ in 2:ncol(BCV_mtx))
 plot(spatstat::CDF(density(BCV_mtx[,JJ])), col=groupCol[colnames(BCV_mtx)[JJ]], lwd=2, add=T, xlim=c(0, 2.0))

 legend('bottomright', legend=names(groupCol), col=groupCol, lwd=2)
  
 BCV_90perc_vec <- apply(BCV_mtx,2,quantile, prob=0.90)
 BCV_50perc_vec <- apply(BCV_mtx,2,quantile, prob=0.50)
 for(JJ in 1:length(BCV_90perc_vec)) 
 rug(BCV_90perc_vec[JJ], lwd=2, col=groupCol[names(BCV_90perc_vec)[JJ]])

