filteredCountsA_voom_tfit <- limma::treat(filteredCountsA_voom_fit, lfc=log2(1.10))
filteredCountsA_voom_tfit_dt <- limma::decideTests(filteredCountsA_voom_tfit)

cat("10% FC Gene Identification Summary - voom, adjust.method = BH, p.value = 0.05:\n")
summary(filteredCountsA_voom_tfit_dt)

# log-fold-change vs ave-expr
limma::plotMD(filteredCountsA_voom_efit,
 ylim = c(-0.5, 0.5),
 column='HCCvsControl',
 status=filteredCountsA_voom_tfit_dt[,'HCCvsControl'],
 hl.pch = 16, hl.col = c("lightblue", "pink"), hl.cex = .5,
 bg.pch = 16, bg.col = "grey", bg.cex = 0.5,
 main = '',
 xlab = paste0(
    "Average log-expression: IQR=",
    paste(round(quantile(filteredCountsA_voom_efit$Amean, prob = c(1, 3) / 4), 2),
      collapse = ", "
    )
  ),
  ylab = paste0(
    "log-fold-change: IQR=",
    paste(round(quantile(filteredCountsA_voom_efit$coefficients[, 'HCCvsControl'], prob = c(1, 3) / 4), 2),
      collapse = ", "
    )
  ),
  legend = F
)
abline(h = 0, col = "black")
rug(quantile(filteredCountsA_voom_efit$coefficients[, 'HCCvsControl'], prob = c(1, 2, 3) / 4),
  col = "purple",
  ticksize = .03, side = 2, lwd = 2
)
rug(quantile(filteredCountsA_voom_efit$Amean, prob = c(1, 2, 3) / 4),
  col = "purple",
  ticksize = .03, side = 1, lwd = 2
)


