# Preprocessing {#preproc}

<!--
 FN <- 'tmp'
 # Shotcuts for knitting and redering while in R session (Invoke interactive R from R/Scripts folder)
 kk <- function(n='') knitr::knit2html(paste("t", n, sep=''), envir=globalenv(),
       output=paste(FN,".html", sep=''))

 rr <- function(n='') rmarkdown::render(paste("t", n, sep=''), envir=globalenv(),
       output_file=paste(FN,".html", sep='')) ##, output_dir='Scripts')

 bb <- function(n='') browseURL(paste(FN,".html", sep=''))

 # The usual shotcuts
 zz <- function(n='') source(paste("t", n, sep=''))
-->



```{r setup-preproc, include=F}
   # file rmarkdown file management options: cache, figures
 figures_DIR <- file.path('Static', 'figures/')
 suppressMessages(dir.create(figures_DIR, recursive=T))
 knitr::opts_chunk$set(fig.path=paste0(figures_DIR))
```

## Load the data

<!-- THIS ENSURES NO EVALUATION TAKES PLACE BY DEFAULT -->
<!-- TO TURN ON, SET eval=T                            -->
```{r chunk-options, include=FALSE, eval=F}
library("knitr")
opts_chunk$set(eval = FALSE)
```

<!-- Add base libraries -->
```{r libraries, include=FALSE, eval=T}
library("magrittr")
```


The data that are available from NCBI GEO
[Series GSE112679](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE112679)
can be conveniently accessed through an R data package.
Attaching the GSE112679 package makes the count data tables 
available as well as a gene annotation table and a sample description table.
See [GSE112679 R Data Package page](https://12379monty.github.io/GSE112679/).
For the Cai et al. [@Cai:2019aa] model fitting and analysis, samples were separated into
`Train` and `Val-1` subsets.  `Val-2` was an external validation set.

```{r loadData, cache=F}

if (!("GSE112679" %in% rownames(installed.packages()))) {
  if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
  }
  devtools::install_github("12379Monty/GSE112679")
}
library(GSE112679)
sampDesc$DxStage <- with(sampDesc, ifelse(outcome=='HCC', 
   paste0(outcome,':',stage), outcome))

with(
  sampDesc %>% dplyr::filter(sampType == "blood"),
  knitr::kable(table(DxStage, trainValGroup, exclude = NULL),
   caption="GSE112679 Samples by Dx Group and Subset") %>% 
  kableExtra::kable_styling(full_width = F)
)

```

For this analysis, we will consider early stage cancer samples
and healthy or benign samples from the `Train` or `Val-1` subsets.
The appropriate outcome variable will be renamed or aliased `group`

```{r subsetSamples, cache=T, cache.vars=c('sampDescA','groupCol'), echo=F}

# Use suffix 'A' for Analysis samples
sampDescA <-
  sampDesc %>%
  dplyr::filter(sampType == "blood" &
    (trainValGroup %in% c("Train", "Val-1")) &
    ((outcome2 == "BenignHealthy") |
      (outcome2 == "HCC" & stage == "Early"))) %>%
  dplyr::rename(group = outcome2) %>%
  dplyr::arrange(group, sampID)
# Recode group
sampDescA$group <- with(
  sampDescA,
  ifelse(group == "BenignHealthy", "Control", group)
)
# set groupCol for later
groupCol <- c("#F3C300", "#875692")
names(groupCol) <- unique(sampDescA$group)

with(sampDescA, 
 knitr::kable(table(group, exclude = NULL),
  caption="Samples used in this analysis") %>%
  kableExtra::kable_styling(full_width = F)
)

```


The features are counts of reads captured by chemical labeling, and indicate
the level of 5-hydroxymethylcytosines within each gene body.  Cai et al. (2019),
Li et al. (2017) and Song et al. (2017) [@Cai:2019aa;@Li:2017aa;@Song:2017aa]
all analyze 5hmC gene body counts using standard RNA-Seq methodologies, and we will
do the same here.  

Note that before conducting any substantive analyses, the data would normally
be very carefully examined for any sign of quality variation between groups
of samples.  This analysis would integrate sample meta data - where and when were
the blood samples collected - as well as library preparation and sequencing metrics
in order to detect any sign of processing artifacts that may be present in the dataset.
This is particularly important when dealing with blood samples as variable
DNA quality degradation is a well known challenge that is encountered when dealing with
such samples [@Huang:2017aa].  Although blood specimen handling protocols can be 
put in place to minimize quality variation [@Permenter:2015aa], variability
can never be completely eradicated, especially in the context of blood samples
collected by different groups, working in different environments.  The problem
of variable DNA quality becomes paricularly pernicuous when it is compounded
with a confounding factor that sneaks in when the control sample collection
events are separated in time and space from the cancer sample collection events;
an all too common occurence.  

As proper data QC requires an intimate familiarity with the details of
data collection and processing, such a task cannot be untertaken here.
We will simply run a *minimal set of QC sanity checks* to make sure that
there are no apparent systematic effects in the data.


```{r getfeatures, cache=T, cache.vars=c('featureCountsA'), echo=F}

featureCountsA <- cbind(
  Train_featureCount,
  Val1_featureCount,
  Val2_featureCount
)[, rownames(sampDescA)]

```

We first look at coverage - make sure there isn't too much disparity of coverage 
across samples. To detect shared variability, samples can be annotated and ordered
according to sample features that may be linked to sample batch processing.  Here we 
the samples have been ordered by group and sample id (an alias of geoAcc).

```{r lcpmBxp, cache=T, fig.height=4, fig.width=10, fig.cap='Sample log2 count boxplots', echo=F}

par(mar = c(1, 3, 2, 1))
boxplot(log2(featureCountsA + 1),
  ylim = c(3, 11), ylab='log2 Count',
  staplewex = 0,       # remove horizontal whisker lines
  staplecol = "white", # just to be totally sure :)
  outline = F,         # remove outlying points
  whisklty = 0,        # remove vertical whisker lines
  las = 2, horizontal = F, xaxt = "n",
  border = groupCol[sampDescA$group]
)
legend("top", legend = names(groupCol), text.col = groupCol, 
  ncol = 2, bty = "n")
# Add reference lines
SampleMedian <- apply(log2(featureCountsA + 1), 2, median)
abline(h = median(SampleMedian), col = "grey")
axis(side = 4, at = round(median(SampleMedian), 1), 
  las = 2, col = "grey", line = -1, tick = F)

```

<!--
Coverage level looks fairly comparable across samples.  It is sometimes helpful to
keep track of the actual coverage which can be adequetely tracked by distribution
quantiles.
-->

```{r quantCountsA, echo=F}

featureCountsA_quant <- apply(featureCountsA, 2, function(CC) {
  c(quantile(CC, prob = c(.15, (1:3) / 4)), totCovM = sum(CC) / 1e6)
})

featureCountsA_quant2 <- apply(featureCountsA_quant, 1, function(RR) {
  quantile(RR, prob = (1:3) / 4)
})

knitr::kable(featureCountsA_quant2,
  digits = 1,
  caption = paste(
    "Coverage Summary - Columns are sample coverage quantiles and total coverage",
    "\nRows are quartiles across samples"
  )
) %>% kableExtra::kable_styling(full_width = F)

```


From this table, we see that 25% of the samples have total coverage exceeding
`r round(featureCountsA_quant2["75%", "totCovM"],1)`M reads, 25% of samples
have a 15 percentile of coverage lower than
`r featureCountsA_quant2["25%", "15%"]`, etc.  


<!-- SKIP
We next look at relative log representation (RLR) (in the context of measuring the density of 
5hmC marks in genes, we refer to `representation` as opposed to `expression`; the
two can be used interchangibly) -
make sure the shapes of the distributions are not widely different.
-->
```{r rlr, cache=T, cache.vars='lcpm_mtx', fig.height=4, fig.width=10, fig.cap='Sample RLR', eval=T, echo=F,include=F}

lcpm_mtx <- edgeR::cpm(featureCountsA, log = T)
median_vec <- apply(lcpm_mtx, 1, median)
RLR_mtx <- sweep(lcpm_mtx, 1, median_vec, "-")

par(mar = c(1, 3, 2, 1))
boxplot(RLR_mtx,
  xlab = "", ylab='Relative Log Representation', ylim = c(-.6, .6),
  staplewex = 0, # remove horizontal whisker lines
  staplecol = "white", # just to be totally sure :)
  outline = F, # remove outlying points
  whisklty = 0, # remove vertical whisker lines
  las = 2, horizontal = F, xaxt = "n",
  border = groupCol[sampDescA$group]
)
legend("top", legend = names(groupCol), 
  text.col = groupCol, ncol = 2, bty = "n")
# Add group Q1, Q3
for (GRP in unique(sampDescA$group)) {
  group_ndx <- which(sampDescA$group == GRP)
  group_Q1Q3_mtx <- apply(RLR_mtx[, group_ndx], 2, 
     quantile, prob = c(.25, .75))
  abline(h = apply(group_Q1Q3_mtx, 1, median), 
     col = groupCol[GRP], lwd = 2)
}


```

<!-- SKIPPED
We note that the HCC samples have slightly more variable coverage distribution.
A few samples are quite different.
-->

## Differential representation analysis {#dra}


In the remainder of this section, we will process the data and
perform differential expression analysis as outlined in 
Law et al. (2018) [@Law:2018aa].   The main analysis steps are: 

* remove lowly expressed genes
* normalize gene expression distributions
* remove heteroscedascity
* fit linear models and examine DE results


It is good practice to perform this differential expression analysis prior to 
fitting models to get an idea of how difficult it will be to discriminate 
between samples belonging to the different subgroups. The pipeline
outlined in Law et al. (2018) [@Law:2018aa] also provides some
basic quality assessment opportunities.


### Remove lowly expressed genes {-}

Genes that are not expressed at a biologically 
meaningful level in any condition should be discarded to reduce the 
subset of genes to those that are of interest, and to reduce the number of tests 
carried out downstream when looking at differential expression.  Carrying
un-informative genes may also be a hindrance to classification and other
downtream analyses.  

To determine a sensible threshold we can begin by examining the shapes of the distributions.

```{r densityLcpm, fig.height=4, fig.width=10, fig.cap='Sample $log_2$ CPM densities', eval=T, echo=F}

par(mar = c(4, 3, 2, 1))
plot(density(lcpm_mtx[, 1]),
  col = groupCol[sampDescA$group[1]],
  lwd = 2, ylim = c(0, .25), las = 2, main = "", xlab = "log2 CPM"
)
abline(v = 0, col = 3)
# After verifying no outliers, can plot a random subset 
for (JJ in sample(2:ncol(lcpm_mtx), size = 100)) {
  den <- density(lcpm_mtx[, JJ])
  lines(den$x, den$y, col = groupCol[sampDescA$group[JJ]], lwd = 2)
} # for(JJ
legend("topright", legend = names(groupCol), 
  text.col = groupCol, bty = "n")

```

`r  LibSizeSum <- summary( colSums(featureCountsA) / 1e6 )`
`r CPM_THR <- 3; SAMP_THR <- 25` 

As is typically the case with RNA-Seq data, we notice many weakly represented genes 
in this dataset.  A cpm value of 1 appears to adequatly separate 
the expressed from the un-expressed genes, but we will be slightly more strict here
and require a CPM threshold of `r CPM_THR` .  Using a nominal CPM value of 
`r CPM_THR`, genes are deeemed to be `represented` if their expression is 
above this threshold, and not represented otherwise. 
For this analysis we will require that genes be `represented` in at least 
`r SAMP_THR` samples across the entire dataset to be retained for downstream analysis.
Here, a CPM value of `r CPM_THR` means that a gene is represented if it 
has at least `r round(CPM_THR*LibSizeSum['Min.'])` reads in the sample with the 
lowest sequencing depth (library size `r round(LibSizeSum['Min.'],1)` million).
Note that the thresholds used here are arbitrary as there are no hard and fast 
rules to set these by.
The voom-plot, which is part of analyses done to remove heteroscedasticity,
can be examined to verify that the filtering performed is adequate.


<!-- or at least 
`r round(CPM_THR*LibSizeSum['Max.'])` counts in the sample with the 
greatest sequencing depth (library size `r round(LibSizeSum['Max.'],1)` million).
-->

Remove weakly represented genes and replot densities.

`r weak_flg <- rowSums(edgeR::cpm(featureCountsA) > CPM_THR) < SAMP_THR `
Removing `r round(100 * mean(weak_flg), 1)`%  of genes...

```{r removeWeak, cache=T, cache.vars=c('featureCountsAF', 'genes_annotAF', 'lcpm_mtx'), echo=F, include=F}

# Use suffix 'F' for Filtered genes
featureCountsAF <- featureCountsA[!weak_flg, ]

genes.ndx <- match(rownames(featureCountsAF), genes_annot$Symbol)
if(sum(is.na(genes.ndx))) stop("featureCountsA/genes_annot mismatch")
genes_annotAF <- genes_annot[genes.ndx,]

lcpm_mtx <- edgeR::cpm(featureCountsAF, log = T)
dim(lcpm_mtx)

rm(featureCountsA)

```


```{r densityLcpm2, fig.height=4, fig.width=10, fig.cap='Sample $log_2$ CPM densities after removing weak genes', eval=T, echo=F}

par(mar = c(4, 3, 2, 1))
plot(density(lcpm_mtx[, 1]),
  col = groupCol[sampDescA$group[1]],
  lwd = 2, ylim = c(0, .25), las = 2, main = "", xlab = "log2 CPM"
)
#abline(v = 0, col = 3)
# After verifying no outliers, can plot a random subset 
for (JJ in sample(2:ncol(lcpm_mtx), size = 100)) {
  den <- density(lcpm_mtx[, JJ])
  lines(den$x, den$y, col = groupCol[sampDescA$group[JJ]], lwd = 2)
} # for(JJ
legend("topright", legend = names(groupCol), 
  text.col = groupCol, bty = "n")

```

<!--
Note that the $log_2(CMP)$ distribution is not quite symmetric.
-->

As another sanity check, we will look at a 
multidimensional scaling plot of distances between gene expression
profiles.  We use `plotMDS` in limma package [@Ritchie:2015aa]),
which plots samples on a two-dimensional scatterplot so that distances on
the plot approximate the typical log2 fold changes between the
samples.   

Before producing the MDS plot we will normalize the distributions.
We will store the data into s `DGEList` object as this is convenient
when running many of the analyses implemented in the edgeR and limma packages.
Call the set 'AF', for set 'A', 'Filtered'.

```{r getDGEL, cache=T, cache.var=c('AF_dgel', 'AF_lcmp_mtx')}

AF_dgel <- edgeR::DGEList(
  counts = featureCountsAF,
  genes = genes_annotAF,
  samples = sampDescA,
  group = sampDescA$group
)
AF_dgel <- edgeR::calcNormFactors(AF_dgel)
AF_lcmp_mtx <- edgeR::cpm(AF_dgel, log = T)

# Save AF_dgel to facilitate restarting
# remove from final version
save(list = "AF_dgel", file = "RData/AF_dgel")
```

Verify that the counts are properly normalized.


```{r normedLcpmBxp, cache=T, fig.height=4, fig.width=10, fig.cap='Sample log2 count boxplots', echo=F}

par(mar = c(1, 3, 2, 1))
boxplot(AF_lcmp_mtx,
  ylim = c(1, 8), ylab='Normalized Log CPM',
  staplewex = 0,       # remove horizontal whisker lines
  staplecol = "white", # just to be totally sure :)
  outline = F,         # remove outlying points
  whisklty = 0,        # remove vertical whisker lines
  las = 2, horizontal = F, xaxt = "n",
  border = groupCol[sampDescA$group]
)
legend("top", legend = names(groupCol), text.col = groupCol,
  ncol = 2, bty = "n")
# Add reference lines
SampleMedian <- apply(AF_lcmp_mtx, 2, median)
abline(h = median(SampleMedian), col = "grey")
axis(side = 4, at = round(median(SampleMedian), 1),
  las = 2, col = "grey", line = -1, tick = F)

```

Proceed with MDS plots.

```{r plotMDS, cache=T, fig.height=5, fig.width=10, fig.cap='MDS plots of log-CPM values', echo=F}

par(mfcol = c(1, 2), mar = c(4, 4, 2, 1), xpd = NA, oma = c(0, 0, 2, 0))

# wo loss of generality, sample 500 samples
# simply a matter of convenience to save time
# remove from final version
set.seed(1)
samp_ndx <- sample(1:ncol(AF_lcmp_mtx), size = 500)
MDS.out <- limma::plotMDS(AF_lcmp_mtx[, samp_ndx],
  col = groupCol[sampDescA$group[samp_ndx]], pch = 1
)
legend("topleft",
  legend = names(groupCol),
  text.col = groupCol, bty = "n"
)

MDS.out <- limma::plotMDS(AF_lcmp_mtx[, samp_ndx],
  col = groupCol[sampDescA$group[samp_ndx]], pch = 1,
  dim.plot = 3:4
)
```

The MDS plot, which is analogous to a PCA plot adapted to gene exression data,
does not indicate strong clustering of samples.  The fanning pattern observed in the
first two dimensions indicates that a few samples are drifting way from the
core set, but in no particular direction.  There is some structure in the
3rd and 4th dimension plot which should be investigated.  
`glMDSPlot` from package `Glimma` provides an interactive MDS 
plot that can extremely usful for exploration

```{r GlMDSplot, echo=T,cache=T, cache.vars='', fig.height=6, fig.width=11,fig.cap="MDS plots of log-CPM values", echo=F}

Glimma::glMDSPlot(AF_dgel[, samp_ndx],
  groups = AF_dgel$samples[
    samp_ndx,
    c("group", "trainValGroup", "sampType", "tissue", "title", "stage")
  ],
  main = paste("MDS plot: filtered counts"), #### , Excluding outlier samples"),
  path = ".", folder = figures_DIR,
  html = paste0("GlMDSplot"), launch = F
)
```

Link to glMDSPlot: 
[Here](`r file.path(figures_DIR, paste0("GlMDSplot.html"))`)  

No obvious factor links the samples in the 3 clusters observed on the
4th MDS dimensions. The percent of variance exaplained by this dimension or 
$\approx$ 4%.   The glMDSPlot indicates further segregation along
the 6th dimension.  The percent of variance exaplained by this dimension or 
$\approx$ 2%.  Tracking down this source of variability may be quite challenging,
especially without having the complete information about the sample attributes 
and provenance.  

Unwanted variability is a well-documented problem in the analysis of RNA-Seq data
(see Peixoto et al. (2015) [@Peixoto:2015aa]), and many procedures have been proposed
to reduce the effect of unwanted variation on RNA-Seq analsys results 
([@Gandolfo:2018aa;@Peixoto:2015aa;@Risso:2014aa]).  There are undoubtedly
some similar sources of systematic variation in the 5hmC data, but it is
beyond the scope of this work to investigate these in this particular dataset.
Given that the clustering of samples occurs in MDS dimensions that explain
a small fraction of variability, and that these is no assocation with the
factor of interest, HCC vs Control, these sources of variability should not
interfere too much with our classification analysis.  It would nonetheless be interesting
to assess whether downstream results can be improved by removing this variability.



### Creating a design matrix and contrasts  {-}

Before proceeding with the statistical modeling used for the 
differential expression analysis, we need to set up a
model design matrix.

```{r DEADesign, cache=F, include=T, echo=F, include=T}

Design_mtx <- model.matrix( ~  -1 + group, data=AF_dgel$samples)
colnames(Design_mtx) <- sub('group', '', colnames(Design_mtx))

cat("colSums(Design_mtx):\n")
colSums(Design_mtx)

Contrasts_mtx <- limma::makeContrasts(
  HCCvsControl = HCC  - Control,
  levels=colnames(Design_mtx))

cat("Contrasts:\n")
Contrasts_mtx

```

```{r printDesign, echo=F, include=F}
 knitr::kable(head(Design_mtx), caption='Design Matrix') %>%
  kableExtra::kable_styling(full_width = F)
```


### Removing heteroscedasticity from the count data {-}


As for RNA-Seq data, for 5hmC count data the variance is not independent of the mean.
In `limma`, the R package we are using for our analyses, 
linear modeling is carried out on the log-CPM values which are assumed to be 
normally distributed and the mean-variance relationship is accommodated using precision 
weights calculated by the voom function.  We apply this transformation next.


```{r Voom1, cache=T, cache.vars=c('filteredCountsAF_voom'), fig.height=6, fig.width=11, fig.cap="Removing heteroscedascity", echo=T}

par(mfrow=c(1,1))
filteredCountsAF_voom <- limma::voom(AF_dgel, Design_mtx, plot=T)

```

Note that the voom-plot provides a visual check on the level of filtering performed upstream.
If filtering of lowly-expressed genes is insufficient, 
a drop in variance levels can be observed at the low end of the 
expression scale due to very small counts. 

<!--
Means (x-axis) and variances (y-axis) of each gene are plotted to show the dependence between the two before voom is applied to the data ( A) and how the trend is removed after voom precision weights are applied to the data ( B). The plot on the left is created within the voom function which extracts residual variances from fitting linear models to log-CPM transformed data. Variances are then re-scaled to quarter-root variances (or square-root of standard deviations) and plotted against the mean expression of each gene. The means are log 2-transformed mean-counts with an offset of 2. The plot on the right is created using plotSA which plots log 2 residual standard deviations against mean log-CPM values. The average log 2 residual standard deviation is marked by a horizontal blue line. In both plots, each black dot represents a gene and a red curve is fitted to these points.
-->

<!--
To get a sense of how this compares with RNA-Seq dsta, we can take a look at Figure 4 in Law et al. [@Law:2018aa]:

`r knitr::include_graphics(
  "Static/images/Law-fig4.gif",  dpi=100)`

We observe that the variability in the 5hmC data is quite a bit lower.  Statistical summaries would give
a better idea.  

-->


### Fit linear models and examine the results {-}

Having properly filtered and  normalized the data,
the linear models can be fitted to each gene and the results
examined to assess differential expression between the two groups
of interest, in our case HCC vs Control.

Table \@ref(tab:lmFit) displays the counts of genes in each DE category:

```{r lmFit, cache=T, echo=T, cache.vars=c('filteredCountsAF_voom_efit','filteredCountsAF_voom_efit_dt'),echo=F}


 filteredCountsAF_voom_fit <- limma::lmFit(filteredCountsAF_voom, Design_mtx)
 colnames(filteredCountsAF_voom_fit$coefficients) <- sub("\\(Intercept\\)", "Intercept",
 colnames(filteredCountsAF_voom_fit$coefficients) )

 filteredCountsAF_voom_fit <- limma::contrasts.fit(
    filteredCountsAF_voom_fit, contrasts=Contrasts_mtx)

 filteredCountsAF_voom_efit <- limma::eBayes(filteredCountsAF_voom_fit)

 filteredCountsAF_voom_efit_dt <-
 limma::decideTests(filteredCountsAF_voom_efit,adjust.method = "BH", p.value = 0.05)
 
 knitr::kable(summary(filteredCountsAF_voom_efit_dt),
  caption="DE Results at FDR = 0.05") %>% 
  kableExtra::kable_styling(full_width = F)

```

### Graphical representations of DE results: MD Plots {-}

To summarise results for all genes visually, mean-difference plots
(aka MA plot), which display log-FCs from the linear model fit against 
the average log-CPM values can be generated using the plotMD function,
with the differentially expressed genes highlighted.

We may also be interested in whether certain gene features are 
related to gene identification.  Gene GC content, for example, might be
of interest.

```{r mdPlotEfit, cache=T, cache.vars='', fig.height=5, fig.width=11, fig.cap="HCC vs Control - Identified Genes at FDR = 0,05", echo=F}

par(mfrow=c(1,3), mar=c(4,4,2,1),oma=c(0,0,2,0))

# log-fold-change vs ave-expr
limma::plotMD(filteredCountsAF_voom_efit,
 ylim = c(-0.5, 0.5),
 column='HCCvsControl',
 status=filteredCountsAF_voom_efit_dt[,'HCCvsControl'],
 hl.pch = 16, hl.col = c("lightblue", "pink"), hl.cex = .5,
 bg.pch = 16, bg.col = "grey", bg.cex = 0.5,
 main = '',
 xlab = paste0(
    "Average log-expression: IQR=",
    paste(round(quantile(filteredCountsAF_voom_efit$Amean, prob = c(1, 3) / 4), 2),
      collapse = ", "
    )
  ),
  ylab = paste0(
    "log-fold-change: IQR=",
    paste(round(quantile(filteredCountsAF_voom_efit$coefficients[, 'HCCvsControl'], prob = c(1, 3) / 4), 2),
      collapse = ", "
    )
  ),
  legend = F
)
abline(h = 0, col = "black")
rug(quantile(filteredCountsAF_voom_efit$coefficients[, 'HCCvsControl'], prob = c(1, 2, 3) / 4),
  col = "purple",
  ticksize = .03, side = 2, lwd = 2
)
rug(quantile(filteredCountsAF_voom_efit$Amean, prob = c(1, 2, 3) / 4),
  col = "purple",
  ticksize = .03, side = 1, lwd = 2
)

# log-fold-change vs identification

gplots::boxplot2(split(
 filteredCountsAF_voom_efit$coefficients[, 'HCCvsControl'],
 filteredCountsAF_voom_efit_dt[,'HCCvsControl']),
 outline=F,
 border=c("pink", "grey", "lightblue"), xaxt='n',
 ylab='log-fold-change', ylim=c(-.5, .5))
 axis(side=1, at=1:3, c('down', 'notDE', 'up'))


# gc vs identification
genes_ndx <- match(rownames(filteredCountsAF_voom_efit), genes_annotAF$Symbol)
if(sum(is.na(genes_ndx))) stop("filteredCountsAF_voom_efit/genes_annotAF: genes mismatch")
GC_vec <- with(genes_annotAF[genes_ndx,],(G+C)/(A+C+G+T))


boxplot(split(
 GC_vec,
 filteredCountsAF_voom_efit_dt[,'HCCvsControl']),
 outline=F,
 border=c("pink", "grey", "lightblue"), xaxt='n',
 ylab='gene-gc')
 axis(side=1, at=1:3, c('down', 'notDE', 'up'))

 mtext(side=3, outer=T, cex=1.25, "Genes identified at adjusted p-value=0.05")
```


```{r quantlogFC,echo=F}

featureCountsAF_logFC_sum <- sapply(
 split(
 filteredCountsAF_voom_efit$coefficients[, 'HCCvsControl'],
 filteredCountsAF_voom_efit_dt[,'HCCvsControl']),
 quantile, prob = (1:3) / 4)

colnames(featureCountsAF_logFC_sum) <- as.character(factor(
 colnames(featureCountsAF_logFC_sum), 
 levels=c("-1", "0", "1"),
 labels=c('down', 'notDE', 'up')
))


knitr::kable(featureCountsAF_logFC_sum,
  digits = 2,
  caption = "log FC quartiles by gene identification") %>%
  kableExtra::kable_styling(full_width = F)

```


We see that while many genes are identified, the effect sizes are quite small.
The GC content of down regulated genes tends to be slightly lower than the
rest of the genes.



### DE genes at 10% fold change {-}

For a stricter definition on significance, one may require log-fold-changes 
(log-FCs) to be above a minimum value. The treat method 
(McCarthy and Smyth 2009 [@McCarthy:2009aa]) can be used to calculate p-values 
from empirical Bayes moderated t-statistics with a minimum log-FC requirement. 
The number of differentially expressed genes are greatly reduced if we 
impose a minimal fold-change requirement of 10%.

```{r mdPlotTfit, cache=T, cache.vars='', fig.height=5, fig.width=11, fig.cap="HCC vs Control - Identified Genes at FDR = 0,05 and logFC > 10%",echo=F}

filteredCountsAF_voom_tfit <- limma::treat(filteredCountsAF_voom_fit, lfc=log2(1.10))
filteredCountsAF_voom_tfit_dt <- limma::decideTests(filteredCountsAF_voom_tfit)

cat("10% FC Gene Identification Summary - voom, adjust.method = BH, p.value = 0.05:\n")
summary(filteredCountsAF_voom_tfit_dt)

# log-fold-change vs ave-expr
limma::plotMD(filteredCountsAF_voom_efit,
 ylim = c(-0.5, 0.5),
 column='HCCvsControl',
 status=filteredCountsAF_voom_tfit_dt[,'HCCvsControl'],
 hl.pch = 16, hl.col = c("blue", "red"), hl.cex = .7,
 bg.pch = 16, bg.col = "grey", bg.cex = 0.5,
 main = '',
 xlab = paste0(
    "Average log-expression: IQR=",
    paste(round(quantile(filteredCountsAF_voom_efit$Amean, prob = c(1, 3) / 4), 2),
      collapse = ", "
    )
  ),
  ylab = paste0(
    "log-fold-change: IQR=",
    paste(round(quantile(filteredCountsAF_voom_efit$coefficients[, 'HCCvsControl'], prob = c(1, 3) / 4), 2),
      collapse = ", "
    )
  ),
  legend = F
)
abline(h = 0, col = "black")
rug(quantile(filteredCountsAF_voom_efit$coefficients[, 'HCCvsControl'], prob = c(1, 2, 3) / 4),
  col = "purple",
  ticksize = .03, side = 2, lwd = 2
)
rug(quantile(filteredCountsAF_voom_efit$Amean, prob = c(1, 2, 3) / 4),
  col = "purple",
  ticksize = .03, side = 1, lwd = 2
)


```


## Analysis of coverage variability

We will use the methods described in Hart et al. (2013) [@Hart:2013aa] 
to characterize coverage variability in these data. 
These methods do not take multiple comparisons into account. 
Other tools for sample size calculation in RNA-Seq studies include 
Bi and Liu (2016) [@Bi:2016aa,], Baccarella (2018) [@Baccarella:2018aa], 
Guo (2014) [@Guo:2014aa], Yu (2017) [@Yu:2017aa], and 
Zhao (2018) [@Zhao:2018ab]. 
Poplawski (2018) [@Poplawski:2017aa] 
evaluated RNA-seq sample size tools identified from a systematic search. They found the six evaluated tools provided widely different answers, which were strongly affected by fold change.


The references listed above aim at providing guidance for RNA-Seq experimental design.
There is much discussion and a wide range of opinion on sample size requirements
to ensure reproducibility in RNA-Seq results.  At one end of the spectrum, 
Ein-Dor et al. (2006) [@Ein-Dor:2006aa] argue that thousands of samples are needed to 
generate a robust gene list for predicting outcome in cancer.  At the other end,
Dobbin et al. (2007, 2008) [@Dobbin:2007aa;@Dobbin:2008aa] claim that sample sizes in 
the range of 20–30 per class may be adequate for building a good predictor in many cases. 
Part of the disparity in sample size requirement recomendation comes from 
differences of opinion in terms of what constitutes `reproducible results`. 
In the context of sample classification, if we focus on  the predicted probabilities
for individual samples, we may find good  reproducibility across studies with
moderate samples sizes.  If, on the other hand, we closely inspect the gene
signatures reported across studies, much greater sample sizes may be required to
achieve concordance.  Kim (2009) [@Kim:2009aa], like Ein-Dor et al., also find
issues in RNA-Seq research in terms of the instability of identified prognostic 
gene signatures, few overlap between independently developed prognostic gene signatures, and
poor inter-study applicability of gene signatures.  Fan et al. (2006) [@Fan:2006aa],
on the other hand, found good concordance among gene-expression–based predictors 
for breast cancer.  We will return to this question when we examine
the relationship between classification model results and sample size in this
dataset later on this paper.  


For two groups comparisons, the basic formula for the required number of samples per group is:

$$ n = 2(z_{1-\alpha/2} + z_{\beta})^2 \frac{(1/\mu + \sigma^2)}{ln(\Delta^2)} $$

* The parameters $\alpha$ and $\beta$ are **size** and **power** of the test.  
* $\Delta$ is the targeted **effect size**.   
* $\mu$ and $\sigma$ are the **mean** and **coefficient of variation** 
of the distribution of measurement, gene representation indices in this case.

These three parameters will be fixed across genes or a given study, and are often dictated by
external requirements. Typical values might be an effect size of $\Delta = 1.5$ 
(a.k.a fold change), corresponding to detection of a 50% change in gene expression 
between the two groups.
$z_{1 - .05/2} = 1.96$, corresponding to a two sided test at size $\alpha = 0.05$;
and $z_{.90}= 1.28$ corresponding to 90% power.
The other two variables will be gene and experiment dependent: the normalized depth 
of coverage $\mu$ of the gene, and the coefficient of variation $\sigma$ in this 
gene between biological replicates.  The technical variation of the 
comparison is inversely proportional to the number of sequenced reads
for the gene and therefore decreases with sequencing depth.
The biological variation is a property of the particular gene/model system/condition under study.
One would expect it to be smaller for uniform systems such as cell culture and/or products that
are under tight regulatory control, and larger for less uniform replicates such 
as human subject samples.  The dataset under study in this report is the first of its
kind to give us an idea of variability levels in 5hmC representation.


<!-- SKIP THIS
Before examining biological variability lets take a look at
the sequencing rates.

-->
```{r sumCount, cache=T, cache.vars='DD.cpmSum.frm', eval=F, echo=F}

countSum_frm <- do.call(
  "rbind",
  lapply(unique(sampDescA$group), function(GRP) {
    GRP_featureCountsAF <- featureCountsAF[, sampDescA$group == GRP]
    N_GRP <- ncol(GRP_featureCountsAF)
    MeanTotCovM <- mean(colSums(featureCountsAF[, sampDescA$group == GRP])) / 1e6
    GRP_countSum_mtx <- do.call("cbind", lapply(
      1:ncol(GRP_featureCountsAF),
      function(CC) {
        CC_cpm_vec <- GRP_featureCountsAF[, CC]
        100 * table(cut(CC_cpm_vec, breaks = c(0, .1, 1, 10, 50, 100, 1000, Inf))) /
          length(CC_cpm_vec)
      }
    ))
    GRP_countSum_vec <- apply(GRP_countSum_mtx, 1, median)
    GRP_Zeros_vec <- colSums(GRP_featureCountsAF == 0)
    GRP_ZerosPct_vec <- 100 * GRP_Zeros_vec / nrow(GRP_featureCountsAF)
    data.frame(
      Group = GRP, N = N_GRP, MeanTotCovM = round(MeanTotCovM, 1),
      Zeros = median(GRP_Zeros_vec), ZeroPct = median(GRP_ZerosPct_vec),
      t(GRP_countSum_vec), check.names = F
    )
  })
)
knitr::kable(countSum_frm) %>%
  kableExtra::kable_styling(full_width = F)

```


As in Hart et al. (2013) [@Hart:2013aa] we estimate the 
biological coefficient of variation (CV) in expression across samples in the data set
using a negative binomial model.

<!--
From [edgeR User's Guide](https://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf)
, the biolgical coefficient of variation (BCV), is the square root of the dispersion
parameter under the negative binomial model. Hence, it is equivalent to estimating the
dispersion(s) of the negative binomial model.
-->

```{r estBCV, cache=T, cache.vars=c('BCV_mtx','BCV_90perc_vec', 'BCV_50perc_vec'),message=F, echo=F,include=F}

BCV_mtx <- do.call("cbind", lapply(
  unique(sampDescA$group),
  function(GRP) {
    GRP_dgel <-
      edgeR::DGEList(counts = featureCountsAF[, sampDescA$group == GRP])
    GRP_dgel <- edgeR::estimateDisp(GRP_dgel)

    sqrt(GRP_dgel$tagwise.dispersion)
  }
))
colnames(BCV_mtx) <- unique(sampDescA$group)

```

```{r plotBCV,  fig.height=5, fig.width=10, fig.cap="Cumulative Distribution of CV - rug = 90th percentile",message=F, echo=F,include=T}

plot(spatstat::CDF(density(BCV_mtx[, 1])),
  col = groupCol[colnames(BCV_mtx)[1]],
  lwd = 2, ylab = "Prob(BCV<x)",
  xlim = c(0, 0.5)
)
for (JJ in 2:ncol(BCV_mtx)) {
  plot(spatstat::CDF(density(BCV_mtx[, JJ])),
    col = groupCol[colnames(BCV_mtx)[JJ]],
    lwd = 2, add = T, xlim = c(0, 2.0)
  )
}
legend("bottomright", legend = names(groupCol), col = groupCol, lwd = 2)
BCV_90perc_vec <- round(apply(BCV_mtx, 2, quantile, prob = 0.90), 2)
BCV_50perc_vec <- round(apply(BCV_mtx, 2, quantile, prob = 0.50), 2)
for (JJ in 1:length(BCV_90perc_vec)) {
  rug(BCV_90perc_vec[JJ],
    lwd = 2, ticksize = 0.05,
    col = groupCol[names(BCV_90perc_vec)[JJ]]
  )
}

```
<br/>

BCV values are quite low:

```{r bspBCV,  fig.height=5, fig.width=10, fig.cap="BCV distributions",message=F, echo=F,include=T}

 boxplot(BCV_mtx, outline=F, ylab='BCV')

```

<br/>

<!-- 
 SampleSize_f <- function(Alpha=.01, Beta=.80, Delta=2, Mu, Disp)
 {temp <- qnorm(1-Alpha/2) + qnorm(Beta)
  vtemp <- 1/Mu + Disp
  2*vtemp*(temp/log(Delta))^2}

 cat("SampleSize_f(Mu=10, Disp=.4^2):\n")
 SampleSize_f(Mu=10, Disp=.4^2)
-->

We can now look at sample size estimates to required to detect 
various effect sizes.  The effect sizes examined here are selected
based on the differential representation analysis in Section \@ref(dra) below.


```{r sampleSizeCurve, cache=T, cache.vars='', fig.height=8, fig.width=11,fig.cap='Sample Size Estimates', echo=F}

par(mfrow = c(2, 1), mar = c(3, 3, 2, 1), oma = c(2, 2, 2, 0))
for (EFFECT in c(1.10, 1.05)) {
  plot(
    x = 10:100, 
    y = RNASeqPower::rnapower(depth = 10:100, cv = max(BCV_90perc_vec), 
         effect = EFFECT, alpha = .05, power = .80),
    lwd = 2, ylab = "", xlab = "", type='l'
  )
  title(paste("Effect =", EFFECT, "cv =", 
     max(BCV_90perc_vec), "Alpha=0.05, Beta=0.80"))
}
mtext(side = 1, outer = T, "Gene Coverage")
mtext(side = 2, outer = T, "Samples Needed")
mtext(side = 3, outer = T, cex = 1.25, 
    "Negative Binomial 2 Group Sample Size Curves")

```

For the filtered reads, coverage looks like this:

```{r quantCountsAF, echo=F}

featureCountsAF_quant <- apply(featureCountsAF, 2, function(CC) {
  c(quantile(CC, prob = c(.15, (1:3) / 4)))###, totCovM = sum(CC) / 1e6)
})

featureCountsAF_quant2 <- apply(featureCountsAF_quant, 1, function(RR) {
  quantile(RR, prob = (1:3) / 4)
})

knitr::kable(featureCountsAF_quant2,
  digits = 1,
  caption = paste(
    "Coverage Summary - Columns are sample coverage quantiles and total coverage",
    "\nRows are quartiles across samples"
  )
) %>%
  kableExtra::kable_styling(full_width = F)

```


From this table, we see that 25% of the samples have an upper quartile of gene coverage
exceeding `r round(featureCountsAF_quant2["75%", "75%"],1)` reads, 25% of samples
have a 15 percentile of coverage lower than
`r featureCountsAF_quant2["25%", "15%"]`, etc.


Note that with these data, moderate sample sizes are adequate 
to detect genes with effect sizes as small as 1.05 (5% fold change).  This is due to the fact
that the biological variability in gene body 5hmC density is quite low.
For human samples, RNA-Seq within group biological variability is
typically in the 0.4-1.0 range [@Hart:2013aa].  

### Another look at BCV {-}

Alternatively, we can extract sigma and signal from the fit objects.

```{r altCV, cache=T, cache.vars='',fig.height=4, fig.width=6, fig.cap="Alternative CV Calculation"}
lib.size <- colSums(AF_dgel$counts)

fit <- filteredCountsAF_voom_efit
sx <- fit$Amean + mean(log2(lib.size + 1)) - log2(1e+06)
sy <- sqrt(fit$sigma)

CV <- sy/sx    

cor(cbind(BCV_mtx, CV)) 
boxplot(cbind(BCV_mtx, CV), outline=F) 

```

CV computed frm sigma and coverage extracted from the fit are slighly higher.
For purposes of gene identification, we really want to look at effect size in relation to residual 
standard deviation.  


```{r plotSNR,  fig.height=5, fig.width=10, fig.cap="Cumulative Distribution of SNR - rug = 25, 50, 75 and 90th percentile",message=F, echo=F,include=T}

Effect <- abs(filteredCountsAF_voom_efit$coefficients[,'HCCvsControl'])
Noise <- filteredCountsAF_voom_efit$sigma
SNR <- Effect/Noise

plot(spatstat::CDF(density(SNR)),
  col = 1, lwd = 2, ylab = "Prob(SNR<x)",
  xlim = c(0, 0.2)
)

SNR_quant <- quantile(SNR, prob=c((1:3)/4,.9))
rug(SNR_quant,
    lwd = 2, ticksize = 0.05, col = 1
  )


knitr::kable(SNR_quant,
  digits = 3,
  caption = paste(
    "SNR Summary") 
) %>% kableExtra::kable_styling(full_width = F)


```


