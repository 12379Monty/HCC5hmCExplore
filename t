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
sampDesc$DxStage <- with(sampDesc, ifelse(outcome=='HCC', paste0(outcome,':',stage), outcome))

with(
  sampDesc %>% dplyr::filter(sampType == "blood"),
  table(DxStage, trainValGroup, exclude = NULL)
)

```

For this analysis, we will consider early stage cancer samples
and healthy or benign samples from the `Train` or `Val-1` subsets.

```{r subsetSamples, cache=T, cache.vars=c('sampDescA','groupCol')}

sampDescA <-
  sampDesc %>%
  dplyr::filter(sampType == "blood" & (trainValGroup %in% c("Train", "Val-1")) &
    ((outcome2 == "BenignHealthy") | 
     (outcome2 == "HCC" & stage == "Early"))) %>%
  dplyr::rename(group = outcome2) %>%
  dplyr::arrange(group, sampID)
# Recode group
sampDescA$group <- with(sampDescA, 
   ifelse(group == "BenignHealthy", "Control", group))
# set groupCol for later
groupCol <- c("#F3C300", "#875692")
names(groupCol) <- unique(sampDescA$group)

with(sampDescA, table(group, exclude = NULL))

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


```{r getfeatures, cache=T, cache.vars=c('featureCountsA')}

  featureCountsA <- cbind(Train_featureCount, 
                          Val1_featureCount, 
                          Val2_featureCount)[,rownames(sampDescA)]

```

We first look at coverage - make sure there isn't too much disparity of coverage 
across samples. To detect shared variability, samples can be annotated and ordered
according to sample features that may be linked to sample batch processing.  Here we 
the samples have been ordered by group and sample id (an alias of geoAcc).

```{r lcpmBxp, cache=T, fig.height=4, fig.width=10, fig.cap='Sample log2 count boxplots'}

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
axis(side = 2, at = round(median(SampleMedian), 1), 
  las = 2, col = "grey", line = -1, tick = F)

```

Coverage level looks fairly comparable across samples.  It is sometimes helpful to
keep track of the actual coverage which can be adequetely tracked by distribution
quantiles.

```{r quantCounts}
 featureCountsA_quant <- apply(featureCountsA, 2, function(CC)
  c(quantile(CC, prob=c(.1,(1:3)/4)), totCovM=sum(CC)/1e6))

 featureCountsA_quant2 <- apply(featureCountsA_quant, 1, function(RR)
   quantile(RR, prob=(1:3)/4))

 knitr::kable(featureCountsA_quant2, digits=1,
 caption=paste("Coverage Summary - Columns are sample coverage quantiles and total coverage",
 "\\nRows are quartiles across samples"))

```

From this table, we see that 25% of the samples have total coverage exceeding
`r round(featureCountsA_quant2["75%", "totCovM"],1)`, 25% of samples
have a 10 percentile of coverage lower than
`r featureCountsA_quant2["25%", "10%"]`, etc.  


We next look at relative log representation (RLR) (in the context of measuring the density of 
5hmC marks in genes, we refer to `representation` as opposed to `expression`; the
two can be used interchangibly) -
make sure the shapes of the distributions are not widely different.

```{r rlr, cache=T, cache.vars='lcpm_mtx', fig.height=4, fig.width=10, fig.cap='Sample RLR', eval=T, echo=T}

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

We note that the HCC samples have slightly more variable coverage distribution.
A few samples are quite different.

## Differential representation analysis {#dra}

   - word on GC content


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

```{r densityLcpm, fig.height=4, fig.width=10, fig.cap='Sample $log_2$ CPM densities', eval=T, echo=T}

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

```{r removeWeak, cache=T, cache.vars=c('featureCountsA', 'genes_annotA', 'lcpm_mtx')}

featureCountsA <- featureCountsA[!weak_flg, ]
genes_annotA <- genes_annot[rownames(featureCountsA), ]
lcpm_mtx <- edgeR::cpm(featureCountsA, log = T)
dim(lcpm_mtx)

```


```{r densityLcpm2, fig.height=4, fig.width=10, fig.cap='Sample $log_2$ CPM densities after removing weak genes', eval=T, echo=T}

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

