# Fitted Model Suite {#model-suite}

We now examine the results of fitting a suite of models to
investigate the effect of sample size on 
various aspects of model performance:  

* assessed accuracy: out-of-fold estimates of precision and variability  and
cv assessed accuracy and bias.

* selected feature profile stability - to what extent does the
feature set implicitly selected by the lasso vary across random
sampling and what is the effect of sample size.

It is hypothesized that below a certain threshold,
sample sizes are too small to provide reliable estimates
of performance or stable selected feature profiles.  

We will attempt to separate variability which is due to
sample size from variability due to sample composition.
To to this we will track sample quality.
Predicted probabilities from fitted model can be transformed into sample
quality scores: $Q_i = p_i^{y_i}(1-p_i)^{1-y_i}$, where $p_i$ is the
estimated probability of HCC for sample i and $y_i$ is 1 for HCC samples and
0 for Controls.  ie. we use the fitted sample contribution to the
likelihood function as a sample quality score.  To derive the quality scores,
we will use the predicted response from a lasso model fitted to the entire data set.

Hard to classify samples will have low quality scores.
In the results that we discuss below, when we look at variability across repeated 
random sampling of different sizes, we can use sample quality scores to investigate 
how much of the variability is due to sample selection.
Note that quality here is not used to say anything about the sample data quality.
Low quality here only means that a sample is different from the 
core of the data set in a way that makes it hard to properly classify.
That could happen if the sample were mislabeled, in which case we could 
think of this sample as being poor quality of course.

## Sample quality scores

To get sample quality scores, fit a lasso model, extracted predicted probabilities and
convert to a quality score.

```{r get-all-data, cache=T, cache.vars=c('all_lcpm_mtx', 'all_group_vec')}

# combine train and test 
all_lcpm_mtx <- rbind(train_lcpm_mtx, test_lcpm_mtx)

# we have to be careful with factors!
# We'll keep as a character and change to factor when needed
all_group_vec <- c(
 as.character(train_group_vec), 
 as.character(test_group_vec)
)
# I suspect adding names to vectors breaks one of the tidy commandments,
# but then again I am sure I have already offended the creed beyond salvation
names(all_group_vec) <- c(
 names(train_group_vec),
 names(test_group_vec)
)

knitr::kable(table(group = all_group_vec),
  caption = "samples by group") %>%
   kableExtra::kable_styling(full_width = F)

```


```{r lasso-fit-all, cache=T, cache.vars=c('cv_lassoAll')}
start_time <-  proc.time()

# since we have not set the foldid, set seed here
set.seed(1)
cv_lassoAll <- glmnet::cv.glmnet(
 x = all_lcpm_mtx,
 y = factor(all_group_vec,levels = c('Control', 'HCC')),
 alpha = 1,
 family = 'binomial',
 type.measure  =  "class",
 keep = F,
 nlambda = 100
)

message("lassoAll time: ", round((proc.time() - start_time)[3],2),"s")

```

We can examine the fit.

```{r look-lassoAll, as.is=T, comments='', cache=T, cache.vars=''}

cv_lassoAll

```

```{r plot-lassoAll, cache=T, cache.vars='', fig.height=5, fig.width=6, fig.cap='lasso model fitted to all samples'}

plot(cv_lassoAll)

```

More features are now selected, but the cv error rate is comparable to the results
of the lasso fit to the training set.  We can use this fit to compute sample
quality scores.  We will use the minimum lambda model to provide
the fitted probabilities.

```{r get-sample-qual, cache=T, cache.vars=c('sample_qual_vec','lassoAll_conf_vec', 'lassoAll_conf_mtx')}

# predicted probs
lassoAll_predResp_vec <- predict(
 cv_lassoAll,
 newx = all_lcpm_mtx,
 s = "lambda.min",
 type = "resp"
 )

# also get predicted class
lassoAll_predClass_vec <- predict( 
 cv_lassoAll,
 newx = all_lcpm_mtx,
 s = "lambda.min",
 type = "class"
 )

# get qual scores

# Note here that factors are a little but of a pain.
y <- as.numeric(all_group_vec == 'HCC')
p <- lassoAll_predResp_vec
sample_qual_vec <- p^y*(1-p)^(1-y)
sample_qual_vec <- sample_qual_vec[,1] # to drop the matrix structure

# Look at quality scores as a finction of classification
lassoAll_conf_vec <- paste(
 y, 
 as.numeric(lassoAll_predClass_vec=='HCC'),
 sep = ':'
)

lassoAll_conf_mtx <- table(
 truth = all_group_vec,
 pred = lassoAll_predClass_vec
)

```

```{r print-conf, fig.cap=''}

knitr::kable(lassoAll_conf_mtx,
  caption = "Rows are truth.  Columns are predicted") %>%
   kableExtra::kable_styling(full_width = F)

```


The model makes few errors (although in most liquid biopsy applications,
this error rate would arguably be too high).  We can examine how quality
varies among classification groups.

```{r plot-qual-conf, cache=T, cache.vars='', fig.height=5, fig.width=5, fig.cap='quality scores by classification - HCC=1'}

gplots::boxplot2(split(sample_qual_vec, lassoAll_conf_vec), ylab = 'Quality Score')
title(sub = 'Classification - Truth:Predicted')

```

## Simulation Design

We are now ready to run the simulations.

```{r simParms, cahce=F}
 SIM <- 30
 SIZE <- c(25, 50, 100, 200, 300)
 CV_REP <- 30

```

Simluation parameters:  

* Number of simulations : SIM = $`r SIM`$

* Sample sizes: SIZE = $`r SIZE`$  

* Number of CV Replicates:  CV_REP = $`r CV_REP`$


We will repeat the simulation process SIM = $`r SIM`$ times.
For each simulation iteration, we will select $`r max(SIZE)`$ Control and 
$`r max(SIZE)`$ HCC samples at random.  Models will be fitted and analyzed
to balanced subsets of sise SIZE = $`r SIZE`$, in a telescopic manner to
emulated a typical sample accrual process.  Note that in this accural process
there is no time effect - the accrual process is completely randomized.  In practice,
there could be significant time effects.  For example, the first 25 HCC samples could come
from Center A, while the next 25 could come from Center B.  In other words,
there is no batch effect or shared variability in our simulation,
while these are almost always present in real data, including 
batch effects that are associated with class labels - controls being in
different batches than affected samples is an all too common occurence,
for example.  One should be especially watchful of potential batch effects
when dealing with blood samples as blood is notoriously finicky in
character [@Huang:2017aa; @Permenter:2015aa;].
Presented with results that look impressively good based on a small data set,
one should definitely be skeptical of the promise of future equally good results.

For a given simulation and a given sample size, we will obtain
CV_REP = $`r CV_REP`$ cross-validated lasso fits.  From these fits,
we can obtain $`r CV_REP`$ out-of-fold assessments of classification accuracy 
to get a sense if its variability. From each cv replicate, we also obtain
an estimated model size and a set of selected features.  We will want
to examine how these stabilize as the sample size increases.

Note that we limit the simulations to a maximum of sample size of 300 in 
order to allow us to have repeated simulations with low overlap.  With 300
randomly selected HCC samples, the expected overlap between two randomly
selected sets of HCC samples is $`r round(100*(300/sum(all_group_vec=='HCC'))^2,1)`$%.
For Controls the expected overlap is $`r round(100*(300/sum(all_group_vec=='Control'))^2,1)`$%. 


## Setup simulation 

To setup the simulation, we only need two master tables: one for the selection of Controls
and one for the selection of HCC samples.

```{r get-all-vec, cache=T, cache.vars=c('all_control_vec', 'all_affected_vec')}

all_control_vec <- names(all_group_vec[all_group_vec=='Control']) 
all_affected_vec <- names(all_group_vec[all_group_vec=='HCC'])  

```

We have $`r length(all_control_vec)`$ IDs of control samples in stored in `all_control_vec`
and $`r length(all_affected_vec)`$ IDs of affected samples in stored in `all_affected_vec`.
To create random samples from these we only need to randomly select indices from
each vector.

  
```{r getSimTable, cache=T, cache.vars=c('sim_control_mtx', 'sim_affected_mtx')}

set.seed(12379)

sim_control_mtx <- sapply(
 1:SIM, 
 function(dummy) 
   sample(1:length(all_control_vec), size =  max(SIZE))
)


sim_affected_mtx <- sapply(
 1:SIM, 
 function(dummy) 
   sample(1:length(all_affected_vec), size =  max(SIZE))
)


```

Each simulation is specified by a given column of the simulation design matrices:
`sim_control_mtx` and `sim_affected_mtx`.  Within each simulation, we can run
the analyses of size $`r SIZE`$.

We can examine how much variability we have in the quality scores of the selected samples.

```{r look-sim-qual, cache=T, cache.vars=c('sim_control_qual_mtx', 'sim_affected_qual_mtx'), fig.height=8, fig.width=10, fig.cap='sample quality by simulation run'}

all_control_qual_vec <- sample_qual_vec[all_control_vec]
sim_control_qual_mtx <- sapply(
  1:ncol(sim_control_mtx), 
  function(CC) all_control_qual_vec[sim_control_mtx[,CC]]
 )

all_affected_qual_vec <- sample_qual_vec[all_affected_vec]
sim_affected_qual_mtx <- sapply(
  1:ncol(sim_affected_mtx),  
  function(CC) all_affected_qual_vec[sim_affected_mtx[,CC]]
 )

# Get stage from SIZE 
stage_vec <- cut(1:nrow(sim_control_qual_mtx), c(0,SIZE), include.lowest = T)

sim_control_qual_byStage_lst <- do.call('c', 
 lapply(1:ncol(sim_control_qual_mtx), 
  function(CC) c(split(sim_control_qual_mtx[,CC], stage_vec),NA)
 )
)

sim_affected_qual_byStage_lst <- do.call('c', 
 lapply(1:ncol(sim_affected_qual_mtx), 
  function(CC) c(split(sim_affected_qual_mtx[,CC], stage_vec),NA)
 )
)

# PLOT
par(mfrow=c(2,1), mar = c(2,5,2,1))
# control
boxplot(
  sim_control_qual_byStage_lst, 
  outline = F, 
  border = 1:6,
  ylab = 'Quality Score',
  xaxt = 'n'
)
legend('bottomright', title = 'Stage', ncol = 2,
 legend = names(sim_control_qual_byStage_lst[1:5]), 
 text.col = 1:5,
 bty = 'n', horiz = F
)
sim_ndx <- which(names(sim_control_qual_byStage_lst) =='')
abline(v = sim_ndx, col = 'grey')
axis(
  side = 1, 
  at = sim_ndx-2, 
  label = 1:length(sim_ndx),
  tick = F, 
  line = -1, las = 2,
  cex.axis = 0.8)
title("Control sample quality by stage and simulation")

# affected
boxplot(
  sim_affected_qual_byStage_lst, 
  outline = F, 
  border = 1:6,
  ylab = 'Quality Score',
  xaxt = 'n'
)
sim_ndx <- which(names(sim_affected_qual_byStage_lst)=='')
abline(v = which(names(sim_affected_qual_byStage_lst)==''), col = 'grey')
axis(
  side=1, 
  at = sim_ndx-2,        
  label = 1:length(sim_ndx),
  tick = F, 
  line = -1, las = 2,
  cex.axis = 0.8)
title("Affected sample quality by stage and simulation")

```

We see some variability in sample quality in the smaller analysis stages.
This may lead observers to be overly optimistic, or overly pessimistic,
in the early accrual stages.  

## Run simulations


As these make take a while to run, 
we will save the results of each similation to a different
object and store to disk.  These can be easily read from disk
when needed for analysis.


The simulation saves results to the file system and
only needs to be run once.  The simulation takes $\approx$ 8 minutes
per iteration, or 4 hours of run time on a laptop.
(Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS Mojave 10.14.6)

```{r runSim, cache=T, cache.vars='', eval=F}

# Get stage from SIZE
stage_vec <- cut(1:nrow(sim_control_qual_mtx), c(0, SIZE), include.lowest = T)

for (SIMno in 1:ncol(sim_control_qual_mtx)) {
  start_time <- proc.time()

  cat("Running simulation ", SIMno, "\n")

  sim_cv_lst <- lapply(1:length(levels(stage_vec)), function(STGno) {
    Stage_rows_vec <- which(stage_vec %in% levels(stage_vec)[1:STGno])
    #cat("Stage ", STGno, "- analyzing", length(Stage_rows_vec), "paired samples.\n")

    sim_stage_samples_vec <- c(
      all_control_vec[sim_control_mtx[Stage_rows_vec, SIMno]],
      all_affected_vec[sim_affected_mtx[Stage_rows_vec, SIMno]]
    )
    sim_stage_lcpm_mtx <- all_lcpm_mtx[sim_stage_samples_vec, ]
    sim_stage_group_vec <- all_group_vec[sim_stage_samples_vec]
    print(table(sim_stage_group_vec))

    sim_stage_cv_lst <- lapply(1:CV_REP, function(CV) {
      cv_fit <- glmnet::cv.glmnet(
        x = sim_stage_lcpm_mtx,
        y = sim_stage_group_vec,
        alpha = 1,
        family = "binomial",
        type.measure = "class",
        keep = T,
        nlambda = 30
      )
      ndx_1se <- which(cv_fit$lambda == cv_fit$lambda.1se)

      nzero_1se <- cv_fit$nzero[ndx_1se]
      cvm_1se <- cv_fit$cvm[ndx_1se]

      # oof error
      oofPred_1se_vec <- ifelse(
        cv_fit$fit.preval[, ndx_1se] > 0.5, "HCC", "Control"
      )
      oofPred_1se_error <- mean(oofPred_1se_vec != sim_stage_group_vec)

      # test error
      sim_stage_test_samples_vec <- setdiff(rownames(all_lcpm_mtx), sim_stage_samples_vec)
      sim_stage_test_lcpm_mtx <- all_lcpm_mtx[sim_stage_test_samples_vec,]
      sim_stage_test_group_vec <- all_group_vec[sim_stage_test_samples_vec]

      test_pred_1se_vec <- predict(
       cv_fit,
       newx=sim_stage_test_lcpm_mtx,
       s="lambda.1se",
       type="class"
      )
      test_1se_error <- mean(test_pred_1se_vec != sim_stage_test_group_vec)

      # genes
      coef_1se <- coef(
        cv_fit,
        s = "lambda.1se"
      )
      genes <- coef_1se@Dimnames[[1]][coef_1se@i[-1]]

      list(p = nzero_1se, cv_error = cvm_1se, oof_error = oofPred_1se_error, 
           test_error=test_1se_error, genes = genes)
    })
    sim_stage_cv_lst
  })

  # save  sim_cv_lst
  fName <- paste0("sim_", SIMno, "_cv_lst")
  assign(fName, sim_cv_lst)
  save(list = fName, file=file.path("RData", fName))

  message("simulation time: ", round((proc.time() - start_time)[3], 2), "s")
}

```


