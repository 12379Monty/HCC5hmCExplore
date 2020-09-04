# The bet on sparsity {#explore-sparsity}

## CV analysis setup {-}

```{r setParameters}

K_FOLD <- 10
trainP <- 0.8
EPS <- 0.02    # Have no idea what "small" epsilon means

```

First we divide the analysis dataset into `train` and `test` in a `r trainP/(1-trainP)`:1 ratio.  

```{r getTrainVal, cache=T, cache.vars=c('train_sampID_vec', 'test_sampID_vec','train_group_vec','test_group_vec','train_lcpm_mtx','test_lcpm_mtx')}

set.seed(1)
train_sampID_vec <- with(AF_dgel$samples,
AF_dgel$samples$sampID[caret::createDataPartition(y=group, p=trainP, list=F)]
)

test_sampID_vec <- with(AF_dgel$samples,
setdiff(sampID, train_sampID_vec)
)

train_group_vec <- AF_dgel$samples[train_sampID_vec, 'group']
test_group_vec <- AF_dgel$samples[test_sampID_vec, 'group']

knitr::kable(table(train_group_vec),
  caption="Train set") %>%
   kableExtra::kable_styling(full_width = F)

knitr::kable(table(test_group_vec),
  caption="Test set") %>%
   kableExtra::kable_styling(full_width = F)

train_lcpm_mtx <- t(lcpm_mtx[,train_sampID_vec])
test_lcpm_mtx <- t(lcpm_mtx[,test_sampID_vec])

```

We explore some glmnet fits and the "bet on sparsity"    
  
* Consider models:  
    - lasso: $\alpha = 1.0$ - sparse model  
    - ridge $\alpha = 0$ - shrunken coefficients model
    - elastic net:  $\alpha = 0.5$  - semi sparse model
    - lassoC: $\alpha = 1-\epsilon =$ `r 1- EPS` - lasso for correlated predictors  
* Does the relaxed lasso improve performance?    
* Does the shrunken relaxed lasso improve performance  
* How sparse is the model undelying best 5hmC classifier for Early HCC vs Control?    
* Is the degree of sparsity, or the size of the model, a stable feature of the problem and data set?  

In this analysis, we will only evaluate models in terms of 
model size, stability and performance.  We leave the question
of significance testing of hypotheses about model parameters
completely out.  See Lockhart et al. (2014) [@Lockhart:2014aa]
and Wassermam (2014) [@Wasserman:2014aa] for a discussion of this topic.


Next we create folds for `r K_FOLD`-fold cross-validation of models fitted to
training data.  We'll use caret::createFolds to assign samples
to folds while keeping the outcome ratios constant across folds.


```{r getTrainFolds, cache=T, cache.vars='train_foldid_vec'}
# This is too variable, both in terms of fold size And composition
#foldid_vec <- sample(1:10, size=length(train_group_vec), replace=T)

set.seed(1)
train_foldid_vec <- caret::createFolds(
 factor(train_group_vec), 
 k=K_FOLD,
 list=F)

knitr::kable(sapply(split(train_group_vec, train_folds_vec), 
  table), caption="training samples fold composition") %>%
   kableExtra::kable_styling(full_width = F)
 
```

Note that the folds identify samples that are left-out of the training
data for each fold fit.


## Fit and compare models {-}


* cross-validated accuracy
* test set accuracy  
* sparsity
    - for lasso, enet and lassoC, examine number of selected variables

Although “the one standard error rule” can produce a model with fewer predictors, it usually results in increased MSE and more biased parameter estimates
(see Engebretsen et al. (2019) [@Engebretsen:2019aa] for example).
We will look at both the minimum cv error and the one standard error rule model 
preformance.

```{r fitModels, cache=T, cache.vars=c('cv_lasso', 'cv_ridge', 'cv_enet', 'cv_lassoC')}

require(doMC)
registerDoMC(cores=14)

start_time <-  proc.time()

cv_lasso <- glmnet::cv.glmnet(
 x=train_lcpm_mtx,
 y=train_group_vec,
 foldid=train_foldid_vec,
 alpha=1,
 family='binomial', 
 type.measure = "class")

message("lasso time: ", round((proc.time() - start_time)[3],2),"s")

start_time <-  proc.time()

cv_ridge <- glmnet::cv.glmnet(
 x=train_lcpm_mtx,
 y=train_group_vec,
 foldid=train_foldid_vec,
 alpha=0,
 family='binomial', 
 type.measure = "class")

message("ridge time: ", round((proc.time() - start_time)[3],2),"s")

start_time <-  proc.time()

cv_enet <- glmnet::cv.glmnet(
 x=train_lcpm_mtx,
 y=train_group_vec,
 foldid=train_foldid_vec,
 alpha=0.5,
 family='binomial',
 type.measure = "class")

message("enet time: ", round((proc.time() - start_time)[3],2),"s")

start_time <-  proc.time()

cv_lassoC <-  glmnet::cv.glmnet(
 x=train_lcpm_mtx,
 y=train_group_vec,
 foldid=train_foldid_vec,
 alpha=1-EPS,
 family='binomial',
 type.measure = "class")

message("lassoC time: ", round((proc.time() - start_time)[3],2),"s")

```


