
plot_cv_f <- function(cv_fit, Nzero=T, ...) {
 
 suppressPackageStartupMessages(require(glmnet))

 # No nonger used
 #lambda.1se_p <- cv_fit$nzero[cv_fit$lambda == cv_fit$lambda.1se]
 #lambda.min_p <- cv_fit$nzero[cv_fit$lambda == cv_fit$lambda.min]
 
 # Get oof error - cv errors produced by extraction method ARE oof!!!
 ndx_1se <- match(cv_fit$lambda.1se,cv_fit$lambda)
 train_oofPred_1se_vec <- ifelse(
  logistic_f(cv_fit$fit.preval[,ndx_1se]) > 0.5, 'HCC', 'Control')
 train_oofPred_1se_error <- mean(train_oofPred_1se_vec != train_group_vec)

 ndx_min <- match(cv_fit$lambda.min,cv_fit$lambda)
 train_oofPred_min_vec <- ifelse(
  logistic_f(cv_fit$fit.preval[,ndx_min]) > 0.5, 'HCC', 'Control')
 train_oofPred_min_error <- mean(train_oofPred_min_vec != train_group_vec)

 # Get test set error
 test_pred_1se_vec <- predict(
  cv_fit, 
  newx=test_geneExpr_mtx, 
  s="lambda.1se",
  type="class"
 )
 test_pred_1se_error <- mean(test_pred_1se_vec != test_group_vec)
 
 test_pred_min_vec <- predict(
  cv_fit, 
  newx=test_geneExpr_mtx, 
  s="lambda.min",
  type="class"
 )
 test_pred_min_error <- mean(test_pred_min_vec != test_group_vec)
 
  
 plot(
  log(cv_fit$lambda),
  cv_fit$cvm,
  pch=16,col="red",
  xlab='',ylab='',
  ...
 )
 abline(v=log(c(cv_fit$lambda.1se, cv_fit$lambda.min)))
 if(Nzero)
 axis(side=3, tick=F, at=log(cv_fit$lambda), 
  labels=cv_fit$nzero, line = -1
 )
 LL <- 2
 #mtext(side=1, outer=F, line = LL, "log(Lambda)")
 #LL <- LL+1
 mtext(side=1, outer=F, line = LL, paste(
  #ifelse(Nzero, paste("1se p =", lambda.1se_p),''),
  "1se: train =", round(100*cv_fit$cvm[cv_fit$lambda == cv_fit$lambda.1se], 1),
  ##"oof =", round(100*train_oofPred_1se_error, 1), ### REDUNDANT
  "test =", round(100*test_pred_1se_error, 1)
 ))
 LL <- LL+1
 mtext(side=1, outer=F, line = LL, paste(
  #ifelse(Nzero, paste("min p =", lambda.min_p),''),
  "min: train =", round(100*cv_fit$cvm[cv_fit$lambda == cv_fit$lambda.min], 1),
  ##"oof =", round(100*train_oofPred_min_error, 1),  ### REDUNDANT
  "test =", round(100*test_pred_min_error, 1)
 ))
 
 tmp <-
 cbind(
  error_1se = c(
   p = cv_fit$nzero[cv_fit$lambda == cv_fit$lambda.1se],
   train  = 100*cv_fit$cvm[cv_fit$lambda == cv_fit$lambda.1se],
   #train_oof = 100*train_oofPred_1se_error,  ### REDUNANT
   test = 100*test_pred_1se_error),
  error_min = c(
   p = cv_fit$nzero[cv_fit$lambda == cv_fit$lambda.min],
   train  = 100*cv_fit$cvm[cv_fit$lambda == cv_fit$lambda.min],
   #train_oof = 100*train_oofPred_min_error, ### REDUNDSANT
   test = 100*test_pred_min_error)
  )
  # Need to fix names  
  rownames(tmp) <- c('p', 'train', 'test')
  tmp 
}

