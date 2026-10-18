# Refit corrected receiver-only HCC models from frozen patient-level features.
# Internal scaling and penalty selection are nested; external models are locked.
run_clinical <- function(d,out) {
internal <- copy(d$hcc_internal); external <- copy(d$hcc_external)
auc_pairwise <- function(y,p)mean(outer(p[y==1],p[y==0],function(a,b)(a>b)+.5*(a==b)))
ridge_predict <- function(train_x, train_y, test_x, lambda) {
  mu <- colMeans(train_x)
  sig <- apply(train_x, 2, sd)
  sig[!is.finite(sig) | sig == 0] <- 1
  tx <- sweep(sweep(train_x, 2, mu, "-"), 2, sig, "/")
  vx <- sweep(sweep(test_x, 2, mu, "-"), 2, sig, "/")
  X <- cbind(1, tx)
  objective <- function(beta) {
    eta <- as.numeric(X %*% beta)
    sum(log1p(exp(-abs(eta))) + pmax(eta, 0) - train_y * eta) +
      lambda * sum(beta[-1]^2) / 2
  }
  gradient <- function(beta) {
    eta <- as.numeric(X %*% beta)
    g <- as.numeric(crossprod(X, plogis(eta) - train_y))
    g[-1] <- g[-1] + lambda * beta[-1]
    g
  }
  fit <- optim(rep(0, ncol(X)), objective, gradient, method = "BFGS",
               control = list(maxit = 500, reltol = 1e-10))
  plogis(cbind(1, vx) %*% fit$par)[, 1]
}

choose_lambda <- function(X, y, grid = c(0.01, 0.1, 1, 10, 100)) {
  losses <- sapply(grid, function(lam) {
    pred <- vapply(seq_along(y), function(i) {
      ridge_predict(X[-i, , drop = FALSE], y[-i], X[i, , drop = FALSE], lam)
    }, numeric(1))
    -sum(y * log(pmax(pred, 1e-8)) + (1 - y) * log(pmax(1 - pred, 1e-8)))
  })
  grid[which.min(losses)]
}

sets <- list('Direct killing only'=c('delta_apoptosis','post_apoptosis_q95'),
  'Pretreatment state'=paste0(c('interferon','cytotoxic','fibrosis','hypoxia'),'_Pre'),
  'SPF geometry only'=paste0('field_slope_',c('interferon','cytotoxic','fibrosis','hypoxia')),
  'Direct + SPF'=c('delta_apoptosis','post_apoptosis_q95',paste0('field_slope_',c('interferon','cytotoxic','fibrosis','hypoxia'))))
pred <- rbindlist(lapply(names(sets),function(model){X<-as.matrix(internal[,sets[[model]],with=FALSE]);y<-internal$outcome;p<-chosen<-numeric(length(y));for(i in seq_along(y)){lam<-choose_lambda(X[-i,,drop=FALSE],y[-i]);p[i]<-ridge_predict(X[-i,,drop=FALSE],y[-i],X[i,,drop=FALSE],lam);chosen[i]<-lam};data.table(patient_id=internal$patient_id,Response=internal$Response,outcome=y,model,probability=p,selected_lambda=chosen)}))
metric <- pred[,.(auc=auc_pairwise(outcome,probability),brier=mean((outcome-probability)^2),log_loss=-mean(outcome*log(pmax(probability,1e-8))+(1-outcome)*log(pmax(1-probability,1e-8))),accuracy=mean((probability>=.5)==outcome)),by=model]
fwrite(pred,file.path(out,'hcc_nested_loo_predictions.tsv'),sep='\t');fwrite(metric,file.path(out,'hcc_nested_loo_metrics.tsv'),sep='\t');print(metric)
ridge_fit <- function(X, y, lambda) {
  mu <- colMeans(X); sig <- apply(X, 2, sd); sig[!is.finite(sig) | sig == 0] <- 1
  Z <- sweep(sweep(X, 2, mu, "-"), 2, sig, "/"); D <- cbind(1, Z)
  fn <- function(beta) {
    eta <- as.numeric(D %*% beta)
    sum(log1p(exp(-abs(eta))) + pmax(eta, 0) - y * eta) + lambda * sum(beta[-1]^2) / 2
  }
  gr <- function(beta) {
    g <- as.numeric(crossprod(D, plogis(D %*% beta) - y)); g[-1] <- g[-1] + lambda * beta[-1]; g
  }
  fit <- optim(rep(0, ncol(D)), fn, gr, method = "BFGS")
  list(beta = fit$par, mu = mu, sig = sig)
}
ridge_pred <- function(fit, X) {
  Z <- sweep(sweep(X, 2, fit$mu, "-"), 2, fit$sig, "/")
  plogis(cbind(1, Z) %*% fit$beta)[, 1]
}
choose_lambda <- function(X, y, grid = c(0.01, 0.1, 1, 10, 100)) {
  loss <- sapply(grid, function(lam) {
    p <- vapply(seq_along(y), function(i) {
      ridge_pred(ridge_fit(X[-i, , drop = FALSE], y[-i], lam), X[i, , drop = FALSE])
    }, numeric(1))
    -sum(y * log(pmax(p, 1e-8)) + (1-y) * log(pmax(1-p, 1e-8)))
  })
  grid[which.min(loss)]
}
auc <- function(y, p) mean(outer(p[y == 1], p[y == 0], function(a,b) (a>b) + 0.5*(a==b)))

external_sets <- list('Post-treatment direct'='post_apoptosis_q95','SPF geometry'=paste0('field_slope_',c('interferon','cytotoxic','fibrosis','hypoxia')),'Direct + SPF'=c('post_apoptosis_q95',paste0('field_slope_',c('interferon','cytotoxic','fibrosis','hypoxia'))))
pred <- rbindlist(lapply(names(external_sets),function(model){f<-external_sets[[model]];Xi<-as.matrix(internal[,..f]);Xe<-as.matrix(external[,..f]);lam<-choose_lambda(Xi,internal$outcome);fit<-ridge_fit(Xi,internal$outcome,lam);p<-ridge_pred(fit,Xe);saveRDS(list(features=f,fit=fit,lambda=lam),file.path(out,paste0('locked_',gsub('[^a-z0-9]+','_',tolower(model)),'.rds')));data.table(sample_id=external$sample_id,Response=external$Response,outcome=external$outcome,model,probability=p,lambda=lam)}))
metric<-pred[,{obs<-auc(outcome,probability);aa<-combn(.N,sum(outcome),simplify=FALSE);null<-vapply(aa,function(ii){yy<-integer(.N);yy[ii]<-1;auc(yy,probability)},numeric(1));list(auc=obs,brier=mean((outcome-probability)^2),accuracy=mean((probability>=.5)==outcome),exact_auc_p=mean(null>=obs-1e-12),exact_assignments=length(null))},by=model]
fwrite(pred,file.path(out,'hcc_external_locked_predictions.tsv'),sep='\t');fwrite(metric,file.path(out,'hcc_external_locked_metrics.tsv'),sep='\t');print(metric)

invisible(metric)
}
