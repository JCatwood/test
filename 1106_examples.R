# When writing codes, please keep the following in mind:
# - Easy to read, including organization, variable naming, documentation
# - Try to use existing codes/packages, also for efficiency, 
#   try to reduce loops in R
# - Organize your codes well in Github and reuse them in the future. 
#   Also this is good for reproducing your results

cond_corr <- function(j,c_i) { # corr( y_n, y_j | y_c(i) )
  ind_set <- c(j, n)
  cond_covmat <- covmat[ind_set, ind_set] - covmat[ind_set, c_i] %*% 
    solve(covmat[c_i, c_i]) %*% covmat[c_i, ind_set]
  cond_corr <- cond_covmat[1, 2] / sqrt(cond_covmat[1, 1] * cond_covmat[2, 2])
  abs(cond_corr)
}

cond_var <- function(c_i) { # var( y_n | y_c(i) )
  covmat[n, n] - covmat[n, c_i] %*% 
    solve(covmat[c_i, c_i]) %*% covmat[c_i, n]
}

fig <- function(loc,k_hat,opt_set,seed){
  plot(loc[,1],loc[,2],cex=2,main=paste("seed =",seed))
  points(loc[n,1],loc[n,2],pch=8,col='red')
  points(loc[k_hat,1],loc[k_hat,2],pch=2,col='blue')
  points(loc[opt_set,1],loc[opt_set,2],pch=17,col='green') # optimal conditioning set
  
  NNOrd <- GPvecchia::order_dist_to_point(locs = loc[-n,],
                                          loc0 = t(as.matrix(loc[n,])))
  for(i in seq(m)){
    points((loc[-n,])[NNOrd[i],1],(loc[-n,])[NNOrd[i],2],cex=3,col='red') # NN
  }
}

m <- 3

### random covmat ####
n <- 5
seeds <- c()
for (s in seq(100) ) {
  set.seed(s)
  tmp_mat <- matrix(runif(n * n), n, n)
  covmat <- tmp_mat %*% t(tmp_mat)
  
  for(j in 1 : m) {
    if (j == 1) {
      k_hat <- which.max(covmat[1:(n - 1), n]/sqrt(diag(covmat)[-n]))
    } else {
      v <- seq(n-1)[-k_hat]
      k <- sapply(v, function(j) cond_corr(j,k_hat))
      k_hat <- c(k_hat,v[which.max(k)])
    }
  }
  
  cond_sets <- utils::combn(seq(n-1),m) # cols are all possible subsets of {1,...,i-1} with size m
  cond_vars <- apply(cond_sets,2,cond_var)
  opt_set <- cond_sets[,which.min(cond_vars)]
  
  if (all(sort(k_hat)==opt_set)) {seeds <- seeds}
  else {seeds <- c(seeds,s)}
}


### matern covmat ####
par(mfrow=c(2,2))
covpar <- c(1,1,2.5,0) #c(var,alpha,smoothness,nugget)
for (seed in c(12, 49, 128, 189)){
  n <- 5
  set.seed(seed)
  loc <- cbind(runif(n^2),runif(n^2))
  
  covmat <- GpGp::matern_isotropic(covparms = covpar, locs = loc)
  
  n <- dim(covmat)[1] #sample size
  
  for(j in 1 : m) {
    if (j == 1) {
      k_hat <- which.max(covmat[1:(n - 1), n]/sqrt(diag(covmat)[-n]))
    } else {
      v <- seq(n-1)[-k_hat]
      k <- sapply(v, function(j) cond_corr(j,k_hat))
      k_hat <- c(k_hat,v[which.max(k)])
    }
  }
  
  cond_sets <- utils::combn(seq(n-1),m) # cols are all possible subsets of {1,...,i-1} with size m
  cond_vars <- apply(cond_sets,2,cond_var)
  opt_set <- cond_sets[,which.min(cond_vars)]
  
  fig(loc,k_hat,opt_set,seed)
}

covpar <- c(1,1,.9,0) #c(var,alpha,smoothness,nugget)
for (seed in c(49, 167, 189, 315)){
  n <- 5
  set.seed(seed)
  loc <- cbind(runif(n^2),runif(n^2))
  
  covmat <- GpGp::matern_isotropic(covparms = covpar, locs = loc)
  
  n <- dim(covmat)[1] #sample size
  
  for(j in 1 : m) {
    if (j == 1) {
      k_hat <- which.max(covmat[1:(n - 1), n]/sqrt(diag(covmat)[-n]))
    } else {
      v <- seq(n-1)[-k_hat]
      k <- sapply(v, function(j) cond_corr(j,k_hat))
      k_hat <- c(k_hat,v[which.max(k)])
    }
  }
  
  cond_sets <- utils::combn(seq(n-1),m) # cols are all possible subsets of {1,...,i-1} with size m
  cond_vars <- apply(cond_sets,2,cond_var)
  opt_set <- cond_sets[,which.min(cond_vars)]
  
  fig(loc,k_hat,opt_set,seed)
}