{
  if ('dynCluster' %in% installed.packages()[,1]) {
    path <- paste0(.libPaths()[1], "/dynCluster", "/lib/dynCluster.so")
  } else {
    path <- "./lib/dynCluster.so"
  }

  print(path)

  ztm <- Rcpp:::sourceCppFunction(function(nsim, maxiter, eps, numModelsM, nu, OP_w_ij_inp, inputYears, numThreads, total_mc_trials, q_filename, mu_filename, sigmasq_filename, pi_filename, zeta_filename, seedOffset) {}, FALSE, dyn.load(path), 'sourceCpp_1_ztm')

  mainRcpp <- Rcpp:::sourceCppFunction(function(configTxt, randomSeed) {}, FALSE, dyn.load(path), 'sourceCpp_1_mainRcpp')

  rm(path)

}
