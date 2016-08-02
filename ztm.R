library(Rcpp)
Sys.setenv("PKG_CXXFLAGS" = "-fopenmp -lboost_regex")
Sys.setenv("PKG_LIBS" = "-fopenmp -lboost_regex")
sourceCpp("Analytics.cpp")
sourceCpp("ChkpBase.cpp")
sourceCpp("ZeroTradeModel.cpp")
sourceCpp("ZeroTradeModelIO.cpp")
sourceCpp("mainRcpp.cpp")

ztmR <- function(nsim = 1,
                 maxiter = 5, 
                 eps = 0.001,
                 numModelsM = 3,
                 nu = 0.0,
                 OP_w_ij_inp = 0.000001,
                 years = c(1962,1972,1982,1992,2002),
                 numThreads = 1,
                 total_mc_trials = 2,
                 q_filename = "CHECKPOINT_Q.txt",
                 mu_filename = "CHECKPOINT_MU.txt",
                 sigmasq_filename = "CHECKPOINT_SIGMASQ.txt",
                 pi_filename = "CHECKPOINT_PI.txt",
                 zeta_filename = "CHECKPOINT_ZETA.txt",
                 seedOffset = 1) {

    # convert to a string of individual years
    ys <- paste(years, collapse=",")

    # make the Rcpp call
    flag = ztm(nsim, maxiter, eps, numModelsM, nu, OP_w_ij_inp, ys, numThreads, total_mc_trials, 
               q_filename, mu_filename, sigmasq_filename, pi_filename, zeta_filename, seedOffset)

    return(flag)
}

# Examples 

ztmR(numThreads = 6)


# The following are the same, as long as config.txt contains the same parameters
ztmR()
mainR("config.txt", 1)
ztm(1, 5, 0.001, 3, 0.0, 0.000001, "1962,1972,1982,1992,2002", 1, 2, "CHECKPOINT_Q.txt", "CHECKPOINT_MU.txt", "CHECKPOINT_SIGMASQ.txt", "CHECKPOINT_PI.txt", "CHECKPOINT_ZETA.txt", 1)

