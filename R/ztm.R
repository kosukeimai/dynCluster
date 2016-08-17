library(Rcpp)
Sys.setenv("PKG_CXXFLAGS" = "-fopenmp -lboost_regex")
Sys.setenv("PKG_LIBS" = "-fopenmp -lboost_regex")
# These will be compiled by the Makefile, no need to use Rcpp compile here as the functions in them are all internal
# sourceCpp("./src/Analytics.cpp")
# sourceCpp("./src/ChkpBase.cpp")
# sourceCpp("./src/ZeroTradeModel.cpp")
# sourceCpp("./src/ZeroTradeModelIO.cpp")
sourceCpp("./src/mainRcpp.cpp")

#
# This will automatically generate the binary objects of mainRcpp and ztm (same as click "Build" in RStudio:
#
#     devtools::document(roclets=c('rd', 'collate', 'namespace'))
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'
#


#' Title This is a function that calls the CPP ZTM code via a list of parameters
#'
#' @param nsim                number of simulation
#' @param maxiter             max number of iteration
#' @param eps                 epsilon
#' @param numModelsM          number of models
#' @param nu                  nu
#' @param OP_w_ij_inp         w_ij
#' @param years               data years (such as 1962,1972,1982,1992,2002). Data file names are like dynamic_1962.csv
#' @param numThreads          number of threads to run parallel (OpenMP)
#' @param total_mc_trials     number of MC trials
#' @param q_filename          Q_cluster file name
#' @param mu_filename         MU_cluster file name
#' @param sigmasq_filename    SIGMA_cluster file name
#' @param pi_filename         PI_cluster file name
#' @param zeta_filename       ZETA file name
#' @param seedOffset          random seed offset
#'
#' @return                    1 for error and 0 for success
#' @export
#'
#' @examples
#' library(Rcpp)
#' Sys.setenv("PKG_CXXFLAGS" = "-fopenmp -lboost_regex")
#' Sys.setenv("PKG_LIBS" = "-fopenmp -lboost_regex")
#' \dontrun{sourceCpp("./src/mainRcpp.cpp")}
#' \dontrun{setwd("./example")}
#' \dontrun{ztm()}
#' \dontrun{setwd("./..")}
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

#' Title This is an example of using ztmR, assuming the required data is in the dataDirectory
#'
#' @param dataDirectory    the working directory contains the data
#' @param nThreads         number of threads
#' @param comeBack         return from the working directory
#'
#' @return                 1 for error and 0 for success
#' @export
#'
#' @examples
#' library(Rcpp)
#' Sys.setenv("PKG_CXXFLAGS" = "-fopenmp -lboost_regex")
#' Sys.setenv("PKG_LIBS" = "-fopenmp -lboost_regex")
#' \dontrun{sourceCpp("./src/mainRcpp.cpp")}
#' \dontrun{testExample("./example", nThreads=4, comeBack=TRUE)}
testExample <- function(dataDirectory = "./", nThreads = 6, comeBack = FALSE) {
  # take care of the extreme input
  if (is.null(dataDirectory)) {
    dataDirectory = "./"
  }
  if (is.na(dataDirectory)) {
    dataDirectory = "./"
  }
  if (dataDirectory == "") {
    dataDirectory = "./"
  }

  # hold the value for the current directory, change directory if needed
  prevDirectory = NA
  if (!dir.exists(dataDirectory)) {
    print(paste("data directory does not exist:", dataDirectory))
    return(1)
  } else {
    if (dataDirectory != "./") {
      prevDirectory = getwd()
      setwd(dataDirectory)
    }

    # the parametes are preset by the R API interface as default values
    flag = ztmR(numThreads = nThreads)

    # The following are doing the same, as long as config.txt contains the same parameters
    # ztmR()
    # mainRcpp("config.txt", 1)
    # ztm(1, 5, 0.001, 3, 0.0, 0.000001, "1962,1972,1982,1992,2002", 1, 2, "CHECKPOINT_Q.txt", "CHECKPOINT_MU.txt", "CHECKPOINT_SIGMASQ.txt", "CHECKPOINT_PI.txt", "CHECKPOINT_ZETA.txt", 1)

    # recover the previous working direcoty
    if (!is.na(prevDirectory) & comeBack) {
      if (prevDirectory != "./") {
        setwd(prevDirectory)
      }
    }
    return(flag)
  }
}

#' Title This is a function that calls the CPP ZTM code via a config.txt which specifies parameters and data
#'
#' @param dataDirectory     the working directory contains the data
#' @param comeBack          return from the working directory
#'
#' @return                  1 for error and 0 for success
#' @export
#'
#' @examples
#' library(Rcpp)
#' Sys.setenv("PKG_CXXFLAGS" = "-fopenmp -lboost_regex")
#' Sys.setenv("PKG_LIBS" = "-fopenmp -lboost_regex")
#' \dontrun{sourceCpp("./src/mainRcpp.cpp")}
#' \dontrun{mainZTM("./example", comeBack=TRUE)}
mainZTM <- function(dataDirectory = "./", comeBack = FALSE) {
  # take care of the extreme input
  if (is.null(dataDirectory)) {
    dataDirectory = "./"
  }
  if (is.na(dataDirectory)) {
    dataDirectory = "./"
  }
  if (dataDirectory == "") {
    dataDirectory = "./"
  }

  # check if the directory is valid
  if (!dir.exists(dataDirectory)) {
    print(paste("data directory does not exist:", dataDirectory))
    return(1)
  } else {
    # set working directory
    prevDirectory = NA
    if (dataDirectory != "./") {
      prevDirectory = getwd()
      setwd(dataDirectory)
    }

    # Note: in this configuration, we use the parameters in config.txt, the same way as the C++ code ZTM
    if (file.exists("config.txt")) {
      flag = mainRcpp("config.txt", 1)
    } else {
      print(paste("config.txt does not exist in", getwd()))
      flag = 1
    }

    # recover the previous working direcoty
    if (!is.na(prevDirectory) & comeBack) {
      if (prevDirectory != "./") {
        setwd(prevDirectory)
      }
    }
    return(flag)
  }
}


