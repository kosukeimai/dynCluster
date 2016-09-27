## dynCluster: dynamic Cluster algorithm

# There are two ways to use the code, i.e. use it as CPP program, or create R functions via Rcpp.

* To compile the CPP exectutable ZTM, at the UNIX/LINUX prompt:

   * make clean
   * make

   Note: The dynCluster.so library is for the functions to be called in R.

# To build, check and the package,

    Dependency:
      System packages need to be installed already: OpenMP and boost. 
      R package need to be installed: Rcpp

    Based on the installation, src/Makefile need to be modified accordingly to specify the compiler, the location of header files, and the libraries.

    It is always a good practice to: 'make clean', and then you do the following 

    R CMD build .
    R CMD check .
    R CMD check --as-cran .
    R CMD INSTALL .
 
# Usage
    
    The CPP functions available are 
 
    (1) mainRcpp - use the parameters specified in the default file config.txt
    (2) ztm - use the parameters provided in the function input 

    Here is the R code to load the CPP functions (assume it is from the installed R library):
      path <- paste0(.libPaths()[1], "/dynCluster", "/libs/dynCluster.so")
      mainRcpp <- Rcpp:::sourceCppFunction(function(configTxt, randomSeed) {}, FALSE, dyn.load(path), 'dynCluster_mainRcpp')
      ztm <- Rcpp:::sourceCppFunction(function(nsim, maxiter, eps, numModelsM, nu, OP_w_ij_inp, inputYears, numThreads, total_mc_trials, q_filename, mu_filename, sigmasq_filename, pi_filename, zeta_filename, seedOffset) {}, FALSE, dyn.load(path), 'dynCluster_ztm')

    The corresponding R wrap functions are

    (1) mainZTM - calls mainRcpp with the directory that contains the config.txt
    (2) testExample - a test function that runs the provided ZTM examplea

# Examples

    Assumes that you are in the directory of dynCluster, and the package is already installed.

    (1) Under UNIX:

        cd example
        ./../src/ZTM ./config.txt 1

        or for SLURM

        cd example
        sbatch skip9.slurm

    (2) In R:

        library(dynCluster)
        testExample("./example", nThreads=4, comeBack=TRUE)

        Or,

        library(dynCluster)
        mainZTM("./example", comeBack=TRUE)
    
    Note: The default value of comeBack is FALSE, meaning it will stay in the directory specified by the first parameter, instead of returning to the current directory after the job is done.


