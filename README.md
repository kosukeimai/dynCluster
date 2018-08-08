# dynCluster: Dynamic Clustering Algorithm [![Build Status](https://travis-ci.org/kosukeimai/dynCluster.svg?branch=master)](https://travis-ci.org/kosukeimai/dynCluster)

Authors: [In Song Kim](http://web.mit.edu/insong/www/index.html), [Steven Liao](https://www.stevenliao.org/), [Kosuke Imai](http://imai.princeton.edu/)

For more details of the method and applications, see our paper:

+ [Measuring Trade Profile with Granular Product-level Trade Data](https://www.stevenliao.org/uploads/2/5/6/9/25699716/bigtrade.pdf)

## Running dynCluster on Amazon Web Services: A Toy Example
1. Install dynCluster on [Amazon Web Services (AWS)](https://aws.amazon.com/). See our [Wiki page](https://github.com/kosukeimai/dynCluster/wiki) for step-by-step instructions.

2. Once dynCluster is installed, the R code below creates a small simulated dataset that follows the data generating process described in our paper. The toy example runs on *t2.micro* AWS instances, which is free.
  
    + Set up: `10` countries (`90` dyads) trading `40` products over `10` time periods. Dyads belong to `3` different clusters.
        ```R
        # set clean slate
        rm(list = ls())
        date()
        
        # load packages
        pkg <- c("foreign", "dplyr", "ggplot2", "stringr", "tidyr", "purrr", 
             "RColorBrewer", "plotrix")
        install.packages(pkg, Sys.getenv("R_LIBS_USER"), repos = "http://cran.case.edu")
        lapply(pkg, require, character.only = TRUE)
        
        # set simulation ID number
        sim.id <- 25
        
        # set file directories
        MAIN_DIR <- "~/dynCluster-master"
        SIM_INPUT <- paste(MAIN_DIR, "/sim-", sim.id, sep = "")
        
        # set seed
        set.seed(1234)
        
        # set constants
        N <- 10 # number of units
        T <- 10 # number of time periods
        M <- 3 # number of clusters
        K <- 40 # number of products
        ndyads <- choose(N, 2)*2 # number of dyads
        
        # set M by M transition matrix 
        # where (i,j)th element specifies the probabilty that a dyad in 
        # cluster i in time t-1 will move to cluster j in time t (each 
        # row of this matrix should sum up to 1)
        P <- diag(.85, M)
        transition <- 0.15/(M-1)
        P[P!=.85] <- transition
        
        # function to simulate Markov transitions
        run.sim <- function(P, num.periods = T, first.cluster = 1) {
          
          # number of possible clusters
          num.clusters <- nrow(P)
          
          # stores the clusters Z_t through time
          clusters <- numeric(num.periods)
          
          # initialize variable for first state 
          clusters[1] <- first.cluster
          
          for(t in 2:num.periods) {
            
            # probability vector to simulate next state X_{t+1}
            p <- P[clusters[t-1], ]
            
            ## draw from multinomial and determine state
            clusters[t] <-  which(rmultinom(n = 1, size = 1, p) == 1)
          }
          return(clusters)
        }
        ```
    + Create dyadic cluster membership data, i.e., the cluster assignment of each dyad for each time period. The ultimate goal of this toy example is to see how well dynCluster can recover the true cluster membership data.
        ```R
        # create dyad-time indices
        data <- expand.grid(cty1 = seq(1:N), 
                            cty2 = seq(1:N), 
                            time = seq(1,T))
        
        data <- data %>% 
          filter(cty1 != cty2) %>% # drop rows with cty1 == cty2
          mutate(dyad = ifelse(as.numeric(cty1) < as.numeric(cty2), # create dyad ID
                               paste(cty1, cty2, sep = "_"), 
                               paste(cty2, cty1, sep = "_"))) %>%
          arrange(time, cty1, cty2) # sort
        
        # extract dyads from the first period and initiate their cluster membership
        d1 <- data %>%
          filter(time == 1) %>%
          distinct(dyad, .keep_all = TRUE) %>%
          mutate(z1 = sample(1:M, ndyads/2, replace = TRUE))
        
        # simulate cluster membership for T transitions
        trans.list <- list()
        for(i in 1:nrow(d1)){
          sim.out <- run.sim(P, T, d1$z1[i])
          trans.list[[i]] <- sim.out # add it to list
        }
        
        # combine into a dataframe
        trans.df <- do.call(rbind, trans.list)
        
        # add time labels
        colnames(trans.df) <- paste("z", seq(1:T), sep = "")
        
        # combine dyads and cluster membership info
        df.1 <- cbind(d1 %>% select(-z1, -time), trans.df)
        
        # create similar df but flipping the direction of dyads
        df.2 <- cbind(df.1$cty2, df.1$cty1, df.1)
        df.2 <- df.2 %>%
          select(-cty1, -cty2) %>%
          rename(cty1 = `df.1$cty2`,
                 cty2 = `df.1$cty1`)
        
        # stack dfs to get the full directed-dyad cluster membership data
        # each row contains one directed-dyad and its cluster membership in each period
        # colnames: cty1, cty2, dyad, z1, z2, ..., z_t
        df <- rbind(df.1, df.2) %>% 
          arrange(dyad)
        
        # view data  
        head(df, 10)
        ##   cty1 cty2 dyad z1 z2 z3 z4 z5 z6 z7 z8 z9 z10
        ## 1    1   10 1_10  2  2  2  2  2  2  2  2  2   2
        ## 2   10    1 1_10  2  2  2  2  2  2  2  2  2   2
        ## 3    1    2  1_2  1  1  1  1  1  1  1  1  1   1
        ## 4    2    1  1_2  1  1  1  1  1  1  1  1  1   1
        ## 5    1    3  1_3  2  2  2  2  2  2  2  2  2   2
        ## 6    3    1  1_3  2  2  2  2  2  2  2  2  2   2
        ## 7    1    4  1_4  2  2  2  2  2  1  1  1  1   1
        ## 8    4    1  1_4  2  2  2  2  2  1  1  1  1   1
        ## 9    1    5  1_5  2  2  2  2  2  2  2  2  2   2
        ## 10   5    1  1_5  2  2  2  2  2  2  2  2  2   2
          
        # write/save data
        if (!dir.exists(SIM_INPUT)){
          print("Directory does not exist! Creating...")
          dir.create(SIM_INPUT)
          write.csv(df, paste(SIM_INPUT, "/truth.csv", sep = ""), row.names = FALSE)
          
        } else {
          write.csv(df, paste(SIM_INPUT, "/truth.csv", sep = ""), row.names = FALSE)
        }
        ```
        
      + Simulate dyadic product-level trade data using the cluster membership data above.
        ```R
        # 1. set q_{kz}
        Q <- matrix(NA, nrow = K, ncol = M)
        Q[,1] <- runif(K, min=0.8, max=1) # cluster 1 is sparse trade
        Q[,2] <- runif(K, min=0.5, max=0.8)
        Q[,3] <- runif(K, min=0, max=0.2)
        
        # 2. set mu_{kz}
        MU <- matrix(NA, nrow = K, ncol = M)
        MU[,1] <- runif(K, min=100, max=1000)
        MU[,2] <- runif(K, min=900, max=700000)
        MU[,3] <- runif(K, min=300000, max=10000000)
        
        # generate product-trade for each directed dyad
        X <- matrix(NA, nrow = nrow(data), ncol = K)
        
        for(i in 1:nrow(data)){
            cty1 <- data$cty1[i]
            cty2 <- data$cty2[i]
            time <- data$time[i]
        
            #for(t in 1:T){
                idx <- which(df$cty1 == cty1 & df$cty2 == cty2)
                z_ijt <- df[idx,time+3] # extract cluster for given period 
            #}
        
            X_ijtk <- rep(NA, length=K)
            for(k in 1:K){
                q_kz <- Q[k,z_ijt]
                D_ijtk <- rbinom(1, 1, prob = q_kz)
                if(D_ijtk == 0){ # positive trade
                    W_ijtk <- max(100, rnorm(1,MU[k,z_ijt], 100))
                    X_ijtk[k] <- W_ijtk
                } else {
                    X_ijtk[k] <- 0
                }
            }
            X[i,] <- X_ijtk
        
        }
        
        # combine dyad, time, and trade data
        data <- cbind(data, X)
        
        # reorder vars
        data.trade <- data %>% 
          select(time, cty1, cty2, everything()) %>%
          select(-dyad)
          #select(time, cty1, cty2, `1`:`20`)
        
        # change var names (package very picky about names)
        prods <- paste("SITC0_", 1:K, sep="")
        colnames(data.trade) <- c("year", "importer_ISO", "exporter_ISO", prods)
        
        # view first 5 rows and 10 columns of data
        data.trade[1:5, 1:10]
        ##   year importer_ISO exporter_ISO SITC0_1 SITC0_2 SITC0_3  SITC0_4 SITC0_5      SITC0_6   SITC0_7
        ## 1    1            1            2       0       0       0      0.0       0     177.0528       0.0
        ## 2    1            1            3       0       0       0 664344.6       0  536712.5031  133149.7
        ## 3    1            1            4       0       0       0      0.0       0  536385.8453       0.0
        ## 4    1            1            5  372390       0       0      0.0       0  536450.8843       0.0
        ## 5    1            1            6 3171746 2797487 4872051 981809.8 2497946 1412988.9059 1075698.3
        
        # write/save data by time period (1961, 1962, ..., 1970)
        for(t in 1:T){
          sub <- data.trade[data.trade$year == t,]
          filename <- paste(SIM_INPUT, "/dynamic_", as.character(1960+t), ".csv", sep="")
          write.table(sub, filename, sep ="\t", row.names = FALSE, quote= FALSE)
        }
        ```
      + Note that users can replace these files (`dynamic_1961.csv`, `dynamic_1962.csv`, ..., `dynamic_1970.csv`) with their own data in real applications. However, file names and format need to be exactly the same in the current version for dynCluster to work. For example, file names need to begin with `dynamic_` and column names need to start with `SITC0_` (4 characters + 1 integer + underscore).

3. Besides the dyadic trade data, dynCluster also requires `config.txt` to run. The file contains the parameters used in the simulation.

    ![](images/config.png)

      + `NSIM`: the number of simulation. Each simulation will use different starting values. The final output will select the simulation that yields the highest log-likelihood.
      + `MAXITER`: the max number of iterations
      + `EPS`: the epsilon for convergence
      + `Z`: the number of clusters
      + `NU`: the scaling parameter for zero-trade probability
      + `OP_w_ij`: the mixture responsibilities
      + `file`: the names of the input files excluding "dynamic_", e.g., 1961, 1962, ..., 1970
      + `threads`: the number of threads to run parallel (OpenMP)
      + `MCITERATIONS`: the number of MC iterations
      + `QFILE`: Q_cluster file to output Bernoulli probability for zero trade
      + `MUFILE`: MU_cluster file to output the mean of the normal distribution
      + `SIGMAQFILE`: SIGMA_cluster file to output variance of the normal distribution
      + `PIEFILE`: PI_cluster file to output mixture probability
      + `ZETAFILE`: ZETA file to output mixture responsibilities
      + `TFFILE`: 
      + `PMATFILE`:

4. To run dynCluster,
    + Copy `config.txt` to the folder containing the simulated data and adjust the parameters as needed.
        ```sh
        cp ~/dynCluster-master/example/config.txt ~/dynCluster-master/sim-25
        ```
    + In R, run the function `mainZTM`. This function wraps and calls C++ functions (e.g., `mainRcpp`) from dynCluster.
        ```R
        ptm <- proc.time() # start the clock
        mainZTM("./sim-25", comeBack=TRUE)
        proc.time() - ptm # stop the clock
        ```
    + Results will be saved to the same folder

## Assessing the Performance of dynCluster using the Toy Example
1. Set up parameters in R
    ```R
    # clean slate
    rm(list = ls())
    date()
          
    # load packages
    pkg <- c("foreign", "dplyr", "ggplot2", "stringr", "tidyr", "purrr", 
             "RColorBrewer", "plotrix")
    lapply(pkg, require, character.only = TRUE)
          
    # set simulation ID number
    sim.id <- 25
          
    # set parameters
    n.cluster <- 3 # the Z in config.txt
          
    # set main directory
    MAIN_DIR <- "~/dynCluster-master"
          
    # set other directories
    SIM_FOLDER <- paste(MAIN_DIR, "/sim-", sim.id, sep = "")
          
    # set file directory
    TRUTH_FILE <- paste(SIM_FOLDER, "/truth.csv", sep = "")
    NEWTF_FILE <- paste(SIM_FOLDER, "/NewTF.csv", sep = "")
    TF_TRUTH_FILE <- paste(SIM_FOLDER, "/TF-truth.csv", sep = "")
          
    # source functions
    source(paste(MAIN_DIR, "R/functions/process-functions.R", sep = "/"))
    source(paste(MAIN_DIR, "R/functions/visualize-functions.R", sep = "/"))

    # list all files in RESULTS folder
    files <- list.files(SIM_FOLDER)
    files
          
    # extract vector of years from raw input data
    raw.data.files <- list.files(SIM_FOLDER)
    raw.data.files <- raw.data.files[grep("^dynamic_", raw.data.files)]
    raw.data.files <- str_replace_all(raw.data.files, "dynamic_", "")
    raw.data.files <- as.numeric(str_replace_all(raw.data.files, "\\.csv", ""))
    raw.data.files
          
    # set parameters
    periods <- seq(0, n_distinct(raw.data.files)-1, by = 1)
    periods
    years <- raw.data.files
    years
    cl.names <- seq(0, n.cluster - 1, by = 1)
    cl.names
    ```
    
2. Read in dyadic cluster membership output from dynCluster (e.g., `DYADS_cluster_T0_Z0.txt`, ...).

    ```R
    # compile output into dataframe format
    dyad.df <- map_df(periods, function(x) {
      # index
      element <- x + 1
    
      # extract specific period index
      cat(paste("Extracting data for Year ", years[[element]], "\n", sep = ""))
      
      # read into df dyads for each cluster 
      out <- map_df(cl.names, function(y) {
        file <- paste(SIM_FOLDER, "/DYADS_cluster_T", x, "_Z", y, ".txt", sep = "")
    
        tryCatch({
         #if(file.exists(file)) {
            data <- read.delim(file, header = TRUE, check.names = FALSE)
            dyads <- data[, 1:2]
            trade <- data[,-c(1:2)]
            dyads <- dyads %>%
              mutate(cluster = y + 1)
            dyads <- cbind(dyads, trade)
          #  }
        }, error = function(e) NULL) # if no dyads are assigned to a cluster-year, no corresponding txt file will be created, assign NULL and skip
      })
      
      # add variable for year
      out <- out %>%
        mutate(year = years[[element]])
    })
    
    # create/subset vars
    dyad.df.sub <- dyad.df %>%
      select(imp, exp, year, cluster) %>%
      mutate(dyad = ifelse(as.numeric(imp) < as.numeric(exp), 
                           paste(imp, exp, sep = "_"), 
                           paste(exp, imp, sep = "_"))) %>%
      arrange(imp, exp, year)
    
    # dataframe contains:
    # 450 rows (45 dyads over 10 years)
    # 5 columns (importer, exporter, year, cluster, dyad id)
    # view first 10 rows of the data
    head(dyad.df.sub, 10)
    ##  imp exp year cluster dyad
    ## 1   1   2 1963       1  1_2
    ## 2   1   2 1967       1  1_2
    ## 3   1   2 1970       1  1_2
    ## 4   1   3 1961       2  1_3
    ## 5   1   3 1962       2  1_3
    ## 6   1   3 1966       2  1_3
    ## 7   1   3 1968       2  1_3
    ## 8   1   3 1969       2  1_3
    ## 9   1   3 1970       2  1_3
    ## 10  1   4 1961       2  1_4
    ```

3. Read in raw input trade data (`dynamic_1962.csv`, ...)
    ```R
    # load raw data for total trade value by dyad-cluster-year-sitc4d
    dyad.trade.sitc.df <- map_df(years, function(x) {
      
      # read in raw trade data by periods
      cat(paste("Extracting data for Year ", x, "\n", sep = ""))
      
      # set file directory
      filename <- paste(SIM_FOLDER, "/dynamic_", x, ".csv", sep="")
      
      # sample data
      sampleData <- read.table(filename, header = TRUE, 
                               check.names = FALSE, nrow = 1)
      ncols <- ncol(sampleData)
      
      # read in all data
      D <- read.table(filename,
                      header = TRUE,
                      check.names = FALSE,
                      colClasses = c(rep("character", 3),
                                     rep("numeric", (ncols - 3))))
    })
    
    # convert to fake years
    dyad.trade.sitc.df <- dyad.trade.sitc.df %>%
      mutate(year = as.numeric(year) + 1960)
    
    # dataframe contains:
    # 900 rows (90 directed-dyads over 10 years)
    # 43 columns (year, importer, exporter, 40 products)
    # view first 5 rows and 10 columns of the data
    dyad.trade.sitc.df[1:5, 1:10]
    ##   year importer_ISO exporter_ISO SITC0_1 SITC0_2 SITC0_3  SITC0_4 SITC0_5      SITC0_6   SITC0_7
    ## 1 1961            1            2       0       0       0      0.0       0     177.0528       0.0
    ## 2 1961            1            3       0       0       0 664344.6       0  536712.5031  133149.7
    ## 3 1961            1            4       0       0       0      0.0       0  536385.8453       0.0
    ## 4 1961            1            5  372390       0       0      0.0       0  536450.8843       0.0
    ## 5 1961            1            6 3171746 2797487 4872051 981809.8 2497946 1412988.9059 1075698.3
    ```
    
4. Create a cluster proportion file using dyadic cluster membership **estimated** by dynCluster. That is, a file that contains the typical composition of product trade for each cluster. 
    ```R
    # calculate total trade volume for each product disaggregated by year
    year.total <- dyad.trade.sitc.df %>%
      select(-importer_ISO, -exporter_ISO) %>%
      #filter(!(SITC0_0011 == -30)) %>%
      group_by(year) %>%
      summarise_all(funs(sumPlusOne)) %>% # plus one in case of zero
      ungroup() 
    
    # calculate trade weighted by yearly total
    weight.trade.year <- map_df(years, function(x){
      
      cat(paste("Calculating for Year ", x, "\n", sep = ""))
      
      # subset total by year
      year.total.sub <- year.total %>%
        filter(year == x) %>%
        select(-year)
      
      # subset raw trade by year
      raw.df.sub <- dyad.trade.sitc.df %>%
        select(-importer_ISO, -exporter_ISO) %>%
        filter(year == x) %>%
        select(-year)
      
      # set negative values to 0
      raw.df.sub[raw.df.sub < 0] <- 0  
      
      # divide raw directed-dyad trade volume by vector of yearly product trade total
      out <- sweep(raw.df.sub, MARGIN = 2, as.matrix(year.total.sub), FUN = "/")
      
      # extract id 
      id <- dyad.trade.sitc.df %>%
        select(year, importer_ISO, exporter_ISO) %>%
        filter(year == x)
      
      # combine with id and year
      out <- cbind(out, id)
      
      return(out)
    })
    
    # coerce to numeric to be consistent
    weight.trade.year <- weight.trade.year %>%
      mutate(year = as.numeric(year),
             importer_ISO = as.character(importer_ISO),
             exporter_ISO = as.character(exporter_ISO))
    
    # estimated proportions
    # extract dyad-year cluster info from estimated
    dyad.cluster.year <- dyad.df %>%
      select(imp, exp, year, cluster) %>%
      mutate(imp = as.character(imp),
             exp = as.character(exp))
    
    # merge trade proportions with cluster info
    weight.trade.year.impexp <- left_join(dyad.cluster.year, weight.trade.year,
                                          by = c("imp" = "importer_ISO",
                                                 "exp" = "exporter_ISO",
                                                 "year" = "year"))
    weight.trade.year.expimp <- left_join(dyad.cluster.year, weight.trade.year,
                                          by = c("imp" = "exporter_ISO",
                                                 "exp" = "importer_ISO",
                                                 "year" = "year"))
    
    # extract preferred ordering of products
    prod.order <- names(weight.trade.year.impexp)[-c(1:4)]
    
    # calculate mean cluster-year proportions: direction 1
    weight.trade.year.impexp.sum <- weight.trade.year.impexp %>%
      select(-imp, -exp) %>%
      arrange(year, cluster) %>%
      group_by(year, cluster) %>%
      summarise_all(funs(mean)) %>%
      ungroup()
    
    # calculate mean cluster proportions: direction 1
    weight.trade.impexp.sum <- weight.trade.year.impexp.sum %>%
      arrange(cluster, year) %>%
      group_by(cluster) %>%
      summarise_all(funs(mean)) %>%
      ungroup() %>%
      select(-year) %>%
      gather(prod, value, -cluster) %>%
      arrange(cluster) %>%
      group_by(cluster) %>%
      mutate(cluster_total = sum(value)) %>%
      ungroup() %>%
      mutate(weighted_value = value/cluster_total) %>%
      select(-value, -cluster_total) %>%
      mutate(prod = factor(prod, levels = prod.order)) %>%
      spread(prod, weighted_value)
    
    # calculate mean cluster-year proportions: direction 2
    weight.trade.year.expimp.sum <- weight.trade.year.expimp %>%
      select(-imp, -exp) %>%
      arrange(year, cluster) %>%
      group_by(year, cluster) %>%
      summarise_all(funs(mean)) %>%
      ungroup()
    
    # calculate mean cluster proportions: direction 2
    weight.trade.expimp.sum <- weight.trade.year.expimp.sum %>%
      arrange(cluster, year) %>%
      group_by(cluster) %>%
      summarise_all(funs(mean)) %>%
      ungroup() %>%
      select(-year) %>%
      gather(prod, value, -cluster) %>%
      arrange(cluster) %>%
      group_by(cluster) %>%
      mutate(cluster_total = sum(value)) %>%
      ungroup() %>%
      mutate(weighted_value = value/cluster_total) %>%
      select(-value, -cluster_total) %>%
      mutate(prod = factor(prod, levels = prod.order)) %>%
      spread(prod, weighted_value)
    
    # combine both directions
    newtf <- cbind(weight.trade.impexp.sum, 
                   weight.trade.expimp.sum %>%
                     select(contains("SITC0")))
    
    newtf$cluster <- NULL
    
    # dataframe contains: 
    # 3 rows (one for each cluster)
    # 80 total columns (40 products*2 directions)
    # note that proportions for each row should sum up to 2 
    # view data first 5 columns of the data
    newtf[,1:5]
    ##      SITC0_1     SITC0_2      SITC0_3    SITC0_4    SITC0_5
    ## 1 0.02104624 0.008787409 0.0007205318 0.00000000 0.09258156
    ## 2 0.02053428 0.039401288 0.0239437045 0.07493571 0.01843604
    ## 3 0.02557064 0.021981360 0.0263553161 0.01657127 0.02597570
    
    # save 
    write.csv(newtf, file = paste(SIM_FOLDER, "/NewTF.csv", sep = ""), 
              row.names = FALSE)
    ```    
    
5. Create a cluster proportion file using the **true** cluster membership data.
    ```R
    # read in truth file
    truth <- read.csv(TRUTH_FILE)
    truth.l <- truth %>%
      gather(year, cluster, -cty1, -cty2, -dyad) %>%
      mutate(year = as.numeric(str_replace_all(year, "z", "")) + 1960) %>%
      arrange(cty1, cty2, year)
      
    # extract dyad-year cluster info from truth
    dyad.cluster.year.truth <- truth.l %>%
      distinct(dyad, year, .keep_all = TRUE) %>%
      select(cty1, cty2, year, cluster) %>%
      rename(imp = cty1, 
             exp = cty2) %>%
      mutate(imp = as.character(imp),
             exp = as.character(exp))
    
    # merge trade proportions with cluster info
    weight.trade.year.impexp.truth <- left_join(dyad.cluster.year.truth, 
                                                weight.trade.year,
                                                by = c("imp" = "importer_ISO",
                                                       "exp" = "exporter_ISO",
                                                       "year" = "year"))
    weight.trade.year.expimp.truth <- left_join(dyad.cluster.year.truth, 
                                                weight.trade.year,
                                                by = c("imp" = "exporter_ISO",
                                                       "exp" = "importer_ISO",
                                                       "year" = "year"))
    
    # extract preferred ordering of products
    prod.order <- names(weight.trade.year.impexp.truth)[-c(1:4)]
    
    # calculate mean cluster-year proportions: direction 1
    weight.trade.year.impexp.sum.truth <- weight.trade.year.impexp.truth %>%
      select(-imp, -exp) %>%
      arrange(year, cluster) %>%
      group_by(year, cluster) %>%
      summarise_all(funs(mean)) %>%
      ungroup()
    
    # calculate mean cluster proportions: direction 1
    weight.trade.impexp.sum.truth <- weight.trade.year.impexp.sum.truth %>%
      arrange(cluster, year) %>%
      group_by(cluster) %>%
      summarise_all(funs(mean)) %>%
      ungroup() %>%
      select(-year) %>%
      gather(prod, value, -cluster) %>%
      arrange(cluster) %>%
      group_by(cluster) %>%
      mutate(cluster_total = sum(value)) %>%
      ungroup() %>%
      mutate(weighted_value = value/cluster_total) %>%
      select(-value, -cluster_total) %>%
      mutate(prod = factor(prod, levels = prod.order)) %>%
      spread(prod, weighted_value)
    
    # calculate mean cluster-year proportions: direction 2
    weight.trade.year.expimp.sum.truth <- weight.trade.year.expimp.truth %>%
      select(-imp, -exp) %>%
      arrange(year, cluster) %>%
      group_by(year, cluster) %>%
      summarise_all(funs(mean)) %>%
      ungroup()
    
    # calculate mean cluster proportions: direction 2
    weight.trade.expimp.sum.truth <- weight.trade.year.expimp.sum.truth %>%
      arrange(cluster, year) %>%
      group_by(cluster) %>%
      summarise_all(funs(mean)) %>%
      ungroup() %>%
      select(-year) %>%
      gather(prod, value, -cluster) %>%
      arrange(cluster) %>%
      group_by(cluster) %>%
      mutate(cluster_total = sum(value)) %>%
      ungroup() %>%
      mutate(weighted_value = value/cluster_total) %>%
      select(-value, -cluster_total) %>%
      mutate(prod = factor(prod, levels = prod.order)) %>%
      spread(prod, weighted_value)
    
    # combine both directions
    tf.truth <- cbind(weight.trade.impexp.sum.truth, 
                      weight.trade.expimp.sum.truth %>%
                        select(contains("SITC0")))
    
    tf.truth$cluster <- NULL
    
    # dataframe contains: 
    # 3 rows (one for each cluster)
    # 80 total columns (40 products*2 directions)
    # note that proportions for each row should sum up to 2 
    # view data first 5 columns of the data
    tf.truth[,1:5]
    ##      SITC0_1    SITC0_2      SITC0_3    SITC0_4    SITC0_5
    ## 1 0.01255561 0.01115383 0.0008767924 0.01262706 0.06585254
    ## 2 0.02036431 0.05316646 0.0229124523 0.10074567 0.01593227
    ## 3 0.02507445 0.02339466 0.0254539034 0.01728950 0.02655064

    # save 
    write.csv(tf.truth, file = paste(SIM_FOLDER, "TF-truth.csv", sep = ""), 
                row.names = FALSE)
    ```

6. Plot the **estimated** cluster proportion data as a heatmap.
    ```R
    # load NewTF file
    newtf <- read.csv(paste(SIM_FOLDER, "/NewTF.csv", sep = ""), 
                      as.is = TRUE) # file contains mean proportions
    
    # set ordering of clusters and label location
    order.est <- c(2, 1, 4, 3, 6, 5)
    legend.label.high.y <- 38
    legend.label.low.y <- 12
    
    # set industry label cut point
    indus.vec <- seq(1, 40, 10)
    
    # set y-axis product description
    descr <- paste("Product", indus.vec, sep = "-")
    
    # set filename for output
    filename.est.out <- paste(SIM_FOLDER, "/TF_heatmap_demeaned_est.pdf", sep = "")
    
    # plot
    product_heatmap(data = newtf,
                    cluster.order = order.est,
                    filename = filename.est.out,
                    product.description = descr,
                    product.txt.cex = 1.2,
                    sitc.1d.dash.shift = -0.5,
                    industry.label.loc = indus.vec,
                    cluster.label.hjust = 0.9,
                    legend.cex = 1.75,
                    legend.label.high.txt = "Higher Trade\nProportion",
                    legend.label.low.txt = "Lower Trade\nProportion",
                    pdf.width = 10,
                    pdf.height = 10,
                    legend.label.high.y = legend.label.high.y,
                    legend.label.low.y = legend.label.low.y)
    ```

7. Plot the **true** cluster proportion data as a heatmap.
    ```R
    # load tf file
    tf.truth <- read.csv(paste(SIM_FOLDER, "TF-truth.csv", sep = ""), 
                         as.is = TRUE) # file contains mean proportions
    
    # set ordering of clusters
    order.truth <- seq(1, 6, 1)
    
    # set label location
    legend.label.high.y <- 38
    legend.label.low.y <- 12
      
    # set filename for output
    filename.truth.out <- paste(SIM_FOLDER, "/TF_heatmap_demeaned_truth.pdf", sep = "")
    
    # plot
    product_heatmap(data = tf.truth,
                    cluster.order = order.truth,
                    filename = filename.truth.out,
                    product.description = descr,
                    product.txt.cex = 1.2,
                    sitc.1d.dash.shift = -0.5,
                    industry.label.loc = indus.vec,
                    cluster.label.hjust = 0.9,
                    legend.cex = 1.75,
                    legend.label.high.txt = "Higher Trade\nProportion",
                    legend.label.low.txt = "Lower Trade\nProportion",
                    pdf.width = 10,
                    pdf.height = 10,
                    legend.label.high.y = legend.label.high.y,
                    legend.label.low.y = legend.label.low.y)
    ```
    
8. A side-by-side comparison of the two heatmaps show that the composition of product trade is very similar. This suggests that dynCluster did well in recovering the original clusters.

    True Product Proportion             |  Estimated Product Proportion
    :-------------------------:|:-------------------------:
    ![](images/TF_heatmap_demeaned_truth.png)  |  ![Estimated](images/TF_heatmap_demeaned_est.png)
  
  
9. Overall, dynCluster correctly recovered **98.4%** of the true dyadic cluster membership.
    ```R
    # merge
    compare.df.m <- left_join(truth.l, dyad.df.sub, 
                              by = c("dyad", "year"))
    
    # clean
    compare.df <- compare.df.m %>%
      rename(cluster_truth = cluster.x,
             cluster_out = cluster.y) %>%
      select(cty1, cty2, dyad, imp, exp, year, cluster_truth, cluster_out) 
    
    # summarize
    compare.df <- compare.df %>%
      mutate(
        n_correct = if_else(cluster_truth == cluster_out, 1, 0),
        per_correct = mean(n_correct),
        c_1_1 = if_else(cluster_truth == 1 & cluster_out == 1, 1, 0),
        c_1_2 = if_else(cluster_truth == 1 & cluster_out == 2, 1, 0),
        c_1_3 = if_else(cluster_truth == 1 & cluster_out == 3, 1, 0),
        c_2_1 = if_else(cluster_truth == 2 & cluster_out == 1, 1, 0),
        c_2_2 = if_else(cluster_truth == 2 & cluster_out == 2, 1, 0),
        c_2_3 = if_else(cluster_truth == 2 & cluster_out == 3, 1, 0),
        c_3_1 = if_else(cluster_truth == 3 & cluster_out == 1, 1, 0),
        c_3_2 = if_else(cluster_truth == 3 & cluster_out == 2, 1, 0),
        c_3_3 = if_else(cluster_truth == 3 & cluster_out == 3, 1, 0))
    
    # check overall correct classification
    compare.df %>% 
      pull(per_correct) %>%
      unique()
      
    ## [1] 0.9844444
    ```
    + Note that, in this example, Cluster 1 (truth) happens to correspond to Cluster 1 (estimated), Cluster 2 (truth) to Cluster 2 (estimated), and Cluster 3 (truth) to Cluster 3 (estimated). 
    + This may not always be the case depending on how dynCluster generates the cluster labels. Thus, before making comparisons, users should always visually inspect the correspondence between the "truth" and the output and reorder cluster labels when necessary.
