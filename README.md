# dynCluster: Dynamic Clustering Algorithm [![Build Status](https://travis-ci.org/kosukeimai/dynCluster.svg?branch=master)](https://travis-ci.org/kosukeimai/dynCluster)

Authors: [In Song Kim](http://web.mit.edu/insong/www/index.html), [Steven Liao](https://www.stevenliao.org/), [Kosuke Imai](http://imai.princeton.edu/)

For more details of the method and applications, see our paper:

+ [Measuring Trade Profile with Granular Product-level Trade Data](https://www.stevenliao.org/uploads/2/5/6/9/25699716/bigtrade.pdf)

## dynCluster on Amazon Web Services: A Toy Example
1. Install dynCluster on [Amazon Web Services (AWS)](https://aws.amazon.com/). See our [Wiki page](https://github.com/kosukeimai/dynCluster/wiki/How-to-install-dynCluster-on-AWS) for step-by-step instructions.

2. Once dynCluster is installed, we created a small simulated dataset following the data generating process described in our paper. For detailed code, see our [Wiki page](https://github.com/kosukeimai/dynCluster/wiki/How-to-run-dynCluster-on-AWS).

+ The data covers `10` countries (`90` directed-dyads) trading `40` products over `10` time periods.
```
##   year importer_ISO exporter_ISO SITC0_1 SITC0_2 SITC0_3  SITC0_4 SITC0_5      SITC0_6   SITC0_7 ...
## 1    1            1            2       0       0       0      0.0       0     177.0528       0.0
## 2    1            1            3       0       0       0 664344.6       0  536712.5031  133149.7
## 3    1            1            4       0       0       0      0.0       0  536385.8453       0.0
## 4    1            1            5  372390       0       0      0.0       0  536450.8843       0.0
## 5    1            1            6 3171746 2797487 4872051 981809.8 2497946 1412988.9059 1075698.3
## ...
```

+ Each dyad belongs to 1 of `3` different clusters in a given time period.
```
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
## ...
```
  
+ The ultimate goal is to see how well dynCluster can recover the three true clusters as well as the dyadic cluster membership.

3. Adjust the simulation parameters in `config.txt` and copy the file to the same folder containing the simulated data.

    ![](https://github.com/kosukeimai/dynCluster/blob/master/images/config.png)

  + `NSIM`: the number of simulation. Each simulation will use different starting values. The final output will select the simulation that yields the highest log-likelihood.
  + `MAXITER`: the max number of iterations
  + `EPS`: the epsilon for convergence
  + `Z`: the number of clusters
  + `NU`: the scaling parameter for zero-trade probability
  + `OP_w_ij`: the mixture responsibilities
  + `file`: the names of the input files excluding "dynamic_", e.g., 1961, 1962, ..., 1970
  + `threads`: the number of threads to run parallel (OpenMP)
  + `MCITERATIONS`: the number of MC iterations
  + `QFILE`: output file for the Bernoulli probabilities for zero trade
  + `MUFILE`: output file for the mean of the normal distribution
  + `SIGMAQFILE`: output file for the variance of the normal distribution
  + `PIEFILE`: output file for mixture probabilities
  + `ZETAFILE`: output file for mixture responsibilities
  + `TFFILE`: output file for cluster trade proportions (non-weighted by time)
  + `PMATFILE`: output file for cluster transition probabilities

4. In R, run the function `mainZTM`. This function wraps and calls C++ functions (e.g., `mainRcpp`) from dynCluster.
    ```R
    # load library
    library(dynCluster)
        
    # run and time dynCluster
    ptm <- proc.time() # start the clock
    mainZTM("./sim-25", comeBack=TRUE)
    proc.time() - ptm # stop the clock
    
    ##   user  system elapsed 
    ## 16.976   0.051  17.100 
    ```
+ Note that this toy example runs on **t2.micro** instances on AWS, which is available as a [free tier](https://aws.amazon.com/free/).

5. To assess the performace of dynCluster, we created product-trade heatmaps based on the "true" cluster membership data above and the estimated dyadic cluster membership from dynCluster. 

+ A side-by-side comparison of the two heatmaps show that the composition of product trade is very similar. This suggests that dynCluster did well in recovering the original clusters.

    True Product Proportion             |  Estimated Product Proportion
    :-------------------------:|:-------------------------:
    ![](images/TF_heatmap_demeaned_truth.png)  |  ![Estimated](images/TF_heatmap_demeaned_est.png)
  
+ Overall, dynCluster correctly recovered **98.4%** of the true dyadic cluster membership.
