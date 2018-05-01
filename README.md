# dynCluster: Dynamic Clustering Algorithm [![Build Status](https://travis-ci.org/kosukeimai/dynCluster.svg?branch=master)](https://travis-ci.org/kosukeimai/dynCluster)

This tutorial demonstrates how to install and run dynCluster on [Amazon Linux AMI](http://aws.amazon.com/amazon-linux-ami/). For more details of the algorithm, see [Kim, Liao, Imai (2018)](https://www.stevenliao.org/uploads/2/5/6/9/25699716/bigtrade.pdf).

Dependencies of this package:

    * boost/1.55.0
    * openmpi
    * Rcpp

## Installing dynCluster on Amazon Linux AMI

1. Create an [Amazon AWS account](https://aws.amazon.com/)

2. Launch an Amazon Machine Image (AMI) Instance

 + Choose instance type *Amazon Linux AMI*
 + Choose instance with sufficient CPU, memory, storage, etc., depending on the size of your dataset

3. Log in AMI from local terminal with assigned key (e.g., `~/PATH/KEY.pem`) and the instance's public DNS (`ec2-user@ec2...compute.amazonaws.com`)
```sh
chmod 400 ~/Dropbox/aws/KEY.pem
ssh -i ~/Dropbox/aws/KEY.pem ec2-user@ec2...compute.amazonaws.com
```

4. Update and install dependencies
```sh
sudo yum update
sudo yum -y install make gcc-c++ openmpi-2.1.1 R-3.4.1
```

5. Download the [dynCluster](https://github.com/kosukeimai/dynCluster/archive/master.zip) package and [boost/1.55.0](http://sourceforge.net/projects/boost/files/boost/1.55.0/boost_1_55_0.tar.bz2)

6. Transfer downloaded files from local to AMI
```sh
scp -i ~/Dropbox/aws/KEY.pem ~/Dropbox/aws/dynCluster-master.zip ec2-user@ec2...compute.amazonaws.com:~
scp -i ~/Dropbox/aws/KEY.pem ~/Dropbox/aws/boost_1_55_0.tar.bz2  ec2-user@ec2...compute.amazonaws.com:~
```

7. Unzip dynCluster
```sh
unzip dynCluster-master.zip
rm dynCluster-master.zip
```

8. Un-tar and build Boost
```sh
JOBS=`grep -c ^processor /proc/cpuinfo`
export BOOST_VERSION="1_55_0"
tar xf boost_${BOOST_VERSION}.tar.bz2
cd boost_${BOOST_VERSION}
./bootstrap.sh
./b2 -d1 -j${JOBS} --with-thread --with-filesystem --with-python --with-regex -sHAVE_ICU=1 --with-program_options --with-system link=shared release toolset=gcc stage
sudo ./b2 -j${JOBS} --with-thread --with-filesystem --with-python --with-regex -sHAVE_ICU=1 --with-program_options --with-system toolset=gcc link=shared release install
cd ../
```

9. Set up support for libraries installed in /usr/local/lib
```sh
sudo bash -c "echo '/usr/local/lib' > /etc/ld.so.conf.d/boost.conf"
sudo ldconfig
```

10. Install R packages
```sh
R
```    
```R    
dir.create(Sys.getenv("R_LIBS_USER"), showWarnings = FALSE, recursive = TRUE)
install.packages(c("Rcpp", "knitr"), Sys.getenv("R_LIBS_USER"), repos = "http://cran.case.edu" )
```
```sh
q()
```

11. Build dynCluster package
```sh
cd ~/dynCluster-master/
R CMD build .
R CMD check .
R CMD check --as-cran .
R CMD INSTALL .
```

## Running dynCluster on Amazon Linux AMI

This demo runs a toy example created based on [Kim, Liao, Imai (2018)](https://www.stevenliao.org/uploads/2/5/6/9/25699716/bigtrade.pdf). The smaller dataset in this example contains 4774 dyads trading 625 products in 5 separate years (1962, 1972, 1982, 1992, 2002).

1. Go to the directory where dynCluster was installed
```sh
cd ~/dynCluster-master/
```

2. Open R and load dynCluster
```sh
R
```
```R
library(dynCluster)
```

3. There are two R functions that wrap and call C++ functions:
 + `mainZTM` calls `mainRcpp` from the directory that contains the `config.txt`
```R
mainZTM("./example", comeBack=TRUE)
```
 + `testExample` is a test function that runs the toy example
```R
ptm <- proc.time() # start the clock
testExample("./example", nThreads = 4, comeBack = TRUE)
proc.time() - ptm # Stop the clock (117.269 seconds)
```
 + Adjust `config.txt`as needed
    - Set the number of clusters (*Z*)
    - Set names of input files (*file*)
    - Increase number of CPU cores (*threads*)
    - Increase max number of iterations (*MAXITER*)

4. Download results from dynCluster for further analyses
```sh
scp -i ~/Dropbox/aws/KEY.pem -r ec2-user@ec2...compute.amazonaws.com:~/dynCluster-master/example ~/Dropbox/aws/output/
```


