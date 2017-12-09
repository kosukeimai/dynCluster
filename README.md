# dynCluster: dynamic Clustering algorithm [![Build Status](https://travis-ci.org/kosukeimai/dynCluster.svg?branch=master)](https://travis-ci.org/kosukeimai/dynCluster)

This tutorial shows how to install and run dynCluster on [Amazon Linux AMI](http://aws.amazon.com/amazon-linux-ami/). For more details of the algorithm, see [Imai, Kim, Liao (2017)](https://www.stevenliao.org/uploads/2/5/6/9/25699716/bigtrade.pdf).

Dependency of this package:

    * boost/1.55.0
    * openmpi
    * Rcpp

## Installing dynCluster on Amazon Linux AMI

1. Create an [Amazon AWS account](https://aws.amazon.com/)

2. Log in AMI from local terminal with assigned key
```sh
chmod 400 ~/Dropbox/aws/key.pem
ssh -i ~/Dropbox/aws/key.pem ec2-user@ec2...us-west-1.compute.amazonaws.com
```

3. Update and install dependencies
```sh
sudo yum update
sudo yum -y install make gcc-c++ openmpi-2.1.1 R-3.4.1
```

4. Download [dynCluster](https://github.com/HJ08003/dynCluster/archive/master.zip) package and [boost/1.55.0](http://sourceforge.net/projects/boost/files/boost/1.55.0/boost_1_55_0.tar.bz2)

5. Transfer downloaded files from local to AMI
```sh
scp -i ~/Dropbox/aws/key.pem ~/Dropbox/aws/dynCluster-master.zip ec2-user@ec2...us-west-1.compute.amazonaws.com:~
scp -i ~/Dropbox/aws/key.pem ~/Dropbox/aws/boost_1_55_0.tar.bz2  ec2-user@ec2...us-west-1.compute.amazonaws.com:~
```

6. Unzip dynCluster
```sh
unzip dynCluster-master.zip
rm dynCluster-master.zip
```

7. Un-tar and build Boost
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

8. Set up support for libraries installed in /usr/local/lib
```sh
sudo bash -c "echo '/usr/local/lib' > /etc/ld.so.conf.d/boost.conf"
sudo ldconfig
```

9. Install R packages
```sh
R
```    
```R    
dir.create(Sys.getenv("R_LIBS_USER"), showWarnings = FALSE, recursive = TRUE)
install.packages(c("Rcpp", "knitr"), Sys.getenv("R_LIBS_USER"), repos = "http://cran.case.edu" )
q()
```

10. Build dynCluster package
```sh
cd ~/dynCluster-master/
R CMD build .
R CMD check .
R CMD check --as-cran .
R CMD INSTALL .
```

## Running dynCluster on Amazon Linux AMI

This demo runs a toy example created based on [Imai, Kim, Liao (2017)](https://www.stevenliao.org/uploads/2/5/6/9/25699716/bigtrade.pdf). The smaller dataset in this example contains 4774 dyads trading 625 products in 5 separate years (1962, 1972, 1982, 1992, 2002).

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
    * `mainZTM` calls `mainRcpp` from the directory that contains the config.txt
```R
mainZTM("./example", comeBack=TRUE)
```
    * `testExample` is a test function that runs the toy example
```R
ptm <- proc.time() # start the clock
testExample("./example", nThreads = 4, comeBack = TRUE)
proc.time() - ptm # Stop the clock (117.269 seconds)
```

3. Download results from dynCluster
```sh
scp -i ~/Dropbox/aws/key.pem -r ec2-user@ec2...us-west-1.compute.amazonaws.com:~/dynCluster-master/example ~/Dropbox/aws/out-raw/
```


