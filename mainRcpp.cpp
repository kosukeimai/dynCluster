/*
The MIT License

Copyright (c) 2013 Radhika S. Saksena radhika.saksena@gmail.com,
                   Kosuke Imai kimai@princeton.edu

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include "Analytics.h"
#include "ZeroTradeModelIO.h"
#include <numeric>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <boost/random.hpp>
#include <time.h>

using namespace std;
using namespace boost;


bool sortFunc(pair<int, long double> p1, pair<int, long double> p2);

// [[Rcpp::export]]                                                                                                                                                                         
int ztm(int nsim, int maxiter, long double eps, int numModelsM, long double nu, long double OP_w_ij_inp, 
         string inputYears, int numThreads, int total_mc_trials, string q_filename, string mu_filename,
	string sigmasq_filename, string pi_filename, string zeta_filename, int seedOffset) {

  // int numCountries, numProductsK, numModelsM;                                                                                                                                           
  // long double nu;
  // long double OP_w_ij_inp;
  //long double nu, bayes_beta, bayes_alpha, bayes_tau;                                                                                                                                   
  // string tmpString;
  vector<int> nemptyIdx;
  // int success=0;
  // int numThreads=8;
  // int numModelsM;
  // int  maxiter;
  // int nsim, total_mc_trials 
  int mc_trials;
  // bool successFlag = 0;
  bool resetFlag = 0;
  //string inputFilename;                                                                                                                                                                 
  // string inputYears;
  ofstream logFile;
  // ifstream configFile;
  long double max_ll_new;
  int max_sim;
  // long double eps;
  bool checkpointFlag = 0;
  vector<pair<int, long double> > simDF;
  bool* isSuccess;
  // string q_filename,mu_filename,sigmasq_filename,pi_filename,zeta_filename;
  // int np 
  int rank=0;
  long double** all_TFProportions;
  // double start_time, end_time;

     cout<<"read in: nsim = "<<nsim<<endl;
     cout<<"read in: maxiter = "<<maxiter<<endl;
     cout<<"read in: eps = "<<eps<<endl; 
     cout<<"read in: Z = "<<numModelsM<<endl;
     cout<<"read in: nu = "<<nu<<endl;
     cout<<"read in: OP_w_ij_inp = "<<OP_w_ij_inp<<endl;
     cout<<"read in: q checkpoint file = "<<q_filename<<endl;
     cout<<"read in: mu checkpoint file = "<<mu_filename<<endl;
     cout<<"read in: sigmasq checkpoint file = "<<sigmasq_filename<<endl;
     cout<<"read in: pi checkpoint file = "<<pi_filename<<endl;
     cout<<"read in: zeta checkpoint file = "<<zeta_filename<<endl;
     //cout<<"read in: input file = "<<inputFilename<<endl;
     cout<<"read in: input years = "<<inputYears<<endl;
     cout<<"read in: number of threads = "<<numThreads<<endl;

     cout<<"input: random seed = "<<seedOffset<<endl;

     isSuccess = (bool *)malloc(sizeof(bool)*nsim); 

     omp_set_num_threads(numThreads);
     cout<<"Before call to ZeroTradeModel constructor."<<endl;
	 //ZeroTradeModelIO trade1(numModelsM,nu,OP_w_ij_inp,inputFilename,numThreads); 
	 ZeroTradeModelIO trade1(numModelsM,nu,OP_w_ij_inp,inputYears,maxiter,numThreads); 

     cout<<"After call to ZeroTradeModel constructor."<<endl;
     logFile.open("ZeroTradeModel.log");

#if EXITONZEROVARIANCE
     cout<<"Algorithm exits when it encounters zero variance."<<endl;
     logFile<<"Algorithm exits when it encounters zero variance."<<endl;
#endif
#if FLIPDYADPAIRS
     cout<<"Dyad pairs may be flipped."<<endl;
     logFile<<"Dyad pairs may be flipped."<<endl;
#endif
#if PRIOR_OP
     cout<<"Prior is applied."<<endl;
     logFile<<"Prior is applied."<<endl;
     logFile<<"Executing: ZeroTradeModel with the following parameters: ";
#endif
#ifndef PRIOR_OP
     logFile<<"Prior is not applied."<<endl;
     logFile<<"Executing: ZeroTradeModel with the following parameters: ";
#endif
     logFile<<"read in: maxiter = "<<maxiter<<endl;
     logFile<<"read in: eps = "<<eps<<endl;
     logFile<<"read in: Z = "<<numModelsM<<endl;
     logFile<<"read in: nu = "<<nu<<endl;
     //logFile<<"read in: input file = "<<inputFilename<<endl;
     logFile<<"read in: input years = "<<inputYears<<endl;
     logFile<<"read in: number of threads = "<<numThreads<<endl;


     max_ll_new = -9999999999;
     max_sim = 0;

     for(int isim = 0; isim  < nsim; isim++){

         //int numThreads;
	 long double ll_new = 0.0;
         // long double ll_old = -9999999999999999;
         // double elapsed_time;
         time_t seed;
         // ofstream* clusterFiles = new ofstream[trade1.numModelsM];
         checkpointFlag = 0;

         cout<<endl<<endl<<"Running Simulation "<<isim<<"...."<<endl<<endl;;
         logFile<<endl<<endl<<"Running Simulation "<<isim<<"...."<<endl<<endl;;

         //seed = time(0);
         seed = isim;
         seed = seed + seedOffset;
         cout<<"Using seed = "<<seed<<" seedOffset = "<<seedOffset<<endl;
         //seed = 1377336905;
         //max for Y90

         //cout<<"Modified."<<endl;
         //logFile<<"Modified."<<endl;

         logFile<<"Using seed = "<<seed<<endl;

         //TODO: uncomment
        trade1.initializeModel(seed);
         //0-th iteration was done above
         //trade1.calcMu();

         for(int iter=1; iter < maxiter; iter++){
             cout<<"Running iteration: "<<iter<<endl;
             trade1.expectationTh(iter,logFile);
             cout<<"iter = "<<iter<<": After expectation."<<endl;
             trade1.maximizationTh(0,resetFlag,logFile);
             cout<<"iter = "<<iter<<": After maximization."<<endl;
             //long double ll_new = trade1.logLikelihoodTh(iter);
             //double start_time = omp_get_wtime();
             //long double ll_new = trade1.obsLogPosterior(iter);
             //double end_time = omp_get_wtime();
             cout<<"After ll = "<<ll_new<<endl;
         }//end of iter loop

        trade1.chkpClusters(isim,logFile);
        trade1.chkpModelParams(isim,logFile);

     }//end of isim loop

     omp_set_num_threads(numThreads);

     //mc_trials = (total_mc_trials/np)*np;
     mc_trials = total_mc_trials;
     cout<<"Before call to Analytics constructor."<<endl;
     Analytics mc1(mc_trials,numThreads,numModelsM,trade1.ntimes,trade1.numProdsK,q_filename,
                     mu_filename,sigmasq_filename,pi_filename,zeta_filename);
     cout<<"After call to Analytics constructor."<<endl;

     mc1.analytics(rank);

#pragma omp parallel for
     for(int ith = 0; ith < numThreads; ith++){
        for(int z = 0; z < numModelsM; z++){
            int tmpSum = 0;
            for(int t = 0; t < mc1.ntimes; t++)
            {
                tmpSum += mc1.clusterOcc[t][z].size();
            }
            if(!tmpSum) continue;
            mc1.doMC(z,rank);
        }
    }

     all_TFProportions = (long double**)malloc(numModelsM*sizeof(long double*));
     for(int z = 0; z < numModelsM; z++){
          all_TFProportions[z] = (long double*)malloc(2*trade1.numProdsK*sizeof(long double));
          for(int k = 0; k < 2*trade1.numProdsK; k++){
              all_TFProportions[z][k] = 0.0;
          }
     }

     for(int z = 0; z < numModelsM; z++){
         for(int k = 0; k < 2*trade1.numProdsK; k++){
             for(int ith = 0; ith < numThreads; ith++){
                all_TFProportions[z][k] += mc1.prodProportion[ith][z][k];
             }
         }
     }

     for(int z = 0; z < numModelsM; z++){
        for(int k = 0; k < 2*trade1.numProdsK; k++){
            all_TFProportions[z][k] /= (long double)numThreads;
        }
     }

    ostringstream outfileName;
    outfileName<<"TF_"<<rand()<<".txt";
    ofstream fout;
    fout.open(outfileName.str().c_str());

    //print header with product codes, where available
    fout<<"\"ClusterID\""<<"\t";
    for(int k = 0; k < 2*trade1.numProdsK-1; k++){
        string currProd = mc1.prodNames[k];
        fout<<"\""<<currProd<<"\""<<"\t";
    }
    string currProd = mc1.prodNames[2*trade1.numProdsK-1];
    fout<<"\""<<currProd<<"\""<<endl;

    for(int z = 0; z < numModelsM; z++){
        // long double *max_element_ptr = max_element(all_TFProportions[z],all_TFProportions[z]+2*trade1.numProdsK);
        int tmpSum = 0;
        for(int t = 0; t < trade1.ntimes; t++)
        {
            tmpSum += mc1.clusterOcc[t][z].size();
            cout<<"writing TF: t = "<<t<<" z = "<<z<<" tmpSum = "<<tmpSum<<endl;
        }
        if(!tmpSum) continue;
        fout<<z<<"\t";
        for(int k = 0; k < 2*trade1.numProdsK-1; k++){
            if(isnan(all_TFProportions[z][k])){
                cout<<"Found nan in TF proportions = "<<all_TFProportions[z][k]<<endl;
                exit(-1);
            }
                fout<<all_TFProportions[z][k]<<"\t";
        }
        if(isnan(all_TFProportions[z][2*trade1.numProdsK-1])){
            cout<<"Found nan in TF proportions = "<<all_TFProportions[z][2*trade1.numProdsK-1]<<endl;
            exit(-1);
        }
        fout<<all_TFProportions[z][2*trade1.numProdsK-1]<<endl;
    }

    fout.close();

    logFile.close();

    // HJ: return 1;
    cout << "It is all done!" << endl;
    return 0;
 }

 bool sortFunc(pair<int, long double> p1, pair<int, long double> p2){
     return(p1.second > p2.second);
 } 

// [[Rcpp::export]]
int mainRcpp(string configTxt, int randomSeed) {

     long double nu;
     long double OP_w_ij_inp;
     string tmpString;
     int numThreads=8;
     int numModelsM;
     int maxiter;
     int nsim, total_mc_trials;
     string inputYears;
     ofstream logFile;
     ifstream configFile;
     long double eps;
     string q_filename,mu_filename,sigmasq_filename,pi_filename,zeta_filename;

     configFile.open(&configTxt[0u]);  
     int seedOffset = randomSeed;

     if(!configFile){
         cout<<"Could not open config file."<<endl;
         exit(-1);
     }
     
    cout<<"configFile is "<< "./config.txt" <<endl; 

     configFile >> tmpString;
     configFile >> tmpString;
     configFile >> nsim;

     configFile >> tmpString;
     configFile >> tmpString;
     configFile >> maxiter;

     configFile >> tmpString;
     configFile >> tmpString;
     configFile >> eps;

     configFile >> tmpString;
     configFile >> tmpString;
     configFile >> numModelsM;

     configFile >> tmpString;
     configFile >> tmpString;
     configFile >> nu;

     configFile >> tmpString;
     configFile >> tmpString;
     configFile >> OP_w_ij_inp;

     configFile >> tmpString;
     configFile >> tmpString;
     configFile >> inputYears;

     configFile >> tmpString;
     configFile >> tmpString;
     configFile >> numThreads;

     configFile >> tmpString;
     configFile >> tmpString;
     configFile >> total_mc_trials;   //total MC trials 

     configFile >> tmpString;
     configFile >> tmpString;
     configFile >> q_filename;

     configFile >> tmpString;
     configFile >> tmpString;
     configFile >> mu_filename;

     configFile >> tmpString;
     configFile >> tmpString;
     configFile >> sigmasq_filename;

     configFile >> tmpString;
     configFile >> tmpString;
     configFile >> pi_filename;

     configFile >> tmpString;
     configFile >> tmpString;
     configFile >> zeta_filename;

     configFile.close();

     return ztm(nsim, maxiter, eps, numModelsM, nu, OP_w_ij_inp, inputYears, numThreads, total_mc_trials, 
                q_filename, mu_filename, sigmasq_filename, pi_filename, zeta_filename, seedOffset);
}

