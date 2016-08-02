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
#include <numeric>
#include <algorithm>
#include <sstream> 
#include <iostream>
#include <boost/random.hpp>
#include <time.h>
#include <mpi.h>
using namespace std;
using namespace boost; 

 bool sortFunc(pair<int, long double> p1, pair<int, long double> p2); 

 int main(int argc,char** argv){

     //int numCountries, numProductsK, numModelsM;
     long double nu, OP_w_ij, bayes_beta, bayes_alpha, bayes_tau; 
     string tmpString;  
     vector<int> nemptyIdx;
     int success=0;
     //int nMC;
     //int numThreads; 
     int numProductsK, numModelsM;
     int nsim, total_mc_trials, mc_trials;
     bool successFlag = 0;
     bool resetFlag = 0;
     string inputFilename;
     ofstream logFile;
     ifstream configFile;
     long double max_ll_new;
     int max_sim;
     //long double eps; 
     bool checkpointFlag = 0;
     vector<pair<int, long double> > simDF;
     string q_filename,mu_filename,sigmasq_filename,pi_filename,zeta_filename;
     int np, rank;
     long double** all_TFProportions; 
     double start_time, end_time;

     if(argc != 2){
         cout<<"Usage is: ./MC configFile ";
     }

     MPI_Init(&argc,&argv);
     MPI_Comm_rank(MPI_COMM_WORLD,&rank);
     MPI_Comm_size(MPI_COMM_WORLD,&np);

     start_time = MPI_Wtime(); 

     cout<<"Rank "<<rank<<" of total "<<np<<" is alive."<<endl;

     srand(time(0)); 
     configFile.open(argv[1]);  

     if(!configFile){
         cout<<"Could not open config file."<<endl;
         exit(-1);
     }
     cout<<"configFile is "<<argv[1]<<endl; 

     configFile >> tmpString;
     configFile >> tmpString;
     configFile >> tmpString;   //nsim 

     configFile >> tmpString;
     configFile >> tmpString;
     configFile >> tmpString;   //maxiter 

     configFile >> tmpString;
     configFile >> tmpString;
     configFile >> tmpString;   //eps 

     configFile >> tmpString;
     configFile >> tmpString;
     configFile >> numProductsK;

     configFile >> tmpString;
     configFile >> tmpString;
     configFile >> numModelsM;

     configFile >> tmpString;
     configFile >> tmpString;
     configFile >> tmpString;   //nu 

     configFile >> tmpString;
     configFile >> tmpString;
     configFile >> OP_w_ij;

     configFile >> tmpString;
     configFile >> tmpString;
     configFile >> tmpString;   //input file 

     configFile >> tmpString;
     configFile >> tmpString;
     configFile >> tmpString;   //threads 

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

     if(rank == 0){
        cout<<"read in: total MC trials = "<<total_mc_trials<<endl;
        cout<<"read in: K = "<<numProductsK<<endl;
        cout<<"read in: Z = "<<numModelsM<<endl;
        cout<<"read in: OP_w_ij = "<<OP_w_ij<<endl;
        cout<<"read in: q checkpoint file = "<<q_filename<<endl;
        cout<<"read in: mu checkpoint file = "<<mu_filename<<endl;
        cout<<"read in: sigmasq checkpoint file = "<<sigmasq_filename<<endl;
        cout<<"read in: pi checkpoint file = "<<pi_filename<<endl;
        cout<<"read in: zeta checkpoint file = "<<zeta_filename<<endl;
    }

     //TODO: MPI_Init, etc
     //read in the checkpoint files for which to do the MC simulation
     //MC mc1(maxiter,numThreads,numModelsM,numProductsK,q_filename,mu_filename,sigmasq_filename,pi_filename,zeta_filename);
     mc_trials = (total_mc_trials/np)*np; 
     Analytics mc1(mc_trials,0,numModelsM,numProductsK,q_filename,mu_filename,sigmasq_filename,pi_filename,zeta_filename);

     mc1.readCodeDescriptions("S2_out.txt");

     cout<<"rank = "<<rank<<" before call to analytics."<<endl; 
     mc1.analytics(rank);
     if(rank == 0){
         for(int z = 0; z < mc1.numModelsM; z++){
             cout<<"Cluster "<<z<<" has "<<mc1.clusterOcc[z].size()<<" clusters."<<endl;
         }
     }


     for(int z = 0; z < numModelsM; z++){
        MPI_Barrier(MPI_COMM_WORLD);
 //       mc1.doMC(z,rank);
//#if 0
        if(mc1.clusterOcc[z].size() != 0){
            mc1.doMC(z,rank);
            //mc1.doMomentAverage(z,rank);
            //mc1.doMomentAverage(z,rank);
            //TODO: do an all reduce on the product proportions and average over np
        }
//#endif
     }

     MPI_Barrier(MPI_COMM_WORLD); 

     all_TFProportions = (long double**)malloc(numModelsM*sizeof(long double*));
     for(int z = 0; z < numModelsM; z++){
          all_TFProportions[z] = (long double*)malloc(2*numProductsK*sizeof(long double));
     }

     for(int z = 0; z < numModelsM; z++){
        MPI_Reduce(mc1.prodProportion[z],all_TFProportions[z],2*numProductsK,MPI_LONG_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
        for(int k = 0; k < 2*numProductsK; k++){
            all_TFProportions[z][k] /= (long double)np;
        }
     }

     for(int z = 0; z < numModelsM; z++){
         for(int k = 0; k < 2*numProductsK; k++){
             all_TFProportions[z][k] = mc1.prodProportion[z][k];
         }
     }

     //MPI_Barrier(MPI_COMM_WORLD); 

     if(rank == 0){
     //   mc1.printTFProportions();
        ostringstream outfileName;
        outfileName<<"TF_PROPORTION_ALLCLUSTER_"<<rand()<<".txt";
        cout<<"writing to output file "<<outfileName.str()<<endl;
        ofstream fout;
        fout.open(outfileName.str().c_str());

        //print header with product codes, where available 
        for(int k = 0; k < 2*numProductsK-1; k++){
            string currProd = mc1.prodNames[k];
            //cout<<"looking for product "<<currProd<<endl;
            SITCMapType::iterator mit = mc1.SITCMap.find(currProd);
            fout<<"\""<<currProd<<"\""<<"\t";
        }
        string currProd = mc1.prodNames[2*numProductsK-1];
        SITCMapType::iterator mit = mc1.SITCMap.find(currProd);
        fout<<"\""<<currProd<<"\""<<endl;

        //print header with product descriptions, where available
        //copy(mc1.prodNames.begin(),mc1.prodNames.end(),ostream_iterator<string>(fout,"\t"));
        for(int k = 0; k < 2*numProductsK-1; k++){
            string currProd = mc1.prodNames[k];
            //cout<<"looking for product "<<currProd<<endl;
            SITCMapType::iterator mit = mc1.SITCMap.find(currProd);
            if(mit != mc1.SITCMap.end()){
                string currProdDesc = mit->second;
                fout<<"\""<<currProdDesc<<"\""<<"\t";
            }
            else{
                //no description available
                fout<<"\""<<currProd<<"\""<<"\t";
            }
        }
        currProd = mc1.prodNames[2*numProductsK-1];
        mit = mc1.SITCMap.find(currProd);
        if(mit != mc1.SITCMap.end()){
            string currProdDesc = mit->second;
            fout<<"\""<<currProdDesc<<"\""<<endl;
        }
        else{
            //no description available
            fout<<"\""<<currProd<<"\""<<endl;
        }

        cout<<"Finished writing header of the TF file."<<endl;

        for(int z = 0; z < numModelsM; z++){
            //TODO: sort the element indices in order of maximum product
            //proportion and print them in terms of the column name in a
            //separate file. 
            long double *max_element_ptr = max_element(all_TFProportions[z],all_TFProportions[z]+2*numProductsK);
            //fout<<(max_element_ptr - all_TFProportions[z])<<" ";
            if(mc1.clusterOcc[z].size() != 0){
                fout<<z<<"\t";
                for(int k = 0; k < 2*numProductsK-1; k++){
                    if(isnan(all_TFProportions[z][k])){
                        cout<<"Found nan in TF proportions = "<<all_TFProportions[z][k]<<endl;
                        MPI_Abort(MPI_COMM_WORLD,-99);
                        exit(-1);
                    }
                        fout<<all_TFProportions[z][k]<<"\t";
                }
                if(isnan(all_TFProportions[z][2*numProductsK-1])){
                    cout<<"Found nan in TF proportions = "<<all_TFProportions[z][2*numProductsK-1]<<endl;
                    MPI_Abort(MPI_COMM_WORLD,-99);
                    exit(-1);
                }
                fout<<all_TFProportions[z][2*numProductsK-1]<<endl;
            }
        }
        fout.close();
        cout<<"ran a total of "<<total_mc_trials<<" MC trials."<<endl;
     }

     cout<<"rank = "<<rank<<" before last barrier."<<endl;

     MPI_Barrier(MPI_COMM_WORLD);

     end_time = MPI_Wtime();

     if(rank == 0){
         cout<<"elapsed time is "<<(end_time - start_time)<<endl;
     }

     MPI_Finalize();

	 return 1;
 }
