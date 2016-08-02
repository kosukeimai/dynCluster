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

#ifndef ZEROTRADEMODEL_H_
#define ZEROTRADEMODEL_H_

#include <iostream>
//#include <vector>
#include <numeric>
#include <map>
#include <vector>
#include <cmath>
#include <cstdio>
#include <cassert>
#include <string>
//#include <omp.h>
#include <fstream>
#include <cstdio>
#include <cstring>
#include <iterator>
#include <boost/random.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/math/distributions/inverse_gamma.hpp>
#include <boost/math/distributions/gamma.hpp>
#include <boost/math/distributions/beta.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <algorithm>
#include <boost/regex.hpp>
#include "ZeroTradeModelTypedefs.h"
//#include "ChkpBase.h"


using namespace std;
using namespace boost;
using boost::math::inverse_gamma_distribution;
using boost::math::beta_distribution;
using boost::math::inverse_gamma;
using boost::math::tgamma;

#define PRIOR_OP 1
#define SIM_CHKP 0
#define DEBUG 0 
#define DEBUG2 1
#define FLIPDYADPAIRS 0
#define EXITONZEROVARIANCE 0

// To turn on/of the code I added to skip the (dyad, year) for those missing or excluded
#define HJSKIP 1

	//Global Constants
	const double C = 0.000001; //small constant to take into account zero trade
	const double PI = 3.14159265359;

class ZeroTradeModel {
public:
    ZeroTradeModel(int tmpM, long double tmpNu, long double tmpOP_w_ij_inp, string tmpYears, int nth);
	virtual ~ZeroTradeModel();
	void allocateMemory();
	int initializeModel(time_t seed);
	int initializeDebug(time_t seed);
    long double obsLogPosterior(int iter);
   // long double localLogLikelihood(int iter, int i, int t, int z);
    long double nextStep(int idyad, int t, int z,
                        vector<vector<vector<long double> > >& lll);
    long double logLikelihoodTh(int iter);
	int expectationTh(int iter, ofstream & logFile);
	int maximizationTh(int iter,bool &resetFlag, ofstream &logFile);
    void analytics(int isim, ofstream &logFile);
    // //void chkpClusters(int isim, ofstream &logFile);
    // //void chkpModelParams(int isim, ofstream &logFile);
    //void readCodeDescriptions(string);

    int ntimes;
	int numProdsK; //no. of commodities
	int numProdsK2; //2*(K-1)
	int numModelsM; //M
	int numModelsM2; //2*M
	int numDyads;
	int nu; 
    int maxFreqProd;
    int maxiter;
    time_t ZTMSeed; 
    long double llVal;
    long double llVal_prev;
    long double llVal_prev_prev;
    long double acc_llVal;
    long double bicVal;
    long double bicVal_prev;
    long double OP_w_ij;
    long double OP_w_ij_inp;
    string inputYears;
    vector<string> inputFiles;
    string prodCodePrefix;
    ProdMapType ProdMap;

    base_generator_type generatorb;

	//make all data members public for now
	vector<long double> pi;						//Eqn (3) of ZeroTradeModel.pdf
    vector<long double> piOld;						//Eqn (3) of ZeroTradeModel.pdf
	vector<long double> piTrue;						//Eqn (3) of ZeroTradeModel.pdf
    vector<vector<long double> > transP;
    vector<vector<long double> > debugTransP;
    //making zeta to be long double doesn't seem to make a difference
	vector<vector<vector<long double> > > zeta;	//Eqn (11) of ZeroTradeModel.pdf
	vector<vector<vector<long double> > > resp;	//Eqn (11) of ZeroTradeModel.pdf
	vector<vector<vector<long double> > > zeta1;	//Eqn (11) of ZeroTradeModel.pdf
    vector<vector<vector<vector<long double> > > > ksi;
    vector<vector<long double> > maxk;
    vector<vector<long double> > mink;
    vector<vector<vector<vector<long double> > > > lllk;

    //making q values to be long double don't make a difference
    vector<string> prodNames;   //column header if we view it as a df dyadpairsxproducts
    vector<vector<long double> > q;
    vector<long double> qDenominator;
    vector<vector<long double> > qOld;
    vector<vector<long double> > qTrue;
    vector<long double> qAll;
    vector<long double> qSigmaSqAll;
    vector<long double> muAll;					//overall mean
	vector<long double> sigmaSqAll;				//mean, Eqn (14) of ZeroTradeModel.pdf
	vector<vector<long double> > mu;				//mean, Eqn (14) of ZeroTradeModel.pdf
    vector<vector<long double> > muOld;				//mean, Eqn (14) of ZeroTradeModel.pdf
	vector<vector<long double> > muTrue;				//mean, Eqn (14) of ZeroTradeModel.pdf
	vector<vector<long double> > muNew;				//mean, Eqn (14) of ZeroTradeModel.pdf
	vector<vector<long double> > sigmaSq;				//mean, Eqn (14) of ZeroTradeModel.pdf
	vector<vector<long double> > sigmaSqOld;				//mean, Eqn (14) of ZeroTradeModel.pdf
	vector<vector<long double> > sigmaSqTrue;				//mean, Eqn (14) of ZeroTradeModel.pdf
	vector<vector<long double> > sigmaSqNew;				//mean, Eqn (14) of ZeroTradeModel.pdf

   /* HJ: In the extension, we assume that any negative tradeVolumn 
	  indicates the dyad does not exist for that year. Thus we have 
	  dExist[i][t] to keep track it by {1: exist; 0: non-exist}.
          But becasue of the product data is in xMap and then populated
	  to vMap which is still directional, both (i,j) and (j,i), 
	  we need to keep a record of all the (x,y) pairs:
	    1. record pairs missing for each year
	    2. reset the negative tradeVolumes to 0
	  Thus, we also need eMap to collect those (i, j, year).
	  At the end, eMap data will be used to construct the 
	  dExist that is 0/1 for (dyad, year).
   */
    vector<vector<int> > dExist;   
    vector<MapType> eMap; 
    vector<int> eBegin;
    vector<int> eEnd;

    //xMap file is of the type (i,j)->[ntimes x numProdsK]
    //vMap file is of the type (i,j)->[ntimes x numProdsK]
	VMapType xMap; 	//raw ijk data from .dta file NxNxK
	VMapType vMap; 	//log-proportion data Nx(N-1)x( K-1), Section 1.2 of ZeroTradeModel.pdf
	MapType dyadMap;    //map for stacked dyad (i,j) -> dyadID 
    RMapType revDyadMap;//reverse map for stacked dyad dyadID -> (i,j)  
	MapType dMap;       //map for raw dyad (i,j) -> dyadID 

    //w, tmpD has dimensions ntimes x num_dyads x numProdsK 
	vector<vector<vector<long double> > > w;     //dyad-specific dependent variable, Section 1.2 of ZeroTradeModel.pdf
    vector<vector<vector<int> > > tmpD;      //dyad-specific variable indicating whether a product is 0 or not. 

    //y, d, dDash has dimensions num_dyads x ntimes x numProdsK2
	vector<vector<vector<long double> > > y;     //dyad-specific dependent variable, Section 1.2 of ZeroTradeModel.pdf
    vector<vector<vector<int> > > d;      //dyad-specific variable indicating whether a product is 0 or not. 
    vector<vector<vector<int> > > dDash;      //dyad-specific variable indicating whether a product is 0 or not. 
                                    //TODO: could use a bool for d
    vector<vector<long double> >  kTmpSum; 
    vector<vector<long double> > wts;
    //vector<vector<int> > clusterOcc;
    vector<int> maxZeta;    //not sure if this is necessary
    vector<vector<long double> > bigN;
    //vector<vector<vector<long double> > > bigB;
    //vector<vector<long double> > expectLogPi;
	//vector<vector< vector<long double> > > expectLogOneMinusQ;
    //vector<vector< vector<long double> > > expectLogQ;
	vector< vector<vector<long double> > > expectLogPhi;
	vector< vector<vector<long double> > > expectLogMuSigmaSq;
	vector< vector<long double> > expectLogSigmaSq;
    //ChkpBase chkpObj;

    //void calcR(int iter,ofstream &logFile);
    //void calcExpectLogPi(int iter,ofstream &logFile);
    //void calcExpectLogQ(int iter,ofstream &logFile);
    //void calcExpectLogOneMinusQ(int iter,ofstream &logFile);
    // //int calcExpectLogPhi(int iter,ofstream &logFile);
    // //int calcExpectLogSigmaSq(vector<vector<long double> >& expLogSigmaSq);
    // //int calcExpectLogMuSigmaSq(vector<vector<vector<long double> > >& expLogMuSigmaSq);
    //void calcBigN(int iter,ofstream& logFile);
    // void calcBigB(int iter);
    //void calcPi(void);
    //void calcQ(void);
    //void calcMu(void);
    //void calcSigmaSq(void);
    void chkpModelParams(int isim, ofstream &logFile);

private:
    int numThreads;
    bool zetaFlag;
    long double logDNorm(long double val,long double mean,long double var); 
    long double inverseGamma(long double var, long double alpha, long double beta); 
    long double betaDistribution(long double var, long double alpha, long double beta);
    long double logBetaDistribution2(long double var, long double alpha, long double beta, long double logBetaVal);
    long double calcDiGamma(long double val);
    vector<vector<vector<long double> > >tmpProd;
    //TODO: uncomment
	void calculateMeanAll();
	void calculateCovarianceAll();
    void countProds();
    void getInputFilenames();
	void readInput();
    void reorganizeProducts(); 
    //void fixZeroTradedProds();
    //void flipDyadPair(int);
};

#endif /* ZEROTRADEMODEL_H_ */
