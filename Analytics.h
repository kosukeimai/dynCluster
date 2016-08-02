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

#ifndef Analytics_H_
#define Analytics_H_

#include <iostream>
//#include <vector>
#include <numeric>
#include <map>
#include <vector>
#include <cmath>
#include <cstdio>
#include <cassert>
#include <cstring>
#include <string>
#include <omp.h>
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
//#include "armadillo"
#include <algorithm>
#include <boost/regex.hpp>


using namespace std;
using namespace boost;
//using namespace arma;
using boost::math::inverse_gamma_distribution;
using boost::math::beta_distribution;
using boost::math::inverse_gamma;
using boost::math::tgamma; 

//#define DEBUG //set debug flag

	//Global Constants
	//const double C = 0.000001; //small constant to take into account zero trade
	const double KAPPA = 0.01;
	//const double PI = 3.14159265359;
	//fails const double LAMBDA = 0.001; //tuning parameter for complete-data log-likelihood
	const double LAMBDA = 1.0; //tuning parameter for complete-data log-likelihood
    typedef boost::minstd_rand base_generator_type;
    typedef map<string,string> CodeMapType;

class Analytics{
public:
    Analytics(int,int,int,int,int,string,string,string,string,string);
	virtual ~Analytics();
    void allocateMemory(vector<long double>& piArg,
                        vector<vector<long double> >& sigmaSqArg,
                        vector<vector<long double> >& muArg,
                        vector<vector<long double> >& qArg,
                        vector<vector<vector<long double> > >& zetaArg,
                        vector<vector<vector<int > > >& clusterOccArg,
                        long double ****prodProportionArg);
    void doMC(int z,int rank);
    void analytics(int rank);
    void readCodeDescriptions(string);

	int numProdsK; //no. of products
    int ntimes;
	int numProdsK2; //2*(K-1)
	int numModelsM; //M
	int numModelsM2; //2*M
	int numDyads;
    int maxFreqProd;
    int mc_trials;
    string prodCodePrefix;
    time_t ZTMSeed;
    string q_filename;
    string mu_filename;
    string sigmasq_filename;
    string pi_filename;
    string zeta_filename;
    string TF_filename;


    CodeMapType CodeMap;

    base_generator_type generatorb;

	//make all data members public for now
    vector<string> prodNames;
	vector<long double> pi;						//Eqn (3) of ZeroTradeModel.pdf
	vector<vector<vector<long double> > > zeta;	//Eqn (11) of ZeroTradeModel.pdf

    //making q values to be long double don't make a difference
    vector<vector<long double> > q;
	vector<vector<long double> > mu;				//mean, Eqn (14) of ZeroTradeModel.pdf
	vector<vector<long double> > sigmaSq;				//mean, Eqn (14) of ZeroTradeModel.pdf
    //vector<vector<long double> > expectedMu;
    long double ***prodProportion;
    vector<vector<vector<int> > > clusterOcc;

    vector<vector<long double> > cosineDistMat;
    vector<vector<long double> > euclideanDistMat;

private:
    int numThreads;
    vector<vector<vector<long double> > >tmpProd;
	//void calculateMeanAll();
	//void calculateCovarianceAll();
	void readInput();
    void readInputQ(string, vector<vector<long double> > &, vector<string>&);
    void readInputMu(string, vector<vector<long double> > &);
    void readInputSigmaSq(string, vector<vector<long double> > &);
    void readInputPi(string, vector<long double>&);
    void readInputZeta(string,vector<vector<long double> >&);
};

#endif /* Analytics_H_ */
