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

#ifndef CHKPBASE_H_
#define CHKPBASE_H_

#include "ZeroTradeModelTypedefs.h"

using namespace std;
using namespace boost;

class ChkpBase{
public:
    //int numModelsM;
    int numProdsK;
    int numDyads;
    int ntimes;
    //int numModelsM2;
    //ChkpBase(int numProdsK, int numModelsMu, long double nu, string inputFilename, int numThreads);
    ChkpBase();
	virtual ~ChkpBase();

    //void chkpClusters(int isim, ofstream &logFile);
    void chkpClusters(int isim,ofstream &logFile,
		      int local_numProdsK,
		      int local_numModelsM,
		      int local_numDyads,
		      int local_ntimes,
		      vector<string>& local_prodNames,
		      string local_prodCodePrefix,
		      vector<vector<long double> >& local_q,
		      vector<vector<long double> >& local_mu,
		      vector<vector<long double> >& local_sigmaSq,
		      vector<vector<vector<long double> > >& local_zeta,
		      vector<long double>& local_pi,
		      RMapType& local_revDyadMap,
		      ProdMapType& local_ProdMap,
		      VMapType& local_vMap,
		      long double& local_llVal,
		      long double& local_bicVal,
		      vector<vector<int> >& local_dExist,
		      vector<int>& local_eBegin,
		      vector<int>& local_eEnd);

    void chkpModelParams(int isim,ofstream &logFile,
			 int local_numProdsK,
			 int local_numModelsM,
			 int local_numDyads,
			 int ntimes,
			 vector<string>& local_prodNames,
			 string local_prodCodePrefix,
			 vector<vector<long double> >& local_q,
			 vector<vector<long double> >& local_mu,
			 vector<vector<long double> >& local_sigmaSq,
			 vector<vector<vector<long double> > >& local_zeta,
			 vector<long double>& local_pi,
			 vector<vector<long double> >& local_transP,
			 RMapType& local_revDyadMap,
			 long double& local_llVal,
			 long double& local_bicVal,
			 vector<vector<int> >& local_dExist,
			 vector<int>& local_eBegin,
			 vector<int>& local_eEnd);

    //void chkpModelParams(int isim, ofstream &logFile);
    void readModelParams(int isim,ofstream &logFile,
                                int local_numProdsK,
                                int local_numModelsM,
                                int local_numDyads,
                                vector<string>& local_prodNames,
                                vector<vector<long double> >& local_q,
                                vector<vector<long double> >& local_mu,
                                vector<vector<long double> >& local_sigmaSq,
                                vector<vector<long double> >& local_zeta,
                                vector<long double>& local_pi);
    void readInputQ(string filename,vector<vector<long double> >& qVecVec,vector<string>& namesVec);
    void readInputMu(string filename,vector<vector<long double> >& muVecVec);
    void readInputSigmaSq(string filename,vector<vector<long double> >& sigmaSqVecVec);
    void readInputPi(string filename,vector<long double>& piVec);
};

#endif /* CHKPBASE_H_ */
