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

#include "ZeroTradeModelIO.h"

ZeroTradeModelIO::ZeroTradeModelIO(int tmpM, long double tmpNu,
                                long double tmpOP_w_ij, string inputYears, int maxiter, int nth):
                                ZeroTradeModel(tmpM,tmpNu,tmpOP_w_ij,inputYears,nth){
}

ZeroTradeModelIO::~ZeroTradeModelIO() {

}

void ZeroTradeModelIO::chkpModelParams(int isim,ofstream &logFile){
    chkpObj.chkpModelParams(isim, logFile,
			    numProdsK,
			    numModelsM,
			    numDyads,
			    ntimes,
			    prodNames,
			    prodCodePrefix,
			    q,
			    mu,
			    sigmaSq,
			    zeta,
			    pi,
			    transP,
			    revDyadMap,
			    llVal,
			    bicVal,
			    dExist,
			    eBegin,
			    eEnd);
}

void ZeroTradeModelIO::chkpClusters(int isim,ofstream &logFile){
    chkpObj.chkpClusters(isim, logFile,
			 numProdsK,
			 numModelsM,
			 numDyads,
			 ntimes,
			 prodNames,
			 prodCodePrefix,
			 q,
			 mu,
			 sigmaSq,
			 zeta,
			 pi,
			 revDyadMap,
			 ProdMap,
			 vMap,
			 llVal,
			 bicVal,
			 dExist, 
			 eBegin, 
			 eEnd);
}
