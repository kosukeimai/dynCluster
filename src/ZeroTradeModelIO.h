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

#ifndef ZEROTRADEMODELIO_H_
#define ZEROTRADEMODELIO_H_

#include "ZeroTradeModel.h"
#include "ChkpBase.h"

using namespace std;
using namespace boost;

class ZeroTradeModelIO : public ZeroTradeModel{
public:
    ZeroTradeModelIO(int numModelsMu, long double nu, long double OP_w_ij, string inputFilename, int maxiter, int numThreads);
	virtual ~ZeroTradeModelIO();
    void chkpClusters(int isim, ofstream &logFile);
    void chkpModelParams(int isim, ofstream &logFile);
private:
    ChkpBase chkpObj;
};

#endif /* ZEROTRADEMODELIO_H_ */
