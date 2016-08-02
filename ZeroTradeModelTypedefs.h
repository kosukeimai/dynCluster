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

#ifndef ZEROTRADEMODELTYPEDEFS_H_
#define ZEROTRADEMODELTYPEDEFS_H_

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


using namespace std;
using namespace boost;
using boost::math::inverse_gamma_distribution;
using boost::math::beta_distribution;
using boost::math::inverse_gamma;
using boost::math::tgamma; 

#define SIM_CHKP 0 

    typedef pair<int,int> MatKey;
    typedef map<MatKey,int> MapType;
    typedef map<int,MatKey> RMapType;
    typedef map<MatKey,vector<vector<long double> > > VMapType;
    typedef boost::minstd_rand base_generator_type; 
    typedef map<string,string> ProdMapType;

#endif /*ZEROTRADEMODELTYPEDEFS_H_*/
