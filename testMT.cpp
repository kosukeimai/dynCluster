#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp> 
#include <iostream>
#include <ctime>

using namespace std;

int main(int argc, char** argv){

    boost::mt19937 rnd;
    boost::uniform_real<> uni(0.0,1.0);
    boost::variate_generator<boost::mt19937&,boost::uniform_real<> > gen(rnd,uni);

    rnd.seed(time(NULL));

    cout<<gen()<<endl;
}
