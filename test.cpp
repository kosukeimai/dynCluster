#include <vector>
#include <iostream>
#include <cstdio>

using namespace std;

int main(int argc, char** argv){
    int ntimes = 3;
    int numProdsK = 3;

    vector<vector<long double> > prodCount(ntimes);
    vector<long double> arr(5,0.001);

    for(int t = 0; t < 5; t++){
        cout<<"arr[t] = "<<arr[t]<<endl;
    }

#if 0
    for(int t = 0; t < ntimes; t++){
        prodCount[t] = vector<long double>(numProdsK,0.001);
    }

    for(int t = 0; t < ntimes; t++){
        for(int k = 0; k < numProdsK; k++){
            cout<<"prodCount["<<t<<"]["<<k<<"] = "<<prodCount[t][k]<<"\t";
            printf("%f.\n",prodCount[t][k]);
        }
        cout<<endl;
    }
#endif

    return 1;
}
