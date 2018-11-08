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

#include "ZeroTradeModel.h"

ZeroTradeModel::ZeroTradeModel(int tmpM, long double tmpNu, long double tmpOP_w_ij_inp, string tmpYears, int nth): 
                                                        numModelsM(tmpM), nu(tmpNu), OP_w_ij_inp(tmpOP_w_ij_inp),
                                                        inputYears(tmpYears), numThreads(nth) {

    //ntimes = 6;
    //ntimes = 5;
    //ntimes = 4;
    //ntimes = 3;
    //ntimes = 2;

    printf("ZeroTradeModel constructor: OP_w_ij_inp = %Lf.\n",OP_w_ij_inp);

	//use values passed to constructor to set some variables
	numModelsM2 = 2*numModelsM;
    maxFreqProd = -999;


    getInputFilenames();
    ntimes = inputFiles.size();

    countProds(); //this sets numProdsK

    //numProdsK = 3;
	numProdsK2 = 2*(numProdsK-1);
    srand(time(NULL));
#if DEBUG2
    cout<<"Before call to allocateMemory()."<<endl;
#endif
	allocateMemory();
#if DEBUG2
    cout<<"After call to allocateMemory()."<<endl;
#endif

#if DEBUG2
    cout<<"Before initialize model/debug.\n"<<endl;
#endif
    //TODO: uncomment
	//initializeModel(time(0));
	//initializeDebug(time(0));
#if DEBUG2
    cout<<"After initialize model/debug.\n"<<endl;
#endif
}

//Allocate memory for all the vectors and matrices we are going to use
void ZeroTradeModel::allocateMemory(){

    cout<<"Entering allocateMemory()."<<endl;

    pi.resize(numModelsM2);
    piOld.resize(numModelsM2);

    cout<<"After pi."<<endl;

    transP.resize(numModelsM2);
    for(int z = 0; z < numModelsM2; z++){
        transP[z].resize(numModelsM2);
    }

    debugTransP.resize(numModelsM2);
    for(int z = 0; z < numModelsM2; z++){
        debugTransP[z].resize(numModelsM2);
    }
    cout<<"After setting debugTransP."<<endl;

/* HJ 
#if 0
    muAll.resize(ntimes);
    sigmaSqAll.resize(ntimes);
    for(int t = 0; t < ntimes; t++){
        muAll[t].resize(numProdsK2);
        sigmaSqAll[t].resize(numProdsK2);
    }
#endif
*/

    muAll.resize(numProdsK2);
    sigmaSqAll.resize(numProdsK2);

    cout<<"After muAll."<<endl;

/* HJ
#if 0
    qAll.resize(ntimes);
    qSigmaSqAll.resize(ntimes);
    for(int t = 0; t < ntimes; t++){
        qAll[t].resize(2*numProdsK);
        qSigmaSqAll[t].resize(2*numProdsK);
    }
#endif
*/


    qAll.resize(2*numProdsK);
    qSigmaSqAll.resize(2*numProdsK);

    cout<<"After q."<<endl;

    prodNames.resize(numProdsK);

    cout<<"After prodNames."<<endl;

/*
#if 0
    sigmaSq.resize(ntimes);
    sigmaSqOld.resize(ntimes);
    sigmaSqNew.resize(ntimes);
    mu.resize(ntimes);
    muOld.resize(ntimes);
    muNew.resize(ntimes);
    kTmpSum.resize(ntimes);
    for(int t = 0; t < ntimes; t++){
        sigmaSq[t].resize(numModelsM2);
        sigmaSqOld[t].resize(numModelsM2);
        sigmaSqNew[t].resize(numModelsM2);
        mu[t].resize(numModelsM2);
        muOld[t].resize(numModelsM2);
        muNew[t].resize(numModelsM2);
        kTmpSum[t].resize(numModelsM2);
        for(int i = 0; i < numModelsM2; i++){
            sigmaSq[t][i].resize(numProdsK2);
            sigmaSqOld[t][i].resize(numProdsK2);
            sigmaSqNew[t][i].resize(numProdsK2);
            mu[t][i].resize(numProdsK2);
            muOld[t][i].resize(numProdsK2);
            muNew[t][i].resize(numProdsK2);
            kTmpSum[t][i].resize(numProdsK2); //used in maximization
        }
    }
#endif
*/

    sigmaSq.resize(numModelsM2);
    sigmaSqOld.resize(numModelsM2);
    sigmaSqNew.resize(numModelsM2);
    mu.resize(numModelsM2);
    muOld.resize(numModelsM2);
    muNew.resize(numModelsM2);
    kTmpSum.resize(numModelsM2);
    for(int i = 0; i < numModelsM2; i++){
        sigmaSq[i].resize(numProdsK2);
        sigmaSqOld[i].resize(numProdsK2);
        sigmaSqNew[i].resize(numProdsK2);
        mu[i].resize(numProdsK2);
        muOld[i].resize(numProdsK2);
        muNew[i].resize(numProdsK2);
        kTmpSum[i].resize(numProdsK2);
    }

    cout<<"after sigmaSq."<<endl;

    q.resize(numModelsM2);
    qDenominator.resize(numModelsM2);
    qOld.resize(numModelsM2);
    qTrue.resize(numModelsM2);
    for(int m = 0; m < numModelsM2; m++){
        q[m].resize(2*numProdsK);
        qOld[m].resize(2*numProdsK);
        qTrue[m].resize(2*numProdsK);
    }

    cout<<"after q."<<endl;

    muTrue.resize(numModelsM2);
    sigmaSqTrue.resize(numModelsM2);
    for(int m = 0; m < numModelsM2; m++){
        muTrue[m].resize(numProdsK2);
        sigmaSqTrue[m].resize(numProdsK2);
    } 

    piTrue.resize(numModelsM2);


    muAll.resize(numProdsK2);

   //expectedMu has the same size as mu

/* HJ
#if 0
   clusterOcc.resize(ntimes);
   for(int t = 0; t < ntimes; t++){
        clusterOcc[t].resize(numModelsM);
   }
#endif

#if 0
   expectLogPi.resize(ntimes);
   for(int t = 0; t < ntimes; t++){
	expectLogPi[t].resize(numModelsM2);
   }

   bigN.resize(ntimes);
   bigB.resize(ntimes);
   for(int t = 0; t < ntimes; t++){
    bigN[t].resize(numModelsM2);
    bigB[t].resize(numModelsM2);
    for(int z = 0; z < numModelsM2; z++){
        bigN[t][z].resize(2*numProdsK);
        bigB[t][z].resize(2*(numProdsK-1));
    }
   }
#endif

#if 0
    cout<<"Before allocating memory for izk."<<endl;
	expectLogPhi.resize(numDyads);
	expectLogMuSigmaSq.resize(numDyads);
	for(int i = 0; i < numDyads; i++){
		expectLogPhi[i].resize(numModelsM2);
		expectLogMuSigmaSq[i].resize(numModelsM2);
		for(int z = 0; z < numModelsM2; z++){
			expectLogPhi[i][z].resize(2*(numProdsK-1));
			expectLogMuSigmaSq[i][z].resize(2*(numProdsK-1));
		}
	}
    cout<<"After allocating memory for izk."<<endl;

	expectLogSigmaSq.resize(numModelsM2);
	for(int z = 0; z < numModelsM2; z++){
		expectLogSigmaSq[z].resize(2*(numProdsK-1));
	}
#endif

#if 0
    expectLogOneMinusQ.resize(ntimes);
    expectLogQ.resize(ntimes);
    for(int t = 0; t < ntimes; t++){
        expectLogOneMinusQ[t].resize(numModelsM2);
        for(int z = 0; z < numModelsM2; z++){
            expectLogOneMinusQ[t][z].resize(2*numProdsK);
        }

        expectLogQ[t].resize(numModelsM2);
        for(int z = 0; z < numModelsM2; z++){
            expectLogQ[t][z].resize(2*numProdsK);
        }
    }
#endif
*/

}

//
// HJ on Mar30,2016: This initializeDebug is not used currently, so I am not changing it with eBegin and eEnd. */
//
int ZeroTradeModel::initializeDebug(time_t seed){
    cout<<"in initializeDebug."<<endl;

	ifstream infile[ntimes];
	//ifstream infile[3];
	//ifstream infile[2];
    seed = time(NULL);
    ZTMSeed = seed;
    generatorb.seed(static_cast<unsigned int>(ZTMSeed));
    boost::uniform_real<> uni_dist(-10.0,10.0);
    boost::uniform_real<> uni_dist2(0.0,1.0);
    boost::variate_generator<base_generator_type&,boost::uniform_real<> > uni(generatorb,uni_dist);
    boost::variate_generator<base_generator_type&,boost::uniform_real<> > uni2(generatorb,uni_dist2);

    //open all the input files
    //infile[0].open(inputFile.c_str());
    infile[0].open("test0.txt");
    infile[1].open("test1.txt");
    //infile[2].open("test2.txt");
    //infile[3].open("test3.txt");
    //infile[4].open("test4.txt");
    //infile[5].open("test5.txt");

    for(int t = 0; t < ntimes; t++){
        if(infile[t].fail()){
            cout<<"error opening input file: infile["<<t<<"] = "<<infile[t]<<endl;
        }
    }

    cout<<"Opened input files for reading."<<endl;

    numDyads = 10000;
    numProdsK = 3;
    numProdsK2 = 2*(numProdsK-1);

    OP_w_ij = OP_w_ij_inp/(2*numDyads);
    printf("OP_w_ij_inp = %Lf.\n",OP_w_ij_inp);
    printf("OP_w_ij = %Lf.\n",OP_w_ij);

    //wts.resize(numDyads);
    zeta.resize(numDyads);
    zeta1.resize(numDyads);
    ksi.resize(numDyads);
    for(int i = 0; i < numDyads; i++){
        zeta[i].resize(ntimes);
        zeta1[i].resize(ntimes);
        ksi[i].resize(ntimes);
        //wts[i].resize(ntimes);
        for(int t = 0; t < ntimes; t++){
            long double tmpDenominator = 0.0;

            zeta[i][t].resize(numModelsM2);
            zeta1[i][t].resize(numModelsM2);
            ksi[i][t].resize(numModelsM2);
            for(int z = 0; z < numModelsM2; z++){
                ksi[i][t][z].resize(numModelsM2);
            }

            for(int k = 0; k < 2*numProdsK; k++){
                tmpDenominator += d[i][t][k];
            }
            //wts[i][t] = 1.0/(exp(nu*tmpDenominator));
        }
    }

    cout<<"Initialized values."<<endl;


    //TODO: allocate space for y
    y.resize(numDyads);
    for(int i = 0; i < numDyads; i++){
        y[i].resize(ntimes);
        for(int t = 0; t < ntimes; t++){
            y[i][t].resize(numProdsK2);
        }
    }


    for(int t = 0; t < ntimes; t++){
        cout<<"Reading in file for time t = "<<t<<endl;

        for(int z = 0; z < numModelsM; z++){
            for(int k = 0; k < 2*numProdsK; k++){
                //ofile<<"q["<<z<<"]["<<k<<"]="<<q[z][k]<<endl;
                infile[t]>>qTrue[z][k];
            }
        }

        for(int z = 0; z < numModelsM; z++){
            for(int k = 0; k < numProdsK; k++){
                qTrue[z+numModelsM][k] = qTrue[z][k+numProdsK];
                qTrue[z+numModelsM][k+numProdsK] = qTrue[z][k];
            }
        }

        cout<<"Read in qTrue values from file t = "<<t<<"."<<endl;

        for(int k = 0; k < numProdsK2; k++){
            muTrue[0][k] = 2*k;
            muTrue[0+numModelsM][k] = 2*k;
        }

        for(int k = 0; k < numProdsK2; k++){
            muTrue[1][k] = 10*k;
            muTrue[1+numModelsM][k] = 10*k;
        }

        for(int k = 0; k < numProdsK2; k++){
            muTrue[2][k] = 15*k;
            muTrue[2+numModelsM][k] = 15*k;
        }

        for(int z = 0; z < numModelsM2; z++){
            for(int k = 0; k < numProdsK2; k++){
                sigmaSqTrue[z][k] = 1.0;
            }
        }

        piTrue[0] = 0.33;
        piTrue[1] = 0.33;
        piTrue[2] = 0.34;
        piTrue[3] = 0.0;
        piTrue[4] = 0.0;
        piTrue[5] = 0.0;

        for(int i =0 ; i < numDyads; i++){
            for(int k = 0; k < numProdsK2; k++){
                infile[t] >> y[i][t][k];
            }
        }

        cout<<"Read in y values from file t = "<<t<<"."<<endl;
    }

    cout<<"Read in true q values and y values."<<endl;

    for(int t = 0; t < ntimes; t++)
        infile[t].close();
    cout<<"Finished reading input files."<<endl;

//TODO: change the order of dDash and d so that the indexing is [k][i] rather
//than [i][k].
    d.resize(numDyads);
    dDash.resize(numDyads);
    for(int i = 0; i < numDyads; i++){
        d[i].resize(ntimes);
        dDash[i].resize(ntimes);
        for(int t = 0; t < ntimes; t++){
            d[i][t].resize(2*numProdsK);
            dDash[i][t].resize(numProdsK2);
            //d[i].resize(2*numProdsK);
            //dDash[i].resize(numProdsK2);
            for(int k = 0; k < numProdsK2; k++){
                dDash[i][t][k] = 0;
            }
        }
    }

    for(int i = 0; i < numDyads; i++){
        for(int t = 0; t < ntimes; t++){
            for(int k = 0; k < numProdsK2; k++){
                if(y[i][t][k] == 0.0){
                    dDash[i][t][k] = 1;
                }
                else{
                    dDash[i][t][k] = 0;
                }
            }
        }
    }

    for(int i = 0; i < numDyads;i++){
        for(int t = 0; t < ntimes; t++){
            for(int k = 0; k < numProdsK-1; k++){
                d[i][t][k] = dDash[i][t][k];
            }

//            if(uni2() < 0.5){
               d[i][t][numProdsK-1] = 0;
//            }
//            else{
//                d[i][t][numProdsK-1] = 1;
//            }

            for(int k = numProdsK-1; k < 2*numProdsK - 2; k++){
                d[i][t][k+1] = dDash[i][t][k];
            }

//            if(uni2() < 0.5){
                d[i][t][2*numProdsK - 1] = 0;
//            }
//            else{
//                d[i][t][2*numProdsK - 1] = 1;
//            }

        }
    }

        //TODO: change zeta dimensions to numModelsM2xnumDyads, E- and M- have
        //different access patterns but M- might be giving a bigger hit
        zeta.resize(numDyads);
        zeta1.resize(numDyads);
        wts.resize(numDyads);
        for(int i = 0; i < numDyads; i++){
            zeta[i].resize(ntimes);
            zeta1[i].resize(ntimes);
            wts[i].resize(ntimes);
            for(int t = 0; t < ntimes; t++){
                long double tmpDenominator = 0.0;

                zeta[i][t].resize(numModelsM2);
                zeta1[i][t].resize(numModelsM2);

                for(int k = 0; k < 2*numProdsK; k++){
                    tmpDenominator += d[i][t][k];
                }
                wts[i][t] = 1.0/(exp(nu*tmpDenominator));
            }
        }


    //TODO: uncomment
	//calculate overall mean and covariance for the whole dataset
	calculateMeanAll(); //sets vector muAll
    //cout<<"Calculated overall mean."<<endl;

	calculateCovarianceAll(); //sets matrix sigmaAll
    cout<<"Calculated overall covariance."<<endl;


    for(int z = 0; z < numModelsM2; z++){
        pi[z] = 1.0/numModelsM2;
        piOld[z] = 1.0/numModelsM2;
    }

    //assign initial values to the transition probability matrix
    for(int z = 0; z < numModelsM2; z++){
        for(int zz = 0; zz < numModelsM2; zz++){
            transP[z][zz] = 1.0/numModelsM2;
            debugTransP[z][zz] = 1.0/numModelsM2;
        }
    }

    for(int z = 0; z < numModelsM2; z++){
        double sum = 0;
        for(int zz = 0; zz < numModelsM2; zz++){
            sum += transP[z][zz];
        }
        assert(sum < 1.0000001 && sum > 0.9999999);
    }

    for(int t = 0; t < ntimes; t++){
        for(int k = 0; k < 2*numProdsK; k++){
            qAll[k] = 0.0;
        }
    }

    for(int i = 0; i < numDyads; i++){
        for(int k = 0; k < 2*numProdsK; k++){
	  //HJ qAll[k] += d[i][0][k];
	  qAll[k] += d[i][eBegin[i]][k];
        }
    }

    for(int k = 0; k < 2*numProdsK; k++){
        qAll[k] /= numDyads;
        if(!(qAll[k] <= 1.0)){
            cout<<"qAll["<<k<<"] = "<<qAll[k]<<endl;
            cout.flush();
            exit(-1);
        }
    }

    for(int i = 0; i < numDyads; i++){
      for(int k = 0; k < 2*numProdsK; k++){
	//HJ qSigmaSqAll[k] += (qAll[k] - d[i][0][k])*(qAll[k] - d[i][0][k]);
	qSigmaSqAll[k] += (qAll[k] - d[i][eBegin[i]][k])*(qAll[k] - d[i][eBegin[i]][k]);
      }
    }

    for(int k = 0; k < 2*numProdsK; k++){
        qSigmaSqAll[k] /= (numDyads - 1);
    }

	for(int z = 0; z < numModelsM2; z++){
        //TODO: initialize q using empirical values
        for(int k = 0; k < 2*numProdsK; k++){
            q[z][k] = qAll[k];
            //q[z][k] = qTrue[z][k];
            if(!(q[z][k] <= 1.0)){
                cout<<"q["<<z<<"]["<<k<<"]="<<q[z][k]<<endl;
                cout.flush();
                exit(-1);
            }

        }
    }

	for(int z = 0; z < numModelsM; z++){
	   for(int k = 0; k < numProdsK2; k++){
          mu[z][k] = muAll[k] + uni();
          //mu[z][k] = muAll[k];

          sigmaSq[z][k] = sigmaSqAll[k];
		}
	}

    for(int z = 0; z < numModelsM; z++){
        for(int k = 0; k < numProdsK-1; k++){
            mu[z+numModelsM][k] = mu[z][k+numProdsK-1];
            sigmaSq[z+numModelsM][k] = sigmaSqAll[k+numProdsK-1];

            mu[z+numModelsM][k+numProdsK-1] = mu[z][k];
            sigmaSq[z+numModelsM][k+numProdsK-1] = sigmaSqAll[k];
        }
    }

    printf("OP_w_ij = %10.15Lf.\n",OP_w_ij);
    cout<<"Leaving initializeDebug."<<endl;

    return 1;
}

int ZeroTradeModel::initializeModel(time_t seed){
    vector<double> tmpVector;
    MapType::iterator it;
    VMapType::iterator vit;
    ZTMSeed = seed;
    generatorb.seed(static_cast<unsigned int>(ZTMSeed));
    boost::uniform_real<> uni_dist(-3.0,3.0);
    boost::uniform_real<> uni_dist2(0.0,1.0);
    boost::variate_generator<base_generator_type&,boost::uniform_real<> > uni(generatorb,uni_dist);
    boost::variate_generator<base_generator_type&,boost::uniform_real<> > uni2(generatorb,uni_dist2);

    //initialize x, w and y

#if DEBUG2
    cout<<"Started reading input file."<<endl;
#endif
    //TODO: radhika, uncomment the following
	readInput();

#if DEBUG2
    cout<<"Finished reading input file."<<endl;
#endif

#if DEBUG
    cout<<"Started reading description file."<<endl;
#endif
    //readCodeDescriptions("description.txt");
#if DEBUG
    cout<<"Finished reading description file."<<endl;
#endif

    //only save those dyads in vMap for which there is a reverse flow in the input
    //file
    int dbgCount = 0;
    int bidirCount = 0;
    vMap.clear();

    //Check that each of the dyads read-in from the input file has it's pair
    //information available too, that is they are all bi-directional (possibly
    //zero) trades. Count the total number of bidirectional trades and store in
    //bidir.
    for(vit = xMap.begin(); vit != xMap.end(); vit++){
        pair<int,int> mapKey(vit->first);
        vector<vector<long double> > &tmpVector = vit->second;
        pair<int,int> swapMapKey = make_pair(mapKey.second,mapKey.first);

        dbgCount++;

        VMapType::iterator vitPair = xMap.find(swapMapKey);

        if(vitPair == xMap.end()){
            //the reverse trade dyad doesn't exist in the dataset
            continue;
        }

        bidirCount++;

        vMap.insert(VMapType::value_type(mapKey,tmpVector));
    }

    //vMap now contains list of dyads, each of which has its reverse in vMap

/* HJ
#if 0
#if DEBUG2
    cout<<"Total number of dyads in inputFile = "<<dbgCount<<endl;
    cout<<"Number of dyads that are bidirectional = "<<bidirCount<<endl;
    cout<<"TO DO: reorganizeProducts()."<<endl;
#endif

    //find the most frequently traded product occurs all dyads
    //make it the (K-1)-th last product.
    //TODO: write out a productmap file that maps original product idx to
    //new idx after re-organization.
    reorganizeProducts();

#if DEBUG2
    cout<<"after call to reorganizeProducts()."<<endl;
#endif
#endif
*/

    w.resize(ntimes);
    tmpD.resize(ntimes);
    for(int t = 0; t < ntimes; t++){
        w[t].resize(bidirCount);
        tmpD[t].resize(bidirCount);
        for(int i = 0; i < bidirCount; i++){
            w[t][i].resize(numProdsK);  //the last value of w is to be ignored
            tmpD[t][i].resize(numProdsK);
        }
    }

    //assert that no. of dyad-pairs in vMap are equal to
    //the count of bidirectional dyadic trades
    assert(vMap.size() == bidirCount);

    int count = 0;
    dMap.clear();
    //initialize w[t][i][j][0....K-1] using log-proportion and C (for zero trade)
    //dMap maps the dyad index in w and tmpD to the {i,j} pair in vMap. 


/* HJ
    We use eMap to record the missing of (i, j, year) and (j, i, year).
    And use the dExist table to note the existance of a dyad for certain years.
    dExist = 1 means the dyad exists, and dExist = 0 means no.

    dExist.resize(bidirCount);
    for (int i = 0; i < bidirCount; i++) {
        dExist[i].resize(ntimes);
    } 
*/
    eMap.resize(ntimes);
    for(int t = 0; t < ntimes; t++) {
      eMap[t].clear();
    }


    for(vit = vMap.begin(); vit != vMap.end(); vit++){
	double tmpSum;
	pair<int,int> tmpPair(make_pair(vit->first.first,vit->first.second));
	vector<vector<long double> > &mapVal = vit->second;  //vector containing trade vols for K products 
	//dMap is the map for stacked dyad (i,j) -> dyadID, no time dim
	dMap.insert(MapType::value_type(tmpPair, count));
	//TODO: uncomment the following
	//assert(mapVal.size() == ntimes);
	for(int t = 0; t < ntimes; t++){
	  //TODO: delete the following where missing dyads are set to 0.0
	  if(mapVal.size() != ntimes)
	    {
	      for(int is = mapVal.size(); is <= ntimes; is++)
		{
		  vector<long double> tmpVector(numProdsK,0.0);
		  mapVal.push_back(tmpVector);
		  // HJ: Warn of missing dyad for some data set. Add 0.0 at the end may not be correct, other than prevent the code to fail.
		  cout << "Warning: dyad (" << vit->first.first << ", " << vit->first.second << ") has been patched at t = " << is << endl;
		}
	    }
	  assert(mapVal[t].size() == numProdsK);
	  
	  //mapVal[numProdsK - 1] += C;  //add a small constant to X_{ijK} 
	  

/* HJ: 
          In the extension, we assume that any negative tradeVolumn 
	  indicates the dyad does not exist for that year. Thus we have 
	  dExist[i][t] to keep track it by {1: exist; 0: non-exist}. 
	  First, we need to collect all the (i, j, year) for the missing
          dyad - which is stored in eMap[year].
*/

	  int nonExist = 0;
	  for(int i = 0; i < mapVal[t].size(); i++){
	    if (mapVal[t][i] < 0) {
	      eMap[t].insert(MapType::value_type(tmpPair, 1));
	      nonExist = 1;
	      break;
	    }
	  }

	  if (nonExist == 1) {
	    // To verify the encoding process
	    // cout << "************************" << endl;
	    // HJ: cout << "count = " << count << ", i = " << vit->first.first << ", j = " << vit->first.second << ", t = " << t << endl;
	    // cout << "************************" << endl;
	    tmpSum = 0.0;
/* HJ:
            Please note that if a particular (dyad, year) is not available, 
	    we set mapVal to be all 0.
*/
	    for (int i = 0; i < mapVal[t].size(); i++) {
	      mapVal[t][i] = 0;
	    }
	  } else {
	    tmpSum = accumulate(mapVal[t].begin(),mapVal[t].end(),0.0);
	    if (tmpSum > 0.0) {
	      for(int i = 0; i < mapVal[t].size(); i++){
		mapVal[t][i] = mapVal[t][i]/tmpSum;
	      } //end of loop over products = mapVal[t].size
	    } 
	  }

/* HJ: code being replaced with the above

		tmpSum = accumulate(mapVal[t].begin(),mapVal[t].end(),0.0);

                //vMap now contains, for each dyad i-j and product k the proportion of tradeVolume: X_ijk/sum(X_ijk)
                for(int i = 0; i < mapVal[t].size(); i++){
                    if(tmpSum > 0.0){
                        mapVal[t][i] = mapVal[t][i]/tmpSum;
                    }
                    else{
                        mapVal[t][i] = 0.0;
                    }
                } //end of loop over products = mapVal[t].size?
*/

                //note: w[count][numProdsK-1] is not going to be used in the
                //algorithm
	  for(int k = 0; k < numProdsK; k++){
                   if(mapVal[t][k] > 0.0){
                        tmpD[t][count][k] = 0;
                        //w[count][k] = log((mapVal[k])/(mapVal[numProdsK-1]));
                        if(mapVal[t][numProdsK-1] == 0){
                            w[t][count][k] = log((mapVal[t][k])/(mapVal[t][numProdsK-1] + C));
                        }
                        else{
                            w[t][count][k] = log((mapVal[t][k])/(mapVal[t][numProdsK-1]));
                        }
                    }
                    else{
                        //for all cases where tmpD[count][k] = 1.0, the w value is
                        //not important
                        tmpD[t][count][k] = 1;
                        w[t][count][k] = 0.0;
                    }
                } //end of loop over products
            } //end of for loop over times for this dyad
            count++;    //increment dyad count and process the next dyad in vMap
    } //end of loop over all dyads

    cout << "count = " << count << endl;

/*
#if 0 
    for(int t = 0; t < ntimes; t++){
        int tmpSum = 0;
        for(int i = 0; i < w[t].size(); i++){
            if(w[t][i][0] != 0.0){
                tmpSum++;
            }
        }
        cout<<"Total number of dyads with non-zero value of 0 for first product are: "<<tmpSum<<endl; 
    }
#endif
*/

   numDyads = 0;
   dyadMap.clear();
   y.clear();
   d.clear();
   dDash.clear();
   revDyadMap.clear();
   MapType::iterator dit;
#if DEBUG
   cout<<"size of dMap = "<<dMap.size()<<endl;
#endif


   //dMap is the map for stacked dyad (i,j) -> dyadID, no time dim
   //loop over all the dyads contained in dMap
   for(dit = dMap.begin(); dit != dMap.end(); dit++){
       pair<int,int> mapKey(dit->first);
       pair<int,int> swapMapKey = make_pair(mapKey.second,mapKey.first);
       pair<int,int> tmpPair;
       int wIdx = dit->second;
       MapType::iterator ditPair = dMap.find(swapMapKey);

       //check the pair of the dyad is in the dMap and equivalently is in the
       //vMap. This should be the case as we explicitly removed one-way trades. 
       assert(ditPair != dMap.end());

       int wPairIdx = ditPair->second;

	   (mapKey.first > mapKey.second) ? tmpPair = mapKey : tmpPair = swapMapKey;

       //does this pair already exist in the dyadMap, i.e. have we already
       //processed it?
       //MapType::iterator it = dyadMap.find(tmpPair);
       MapType::iterator it = dyadMap.find(mapKey);
       if(it != dyadMap.end()){
                continue;
       }
       MapType::iterator it2 = dyadMap.find(swapMapKey);
       if(it2 != dyadMap.end()){
                continue;
       }

	   dyadMap.insert(MapType::value_type(mapKey,numDyads));
       revDyadMap.insert(RMapType::value_type(numDyads,mapKey)); //TODO: only for debugging 

       //Stack the D's, W's and dDash's for all time t.
       //dDash is just a trunctaed version of D
       //with 2*(K-1) elements.
	   vector<vector<long double> > pairVector(ntimes,vector<long double>(numProdsK2,0.0));
	   vector<vector<int> > dDashVector(ntimes,vector<int>(numProdsK2,0.0));
	   vector<vector<int> > dVector(ntimes,vector<int>(2*numProdsK,0.0));

       for(int t = 0; t < ntimes; t++){
        for(int i = 0; i < numProdsK - 1; i++){
                pairVector[t][i] = w[t][wIdx][i];
                pairVector[t][i + numProdsK - 1] = w[t][wPairIdx][i];
                dDashVector[t][i] = tmpD[t][wIdx][i];
                dDashVector[t][i + numProdsK - 1] = tmpD[t][wPairIdx][i];
            }
       }

	    y.push_back(pairVector);
        dDash.push_back(dDashVector);

        for(int t = 0; t < ntimes; t++){
            for(int i = 0; i < numProdsK; i++){
                dVector[t][i] = tmpD[t][wIdx][i];
                dVector[t][i + numProdsK] = tmpD[t][wPairIdx][i];
            }
        }
        d.push_back(dVector);

        for(int t = 0; t < ntimes; t++){
            assert(y[numDyads][t].size() == numProdsK2);
            assert(dDash[numDyads][t].size() == numProdsK2);
            assert(d[numDyads][t].size() == 2.0*numProdsK);
        }

	    numDyads++;     //increment the number of dyad-pairs and process the next dyad-pair 
   } //loop over all the dyads in dMap
   assert(y.size() == numDyads);
   assert(dDash.size() == numDyads);
   assert(d.size() == numDyads);

   cout << "numDyads = " << numDyads << endl;

/*
    HJ: Initialize the dExist table, make all elements to be 1, meaning dyad exists for that year.
        Them we compare to the records in eMap, i.e. (i, j, year). 
	
	It is assumed that in the imput files, (i < j) should always hold, but to be safe, 
	we include checking (j, i, year) as well. However, limited items are in eMap anyway.
*/
    dExist.resize(numDyads);
    for (int i = 0; i < numDyads; i++) {
        dExist[i].resize(ntimes);
        for (int j = 0; j < ntimes; j++) {
	dExist[i][j] = 1;
        }
    } 

    for (int t = 0; t < ntimes; t++) {
      MapType::iterator eit;
      for (eit = eMap[t].begin(); eit != eMap[t].end(); eit++) {
	pair<int,int> mapKey(eit->first);
	pair<int,int> swapMapKey = make_pair(mapKey.second, mapKey.first);
	pair<int,int> tmpPair;
	int flag = eit->second;

	// HJ: Check the dyad import/export pair
	// cout << "mapKey = " << mapKey.first << ", " << mapKey.second << endl;
	
	MapType::iterator excluding = dyadMap.find(mapKey);
	if (excluding != dyadMap.end()) {
	  // cout << "Before: dExist[dyad][t] = " << dExist[excluding->second][t] << ";    ";
	  dExist[excluding->second][t] = 0;
	  // cout << "Excluding i = " << excluding->first.first << ", j = " << excluding->first.second;
	  // cout << ", dyad = " << excluding->second << ", t = " << t;
	  // cout << "   ==> dExist[dyad][t] = " << dExist[excluding->second][t] << endl;
	} else {
	  MapType::iterator excluding2 = dyadMap.find(swapMapKey); 
	  if (excluding2 !=  dyadMap.end()) {
	    // cout << "Before: dExist[dyad][t] = " << dExist[excluding2->second][t] << ";    ";
	    dExist[excluding2->second][t] = 0;
	    // cout << "Excluding2 i = " << excluding2->first.first << ", j = " << excluding2->first.second;
	    // cout << ", dyad = " << excluding2->second << ", t = " << t;
	    // cout << "   ==> dExist[dyad][t] = " << dExist[excluding2->second][t] << endl;
	  }
	}
      } // End of loop over all eMap[t] records to fill the 0s in the dExist table
    } // End of all eMap

    eBegin.resize(numDyads);
    eEnd.resize(numDyads);
    for (int i = 0; i < numDyads; i++) {
      eBegin[i] = 0;
      eEnd[i] = ntimes - 1;
      int t = 0;
      while (t < ntimes) {
	if (dExist[i][t] == 0) {
	  t++;
	  eBegin[i] = t;
	} else {
	  t = ntimes;
	}
      }
      t = ntimes - 1;
      while (t >= 0) {
	if (dExist[i][t] == 0) {
	  t--;
	  eEnd[i] = t;
	} else {
	  t = -1;
	}
      }
      if ((eEnd[i] != ntimes -1) | (eBegin[i] != 0)) {
	int impID = (revDyadMap.find(i)->second).first;
	int expID = (revDyadMap.find(i)->second).second;
	cout << "For dyad = (" << impID << ", " << expID << ")     eBegin = " << eBegin[i] << "   eEnd = " << eEnd[i] << endl;
      } 
    }	  

    //TODO: radhika uncomment the following

/* HJ
#if 0
#if FLIPDYADPAIRS
   fixZeroTradedProds();
#endif
#endif
*/

    //Do some more memory allocation once we know numDyads
    maxZeta.resize(numDyads); //used in analytics

    //TODO: change zeta dimensions to numModelsM2xnumDyads, E- and M- have
    //different access patterns but M- might be giving a bigger hit
    zeta.resize(numDyads);
    zeta1.resize(numDyads);
    ksi.resize(numDyads);
    //resp.resize(numDyads);
    wts.resize(numDyads);
    for(int i = 0; i < numDyads; i++){
        zeta[i].resize(ntimes);
        zeta1[i].resize(ntimes);
        //resp[i].resize(ntimes);
        ksi[i].resize(ntimes);
        wts[i].resize(ntimes);
        for(int t = 0; t < ntimes; t++){
            long double tmpDenominator = 0.0;
            zeta[i][t].resize(numModelsM2);
            zeta1[i][t].resize(numModelsM2);
            //resp[i][t].resize(numModelsM2);
            ksi[i][t].resize(numModelsM2);
            for(int z = 0; z < numModelsM2; z++){
                ksi[i][t][z].resize(numModelsM2);
            }
#if HJSKIP
	    // HJ: to by-pass the missing data
	    if ((eBegin[i] <= t) & (t <= eEnd[i])) { 
#endif
	      for(int k = 0; k < 2*numProdsK; k++){
                tmpDenominator += d[i][t][k];
	      }
	      wts[i][t] = 1.0/(exp(nu*tmpDenominator));
#if HJSKIP
	    } else {
	      wts[i][t] = 0.0;
	    }
#endif
	}
    }

#if DEBUG2
    cout<<"ksi[0][0][0][0] = "<<ksi[0][0][0][0]<<endl;
    cout<<"initializeModel: read in OP_w_ij = "<<OP_w_ij<<endl;
#endif
    OP_w_ij = OP_w_ij_inp/(2*numDyads);
    cout<<"Setting OP_w_ij = "<<OP_w_ij<<endl;
#if DEBUG2
    cout<<"initializeModel: using OP_w_ij = "<<OP_w_ij<<" OP_w_ij*2*numDyads = "<<OP_w_ij*2*numDyads<<endl;
#endif

	//calculate overall mean and covariance for the whole dataset
	calculateMeanAll(); //sets vector muAll
    //cout<<"Calculated overall mean."<<endl;

    //debug
#if DEBUG
    int tmpNonZeroDyads = 0;
    for(int i = 0; i < numDyads; i++){
        //if(y[i][0] != 0.0){
        if(dDash[i][148] != 1){
            tmpNonZeroDyads++;
        }
    }
    cout<<"Total number of dyad pairs with non-zero values for 148-th product "<<tmpNonZeroDyads<<endl; 
#endif

	calculateCovarianceAll(); //sets matrix sigmaAll
    //cout<<"Calculated overall covariance."<<endl;
    //mat sigmaSqMat;

    //assign initial values to pi, q, mu, simgmaSq
    for(int z = 0; z < numModelsM2; z++){
        pi[z] = 1.0/numModelsM2;
        piOld[z] = 1.0/numModelsM2;
        //TODO:use kmeans initialization of pi[z]
        //expectLogPi[t][z] = log(pi[t][z]);
    }

    //assign initial values to the transition probability matrix
    for(int z = 0; z < numModelsM2; z++){
        for(int zz = 0; zz < numModelsM2; zz++){
            transP[z][zz] = 1.0/numModelsM2;
        }
    }

    for(int z = 0; z < numModelsM2; z++){
        double sum = 0;
        for(int zz = 0; zz < numModelsM2; zz++){
            sum += transP[z][zz];
        }
        assert(sum < 1.0000001 && sum > 0.9999999);
    }

    //for(int t = 0; t < ntimes; t++){
        for(int k = 0; k < 2*numProdsK; k++){
            qAll[k] = 0.0;
        }
    //}

    for(int t = 0; t < ntimes; t++){
        for(int i = 0; i < numDyads; i++){
#if HJSKIP
          // HJ: to by-pass the missing data                                                                                                                   
          if ((eBegin[i] <= t) & (t <= eEnd[i])) {
#endif
	    for(int k = 0; k < 2*numProdsK; k++){
	      qAll[k] += d[i][t][k];
            }
#if HJSKIP
	  }
#endif
        }
    }

    //for(int t = 0; t < ntimes; t++){
        for(int k = 0; k < 2*numProdsK; k++){
            qAll[k] /= (ntimes*numDyads);
            //cout<<"qAll[k] = "<<qAll[k]<<endl;
            assert(qAll[k] <= 1.0);
        }
    //}

    //for(int t = 0; t < ntimes; t++){
        for(int z = 0; z < numModelsM2; z++){
            for(int k = 0; k < 2*numProdsK; k++){
                q[z][k] = qAll[k];
                //cout<<"Initial q["<<z<<"]["<<k<<"] = "<<q[z][k]<<endl;
                //assert(q[z][k] < 1.0 && q[z][k] > 0.0);
                //expectLogQ[t][z][k] = log(q[t][z][k]);
                //expectLogOneMinusQ[t][z][k] = log(1.0 - q[t][z][k]);
                //TODO: use k-means initialization of q's
                //q[z][k] = uni2();   //this causes empty clusters
                assert(qAll[k] <= 1.0);
            }
        }
    //}

    //for(int t = 0; t < ntimes; t++){
        for(int z = 0; z < numModelsM; z++){
            for(int k = 0; k < numProdsK2; k++){
                mu[z][k] = muAll[k] + uni();
                //mu[z][k] = uni();       //makes LL not decrease

                sigmaSq[z][k] = sigmaSqAll[k];
                //bigB[t][z][k] = sigmaSq[t][z][k];
            }
        }
    //}

    //for(int t = 0; t < ntimes; t++){
        for(int z = 0; z < numModelsM; z++){
            for(int k = 0; k < numProdsK-1; k++){
                mu[z+numModelsM][k] = mu[z][k+numProdsK-1];
                sigmaSq[z+numModelsM][k] = sigmaSqAll[k+numProdsK-1];
                //bigB[t][z+numModelsM][k] = bigB[t][z][k+numProdsK-1];

                mu[z+numModelsM][k+numProdsK-1] = mu[z][k];
                sigmaSq[z+numModelsM][k+numProdsK-1] = sigmaSqAll[k];
                //bigB[t][z+numModelsM][k+numProdsK-1] = bigB[t][z][k];
            }
        }
    //}

/*
#if 0
    for(int t = 0; t < ntimes; t++){
        for(int z = 0; z < numModelsM2; z++){
            for(int k = 0; k < 2*numProdsK; k++){
                bigN[t][z][k] = (1.0 - qAll[t][k]);
            }
        }
    }
#endif

#if 0
    //initialize R
    for(int t = 0; t < ntimes; t++){
        for(int i = 0; i < numDyads; i++){
            for(int z = 0; z < numModelsM2; z++){
                resp[t][i][z] = 1.0/numModelsM2;
            }
        }
    }
#endif

#if 0
    ofstream logFile;
    int isim=-99;
    logFile.open("tmp.txt");
    chkpObj.readModelParams(isim, logFile,
                                numProdsK,
                                numModelsM,
                                numDyads,
                                prodNames,
                                q,
                                mu,
                                sigmaSq,
                                zeta,
                                pi);
    logFile.close();
#endif
*/

#if DEBUG2
    cout<<"Returning from ZeroTradeModel::initializeModel."<<endl;
#endif

    return 1;
} //end of method initializeModel

/* HJ
#if 0
int ZeroTradeModel::initializeModel(time_t seed){
    vector<double> tmpVector;
    MapType::iterator it;
    VMapType::iterator vit;
    ZTMSeed = seed;
    generatorb.seed(static_cast<unsigned int>(ZTMSeed)); 
    boost::uniform_real<> uni_dist(-3.0,3.0);
    boost::uniform_real<> uni_dist2(0.0,1.0);
    boost::variate_generator<base_generator_type&,boost::uniform_real<> > uni(generatorb,uni_dist); 
    boost::variate_generator<base_generator_type&,boost::uniform_real<> > uni2(generatorb,uni_dist2); 
    
    //initialize x, w and y

#if DEBUG
    cout<<"Started reading input file."<<endl;
#endif
	readInput();
#if DEBUG
    cout<<"Finished reading input file."<<endl;
#endif

#if DEBUG
    cout<<"Started reading description file."<<endl;
#endif
    //readCodeDescriptions("description.txt");
#if DEBUG
    cout<<"Finished reading description file."<<endl;
#endif

    //only save those dyads in vMap for which there is a reverse flow in the input
    //file
    int dbgCount = 0;
    int bidirCount = 0;
    vMap.clear();

    //Check that each of the dyads read-in from the input file has it's pair
    //information available too, that is they are all bi-directional (possibly
    //zero) trades. Count the total number of bidirectional trades and store in
    //bidir. 
    for(vit = xMap.begin(); vit != xMap.end(); vit++){
        pair<int,int> mapKey(vit->first);
        vector<vector<long double> > &tmpVector = vit->second;
        pair<int,int> swapMapKey = make_pair(mapKey.second,mapKey.first);

        dbgCount++;

        VMapType::iterator vitPair = xMap.find(swapMapKey);

        if(vitPair == xMap.end()){
            //the reverse trade dyad doesn't exist in the dataset 
            continue;
        }

        bidirCount++; 

        vMap.insert(VMapType::value_type(mapKey,tmpVector));
    } 

    //vMap now contains list of dyads, each of which has its reverse in vMap
#if DEBUG
    cout<<"Total number of dyads in inputFile = "<<dbgCount<<endl;
    cout<<"Number of dyads that are bidirectional = "<<bidirCount<<endl; 
#endif

    //find the most frequently traded product occurs all dyads
    //make it the (K-1)-th last product.
    //TODO: write out a productmap file that maps original product idx to
    //new idx after re-organization.
    reorganizeProducts();

    w.resize(bidirCount);
    tmpD.resize(bidirCount);
    for(int i = 0; i < bidirCount; i++){
        for(int t = 0; t < ntimes; t++){
            w[i][t].resize(numProdsK);  //the last value of w is to be ignored
            tmpD[i][t].resize(numProdsK);
        }
    }

    assert(vMap.size() == bidirCount);

    int count = 0; 
    dMap.clear();
	//initialize w[i][j][0....K-1] using log-proportion and C (for zero trade)
    //dMap maps the dyad index in w and tmpD to the {i,j} pair in vMap. 
	for(vit = vMap.begin(); vit != vMap.end(); vit++){
	        double tmpSum;
            pair<int,int> tmpPair(make_pair(vit->first.first,vit->first.second));
	        vector<long double> &mapVal = vit->second;  //vector containing trade vols for K products 

	        assert(mapVal.size() == numProdsK);

            //mapVal[numProdsK - 1] += C;  //add a small constant to X_{ijK} 
            
	        tmpSum = accumulate(mapVal.begin(),mapVal.end(),0.0);

            //vMap now contains, for each dyad i-j and product k the proportion of tradeVolume: X_ijk/sum(X_ijk)
	        for(int i = 0; i < mapVal.size(); i++){
                if(tmpSum > 0.0){
	                mapVal[i] = mapVal[i]/tmpSum;
                }
                else{
	                mapVal[i] = 0.0;
                }
	        }

            dMap.insert(MapType::value_type(tmpPair, count));

            //note: w[count][numProdsK-1] is not going to be used in the
            //algorithm
	        for(int k = 0; k < numProdsK; k++){
               if(mapVal[k] > 0.0){
                   tmpD[count][k] = 0;
	               //w[count][k] = log((mapVal[k])/(mapVal[numProdsK-1]));
                   if(mapVal[numProdsK-1] == 0){
	                   w[count][k] = log((mapVal[k])/(mapVal[numProdsK-1] + C));
                   }
                   else{
	                   w[count][k] = log((mapVal[k])/(mapVal[numProdsK-1]));
                   }
               }
               else{
                   //for all cases where tmpD[count][k] = 1.0, the w value is
                   //not important 
                   tmpD[count][k] = 1;
                   w[count][k] = 0.0;
               }
	        }

            count++;    //increment dyad count and process the next dyad in vMap 
   }

  int tmpSum = 0;
  for(int i = 0; i < w.size(); i++){
      if(w[i][0] != 0.0){
          tmpSum++;
      }
  }
#if DEBUG
  cout<<"Total number of dyad with non-zero value of 0 for first product are: "<<tmpSum<<endl; 
#endif

   numDyads = 0;
   dyadMap.clear(); 
   y.clear();
   d.clear();
   dDash.clear();
   revDyadMap.clear();
   MapType::iterator dit;
#if DEBUG
   cout<<"size of dMap = "<<dMap.size()<<endl; 
#endif
   for(dit = dMap.begin(); dit != dMap.end(); dit++){
       pair<int,int> mapKey(dit->first);
	   pair<int,int> swapMapKey = make_pair(mapKey.second,mapKey.first);
       pair<int,int> tmpPair;
       int wIdx = dit->second;
       MapType::iterator ditPair = dMap.find(swapMapKey);

       //check the pair of the dyad is in the dMap and equivalently is in the
       //vMap. This should be the case as we explicitly removed one-way trades. 
       assert(ditPair != dMap.end()); 

       int wPairIdx = ditPair->second;

	   (mapKey.first > mapKey.second) ? tmpPair = mapKey : tmpPair = swapMapKey;

       //does this pair already exist in the dyadMap, i.e. have we already
       //processed it?
       //MapType::iterator it = dyadMap.find(tmpPair);
       MapType::iterator it = dyadMap.find(mapKey);
       if(it != dyadMap.end()){
                continue;
       }
       MapType::iterator it2 = dyadMap.find(swapMapKey);
       if(it2 != dyadMap.end()){
                continue;
       }

	   dyadMap.insert(MapType::value_type(mapKey,numDyads));
       revDyadMap.insert(RMapType::value_type(numDyads,mapKey)); //TODO: only for debugging 

       //Stack the D's, W's and dDash's.
       //dDash is just a trunctaed version of D
       //with 2*(K-1) elements.
	   vector<long double> pairVector(numProdsK2,0.0);
	   vector<int> dDashVector(numProdsK2,0.0);
	   vector<int> dVector(2*numProdsK,0.0);

	   for(int i = 0; i < numProdsK - 1; i++){
	        pairVector[i] = w[wIdx][i];
	        pairVector[i + numProdsK - 1] = w[wPairIdx][i];
            dDashVector[i] = tmpD[wIdx][i];
            dDashVector[i + numProdsK - 1] = tmpD[wPairIdx][i];
	    }
       
	    y.push_back(pairVector);
        dDash.push_back(dDashVector);

        for(int i = 0; i < numProdsK; i++){
            dVector[i] = tmpD[wIdx][i];
            dVector[i + numProdsK] = tmpD[wPairIdx][i];
        }
        d.push_back(dVector);

	    assert(y[numDyads].size() == numProdsK2);
	    assert(dDash[numDyads].size() == numProdsK2);
	    assert(d[numDyads].size() == 2.0*numProdsK);

	    numDyads++;     //increment the number of dyad-pairs and process the next dyad-pair 
    }
    assert(y.size() == numDyads);
    assert(dDash.size() == numDyads);
    assert(d.size() == numDyads); 


#if 0
#if DEBUG
    cout<<"numDyads = "<<numDyads<<endl; 
    ofstream outfile;
    outfile.open("ypair.txt");
    for(int i = 0; i < numDyads; i++){
        for(int k = 0; k < 2*(numProdsK-1); k++){
            outfile<<" "<<y[i][k]<<" ";
        }
        outfile<<endl;
    }
    outfile.close();
#endif
#endif

   //after stacking, if it turns out that some products have zero-trade
   //over all dyad-pairs then the following function will fix this
#if FLIPDYADPAIRS
   fixZeroTradedProds();    
#endif

    //Do some more memory allocation once we know numDyads
    maxZeta.resize(numDyads); //used in analytics 

    //TODO: change zeta dimensions to numModelsM2xnumDyads, E- and M- have
    //different access patterns but M- might be giving a bigger hit
    zeta.resize(numDyads);
    zeta1.resize(numDyads);
    r.resize(numDyads);
    wts.resize(numDyads);
    for(int i = 0; i < numDyads; i++){
        long double tmpDenominator = 0.0;

        zeta[i].resize(numModelsM2);
        zeta1[i].resize(numModelsM2);
        r[i].resize(numModelsM2);

        for(int k = 0; k < 2*numProdsK; k++){
            tmpDenominator += d[i][k];
        }
        wts[i] = 1.0/(exp(nu*tmpDenominator));

    }

    //Initialize constant for the PO(Prior Observations) method
    //OP_w_ij = 1.0/(2*numDyads);
    //TODO: radhika, change back ORIG: OP_w_ij = 0.00000001/(2*numDyads);
    //OP_w_ij = 0.00000001/(2*numDyads);
    //OP_w_ij = 50.0/(2*numDyads);
    //OP_w_ij = 0.0001/(2*numDyads);
    //OP_w_ij = 1.0/(2*numDyads);
    //cout<<"initializeModel: Using OP_w_ij = "<<OP_w_ij<<endl;
    //exit(-1);
#if DEBUG
    cout<<"initializeModel: read in OP_w_ij = "<<OP_w_ij<<endl;
#endif
    OP_w_ij = OP_w_ij_inp/(2*numDyads);
#if DEBUG
    cout<<"initializeModel: using OP_w_ij = "<<OP_w_ij<<" OP_w_ij*2*numDyads = "<<OP_w_ij*2*numDyads<<endl;
#endif

	//calculate overall mean and covariance for the whole dataset
	calculateMeanAll(); //sets vector muAll
    //cout<<"Calculated overall mean."<<endl;

    //debug
#if DEBUG
    int tmpNonZeroDyads = 0;
    for(int i = 0; i < numDyads; i++){
        //if(y[i][0] != 0.0){
        if(dDash[i][148] != 1){
            tmpNonZeroDyads++;
        }
    }
    cout<<"Total number of dyad pairs with non-zero values for 148-th product "<<tmpNonZeroDyads<<endl; 
#endif
	calculateCovarianceAll(); //sets matrix sigmaAll
    //cout<<"Calculated overall covariance."<<endl;
    //mat sigmaSqMat;


    //assign initial values to pi, q, mu, simgmaSq 
    for(int z = 0; z < numModelsM2; z++){
        pi[z] = 1.0/numModelsM2;
        //TODO:use kmeans initialization of pi[z]
        expectLogPi[z] = log(pi[z]);
    }

    for(int k = 0; k < 2*numProdsK; k++){
        qAll[k] = 0.0;
    }

    for(int i = 0; i < numDyads; i++){
        for(int k = 0; k < 2*numProdsK; k++){
            qAll[k] += d[i][k];
        }
    }

    for(int k = 0; k < 2*numProdsK; k++){
        qAll[k] /= numDyads;
        assert(qAll[k] <= 1.0);
    }

	for(int z = 0; z < numModelsM2; z++){
        for(int k = 0; k < 2*numProdsK; k++){
            q[z][k] = qAll[k];
            //assert(q[z][k] < 1.0 && q[z][k] > 0.0);
            expectLogQ[z][k] = log(q[z][k]);
            expectLogOneMinusQ[z][k] = log(1.0 - q[z][k]);
            //TODO: use k-means initialization of q's
            //q[z][k] = uni2();   //this causes empty clusters
            assert(qAll[k] <= 1.0);
        }
    }

	for(int z = 0; z < numModelsM; z++){
	   for(int k = 0; k < numProdsK2; k++){
          mu[z][k] = muAll[k] + uni();
          //mu[z][k] = uni();       //makes LL not decrease

          sigmaSq[z][k] = sigmaSqAll[k];
          bigB[z][k] = sigmaSq[z][k];
          //bigB[z][k] = 0.1;
       //   sigmaSq[z][k] = sigmaSqAll[k] + uni();  //makes likelihood not
       //   decrease;
		}
	}

    for(int z = 0; z < numModelsM; z++){
        for(int k = 0; k < numProdsK-1; k++){
            mu[z+numModelsM][k] = mu[z][k+numProdsK-1];
            sigmaSq[z+numModelsM][k] = sigmaSqAll[k+numProdsK-1];
            bigB[z+numModelsM][k] = bigB[z][k+numProdsK-1];

            mu[z+numModelsM][k+numProdsK-1] = mu[z][k];
            sigmaSq[z+numModelsM][k+numProdsK-1] = sigmaSqAll[k];
            bigB[z+numModelsM][k+numProdsK-1] = bigB[z][k];
        }
    }

    for(int z = 0; z < numModelsM2; z++){
        for(int k = 0; k < 2*numProdsK; k++){
            bigN[z][k] = (1.0 - qAll[k]);
        }
    }

    //initialize R
    for(int i = 0; i < numDyads; i++){
        for(int z = 0; z < numModelsM2; z++){
            r[i][z] = 1.0/numModelsM2;
        }
    }

#if 0
    ofstream logFile;
    int isim=-99;
    logFile.open("tmp.txt");
    chkpObj.readModelParams(isim, logFile,
                                numProdsK,
                                numModelsM,
                                numDyads,
                                prodNames,
                                q,
                                mu,
                                sigmaSq,
                                zeta,
                                pi);
    logFile.close();
#endif

#if DEBUG
    cout<<"Returning from ZeroTradeModel::initializeModel."<<endl;
#endif

    return 1;
}
#endif
*/

ZeroTradeModel::~ZeroTradeModel() {

}

/* HJ
#if 0
long double ZeroTradeModel::obsLogPosterior(int iter){
    llVal_prev = llVal;

    llVal = 0.0;
    //long double tmpLLVal = 0.0;

    vector<vector<vector<long double> > > lll(numDyads);
    vector<vector<vector<long double> > > maxk(numDyads);
    vector<vector<vector<long double> > > mink(numDyads);
    vector<vector<vector<vector<long double> > > > lllk(numDyads);
    for(int i = 0; i < numDyads; i++)
    {
        lll[i].resize(ntimes);
        lllk[i].resize(ntimes);
        maxk[i].resize(ntimes);
        mink[i].resize(ntimes);
        for(int t = 0; t < ntimes; t++)
        {
            lll[i][t].resize(numModelsM2);
            lllk[i][t].resize(numModelsM2);
            maxk[i][t].resize(numModelsM2);
            mink[i][t].resize(numModelsM2);
            for(int z = 0; z < numModelsM2; z++)
            {
                lllk[i][t][z].resize(numProdsK2);
            }
        }
    }

    long double maxNormFactor = -9999999999999999;
    for(int i = 0; i < numDyads; i++)
    {
        for(int t = 0; t < ntimes; t++)
        {
            for(int z = 0; z < numModelsM2; z++)
            {
                    //long double tmpDouble = localLogLikelihood(iter,i,t,z,lllk,mink,maxk);
                    //lll[i][t][z] = tmpDouble;
                    //if(maxNormFactor < tmpDouble) maxNormFactor = tmpDouble;
                    //llVal += OP_w_ij*tmpDouble;
            }// end of loop over z
        }//end of loop over t
    }//end of loop over i

#if 0 //DONT USE
    for(int i = 0; i < numDyads; i++)
    {
        for(int t = 0; t < ntimes; t++)
        {
            for(int z = 0; z < numModelsM2; z++)
            {
                lll[i][t][z] /= maxNormFactor;
                llVal += OP_w_ij*lll[i][t][z];
            }
        }
    }
#endif

    for(int i = 0; i < numDyads; i++)
    {
        for(int t = 0; t < ntimes; t++)
        {
            for(int z = 0; z < numModelsM2; z++)
            {
                lll[i][t][z] = 0.0;
                for(int k = 0; k < numProdsK2; k++)
                {
                    lll[i][t][z] += (lllk[i][t][z][k] - (mink[i][t][z] + maxk[i][t][z])/2.0);
#if 0
                    if(isinf(log(exp(lll[i][t][z]))))
                    {
                        cout<<"log(exp(lll[i][t][z])) is inf."<<endl;
                        cout<<"exp("<<lll[i][t][z]<<") = "<<exp(lll[i][t][z])<<endl;
                        cout<<"log(exp(lll[i][t][z])) = "<<log(exp(lll[i][t][z]))<<endl;
                        cout<<"lll["<<i<<"]["<<t<<"]["<<z<<"] = "<<lll[i][t][z]<<endl;
                        cout<<"lllk["<<i<<"]["<<t<<"]["<<z<<"]["<<k<<"] = "<<lllk[i][t][z][k]<<endl;
                        cout<<"maxk["<<i<<"]["<<t<<"]["<<z<<"] = "<<maxk[i][t][z]<<endl;
                        cout<<"mink["<<i<<"]["<<t<<"]["<<z<<"] = "<<mink[i][t][z]<<endl;
                        exit(-1);
                    }
#endif
                }
            }
        }
    }

    // long double maxZVal = -100000000000000;
    // long double maxZVal = 10000000000000;
    vector<vector<long double> > llz(numDyads,vector<long double>(numModelsM2));
    for(int i = 0; i < numDyads; i++)
    {
        for(int z = 0; z < numModelsM2; z++)
        {
            long double tmpDouble = log(pi[z]);
            tmpDouble += lll[i][0][z];
            tmpDouble += nextStep(i,1,z,lll);

            //dont use: if(maxZVal > tmpDouble) maxZVal = tmpDouble;

            llz[i][z] = tmpDouble;
#if 0
            if(isinf(log(exp(llz[i][z]))))
            {
                cout<<"encountered inf after nextStep."<<endl;
                cout<<"llz["<<i<<"]["<<z<<"] = "<<llz[i][z]<<endl;
                cout<<"exp("<<llz[i][z]<<") = "<<exp(llz[i][z])<<endl;
                cout<<"log(exp("<<llz[i][z]<<")) = "<<log(exp(llz[i][z]))<<endl;
                exit(-1);
            }
#endif
        }
    }

    long double tmpISum = 0.0;
    for(int i = 0; i < numDyads; i++)
    {
        long double maxZVal = -999999999999;
        long double tmpZSum = 0.0;

        for(int z = 0; z < numModelsM2; z++)
        {
            if(maxZVal < llz[i][z])
            {
                maxZVal = llz[i][z];
            }
        }

        for(int z = 0; z < numModelsM2; z++)
        {
            tmpZSum += exp(llz[i][z] - maxZVal);
            //tmpZSum += exp(llz[i][z]);
        }

        tmpISum += log(tmpZSum);
        if(isinf(tmpISum))
        {
            cout<<"Found inf tmpISum."<<endl;
            cout<<"tmpZSum = "<<tmpZSum<<endl;
            exit(-1);
        }
    }

    //cout<<"after ending recursion."<<endl;

    llVal += tmpISum;

    return llVal;
}
#endif
*/

/* HJ
#if 0
long double ZeroTradeModel::localLogLikelihood(int iter, int i, int t, int z,
                        vector<vector<vector<vector<long double> > > >& lllk,
                        vector<vector<vector<long double> > >& mink,
                        vector<vector<vector<long double> > >& maxk){
                    //long double tmpDouble = localLogLikelihood(iter,i,t,z,lllk,maxk);

        long double tmpDouble =0.0;
        mink[i][t][z] = 999999999999999999;
        maxk[i][t][z] = -999999999999999999;

        for(int k = 0; k < numProdsK-1; k++){
            long double localTmpDouble;
            if(dDash[i][t][k] == 1.0){
                //assert((double)q[z][k] != 0.0);
                if(q[z][k] > 0.0){
                    localTmpDouble = log(q[z][k]);
                    tmpDouble += localTmpDouble;;
                }
                else{
                    cout<<"LLL exit 1."<<endl;
                    exit(-1);
                }
            }
            else if(dDash[i][t][k] == 0.0){
                if((double)q[z][k] < 1.0){
                    localTmpDouble = log(1-q[z][k]);
                    localTmpDouble += (logDNorm(y[i][t][k],mu[z][k],sigmaSq[z][k]));
                    tmpDouble += localTmpDouble;
                }
                else{
                    cout<<"LLL exit 2."<<endl;
                    exit(-1);
                }
            }
            if(isnan(tmpDouble) || isnan(tmpDouble)){
                cout<<"tmpVal is nan: 1."<<endl;
                exit(-1);
            }
            lllk[i][t][z][k] = localTmpDouble;
            if(mink[i][t][z] > lllk[i][t][z][k])
            {
                mink[i][t][z] = lllk[i][t][z][k];
            }
            if(maxk[i][t][z] < lllk[i][t][z][k])
            {
                maxk[i][t][z] = lllk[i][t][z][k];
            }
        }//end of k loop over numProdsK

        for(int k = numProdsK-1; k < 2*(numProdsK-1); k++){
            long double localTmpDouble;
            if(dDash[i][t][k] == 1.0){
                if((double)q[z][k+1] > 0.0){
                    localTmpDouble = log(q[z][k+1]);
                    tmpDouble += localTmpDouble;
                }
                else{
                    cout<<"LLL exit 3."<<endl;
                    exit(-1);
                }
            }
            if(dDash[i][t][k] == 0.0){
                if((double)q[z][k+1] < 1.0){
                    localTmpDouble = log(1-q[z][k+1]);
                    localTmpDouble += (logDNorm(y[i][t][k],mu[z][k],sigmaSq[z][k]));
                    tmpDouble += localTmpDouble;
                }
                else{
                    cout<<"LLL exit 4."<<endl;
                    exit(-1);
                }
            }
            lllk[i][t][z][k] = localTmpDouble;
            if(mink[i][t][z] > lllk[i][t][z][k])
            {
                mink[i][t][z] = lllk[i][t][z][k];
            }
            else if(maxk[i][t][z] < lllk[i][t][z][k])
            {
                maxk[i][t][z] = lllk[i][t][z][k];
            }
        }//end of k loop from numProdsK to 2*numProdsK

        int kidx=numProdsK-1;
        if(d[i][t][kidx] == 1){
            assert(q[z][kidx] < 1.0000001);
            tmpDouble += log(q[z][kidx]);
            lllk[i][t][z][kidx] = log(q[z][kidx]);
        }
        else{ //else if (d[i][kidx] == 0){
            assert(q[z][kidx] < 1.0000001);
            tmpDouble += log(1-q[z][kidx]);
            lllk[i][t][z][kidx] = log(1-q[z][kidx]);
        }
        if(isnan(tmpDouble) || isinf(tmpDouble)){
            cout<<"tmpDouble is nan or inf: 3."<<endl;
            exit(-1);
        }
        if(mink[i][t][z] > lllk[i][t][z][kidx])
        {
            mink[i][t][z] = lllk[i][t][z][kidx];
        }
        if(maxk[i][t][z] < lllk[i][t][z][kidx])
        {
            maxk[i][t][z] = lllk[i][t][z][kidx];
        }

        kidx = 2*numProdsK-1;
        if(d[i][t][kidx] == 1){
            assert(q[z][kidx] < 1.0000001);
            tmpDouble += log(q[z][kidx]);
            lllk[i][t][z][kidx] = log(q[z][kidx]);
        }
        else{ //else if (d[i][kidx] == 0){
            assert(q[z][kidx] < 1.0000001);
            tmpDouble += log(1-q[z][kidx]);
            lllk[i][t][z][kidx] = log(1-q[z][kidx]);
        }
        if(isnan(tmpDouble) || isinf(tmpDouble)){
            cout<<"tmpDouble is nan: 4."<<endl;
            exit(-1);
        }
        if(mink[i][t][z] > lllk[i][t][z][kidx])
        {
            mink[i][t][z] = lllk[i][t][z][kidx];
        }
        if(maxk[i][t][z] < lllk[i][t][z][kidx])
        {
            maxk[i][t][z] = lllk[i][t][z][kidx];
        }

    return tmpDouble;
}
#endif
*/

long double ZeroTradeModel::nextStep(int idyad, int t, int z,
                        vector<vector<vector<long double> > >& lll)
{
    // cout<<"Entered nextStep with t = "<<t<<endl;
    long double llVal = 0.0;
    //if(t >= ntimes) return llVal;
    if(t >= 2) return llVal;

    for(int zdash = 0; zdash < numModelsM2; zdash++)
    {
        llVal += transP[z][zdash];
        llVal += lll[idyad][t][zdash];
        llVal += nextStep(idyad,t+1,zdash,lll);
    }
    // cout<<"Leaving nextStep with t = "<<t<<endl;
    return llVal;
}

/* HJ
#if 0
long double ZeroTradeModel::logLikelihoodTh(int iter){
    llVal_prev = llVal;

    llVal = 0.0;
    long double tmpLLVal = 0.0;

    for(int i = 0; i < numDyads; i++){
        for(int z = 0; z < numModelsM2; z++){
            llVal += (1 + OP_w_ij)*log(pi[z]);
        }
    }

    if(isnan(llVal)){
        cout<<"llval isnan after first term."<<endl;
        exit(-1);
    }

    //second term will require storing of zeta across all times
    for(int i = 0; i < numDyads; i++){
        for(int t = 1; t < ntimes; t++){
            for(int z = 0; z < numModelsM2; z++){
                for(int zdash = 0; zdash < numModelsM2; zdash++){
                        llVal += log(transP[zdash][z]);
                }
            }
        }
    }

    if(isnan(llVal)){
        cout<<"llval isnan after second term."<<endl;
        exit(-1);
    }

    //third term
    for(int i = 0; i < numDyads; i++){
        for(int t = 0; t < ntimes; t++){
            for(int k = 0; k < numProdsK-1; k++){
                long double tmpVal = 0;
                for(int z = 0; z < numModelsM2; z++){
                    long double tmpFactor;
                    long double tmpDouble;
                    bool zeroFlag = 0; 
                    if(zeta[i][t][z]){
                        tmpDouble = (1.0 + OP_w_ij);
                    }
                    else{
                        tmpDouble = OP_w_ij;
                    }
                    if(dDash[i][t][k] == 1.0){
                        //assert((double)q[z][k] != 0.0);
                        //tmpDouble += dDash[i][k]*log(q[z][k]);
                        if(q[z][k] > 0.0){
                            tmpDouble += log(q[z][k]);
                        }
                        else{
#ifdef PRIOR_OP
                            cout<<"LL exit 1."<<endl;
                            exit(-1);
#endif
                            tmpDouble = -10000;
                            zeroFlag = 1;
                            continue;
                        }
                    }
                    //TODO: do I need the q[z][k] != 1.0 condition???
                    if(dDash[i][t][k] == 0.0){
                        if((double)q[z][k] < 1.0){
                            tmpDouble += log(1-q[z][k]);
                            tmpDouble += (logDNorm(y[i][t][k],mu[z][k],sigmaSq[z][k]));
                        }
                        else{
#ifdef PRIOR_OP
                            cout<<"LL exit 2."<<endl;
                            exit(-1);
#endif
                            //tmpDouble -= 10000;
                            tmpDouble = -10000;
                            zeroFlag = 1;
                            continue;
                        }
                    }
                    if(zeroFlag){
                        tmpVal += 0.0;
                        zeroFlag = 0;
                    }
                    else{
                        tmpVal += tmpFactor*tmpDouble;
                    }
                }//end of z loop over numModelsM2
                //assert((double)tmpVal > 0.0);
                llVal += tmpVal;
                if(isnan(llVal)){
                    cout<<"llVal is nan: 1."<<endl;
                    exit(-1);
                }
            }//end of k loop over numProdsK

            for(int k = numProdsK-1; k < 2*(numProdsK-1); k++){
                long double tmpVal = 0;
                for(int z = 0; z < numModelsM2; z++){
                    long double tmpFactor;
                    long double tmpDouble;
                    if(zeta[i][t][z]){
                        tmpDouble = (1.0 + OP_w_ij);
                    }
                    else{
                        tmpDouble = OP_w_ij;
                    }
                    //do some processing for third term
                    bool zeroFlag = 0;
                    if(dDash[i][t][k] == 1.0){
                        if((double)q[z][k+1] > 0.0){
                            tmpDouble += log(q[z][k+1]);
                        }
                        else{
#ifdef PRIOR_OP
                            cout<<"LL exit 3."<<endl;
                            exit(-1);
#endif
                            //tmpDouble -= 10000;
                            tmpDouble = -10000;
                            zeroFlag = 1;
                            continue;
                        }
                    }
                    //TODO: do I need the q[z][k] != 1.0 condition???
                    if(dDash[i][t][k] == 0.0){
                        if((double)q[z][k+1] < 1.0){
                            tmpDouble += log(1-q[z][k+1]);
                            tmpDouble += (logDNorm(y[i][t][k],mu[z][k],sigmaSq[z][k]));
                        }
                        else{
#ifdef PRIOR_OP
                            cout<<"LL exit 4."<<endl;
                            exit(-1);
#endif
                            //tmpDouble -= 10000;
                            tmpDouble = -10000;
                            zeroFlag = 1;
                            continue;
                        }
                    }
                    if(zeroFlag){
                        tmpVal += 0.0;
                        zeroFlag = 0;
                    }
                    else{
                        tmpVal += tmpFactor*tmpDouble;
                    }
                }//end of loop over z
                llVal += tmpVal;
                if(isnan(llVal)){
                    cout<<"llVal is nan: 2."<<endl;
                    exit(-1);
                }
            }//end of k loop from numProdsK to 2*numProdsK

            int kidx=numProdsK-1;
            long double tmpVal = 0.0;
            for(int z = 0; z < numModelsM2; z++){
                long double tmpFactor;
                long double tmpDouble=0.0;
                if(zeta[i][t][z]){
                   tmpFactor = 1.0 + OP_w_ij;
                }
                else{
                    tmpFactor = OP_w_ij;
                }
                if(d[i][t][kidx] == 1){
                    assert(q[z][kidx] < 1.0000001);
                    tmpDouble += log(q[z][kidx]);
                }
                else{ //else if (d[i][kidx] == 0){
                    assert(q[z][kidx] < 1.0000001);
                    tmpDouble += log(1-q[z][kidx]);
                }
                tmpVal += tmpFactor*tmpDouble;
            }
            llVal += tmpVal;
            if(isnan(llVal)){
                cout<<"llVal is nan: 3."<<endl;
                exit(-1);
            }

        kidx = 2*numProdsK-1; 
        tmpVal = 0.0;
        for(int z = 0; z < numModelsM2; z++){
            long double tmpFactor;
            long double tmpDouble=0.0;
            if(zeta[i][t][z]){
                tmpFactor = 1.0 + OP_w_ij;
            }
            else{
                tmpFactor = OP_w_ij;
            }
            if(d[i][t][kidx] == 1){
                //tmpDouble += d[i][kidx]*log(q[z][kidx]);
                assert(q[z][kidx] < 1.0000001);
                tmpDouble += log(q[z][kidx]);
            }
            else{ //else if (d[i][kidx] == 0){
            //TODO: check if we ever reach here
                //tmpDouble += (1-d[i][kidx])*log(1-q[z][kidx]);
                assert(q[z][kidx] < 1.0000001);
                tmpDouble += log(1-q[z][kidx]);
            }
            tmpVal += tmpFactor*tmpDouble;
        }
        llVal += tmpVal;
        if(isnan(llVal)){
            cout<<"llVal is nan: 4."<<endl;
            exit(-1);
        }

        }//end of t loop over ntimes
    }//end of i loop over numDyads

    if(isnan(llVal)){
        cout<<"llval isnan after third term."<<endl;
        exit(-1);
    }

    return llVal;
}
#endif
*/

/* HJ
#if 0
//Compute log-likelihood which is the stopping criteria and BIC
long double ZeroTradeModel::logLikelihoodTh(int iter){
    int df = 0; 
    long double tmpLLVal = 0.0;
    long double tmpLLVal2 = 0.0; 

    llVal_prev_prev = llVal_prev;
    llVal_prev = llVal;

    bicVal_prev = bicVal;

    llVal = 0.0;
    bicVal = 0.0;


tmpLLVal = 0.0; 
#pragma omp parallel for reduction(+:tmpLLVal)
	for(int i = 0; i < numDyads; i++){
        for(int k = 0; k < numProdsK-1; k++){
	        long double tmpVal = 0;
			for(int z = 0; z < numModelsM2; z++){
                bool zeroFlag = 0;
                long double tmpDouble = log(pi[z]);
                if(dDash[i][k] == 1.0){
                    //assert((double)q[z][k] != 0.0);
                    //tmpDouble += dDash[i][k]*log(q[z][k]);
                    if(q[z][k] > 0.0){
                        tmpDouble += wts[i]*log(q[z][k]);
                    }
                    else{
#ifdef PRIOR_OP
                        cout<<"LL exit 1."<<endl;
                        exit(-1);
#endif
                        //tmpDouble -= 10000;
                        tmpDouble = -10000;
                        zeroFlag = 1;
                        continue;
                    }
                }
                //TODO: do I need the q[z][k] != 1.0 condition???
                if(dDash[i][k] == 0.0){
                    if((double)q[z][k] < 1.0){
                        tmpDouble += wts[i]*log(1-q[z][k]);
                        //tmpDouble += log(1 - q[z][k]);
                        tmpDouble += (logDNorm(y[i][k],mu[z][k],sigmaSq[z][k]));
                    }
                    else{
#ifdef PRIOR_OP
                        cout<<"LL exit 2."<<endl;
                        exit(-1);
#endif
                        //tmpDouble -= 10000;
                        tmpDouble = -10000;
                        zeroFlag = 1;
                        continue;
                    }
                }
                if(zeroFlag){
                    tmpVal += 0.0;
                    zeroFlag = 0; 
                }
                else{
                    tmpVal += exp(tmpDouble);
                }
            }
            //assert((double)tmpVal > 0.0);
            tmpLLVal += log(tmpVal);
        }

        for(int k = numProdsK-1; k < 2*(numProdsK-1); k++){
            long double tmpVal = 0;
			for(int z = 0; z < numModelsM2; z++){
                long double tmpDouble = log(pi[z]);
                bool zeroFlag = 0; 
                if(dDash[i][k] == 1.0){
                    if((double)q[z][k+1] > 0.0){
                        tmpDouble += wts[i]*log(q[z][k+1]);
                    }
                    else{
#ifdef PRIOR_OP
                        cout<<"LL exit 3."<<endl;
                        exit(-1);
#endif
                        //tmpDouble -= 10000;
                        tmpDouble = -10000;
                        zeroFlag = 1;
                        continue;
                    }
                }
                //TODO: do I need the q[z][k] != 1.0 condition???
                //if(dDash[i][k] == 0 && q[z][k+1] != 1.0){
                if(dDash[i][k] == 0.0){
                    if((double)q[z][k+1] < 1.0){
                        tmpDouble += wts[i]*log(1-q[z][k+1]);
                        tmpDouble += (logDNorm(y[i][k],mu[z][k],sigmaSq[z][k]));
                    }
                    else{
#ifdef PRIOR_OP
                        cout<<"LL exit 4."<<endl;
                        exit(-1);
#endif
                        //tmpDouble -= 10000;
                        tmpDouble = -10000;
                        zeroFlag = 1;
                        continue;
                    }
                }
                if(zeroFlag){
                    tmpVal += 0.0;
                    zeroFlag = 0; 
                }
                else{
                    tmpVal += exp(tmpDouble);
                }
            }
             tmpLLVal += log(tmpVal);
        }
        
        int kidx=numProdsK-1;
        long double tmpVal = 0.0;
        for(int z = 0; z < numModelsM2; z++){
            long double tmpDouble = log(pi[z]);
            if(d[i][kidx] == 1){
                //tmpDouble += d[i][kidx]*log(q[z][kidx]);
                assert(q[z][kidx] != 0.0);
                tmpDouble += wts[i]*log(q[z][kidx]);
            }
            else{ //else if (d[i][kidx] == 0){
            //TODO: check if we ever reach here
                //tmpDouble += (1-d[i][kidx])*log(1-q[z][kidx]);
                assert(q[z][kidx] != 1.0);
                tmpDouble += wts[i]*log(1-q[z][kidx]);
            }
            tmpVal += exp(tmpDouble);
        }
        tmpLLVal += log(tmpVal);

        kidx = 2*numProdsK-1; 
        tmpVal = 0.0;
        for(int z = 0; z < numModelsM2; z++){
            long double tmpDouble = log(pi[z]);
            if(d[i][kidx] == 1){
                //tmpDouble += d[i][kidx]*log(q[z][kidx]);
                assert(q[z][kidx] != 0.0);
                tmpDouble += wts[i]*log(q[z][kidx]);
            }
            else{ //else if (d[i][kidx] == 0){
            //TODO: check if we ever reach here
                //tmpDouble += (1-d[i][kidx])*log(1-q[z][kidx]);
                assert(q[z][kidx] != 1.0);
                tmpDouble += wts[i]*log(1-q[z][kidx]);
            }
            tmpVal += exp(tmpDouble);
        }
        tmpLLVal += log(tmpVal);

   } //end of first loop over dyads 

   //cout<<"llVal before addign prior contribution : "<<tmpLLVal<<endl;
//TODO: Thread parallelize

tmpLLVal2 = 0.0; 
#pragma omp parallel for reduction(+:tmpLLVal2)
for(int i = 0; i < numDyads; i++){
   for(int k = 0; k < (numProdsK - 1); k++){
       //TODO: should this be up to M or 2N
       for(int z = 0; z < numModelsM; z++){
           if(dDash[i][k] == 1){
               tmpLLVal2 += 2.0*OP_w_ij*log(q[z][k]);
           }
           else{ //if d[i][k] == 0 
               tmpLLVal2 += 2.0*OP_w_ij*log(1-q[z][k]);
               tmpLLVal2 += 2.0*OP_w_ij*logDNorm(y[i][k],mu[z][k],sigmaSq[z][k]);
           }
                
       }
   }
   for(int k = numProdsK-1; k < 2*(numProdsK - 1); k++){
        for(int z = 0; z < numModelsM; z++){
           if(dDash[i][k] == 1){
               tmpLLVal2 += 2.0*OP_w_ij*log(q[z][k+1]);
           }
           else{ //if d[i][k] == 0 
               tmpLLVal2 += 2.0*OP_w_ij*log(1-q[z][k+1]);
               tmpLLVal2 += 2.0*OP_w_ij*logDNorm(y[i][k],mu[z][k],sigmaSq[z][k]);
           }
       }
   }

   int kidx = numProdsK-1;
   for(int z = 0; z < numModelsM; z++){
        if(d[i][kidx] == 1){
           tmpLLVal2 += 2.0*OP_w_ij*log(q[z][kidx]);
        }
        else{
           tmpLLVal2 += 2.0*OP_w_ij*log(1-q[z][kidx]);
        }
   }

   kidx = 2*numProdsK-1;
   for(int z = 0; z < numModelsM; z++){
        if(d[i][kidx] == 1){
           tmpLLVal2 += 2.0*OP_w_ij*log(q[z][kidx]);
        }
        else{
           tmpLLVal2 += 2.0*OP_w_ij*log(1-q[z][kidx]);
        }
   }
} //end of second loop over dyads

   llVal = tmpLLVal + tmpLLVal2; 

   if(iter < 2){
       acc_llVal = -9999999999;
   }
   else{
       acc_llVal = (llVal - llVal_prev)/(llVal_prev - llVal_prev_prev);
   }

   //cout<<"llVal after addign prior contribution : "<<tmpLLVal<<endl;

    // Compute BIC
    df = numModelsM*2*numProdsK;     //number of q estimated 
    df += numModelsM*2*(numProdsK-1); //number of mu estimated
    df += numModelsM*2*(numProdsK-1); //number of sigmaSq estimated
    df += numModelsM; //number of pi values 

    bicVal = 2.0*llVal - df*log(numDyads);

    cout<<"df = "<<df<<endl;
    cout<<"numDyads = "<<numDyads<<" log(numDyads) = "<<log(numDyads)<<endl;
    cout<<"llVal = "<<llVal<<endl;
    cout<<"bicVal = "<<bicVal<<endl;

    return llVal;
}
#endif
*/

/* HJ
#if 0
void ZeroTradeModel::calcBigN(int iter, ofstream& logFile){

//calc values for z = 0 to numModelsM - 1
    for(int z = 0; z < numModelsM; z++){
        for(int k = 0; k < numProdsK; k++){
            double localSumKFirstHalf = 0.0;
            double localSumKSecondHalf = 0.0;
            for(int i = 0; i < numDyads; i++){
                for(int t = 0; t < ntimes; t++){
                    localSumKFirstHalf += (zeta[i][t][z]+OP_w_ij)*(1.0 - d[i][t][k]) + 
                                                (zeta[i][t][z+numModelsM]+OP_w_ij)*(1.0 - d[i][t][k+numProdsK]);
                    localSumKSecondHalf += (zeta[i][t][z]+OP_w_ij)*(1.0 - d[i][t][k+numProdsK]) +
                                                (zeta[i][t][z+numModelsM]+OP_w_ij)*(1.0 - d[i][t][k]);
                }
            }
            bigN[z][k] = localSumKFirstHalf;
            bigN[z][k+numProdsK] = localSumKSecondHalf;

        }
    }

    //set values for z = numModelsM to 2*numModelsM - 1
    for(int z = 0; z < numModelsM; z++){
        for(int k = 0; k < numProdsK; k++){
            bigN[z+numModelsM][k] = bigN[z][k+numProdsK];
            bigN[z+numModelsM][k+numProdsK] = bigN[z][k];
        }
    }
}
#endif
*/

long double ZeroTradeModel::calcDiGamma(long double val){
    if(val == 0.0 || val < 0.0){
        cout<<"calcDiGamma: val = "<<val<<endl;
        cout.flush();
    }
	return boost::math::digamma<long double>(val);
}

//Equation (12) of ZeroTradeModel.pdf
int ZeroTradeModel::expectationTh(int iter,ofstream &logFile){

#if DEBUG2
    cout<<"Entering expectationTh()."<<endl;
#endif

    vector<vector<vector<long double> > > alpha;
    vector<vector<vector<long double> > > alphaN; //normalized alpha
    vector<vector<vector<long double> > > beta;
    vector<vector<vector<long double> > > betaN; //normalized beta
    vector<vector<long double> > coeff; //coefficient used for normalization

    alpha.resize(numDyads);
    beta.resize(numDyads);
    alphaN.resize(numDyads);
    betaN.resize(numDyads);
    coeff.resize(numDyads);
    for(int i = 0; i < numDyads; i++){
        alpha[i].resize(ntimes);
        beta[i].resize(ntimes);
        alphaN[i].resize(ntimes);
        betaN[i].resize(ntimes);
        coeff[i].resize(ntimes,0.0);
        for(int t = 0; t < ntimes; t++){
            alpha[i][t].resize(numModelsM2);
            beta[i][t].resize(numModelsM2);
            alphaN[i][t].resize(numModelsM2);
            betaN[i][t].resize(numModelsM2);
        }
    }

    vector<vector<long double> > maxNormFactor(numDyads, vector<long double>(ntimes,-999999999999999));

#if DEBUG2
    cout<<"After allocating local memory."<<endl;
#endif

  // HJ: I thought ceiling already rounded up, but it is not the case !!!

  // First to make sure each block has about 100 dyads
  int numBlocks = numDyads/100+1;
  // Then to reduce it to the number of threads
  if (numBlocks > numThreads) numBlocks = numThreads;
  // Then to calculate the size of each block 
  int blockSize = ceil(numDyads/numBlocks) + 1;
  // If the round up error causes problems, adjust it
  while (blockSize * (numBlocks - 1) >= numDyads) numBlocks--;

  cout << "Total of " << numDyads << " dyads, and they are splitted into " << numBlocks << " blocks of " << blockSize << " dyads each."<< endl;
  for (int section = 0; section < numBlocks; section++) {
    cout << "Block " << section << ": from dyads " << section*blockSize << " to " << std::min(section*blockSize+blockSize, numDyads) << endl;
  }

#pragma omp parallel for
  for (int section = 0; section < numBlocks; section++) {
    // cout << "Working on section " << section << ": from dyads " << section*blockSize << " to " << std::min(section*blockSize+blockSize, numDyads) << endl;
    for (int i = section*blockSize; i < std::min(section*blockSize+blockSize, numDyads); i++) {
    //for(int i = 0; i < numDyads; i++){
        for(int z = 0; z < numModelsM2; z++){
            //don't use log scale for beta just yet

/* HJ: replace t=ntimes-1 by eEnd[i], and t=0 by eBegin[i]
*/
            //HJ beta[i][ntimes-1][z] = 1.0;
            //HJ betaN[i][ntimes-1][z] = 1.0; //TODO: check

            beta[i][eEnd[i]][z] = 1.0;                                                                                                                     
            betaN[i][eEnd[i]][z] = 1.0; //TODO: check                                                                                                     
 
            double tmpDouble = 0.0;

            for(int k = 0; k < (numProdsK-1); k++){
	        //HJ if(dDash[i][0][k] == 1){
                if(dDash[i][eBegin[i]][k] == 1){
                    if(q[z][k] > 0.0){
                        //tmpDouble += wts[i][0]*log(q[z][k]);
                        tmpDouble += log(q[z][k]);
                        //cout<<"log(q[z][k]) = "<<log(q[z][k])<<endl;
                    }
                    else{
                        cout<<"weird expectation 1."<<endl;
                        exit(-1);
                    }
                }
                //HJ if(dDash[i][0][k] == 0){
                if(dDash[i][eBegin[i]][k] == 0){
                    if(q[z][k] < 1.0){
                        //tmpDouble += wts[i][0]*log(1 - q[z][k]);
                        tmpDouble += log(1 - q[z][k]);
                        //HJ tmpDouble += logDNorm(y[i][0][k],mu[z][k],sigmaSq[z][k]);
                        tmpDouble += logDNorm(y[i][eBegin[i]][k],mu[z][k],sigmaSq[z][k]);
                        //cout<<"log(1-q[z][k]) = "<<log(1-q[z][k])<<endl;
                        //cout<<"logDNorm = "<<logDNorm(y[i][0][k],mu[z][k],sigmaSq[z][k])<<endl;
                    }
                    else{
                        cout<<"weird expectation 2."<<endl;
                        exit(-1);
                    }
                }
            }

            for(int k = numProdsK-1; k < 2*(numProdsK-1); k++){
	        //HJ if(dDash[i][0][k] == 1){
                if(dDash[i][eBegin[i]][k] == 1){
                    if(q[z][k+1] > 0.0){
                        //tmpDouble += wts[i][0]*log(q[z][k+1]);
                        tmpDouble += log(q[z][k+1]);
                        //cout<<"log(q[z][k+1] = "<<log(q[z][k+1])<<endl;
                    }
                    else{
                        cout<<"weird expectation 3."<<endl;
                        exit(-1);
                    }
                }
                //HJ if(dDash[i][0][k] == 0){
                if(dDash[i][eBegin[i]][k] == 0){
                    if(q[z][k+1] < 1.0){
                        //tmpDouble += wts[i][0]*log(1 - q[z][k+1]);
                        tmpDouble += log(1 - q[z][k+1]);
                        //HJ tmpDouble += logDNorm(y[i][0][k],mu[z][k],sigmaSq[z][k]);
                        tmpDouble += logDNorm(y[i][eBegin[i]][k],mu[z][k],sigmaSq[z][k]);
                        //cout<<"log(1-q[z][k+1]) = "<<log(1-q[z][k+1])<<endl;
                        //cout<<"logDNorm = "<<logDNorm(y[i][0][k],mu[z][k],sigmaSq[z][k])<<endl;
                    }
                    else{
                        cout<<"weird expectation 4."<<endl;
                        exit(-1);
                    }
                }
            }

            int kidx=numProdsK-1;
            //HJ if(d[i][0][kidx]){
            if(d[i][eBegin[i]][kidx]){
                assert(q[z][kidx] != 0.0); //this should neer happen:
                //tmpDouble += wts[i][0]*log(q[z][kidx]);
                tmpDouble += log(q[z][kidx]);
                //cout<<"log(q[z][kidx] = "<<log(q[z][kidx])<<endl;
            }
            else{
                assert(q[z][kidx] != 1.0);  //this should never happen
                //tmpDouble += wts[i][0]*log(1 - q[z][kidx]);
                tmpDouble += log(1 - q[z][kidx]);
                //cout<<"log(1-q[z][kidx] = "<<log(1 - q[z][kidx])<<endl;
            }

            kidx = 2*numProdsK - 1;
            //HJ if(d[i][0][kidx]){
            if(d[i][eBegin[i]][kidx]){
                assert(q[z][kidx] != 0.0); //this should never happen
                //tmpDouble += wts[i][0]*log(q[z][kidx]);
                tmpDouble += log(q[z][kidx]);
                //cout<<"log(q[z][kidx] = "<<log(q[z][kidx])<<endl;
            }
            else{
                assert(q[z][kidx] != 1.0); //this should never happen 
                //tmpDouble += wts[i][0]*log(1 - q[z][kidx]);
                tmpDouble += log(1 - q[z][kidx]);
                //cout<<"log(1-q[z][kidx] = "<<log(1-q[z][kidx])<<endl;
            }

	    /* HJ: Mar23 revision ends here */

            //alpha[i][0][z] = pi[0][z]*exp(tmpDouble);
            //HJ if(maxNormFactor[i][0] < tmpDouble) maxNormFactor[i][0] = tmpDouble;
            if(maxNormFactor[i][eBegin[i]] < tmpDouble) maxNormFactor[i][eBegin[i]] = tmpDouble;

            //HJ alpha[i][0][z] = tmpDouble;
            //HJ alphaN[i][0][z] = tmpDouble;
            alpha[i][eBegin[i]][z] = tmpDouble;
            alphaN[i][eBegin[i]][z] = tmpDouble;

        } //end of loop over numModelsM2

        for(int z = 0; z < numModelsM2; z++){
	    //HJ alpha[i][0][z] = log(pi[z]) + (alpha[i][0][z] - maxNormFactor[i][0]);
            //HJ alphaN[i][0][z] = log(pi[z]) + (alphaN[i][0][z] - maxNormFactor[i][eBegin[i]]);
            alpha[i][eBegin[i]][z] = log(pi[z]) + (alpha[i][eBegin[i]][z] - maxNormFactor[i][eBegin[i]]);
            alphaN[i][eBegin[i]][z] = log(pi[z]) + (alphaN[i][eBegin[i]][z] - maxNormFactor[i][eBegin[i]]);
            //alpha[i][0][z] = log(pi[z]) + alpha[i][0][z];
            //alphaN[i][0][z] = log(pi[z]) + alphaN[i][0][z];
            //HJ coeff[i][0] += exp(alpha[i][0][z]);
            coeff[i][eBegin[i]] += exp(alpha[i][eBegin[i]][z]);
        }

        //HJ coeff[i][0] = log(coeff[i][0]);
        coeff[i][eBegin[i]] = log(coeff[i][eBegin[i]]);

        for(int z = 0; z < numModelsM2; z++){
	   //HJ alphaN[i][0][z] = alphaN[i][0][z] - coeff[i][0];
            alphaN[i][eBegin[i]][z] = alphaN[i][eBegin[i]][z] - coeff[i][eBegin[i]];
        }
    } //end of loop over dyads
  }
  
/* HJ
#if 0 //DEBUG
    for(int i = 0; i < numDyads; i++){
        for(int t = 0; t < ntimes; t++){
            for(int z = 0; z < numModelsM2; z++){
                cout<<"alpha["<<i<<"]["<<t<<"]["<<z<<"] = "<<alpha[i][t][z]<<endl;
            }
        }
    }
#endif
*/

/* HJ
#if 0 //DEBUG
    for(int z = 0; z < numModelsM2; z++){
        int tmpSum = 0.0;
        for(int i = 0; i < numDyads; i++){
            tmpSum += alpha[i][0][z];
        }
        cout<<"pi[0][z] = "<<pi[0][z]<<" sigma_i alpha[i][0][z] = "<<tmpSum<<endl;
    }
#endif
*/

#if DEBUG2
    cout<<"After first loop over dyads."<<endl;
#endif

#pragma omp parallel for
  for (int section = 0; section < numBlocks; section++) {
    // cout << "Working on section " << section << ": from dyads " << section*blockSize << " to " << std::min(section*blockSize+blockSize, numDyads) << endl;
    for (int i = section*blockSize; i < std::min(section*blockSize+blockSize, numDyads); i++) {
      // for(int i = 0; i < numDyads; i++){
      for (int t = eBegin[i]+1; t < eEnd[i]+1; t++) {
	  // for(int t = 1; t < ntimes; t++){

	  for(int z = 0; z < numModelsM2; z++){
                //long double tmpSum = 0.0;
                //long double tmpSumN = 0.0;
                long double tmpLogDouble = 0.0;
		
/* HJ
#if 0
                for(int zdash = 0; zdash < numModelsM2; zdash++){
                        tmpSum += exp(alpha[i][t-1][zdash])*transP[zdash][z];
                }
                for(int zdash = 0; zdash < numModelsM2; zdash++){
                        tmpSumN += exp(alphaN[i][t-1][zdash])*transP[zdash][z];
                }
#endif
*/

                for(int k = 0; k < (numProdsK-1); k++){
                    if(dDash[i][t][k] == 1){
                        if(q[z][k] > 0.0){
                            //tmpLogDouble += wts[i][t]*log(q[z][k]);
                            tmpLogDouble += log(q[z][k]);
                        }
                        else{
                            cout<<"weird expectation 5."<<endl;
                            exit(-1);
                        }
                    }
                    if(dDash[i][t][k] == 0){
                        if(q[z][k] < 1.0){
                            //tmpLogDouble += wts[i][t]*log(1 - q[z][k]);
                            tmpLogDouble += log(1 - q[z][k]);
                            tmpLogDouble += logDNorm(y[i][t][k],mu[z][k],sigmaSq[z][k]);
                        }
                        else{
                            cout<<"weird expectation 6."<<endl;
                            exit(-1);
                        }
                    }
                }

                for(int k = numProdsK-1; k < 2*(numProdsK-1); k++){
                    if(dDash[i][t][k] == 1){
                        if(q[z][k+1] > 0.0){
                            //tmpLogDouble += wts[i][t]*log(q[z][k+1]);
                            tmpLogDouble += log(q[z][k+1]);
                        }
                        else{
                            cout<<"weird expectation 7."<<endl;
                            exit(-1);
                        }
                    }
                    if(dDash[i][t][k] == 0){
                        if(q[z][k+1] < 1.0){
                            //tmpLogDouble += wts[i][t]*log(1 - q[z][k+1]);
                            tmpLogDouble += log(1 - q[z][k+1]);
                            tmpLogDouble += logDNorm(y[i][t][k],mu[z][k],sigmaSq[z][k]);
                        }
                        else{
                            cout<<"weird expectation 8."<<endl;
                            exit(-1);
                        }
                    }
                }

                int kidx=numProdsK-1;
                if(d[i][t][kidx]){
                    assert(q[z][kidx] != 0.0); //this should neer happen:
                    //tmpLogDouble += wts[i][t]*log(q[z][kidx]);
                    tmpLogDouble += log(q[z][kidx]);
                }
                else{
                    assert(q[z][kidx] != 1.0);  //this should never happen
                    //tmpLogDouble += wts[i][t]*log(1 - q[z][kidx]);
                    tmpLogDouble += log(1 - q[z][kidx]);
                }

                kidx = 2*numProdsK - 1;
                if(d[i][t][kidx]){
                    assert(q[z][kidx] != 0.0); //this should never happen
                    //tmpLogDouble += wts[i][t]*log(q[z][kidx]);
                    tmpLogDouble += log(q[z][kidx]);
                }
                else{
                    assert(q[z][kidx] != 1.0); //this should never happen
                    //tmpLogDouble += wts[i][t]*log(1 - q[z][kidx]);
                    tmpLogDouble += log(1 - q[z][kidx]);
                }

                if(maxNormFactor[i][t] < tmpLogDouble) maxNormFactor[i][t] = tmpLogDouble;

                alpha[i][t][z] = tmpLogDouble;
                alphaN[i][t][z] = tmpLogDouble;
            } //end of loop over z

            for(int z = 0; z < numModelsM2; z++)
            {
                long double tmpSum = 0.0;
                long double tmpSumN = 0.0;
                for(int zdash = 0; zdash < numModelsM2; zdash++){
                        tmpSum += exp(alpha[i][t-1][zdash])*transP[zdash][z];
                }
                for(int zdash = 0; zdash < numModelsM2; zdash++){
                        tmpSumN += exp(alphaN[i][t-1][zdash])*transP[zdash][z];
                }

                alpha[i][t][z] = (alpha[i][t][z] - maxNormFactor[i][t]) + log(tmpSum);
                alphaN[i][t][z] = (alphaN[i][t][z] - maxNormFactor[i][t]) + log(tmpSumN);
                //alpha[i][t][z] = alpha[i][t][z] + log(tmpSum);
                //alphaN[i][t][z] = alphaN[i][t][z] + log(tmpSumN);

                coeff[i][t] += exp(alpha[i][t][z]);
                if(coeff[i][t] < 0.0){
                    cout<<"coeff["<<i<<"]["<<t<<"] = "<<coeff[i][t]<<endl;
                    exit(-1);
                }
            } //end of z loop over numModelsM2

            if(isnan(log(coeff[i][t])) || isinf((coeff[i][t])))
            {
                cout<<"Found nan/inf coeff[i][t] 999."<<endl;
                cout<<"coeff["<<i<<"]["<<t<<"] = "<<coeff[i][t]<<endl;
                exit(-1);
            }
            coeff[i][t] = log(coeff[i][t]);

            // complete calculation of coeff[i][t]
            long double tmpDouble = 0.0;
            for(int tdash = 0; tdash < t - 1; tdash++){
                tmpDouble += coeff[i][tdash];
                if(isnan(tmpDouble)){
                    cout<<"Found nan tmpDouble."<<endl;
                    cout<<"tmpDouble = "<<tmpDouble<<endl;
                    cout<<"coeff["<<i<<"]["<<tdash<<"] = "<<coeff[i][tdash]<<endl;
                    exit(-1);
                }
            }

            coeff[i][t] = coeff[i][t] - tmpDouble;

            if(isnan(coeff[i][t]) || isinf(coeff[i][t])){
                cout<<"Found nan/inf coeff[i][t]."<<endl;
                cout<<"coeff["<<i<<"]["<<t<<"] = "<<coeff[i][t]<<endl;
                cout<<"tmpDouble = "<<tmpDouble<<endl;
                exit(-1);
            }

            // complete calculation of alphaN
            for(int z = 0 ; z < numModelsM2; z++){
                alphaN[i][t][z] -= coeff[i][t];
            } //end of second z loop over numModelsM2

        } //end of t loop over ntimes
    } //end of i loop over dyads
  }

#if DEBUG2
    cout<<"After second loop over dyads."<<endl;
#endif

    //set the rest of beta values

#pragma omp parallel for
  for (int section = 0; section < numBlocks; section++) {
    // cout << "Working on section " << section << ": from dyads " << section*blockSize << " to " << std::min(section*blockSize+blockSize, numDyads) << endl;
    for (int i = section*blockSize; i < std::min(section*blockSize+blockSize, numDyads); i++) {
    // for(int i = 0; i < numDyads; i++){
      for(int t = eEnd[i] - 1; t >= eBegin[i]; t--){
      // for(int t = ntimes - 2; t >= 0; t--){

            for(int z = 0; z < numModelsM2; z++){
                beta[i][t][z] = 0.0;
                betaN[i][t][z] = 0.0;
                for(int zdash = 0; zdash < numModelsM2; zdash++){
                    long double tmpSum = log(beta[i][t+1][zdash]) + log(transP[z][zdash]);
                    long double tmpSumN = log(betaN[i][t+1][zdash]) + log(transP[z][zdash]);
                    long double tmpLogDouble = 0.0;
		    
/* HJ
#if 0
                    if(isnan(tmpSum)){
                        cout<<"i = "<<i<<" t+1 = "<<t+1<<" zdash = "<<z<<endl;
                        cout<<"beta[i][t+1][zdash] = "<<beta[i][t+1][zdash]<<endl;
                        cout<<"z = "<<z<<" zdash = "<<zdash<<endl;
                        cout<<"transP[z][zdash] = "<<transP[z][zdash]<<endl;
                        exit(-1);
                    }

                    if(tmpSumN == 0){
                        cout<<"tmpSumN = "<<tmpSumN<<endl;
                        cout<<"i = "<<i<<endl;
                        cout<<"t + 1 = "<<t+1<<endl;
                        for(int zdash = 0; zdash < numModelsM2; zdash++){
                            cout<<"exp(betaN["<<i<<"]["<<t+1<<"]["<<zdash<<"]) = "<<exp(betaN[i][t+1][zdash])<<endl;
                            cout<<"transP["<<z<<"]["<<zdash<<"] = "<<transP[z][zdash]<<endl;
                        }
                        cout<<"log(tmpSumN) = "<<log(tmpSumN)<<endl;
                        exit(-1);
                    }
#endif
*/

                    for(int k = 0; k < (numProdsK-1); k++){
                        if(dDash[i][t+1][k] == 1){
                            if(q[zdash][k] > 0.0){
                                //tmpLogDouble += wts[i][t+1]*log(q[zdash][k]);
                                tmpLogDouble += log(q[zdash][k]);
                            }
                            else{
                                cout<<"weird expectation 9."<<endl;
                                exit(-1);
                            }
                        }
                        if(dDash[i][t+1][k] == 0){
                            if(q[zdash][k] < 1.0){
                                //tmpLogDouble += wts[i][t+1]*log(1 - q[zdash][k]);
                                tmpLogDouble += log(1 - q[zdash][k]);
                                tmpLogDouble += logDNorm(y[i][t+1][k],mu[zdash][k],sigmaSq[zdash][k]);
                            }
                            else{
                                cout<<"weird expectation 10."<<endl;
                                exit(-1);
                            }
                        }
                    }

                    for(int k = numProdsK-1; k < 2*(numProdsK-1); k++){
                        if(dDash[i][t+1][k] == 1){
                            if(q[zdash][k+1] > 0.0){
                                //tmpLogDouble += wts[i][t+1]*log(q[zdash][k+1]);
                                tmpLogDouble += log(q[zdash][k+1]);
                            }
                            else{
                                cout<<"weird expectation 11."<<endl;
                                exit(-1);
                            }
                        }
                        if(dDash[i][t+1][k] == 0){
                            if(q[zdash][k+1] < 1.0){
                                //tmpLogDouble += wts[i][t+1]*log(1 - q[zdash][k+1]);
                                tmpLogDouble += log(1 - q[zdash][k+1]);
                                tmpLogDouble += logDNorm(y[i][t+1][k],mu[zdash][k],sigmaSq[zdash][k]);
                            }
                            else{
                                cout<<"weird expectation 12."<<endl;
                                exit(-1);
                            }
                        }
                    }

                    int kidx=numProdsK-1;
                    if(d[i][t+1][kidx]){
                        assert(q[zdash][kidx] != 0.0); //this should neer happen:
                        //tmpLogDouble += wts[i][t+1]*log(q[zdash][kidx]);
                        tmpLogDouble += log(q[zdash][kidx]);
                    }
                    else{
                        assert(q[zdash][kidx] != 1.0);  //this should never happen
                        //tmpLogDouble += wts[i][t+1]*log(1 - q[zdash][kidx]);
                        tmpLogDouble += log(1 - q[zdash][kidx]);
                    }

                    kidx = 2*numProdsK - 1;
                    if(d[i][t+1][kidx]){
                        assert(q[zdash][kidx] != 0.0); //this should never happen
                        //tmpLogDouble += wts[i][t+1]*log(q[zdash][kidx]);
                        tmpLogDouble += log(q[zdash][kidx]);
                    }
                    else{
                        assert(q[zdash][kidx] != 1.0); //this should never happen
                        //tmpLogDouble += wts[i][t+1]*log(1 - q[zdash][kidx]);
                        tmpLogDouble += log(1 - q[zdash][kidx]);
                    }
                    //WRONG: beta[i][t][z] += beta[i][t+1][zdash]*transP[z][zdash]*exp(tmpLogDouble);
                    //if(tmpSum != 0.0){
                    beta[i][t][z] += exp(tmpSum + (tmpLogDouble - maxNormFactor[i][t+1]));
                    betaN[i][t][z] += exp(tmpSumN + (tmpLogDouble - maxNormFactor[i][t+1]));
                    //beta[i][t][z] += exp(tmpSum + tmpLogDouble);
                    //betaN[i][t][z] += exp(tmpSumN + tmpLogDouble);

                } //end of zdash loop over numModelsM2

                //complete normalization
                betaN[i][t][z] /= exp(coeff[i][t+1]);

                if(isnan(betaN[i][t][z]) || isinf(betaN[i][t][z])){
                    cout<<"betaN["<<i<<"]["<<t<<"]["<<z<<"] = "<<betaN[i][t][z]<<endl;
                    //cout<<"tmpSumN = "<<tmpSumN<<endl;
                    //cout<<"tmpLogDouble = "<<tmpLogDouble<<endl;
                    cout<<"coeff["<<i<<"]["<<t+1<<"] = "<<coeff[i][t+1]<<endl;
            //        exit(-1);
                }
            } //end of z loop over numModelsM2
	    
/* HJ
#if 0
            bool flag = 0;
            for(int z = 0; z < numModelsM2; z++){
                if(exp(betaN[i][t][z]) == 0.0){
                    flag = 1;
                }
            }

            if(flag){
                for(int z = 0; z < numModelsM2; z++){
                    cout<<"betaN["<<i<<"]["<<t<<"]["<<z<<"] = "<<betaN[i][t][z]<<endl;
                }
                exit(-1);
            }
#endif
*/

        } //end of t loop over ntimes
    } // end of i loop over numDyads
  }


        //convert beta to log scale
#pragma omp parallel for
  for (int section = 0; section < numBlocks; section++) {
    // cout << "Working on section " << section << ": from dyads " << section*blockSize << " to " << std::min(section*blockSize+blockSize, numDyads) << endl;
        for (int i = section*blockSize; i < std::min(section*blockSize+blockSize, numDyads); i++) {
	// for(int i = 0; i < numDyads; i++){
	  for (int t = eBegin[i]; t < eEnd[i]+1; t++) {
	  // for(int t = 0; t < ntimes; t++){
	      for(int z = 0; z < numModelsM2; z++){
                    beta[i][t][z] = log(beta[i][t][z]);
                    betaN[i][t][z] = log(betaN[i][t][z]);
                }
            }
        }
  }

#if DEBUG2
    cout<<"After third loop over dyads."<<endl;
#endif

    vector<vector<long double> > maxZeta;
    maxZeta.resize(numDyads);
    for(int i = 0; i < numDyads; i++){
        maxZeta[i].resize(ntimes,-99999999999.0);
    }

#pragma omp parallel for
  for (int section = 0; section < numBlocks; section++) {
    // cout << "Working on section " << section << ": from dyads " << section*blockSize << " to " << std::min(section*blockSize+blockSize, numDyads) << endl;
    for (int i = section*blockSize; i < std::min(section*blockSize+blockSize, numDyads); i++) {
    // for(int i = 0; i < numDyads; i++){
      for (int t = eBegin[i]; t < eEnd[i]+1; t++) {
      // for(int t = 0; t < ntimes; t++){
	  for(int z = 0; z < numModelsM2; z++){
                if(isnan(alpha[i][t][z]) || isnan(beta[i][t][z])){
                        cout<<"found nan."<<endl;
                        cout<<"alpha[i][t][z] = "<<alpha[i][t][z]<<endl;
                        cout<<"beta[i][t][z] = "<<beta[i][t][z]<<endl;
                        exit(-1);
                }
                //zeta[i][t][z] = alpha[i][t][z] + beta[i][t][z];
                zeta[i][t][z] = alphaN[i][t][z] + betaN[i][t][z];
                if(isnan(zeta[i][t][z]) || isinf(zeta[i][t][z])){
                    cout<<"Found zeta nan."<<endl;
                    cout<<"zeta["<<i<<"]["<<t<<"]["<<z<<"] = "<<zeta[i][t][z]<<endl;
                    cout<<"alpha["<<i<<"]["<<t<<"]["<<z<<"] = "<<alpha[i][t][z]<<endl;
                    cout<<"beta["<<i<<"]["<<t<<"]["<<z<<"] = "<<beta[i][t][z]<<endl;
                    exit(-1);
                }
                if(maxZeta[i][t] < zeta[i][t][z]){
                        maxZeta[i][t] = zeta[i][t][z];
                }
                //tmpSum[t][i] += zeta[i][t][z];
            }
        }
    }
  }


    // for(int t = 0; t < ntimes; t++){
#pragma omp parallel for
  for (int section = 0; section < numBlocks; section++) {
    // cout << "Working on section " << section << ": from dyads " << section*blockSize << " to " << std::min(section*blockSize+blockSize, numDyads) << endl;
    for (int i = section*blockSize; i < std::min(section*blockSize+blockSize, numDyads); i++) {
    // for(int i = 0; i < numDyads; i++){
        for(int t = eBegin[i]; t < eEnd[i]+1; t++){
            for(int z = 0; z < numModelsM2; z++){
                //zeta[i][t][z] = exp(zeta[i][t][z] - maxZeta[i][t]);
                zeta[i][t][z] = exp(zeta[i][t][z]);
            }
        }
    }
  }

    vector<vector<long double> > dyadSum(numDyads,vector<long double>(ntimes,0.0));
    
/* HJ
#if 0
    for(int i = 0; i < numDyads; i++){
        for(int t = 0; t < ntimes; t++){
            for(int z = 0; z < numModelsM2; z++){
                dyadSum[i][t] += zeta[i][t][z];
            }
        }
    }
#endif
*/

//#if 0
#pragma omp parallel for
  for (int section = 0; section < numBlocks; section++) {
    // cout << "Working on section " << section << ": from dyads " << section*blockSize << " to " << std::min(section*blockSize+blockSize, numDyads) << endl;
    for (int i = section*blockSize; i < std::min(section*blockSize+blockSize, numDyads); i++) {
    // for(int i = 0; i < numDyads; i++){
        for(int t = eBegin[i]; t < eEnd[i]+1; t++){
	// for(int t = 0; t < ntimes; t++){
	  for(int z = 0; z < numModelsM2; z++){
	    //HJ dyadSum[i][t] += zeta[i][ntimes - 1][z];
	    dyadSum[i][t] += zeta[i][eEnd[i]][z];
	  }
        }
    }
  }
//#endif

#pragma omp parallel for
    for(int z = 0; z < numModelsM2; z++){
      // for(int t = 0; t < ntimes; t++){
      for(int i = 0; i < numDyads; i++){
	for (int t = eBegin[i]; t < eEnd[i]+1; t++) {
	  
/* HJ
#if 0
                    //zeta[i][t][z] /= tmpSum[i][t];
                    //zeta[i][t][z] /= tmpSum[i][t];
                    long double tmpDenominator = 0.0;
                    for(int zdash = 0; zdash < numModelsM2; zdash++){
                        tmpDenominator += zeta[i][t][zdash];
                        //tmpDenominator += exp(alpha[i][ntimes-1][zdash]);
                        //tmpDenominator += exp(alpha[i][t][zdash] + beta[i][t][zdash]);
                    }
                    //zeta[i][t][z] -= log(tmpDenominator);
                    if(tmpDenominator == 0){
                        cout<<"tmpDenominator = "<<tmpDenominator<<endl;
                        cout<<"i = "<<i<<" t = "<<t<<" z = "<<z<<endl;
                    }
                    zeta[i][t][z] /= tmpDenominator;
#endif
*/

/* HJ
#if 0
                    if(dyadSum[i][t] == 0.0){
                        for(int zdash = 0; zdash < numModelsM2; zdash++){
                            cout<<"zeta["<<i<<"]["<<t<<"]["<<zdash<<"] = "<<zeta[i][t][zdash]<<endl;
                            cout<<"alpha["<<i<<"]["<<t<<"]["<<zdash<<"] = "<<alpha[i][t][zdash]<<endl;
                            cout<<"beta["<<i<<"]["<<t<<"]["<<zdash<<"] = "<<beta[i][t][zdash]<<endl;
                        }
                    }
#endif
*/

                    zeta[i][t][z] /= dyadSum[i][t];
                    if(isnan(zeta[i][t][z])){
                            cout<<"Discovered nans."<<endl;
                            cout<<"zeta["<<i<<"]["<<t<<"]["<<z<<"] = "<<zeta[i][t][z]<<endl;
                            cout<<"alpha["<<i<<"]["<<t<<"]["<<z<<"] = "<<alpha[i][t][z]<<endl;
                            cout<<" beta["<<i<<"]["<<t<<"]["<<z<<"] = "<<beta[i][t][z]<<endl;
                            cout<<" dyadSum["<<i<<"]["<<t<<"] = "<<dyadSum[i][t]<<endl;
                            exit(-1);
                    }
            }
        }
    }

#pragma omp parallel for
  for (int section = 0; section < numBlocks; section++) {
    // cout << "Working on section " << section << ": from dyads " << section*blockSize << " to " << std::min(section*blockSize+blockSize, numDyads) << endl;
    for (int i = section*blockSize; i < std::min(section*blockSize+blockSize, numDyads); i++) {
    // for(int i = 0; i < numDyads; i++){
        for(int t = eBegin[i]; t < eEnd[i]+1; t++) {
	// for(int t = 0; t < ntimes; t++){
            if(accumulate(zeta[i][t].begin(),zeta[i][t].end(),0.0) > 1.0000001){
               cout<<"greater than 1.0: "<<accumulate(zeta[i][t].begin(),zeta[i][t].end(),0.0)<<endl;
            }
            else if(accumulate(zeta[i][t].begin(),zeta[i][t].end(),0.0) < 0.9999999){
               cout<<"greater than 1.0: "<<accumulate(zeta[i][t].begin(),zeta[i][t].end(),0.0)<<endl;
            }

            assert(accumulate(zeta[i][t].begin(),zeta[i][t].end(),0.0) < 1.0000001 && 
                        accumulate(zeta[i][t].begin(),zeta[i][t].end(),0.0) > 0.9999999);
        }
    }
  }

#if DEBUG2
    cout<<"After fourth loop over dyads."<<endl;
#endif

    vector<vector<vector<long double> > >maxKsi;

    maxKsi.resize(numDyads);
#pragma omp parallel for
  for (int section = 0; section < numBlocks; section++) {
    // cout << "Working on section " << section << ": from dyads " << section*blockSize << " to " << std::min(section*blockSize+blockSize, numDyads) << endl;
    for (int i = section*blockSize; i < std::min(section*blockSize+blockSize, numDyads); i++) {
    // for(int i = 0; i < numDyads; i++){
        maxKsi[i].resize(ntimes);
        for(int t = 0; t < ntimes; t++){
            maxKsi[i][t].resize(numModelsM2,-999999999.0);
        }
    }
  }
  
/* HJ
#if 0 //not used
    //set ksi for t = 0
    for(int i = 0; i < numDyads; i++){
        //for(int t = 0; t < ntimes; t++){
            for(int zdash = 0; zdash < numModelsM2; zdash++){
                for(int z = 0; z < numModelsM2; z++){
                    //ksi[i][0][zdash][z] = 0.001;
                    int t = 0;
                    //orig: ksi[i][t][zdash][z] = log(transP[zdash][z]) + beta[i][t][z];
                    //ksi[i][t][zdash][z] = log(transP[zdash][z]) + beta[i][t][z];
                    ksi[i][t][zdash][z] = log(transP[zdash][z]) + betaN[i][t][z];
                    long double tmpLogDouble = 0.0;
                    for(int k = 0; k < (numProdsK-1); k++){
                        if(dDash[i][t][k] == 1){
                            if(q[z][k] > 0.0){
                                //tmpLogDouble += wts[i][t]*log(q[z][k]);
                                tmpLogDouble += log(q[z][k]);
                            }
                            else{
                                cout<<"weird expectation."<<endl;
                                exit(-1);
                            }
                        }
                        if(dDash[i][t][k] == 0){
                            if(q[z][k] < 1.0){
                                //tmpLogDouble += wts[i][t]*log(1 - q[z][k]);
                                tmpLogDouble += log(1 - q[z][k]);
                                tmpLogDouble += logDNorm(y[i][t][k],mu[z][k],sigmaSq[z][k]);
                            }
                            else{
                                cout<<"weird expectation."<<endl;
                                exit(-1);
                            }
                        }
                    }

                    for(int k = numProdsK-1; k < 2*(numProdsK-1); k++){
                        if(dDash[i][t][k] == 1){
                            if(q[z][k+1] > 0.0){
                                //tmpLogDouble += wts[i][t]*log(q[z][k+1]);
                                tmpLogDouble += log(q[z][k+1]);
                            }
                            else{
                                cout<<"weird expectation."<<endl;
                                exit(-1);
                            }
                        }
                        if(dDash[i][t][k] == 0){
                            if(q[z][k+1] < 1.0){
                                //tmpLogDouble += wts[i][t]*log(1 - q[z][k+1]);
                                tmpLogDouble += log(1 - q[z][k+1]);
                                tmpLogDouble += logDNorm(y[i][t][k],mu[z][k],sigmaSq[z][k]);
                            }
                            else{
                                cout<<"weird expectation."<<endl;
                                exit(-1);
                            }
                        }
                    }

                    int kidx=numProdsK-1;
                    if(d[i][t][kidx]){
                        assert(q[z][kidx] != 0.0); //this should neer happen:
                        //tmpLogDouble += wts[i][t]*log(q[z][kidx]);
                        tmpLogDouble += log(q[z][kidx]);
                    }
                    else{
                        assert(q[z][kidx] != 1.0);  //this should never happen
                        //tmpLogDouble += wts[i][t]*log(1 - q[z][kidx]);
                        tmpLogDouble += log(1 - q[z][kidx]);
                    }

                    kidx = 2*numProdsK - 1;
                    if(d[i][t][kidx]){
                        assert(q[z][kidx] != 0.0); //this should never happen
                        //tmpLogDouble += wts[i][t]*log(q[z][kidx]);
                        tmpLogDouble += log(q[z][kidx]);
                    }
                    else{
                        assert(q[z][kidx] != 1.0); //this should never happen
                        //tmpLogDouble += wts[i][t]*log(1 - q[z][kidx]);
                        tmpLogDouble += log(1 - q[z][kidx]);
                    }

                    ksi[i][t][zdash][z] += tmpLogDouble;
                    //ksi[i][t][zdash][z] += coeff[i][t];
                    //if(maxKsi[i][t][zdash] < ksi[i][t][zdash][z]){
                    //    maxKsi[i][t][zdash] = ksi[i][t][zdash][z];
                    //}
                    if(isnan(ksi[i][t][zdash][z])){
                        printf("ksi[%d][%d][%d][%d] = %Lf.\n",i,t,zdash,z,ksi[i][t][zdash][z]);
                        printf("log(transP[%d][%d] = %Lf.\n",zdash,z,transP[zdash][z]);
                        printf("beta[%d][%d][%d] = %Lf.\n",i,t,z,beta[i][t][z]);
                        printf("%Lf.\n",tmpLogDouble);
                        exit(-1);
                    }
                }
            }
        //}
    }
#endif
*/

#pragma omp parallel for
  for (int section = 0; section < numBlocks; section++) {
    // cout << "Working on section " << section << ": from dyads " << section*blockSize << " to " << std::min(section*blockSize+blockSize, numDyads) << endl;
    for (int i = section*blockSize; i < std::min(section*blockSize+blockSize, numDyads); i++) {
    // for(int i = 0; i < numDyads; i++){
      
/* HJ
#if 0//this is wrong, most probably
       long double tmpDenominator = 0.0;
       for(int zdash = 0; zdash < numModelsM2; zdash++){
           tmpDenominator += alpha[i][ntimes-1][zdash];
       }
#endif
*/

       for(int t = eBegin[i]+1; t < eEnd[i]+1; t++){
       // for(int t = 1; t < ntimes; t++){
            for(int zdash = 0; zdash < numModelsM2; zdash++){
                for(int z = 0; z < numModelsM2; z++){
                   long double tmpLogDouble = 0;

                    for(int k = 0; k < (numProdsK-1); k++){
                        if(dDash[i][t][k] == 1){
                            if(q[z][k] > 0.0){
                                //tmpLogDouble += wts[i][t]*log(q[z][k]);
                                tmpLogDouble += log(q[z][k]);
                            }
                            else{
                                cout<<"weird expectation 13."<<endl;
                                exit(-1);
                            }
                        }
                        if(dDash[i][t][k] == 0){
                            if(q[z][k] < 1.0){
                                //tmpLogDouble += wts[i][t]*log(1 - q[z][k]);
                                tmpLogDouble += log(1 - q[z][k]);
                                tmpLogDouble += logDNorm(y[i][t][k],mu[z][k],sigmaSq[z][k]);
                            }
                            else{
                                cout<<"weird expectation 14."<<endl;
                                exit(-1);
                            }
                        }
                    }

                    for(int k = numProdsK-1; k < 2*(numProdsK-1); k++){
                        if(dDash[i][t][k] == 1){
                            if(q[z][k+1] > 0.0){
                                //tmpLogDouble += wts[i][t]*log(q[z][k+1]);
                                tmpLogDouble += log(q[z][k+1]);
                            }
                            else{
                                cout<<"weird expectation 15."<<endl;
                                exit(-1);
                            }
                        }
                        if(dDash[i][t][k] == 0){
                            if(q[z][k+1] < 1.0){
                                //tmpLogDouble += wts[i][t]*log(1 - q[z][k+1]);
                                tmpLogDouble += log(1 - q[z][k+1]);
                                tmpLogDouble += logDNorm(y[i][t][k],mu[z][k],sigmaSq[z][k]);
                            }
                            else{
                                cout<<"weird expectation 16."<<endl;
                                exit(-1);
                            }
                        }
                    }

                    int kidx=numProdsK-1;
                    if(d[i][t][kidx]){
                        assert(q[z][kidx] != 0.0); //this should neer happen:
                        //tmpLogDouble += wts[i][t]*log(q[z][kidx]);
                        tmpLogDouble += log(q[z][kidx]);
                    }
                    else{
                        assert(q[z][kidx] != 1.0);  //this should never happen
                        //tmpLogDouble += wts[i][t]*log(1 - q[z][kidx]);
                        tmpLogDouble += log(1 - q[z][kidx]);
                    }

                    kidx = 2*numProdsK - 1;
                    if(d[i][t][kidx]){
                        assert(q[z][kidx] != 0.0); //this should never happen
                        //tmpLogDouble += wts[i][t]*log(q[z][kidx]);
                        tmpLogDouble += log(q[z][kidx]);
                    }
                    else{
                        assert(q[z][kidx] != 1.0); //this should never happen
                        //tmpLogDouble += wts[i][t]*log(1 - q[z][kidx]);
                        tmpLogDouble += log(1 - q[z][kidx]);
                    }
		    
/* HJ
#if 0
                    ksi[i][t-1][zdash][z] = alpha[i][t-1][zdash];
                    ksi[i][t-1][zdash][z] += log(transP[zdash][z]);
                    ksi[i][t-1][zdash][z] += beta[i][t][z];
                    ksi[i][t-1][zdash][z] += tmpLogDouble;
#endif
*/

//#if 0

                    ksi[i][t-1][zdash][z] = coeff[i][t];
                    ksi[i][t-1][zdash][z] += alphaN[i][t-1][zdash];
                    ksi[i][t-1][zdash][z] += log(transP[zdash][z]);
                    ksi[i][t-1][zdash][z] += betaN[i][t][z];
                    ksi[i][t-1][zdash][z] += tmpLogDouble;
//#endif

                    if(isnan(exp(ksi[i][t-1][zdash][z])) || isinf(exp(ksi[i][t-1][zdash][z]))){
                        printf("ksi[%d][%d][%d][%d] = %Lf.\n",i,t,zdash,z,ksi[i][t-1][zdash][z]);
                        printf("exp ksi[%d][%d][%d][%d] = %Lf.\n",i,t,zdash,z,exp(ksi[i][t-1][zdash][z]));
                        printf("log(transP[%d][%d] = %Lf.\n",zdash,z,log(transP[zdash][z]));
                        printf("alpha[i][t-1][zdash] = %Lf.\n",alphaN[i][t-1][zdash]);
                        printf("beta[%d][%d][%d] = %Lf.\n",i,t,z,betaN[i][t][z]);
                        printf("%Lf.\n",tmpLogDouble);
                        exit(-1);
                    }
		    
/* HJ
#if 0
                    if(maxKsi[i][t][zdash] < ksi[i][t][zdash][z]){
                        maxKsi[i][t][zdash] = ksi[i][t][zdash][z];
                    }
#endif
*/

                } //end of z loop over numModelsM2
            } //end of zdash loop over numModelsM2
       }//end of t loop over ntimes
    }//end of i loop over numDyads
  }

//#if 0
    vector<vector<long double> > tmpKsiDen;
    tmpKsiDen.resize(numDyads);

#pragma omp parallel for
  for (int section = 0; section < numBlocks; section++) {
    // cout << "Working on section " << section << ": from dyads " << section*blockSize << " to " << std::min(section*blockSize+blockSize, numDyads) << endl;
    for (int i = section*blockSize; i < std::min(section*blockSize+blockSize, numDyads); i++) {
    // for(int i = 0; i < numDyads; i++){
        tmpKsiDen[i].resize(ntimes);
        for(int t = 0; t < ntimes; t++){
                tmpKsiDen[i][t] = 0.0;
        }
    }
  }

    //compute denominator for ksi
#pragma omp parallel for
  for (int section = 0; section < numBlocks; section++) {
    // cout << "Working on section " << section << ": from dyads " << section*blockSize << " to " << std::min(section*blockSize+blockSize, numDyads) << endl;
    for (int i = section*blockSize; i < std::min(section*blockSize+blockSize, numDyads); i++) {
    // for(int i = 0; i < numDyads; i++){
        for(int t  = eBegin[i]+1; t < eEnd[i]+1; t++){
	// for(int t  = 1; t < ntimes; t++){
            for(int z = 0; z < numModelsM2; z++){
                for(int zdash = 0; zdash < numModelsM2; zdash++){
                    //tmpKsiDen[i][t][z] += exp(alpha[i][t][z] + beta[i][t][z] - maxKsi[i][t][z]);
                    ///tmpKsiDen[i][t] += exp(alpha[i][t][z] + beta[i][t][z]);
                    //tmpKsiDen[i][t-1] += exp(alpha[i][ntimes-1][z]);
                    //tmpKsiDen[i][t-1] += exp(alphaN[i][ntimes-1][z]);
                    tmpKsiDen[i][t-1] += exp(ksi[i][t-1][zdash][z]);
                }
            }
	    
/* HJ
#if 0
            if(tmpKsiDen[i][t-1] == 0.0)
            {
                cout<<"Found tmpKsiDen[i][t-1] == 0."<<endl;
                cout<<"exp(coeff[)"<<i<<"]["<<t<<"] = "<<exp(coeff[i][t])<<endl;
                for(int zdash = 0; zdash < numModelsM2; zdash++)
                {
                    cout<<"exp(alphaN[)"<<i<<"]["<<t-1<<"]["<<zdash<<"]) = "<<exp(alphaN[i][t-1][zdash])<<endl;
                    for(int z = 0; z < numModelsM2; z++)
                    {
                        cout<<"transP["<<zdash<<"]["<<z<<"] = "<<transP[zdash][z]<<endl;
                        cout<<"exp(betaN["<<i<<"]["<<t<<"]["<<z<<"]) = "<<exp(betaN[i][t][z])<<endl;;
                        cout<<"exp(beta["<<i<<"]["<<t<<"]["<<z<<"]) = "<<exp(beta[i][t][z])<<endl;
                        for(int zdd = 0; zdd < numModelsM2; zdd++)
                        {
                            cout<<"\t beta["<<i<<"]["<<t+1<<"]["<<zdd<<"] = "<<beta[i][t+1][zdd]<<endl;
                            cout<<"\t transP["<<z<<"]["<<zdd<<"] = "<<transP[z][zdd]<<endl;
                            cout<<"\t exp(alpha["<<i<<"]["<<t<<"]["<<zdd<<"]) = "<<exp(alpha[i][t][zdd])<<endl;
                        }
                    }
                }
                exit(-1);
            }
#endif
*/

        }
    }
  }

//#endif

    //compute final ksi
#pragma omp parallel for
  for (int section = 0; section < numBlocks; section++) {
    // cout << "Working on section " << section << ": from dyads " << section*blockSize << " to " << std::min(section*blockSize+blockSize, numDyads) << endl;
    for (int i = section*blockSize; i < std::min(section*blockSize+blockSize, numDyads); i++) {
    // for(int i = 0; i < numDyads; i++){
        for(int t = eBegin[i]+1; t < eEnd[i]+1; t++){
	// for(int t = 1; t < ntimes; t++){
            for(int zdash = 0; zdash < numModelsM2; zdash++){
                for(int z = 0; z < numModelsM2; z++){
                    ksi[i][t-1][zdash][z] -= log(tmpKsiDen[i][t-1]);
                    ksi[i][t-1][zdash][z] = exp(ksi[i][t-1][zdash][z]);
                    if(isnan(ksi[i][t-1][zdash][z]) || isinf(ksi[i][t-1][zdash][z]))
                    {
                        //cout<<"tmpKsiDen["<<i<<"]["<<t-1<<"] = "<<tmpKsiDen[i][t-1]<<endl;
                        cout<<"ksi is nan or inf."<<endl;
                        cout<<"log(tmpKsiDen[i][t-1] = "<<tmpKsiDen[i][t-1]<<endl;
                        exit(-1);
                    }
                }
            }
        }
    }
  }

#if DEBUG2
    cout<<"After fifth loop over dyads."<<endl;
#endif

    return 1;
}

//TODO: compute univariate probability density of val given mean and var;
long double ZeroTradeModel::logDNorm(long double val,long double mean,long double var){
    long double result;

    if(isinf(1.0/var) || isnan(1.0/var) || fabs(var) <= 0.0){
        cout<<"ZeroTradeModel::logDNorm : encountered zero variance"<<endl;
        cout<<"var = "<<var<<endl;
        exit(-1);
    }

    result = -0.5*log(2.0*PI*var);
    result -= ((val - mean)*(val - mean))/(2.0*var);

    return result;
}

//TODO: implement inverseGamma function
long double ZeroTradeModel::inverseGamma(long double var, long double shape, long double scale){
    long double retVal = 0.0;

    inverse_gamma_distribution<long double> inverse_gamma(shape,scale);
    retVal = pdf(inverse_gamma, var);

    //cout<<"inverseGamma: retVal = "<<retVal<<endl;
    return retVal;
}

long double ZeroTradeModel::betaDistribution(long double var, long double alpha, long double beta){
    long double retVal = 0.0;


    retVal = pow(var,alpha - 1.0)*pow(1.0 - var, beta - 1.0);
    try{
        retVal /= (boost::math::beta(alpha,beta,1.0));
    }
    catch(int e){
        cout<<"Boost encountered exception "<<e<<endl;
        cout<<"alpha = "<<alpha<<endl;
        cout<<"beta = "<<beta<<endl;
    }

    return retVal;    
}

long double ZeroTradeModel::logBetaDistribution2(long double var, long double alpha, long double beta, long double logBetaVal){
    long double retVal = 0.0;

    if(alpha == 1.0 && beta == 1.0){
        return 0.0;
    }

    retVal = (alpha-1.0)*log(var);
    retVal += (beta - 1.0)*log(1.0 - var);
    retVal -= logBetaVal;
    if(isnan(retVal)){
        cout<<"retVal is nan."<<endl;
        exit(-1);   
    }

    return retVal;    
}

//Calculates overall mean for the dataset x[i][j][k] and stores it in the 2*(K-1) vector muAll[k]
void ZeroTradeModel::calculateMeanAll(){

    //vector<vector<int> > count(ntimes,vector<int>(numProdsK2,0));
    vector<int> count(numProdsK2,0);

    //for(int t = 0; t < ntimes; t++){
    for(int k = 0; k < numProdsK2; k++){
        muAll[k] = 0.0;
    }
    //}

    for(int i = 0; i < numDyads; i++){
      // for(int t = 0; t < ntimes; t++){
      for(int t = eBegin[i]; t < eEnd[i]+1; t++){
	for(int k = 0; k < numProdsK2; k++){
	  if(dDash[i][t][k] == 0.0){
	    count[k]++;
	    muAll[k] += y[i][t][k];
	  }
	}
      }
    }

    //for(int t = 0; t < ntimes; t++){
        for(int k = 0; k < 2*(numProdsK-1); k++){
            //muAll[k] /= numDyads;
            if(count[k] == 0){
                muAll[k] = 0.00001;
            }
            else{
                muAll[k] /= (ntimes*count[k]);
            }
        }
    //}
}

/* HJ
#if 0
//Calculates overall covariance for the dataset y[i][k] and stores it in the matrix sigmaAll
void ZeroTradeModel::fixZeroTradedProds(){

    vector<int> count(numProdsK2,0);
    vector<int> untradedProds; 

    for(int i = 0; i < numDyads; i++){
        for(int k = 0; k < numProdsK2; k++){
            if(dDash[i][k] == 0){
                count[k]++;
            }
        }
    }

    for(int k = 0; k < numProdsK2; k++){
        if(count[k] == 0){
            untradedProds.push_back(k);
            cout<<"Product "<<k<<" is not traded, flipping dyad-pairs to fix this."<<endl;
        }
    }

    //The way we have pre-processed, each product will have at least 2 trades
    //among the dyads, so it should show up once on each side of the stacking.
    // !!!!!DONT THREAD!!!!!
    for(int s = 0; s < untradedProds.size(); s++){
         for(int i = 0; i < numDyads; i++){
             int kidx=untradedProds[s];;
            (kidx < numProdsK - 1) ? kidx = kidx + numProdsK - 1 : kidx = kidx - (numProdsK - 1);
            //if(y[i][kidx] != 0.0){
            if(!dDash[i][kidx]){
                cout<<"In order to fix untraded product "<<untradedProds[s]<<" flipping dyad pair ID "<<i<<"."<<endl; 
                flipDyadPair(i);
                break;
            }
        }
    }

   cout<<"Returning from fixZeroTradeProds."<<endl;
}
#endif
*/

/* HJ
#if 0
void ZeroTradeModel::flipDyadPair(int didx){

    cout<<"In flipDyadPair, flipping dyad "<<didx<<"."<<endl;
    RMapType::iterator revit = revDyadMap.find(didx);
    int impID = revit->second.first; 
    int expID = revit->second.second;   
    pair<int,int> IDPair = make_pair<int,int>(impID,expID);
    pair<int,int> swapIDPair = make_pair<int,int>(expID,impID);

    //reverse the entry in revDyadMap
    revit->second = swapIDPair;

    MapType::iterator it = dyadMap.find(IDPair);
    //reverse the entry in dyadMap, first delete current entry and insert a new
    //one
    dyadMap.erase(it);
    dyadMap.insert(MapType::value_type(swapIDPair,didx));

    //Reverse stack y
    for(int k = 0; k < numProdsK - 1; k++){
        long double tmpVal = y[didx][k];
        y[didx][k] = y[didx][k+numProdsK-1];
        y[didx][k+numProdsK-1] = tmpVal;
    }

    //Reverse stack dDash
    for(int k = 0; k < numProdsK - 1; k++){
        long double tmpVal = dDash[didx][k];
        dDash[didx][k] = dDash[didx][k+numProdsK-1];
        dDash[didx][k+numProdsK-1] = tmpVal;
    }

    //Reverse stack d
    // for(int k = 0; k < numProdsK; k++){
    //     long double tmpVal = d[didx][k];
    //     d[didx][k] = d[didx][k+numProdsK-1];
    //     d[didx][k+numProdsK-1] = tmpVal;
    // }
    
    for(int k = 0; k < numProdsK; k++){
        long double tmpVal = d[didx][k];
        d[didx][k] = d[didx][k+numProdsK];
        d[didx][k+numProdsK] = tmpVal;
    }

   cout<<"Returning from flipDyadPair."<<endl;
}
#endif
*/

//Calculates overall covariance for the dataset y[i][k] and stores it in the matrix sigmaAll
void ZeroTradeModel::calculateCovarianceAll(){
    //vector<vector<int> > count(ntimes,vector<int>(numProdsK2,0));
    vector<int> count(numProdsK2,0);

    //for(int t = 0; t < ntimes; t++){
    for(int k = 0; k < numProdsK2; k++){
        sigmaSqAll[k] = 0.0;
    }
    //}

    for(int i = 0; i < numDyads; i++){
        // for(int t = 0; t < ntimes; t++){
        for(int t = eBegin[i]; t < eEnd[i]+1; t++){
            for(int k = 0; k < numProdsK2; k++){
                if(dDash[i][t][k] == 0){
                //if(w[i][k] != 0){
                    count[k]++;
                    sigmaSqAll[k] += (y[i][t][k] - muAll[k])*(y[i][t][k] - muAll[k]);
                }
            }
        }
    }

#if DEBUG
    for(int t = 0; t < ntimes; t++){
        cout<<"calculateCovarianceAll: count["<<t<<"][148] = "<<count[t][148]<<endl;
    }
#endif

    //for(int t = 0; t < ntimes; t++){
        for(int k = 0; k < numProdsK2; k++){
            if(count[k] == 0){
                cout<<"Product k = "<<k<<" has no non-zero values, count["<<k<<"] = "<<count[k]<<endl;
                sigmaSqAll[k] = 0.0000001;
#if EXITONZEROVARIANCE
                exit(-1);  //prior should handle 0 sigma
#endif
            }
            else if(count[k] == 1){
                //sigmaSqAll[k] /= (count[k]);
                sigmaSqAll[k] = 0.0001;
            }
            else{
                sigmaSqAll[k] /= (count[k]-1);
            }
            //sigmaSqAll[k] /= numDyads;
            if(isinf(1.0/sigmaSqAll[k]) || isnan(1.0/sigmaSqAll[k] ) || sigmaSqAll[k] == 0.0){
                cout<<"sigmaSqAll["<<k<<"] = "<<sigmaSqAll[k]<<" count = "<<count[k]<<endl;
                cout<<"muAll[k] = "<<muAll[k]<<endl;
                for(int i = 0; i < numDyads; i++){
                    for(int t = 0; t < ntimes; t++){
#if HJSKIP
		      if ((eBegin[i] <= t) && (t <= eEnd[i])) {
#endif
                        cout<<"y[i][t][k] = "<<y[i][t][k]<<endl;
#if HJSKIP
		      } else {
			cout<<"This (dyad, t) does not exist: y[i][t][k] = "<<y[i][t][k]<<endl;
		      }
#endif
                    }
                }
                exit(-1);
            }
        }
    //}

    for(int i = 0; i < numDyads; i++){
        for(int t = 0; t < ntimes; t++){
	    // HJ: This is for print out only, so I am not skipping the ones not in [eBegin[i], eEnd[i]]
	    //     I will skip it if it causes problems
            for(int k = 0; k < numProdsK2; k++){
                if(dDash[i][t][k] != 0 && sigmaSqAll[k] < 0.0000000001){
#if HJSKIP
		  if ((eBegin[i] <= t) || (t <= eEnd[i])) {
		    cout << "This (dyad, t) does not exist dyad = " << i << " t = " << t << endl;
		  }
#endif
                    cout<<"Problem in covariance calculation, variance of product "<<k<<" is close to zero!"<<endl;
                    cout<<"sigmaSqAll[k] = "<<sigmaSqAll[k]<<endl;
                    exit(-1);
                }
            }
        }
    }

#if DEBUG
   cout<<"Returning from calculateCovarianceAll "<<endl;
#endif
   //exit(-1);
}

// count the number of products in the input file
void ZeroTradeModel::countProds(){
    //ifstream infile(inputFile.c_str());
    ifstream infile(inputFiles[0].c_str());
    string tmpString;
    char* colName;

    if(infile.fail()){
        cout<<"error opening input file"<<endl;
        exit(-1);
    }

    getline(infile,tmpString);

    numProdsK = 0;
    //regex regex_base("([[:alpha:]]+[[:digit:]]*)_([[:digit:]]+)");
    regex regex_base("([[:alpha:]]+[[:digit:]]+)_([[:alnum:]]+)");
    match_results<string::const_iterator> sm;

    colName = strtok(const_cast<char*>(tmpString.c_str()),"\t");
    string colNameStr(colName);

    while(colName != NULL){
        string colNameStr(colName);

	// HJ
	// cout<<colNameStr<<endl;

        if(regex_match(colNameStr,sm,regex_base)){
            numProdsK++;
        }
        else{
	  cout<<"No match for "<<colNameStr<<endl;
        }
        colName = strtok(NULL,"\t");
    }

    cout<<"countProds: numProdsK = "<<numProdsK<<endl;
}

// name input files
void ZeroTradeModel::getInputFilenames(){
    //ifstream infile("test0.txt");

    char* colName;

    regex regex_base("([[:digit:]]+)");
    match_results<string::const_iterator> sm;

    colName = strtok(const_cast<char*>(inputYears.c_str()),",");
    string colNameStr(colName);

    while(colName != NULL){
        string colNameStr(colName);

        // HJ
	// cout<<colNameStr<<endl;

        if(regex_match(colNameStr,sm,regex_base)){
            ostringstream tmpStr;
            tmpStr<<"dynamic_"<<sm[0]<<".csv";
            //inputFiles.push_back(sprintf("dynamic_%s",sm[0].c_str()));
            inputFiles.push_back(tmpStr.str());
        }
        else{
	    cout<<"No match for "<<colNameStr<<endl;
        }
        colName = strtok(NULL,",");
    }

    cout<<"nameInputFiles: found "<<inputFiles.size()<<endl;
    copy(inputFiles.begin(),inputFiles.end(),ostream_iterator<string>(cout,"\n"));
}

//TODO: Reads .dta file and populates x[i][j][k]
void ZeroTradeModel::readInput(){

    ifstream infile[ntimes];
    string tmpString;
    vector<int> tmpProdCount;
    int tmpI, tmpJ, count=0;
    int tmpK;
    float tmpFloat;
    int idx = 0.0;
    char* colName;

    // HJ: use this to track the countries
    std::set<int> cSet;

    for(int t = 0; t < ntimes; t++){
        infile[t].open(inputFiles[t].c_str());
    }

    for(int i = 0; i < ntimes; i++){
        if(infile[i].fail()){
            cout<<"error opening input file"<<endl;
            exit(-1);
        }
    }

    xMap.clear();

    cout<<"ntimes = "<<ntimes<<endl;

    for(int t = 0; t < ntimes; t++){

        count = 0;

        //read in the header line
        //TODO: parse header to save a map of productindex->productname
        getline(infile[t],tmpString);
	
/* HJ
#if 0
        cout<<"Read in the first line "<<tmpString<<endl;
#endif
*/

        int kidx = 0;
        //regex regex_base("([[:alpha:]]+[[:digit:]]*)_([[:digit:]]+)");
        regex regex_base("([[:alpha:]]+[[:digit:]]+)_([[:alnum:]]+)");
        match_results<string::const_iterator> sm;

        //parse the first line in the input file and store product names
        colName = strtok(const_cast<char*>(tmpString.c_str()),"\t");
        //set the string for product code prefix - HS/SITC
        //cout<<colName<<endl;
        string colNameStr(colName);
        //cout<<colNameStr<<endl;

        while(colName != NULL){
            string colNameStr(colName);

            // HJ	  
            // cout<<colNameStr<<endl;

            if(regex_match(colNameStr,sm,regex_base)){
                prodNames[kidx] = colNameStr;
                //prodNames.push_back(colNameStr);
                kidx++;
                // HJ    cout<<"Pushed product name."<<colNameStr<<endl;
		
/* HJ
#if 0
                cout<<"match size = "<<sm.size()<<endl;
                cout<<"match[0] = "<<sm[0]<<endl;
                cout<<"match[1] = "<<sm[1]<<endl;
                cout<<"match[2] = "<<sm[2]<<endl;
#endif
*/

            }
            else{
	        cout<<"No match for "<<colNameStr<<endl;
            }
	    
/* HJ
#if 0
            cout<<"Read in colname "<<colName<<endl;
#endif
*/

            colName = strtok(NULL,"\t");
        }

        prodCodePrefix = sm[0];
	
/* HJ
#if 0
        cout<<"Setting prodCodePrefix to "<<prodCodePrefix<<endl;
#endif
*/

        //cout<<"prodNames.size = "<<prodNames.size()<<endl;
        cout<<"t = "<<t<<endl;
        cout<<"kidx = "<<kidx<<endl;
        cout<<"numProdsK = "<<numProdsK<<endl;
        assert(kidx == numProdsK);
        tmpProdCount.resize(numProdsK,0);

        //cout<<"prodNames[0] = "<<prodNames[0]<<endl;
        //cout<<tmpString<<endl;
        while(!infile[t].eof()){
            //initialize temporary vector to store products
            vector<long double> tmpVector;

            //getline(infile,tmpString);
            //cout<<tmpString<<endl;
            infile[t] >> tmpK;
            if(infile[t].eof()) break;

            infile[t] >> tmpI;

            //TODO: original format, uncomment the following/
            infile[t] >> tmpJ;

	    // HJ: add the country to the country set
	    if (cSet.find(tmpI) == cSet.end()) cSet.insert(tmpI);
	    if (cSet.find(tmpJ) == cSet.end()) cSet.insert(tmpJ);
	 
            // cout<<"count = "<<count<<" tmpI = "<<tmpI<<" tmpJ = "<<tmpJ<<endl;

            //read the next double and exit while loop if eof has been reached
            infile[t] >> tmpFloat;
            tmpVector.push_back((long double)tmpFloat);

            if(infile[t].eof()) break;
            if(tmpFloat != 0.0){
                tmpProdCount[0]++;
            }

            for(int k = 1; k < numProdsK; k++){
                infile[t] >> tmpFloat;
                tmpVector.push_back((long double)tmpFloat);
                if(tmpFloat != 0.0){
                    tmpProdCount[k]++;
                }
            }

            pair<int,int> tmpPair(make_pair(tmpI,tmpJ));
            //find this pair in xMap
            VMapType::iterator xit = xMap.find(tmpPair);
            if(xit != xMap.end()){
                //if pair is found:
                vector<vector<long double> > &mapVal = xit->second;
                mapVal.push_back(tmpVector);
                xit->second = mapVal;
            }
            else{
                //else if pair is not found, then insert into map first time
                vector<vector<long double> > tmpVecVec(1,tmpVector);
		
/* HJ
#if 0
                cout<<"Inserting the first vector."<<endl;
                cout<<"tmpI = "<<tmpI<<" tmpJ = "<<tmpJ<<endl;
                for(int vv = 0; vv < tmpVecVec.size(); vv++){
                    for(int v = 0; v < tmpVecVec[vv].size(); v++){
                        cout<<"tmpVecVec["<<vv<<"]["<<v<<"] = "<<tmpVecVec[vv][v]<<"\t";
                    }
                    cout<<endl;
                }
#endif
*/

                xMap.insert(VMapType::value_type(MatKey(tmpI,tmpJ),tmpVecVec));
            }

            idx++;
            count++;
            //cout<<"Read in count = "<<count<<endl;
        }

        infile[t].close();
        //exit(-1);
    }
    
/* HJ
#if 0
	// debug: print data stored in xMap
	for(VMapType::iterator xit = xMap.begin(); xit != xMap.end(); xit++){
		pair<int,int> mapKey = xit->first;
		cout<<xit->first.first<<"\t"<<xit->first.second<<endl;
        //For this dyad pair, xit->second stores the product vectors over all times
		vector<vector<long double> > &mapVal = xit->second;
		for(int i = 0; i < ntimes; i++){
			for(int j = 0; j < numProdsK; j++){
				cout<<mapVal[i][j]<<"\t";
			}
            cout<<endl;
		}
	}
#endif
*/

    // HJ: if a dyad does not exist for all the years, remove it
    std::set<int>::iterator iterI;
    std::set<int>::iterator iterJ;
    int oneMin;
    int allMax;
    for (iterI=cSet.begin(); iterI!=cSet.end(); ++iterI) {
      for (iterJ=cSet.begin(); iterJ!=cSet.end(); ++iterJ) {
	if (*iterI != *iterJ) {
	  pair<int,int> tmpPair(make_pair(*iterI,*iterJ));
	  VMapType::iterator xit = xMap.find(tmpPair);
	  oneMin = 0;
	  allMax = -1;
	  if(xit != xMap.end()) {
	    // HJ: if pair is found:  
	    vector<vector<long double> > &mapVal = xit->second;
	    for (int i=0; i<mapVal.size(); i++) {
	      oneMin = 0;
	      for (int j=0; j<mapVal[i].size(); j++) {
		if (mapVal[i][j] < 0.0) {
		  oneMin = -1;
		  break;
		}
	      }
	      if (oneMin == 0) {
		allMax = 0;
		break;
	      }
	    }	  
	    // HJ: if all tmpPair[i] rows has some columns [i][j] that are negative, remove it from xMap 
	    if (allMax < 0) {
	      xMap.erase(xit);
	      cout << "This dyad is removed (" << *iterI << ", " << *iterJ << ")" << endl;
	    }
	  // } else {
	  //  // HJ: else if pair is not found, say it
	  //  cout << "This dyad does not exist (" << *iterI << ", " << *iterJ << ")" << endl;
	  }
	}
      }
    }

    cout<<"Returning from readInput."<<endl;
}

int ZeroTradeModel::maximizationTh(int iter,bool &resetFlag,ofstream &logFile){
    int abortFlag = 0;
    //calculate pi
    fill(pi.begin(),pi.end(),0.0);
    long double piDenominator = 0.0;

    for(int z = 0; z < numModelsM2; z++){
        for(int i = 0; i < numDyads; i++){
	  //HJ pi[z] += zeta[i][0][z] + OP_w_ij;
	  pi[z] += zeta[i][eBegin[i]][z] + OP_w_ij;
        } //end of for loop over i
        piDenominator += pi[z];
    }

    cout<<"numDyads = "<<numDyads<<endl;
    cout<<"piDenominator = "<<piDenominator<<endl;
    cout<<"OP_w_ij*numDyads*numModelsM2 = "<<OP_w_ij*numDyads*numModelsM2<<endl;

    for(int z = 0; z < numModelsM2; z++){
        pi[z] /= piDenominator;
        cout<<"maximizationTh: pi["<<z<<"] = "<<pi[z]<<endl; 
    }

/* HJ
#if 0
    //debug - start - set pi to true P
    pi[0] = 0.333;
    pi[1] = 0.333;
    pi[2] = 0.333;
    pi[3] = 0.0001;
    pi[4] = 0.0001;
    pi[5] = 0.0001;
#endif
*/

    //calculate transP

#pragma omp parallel for
    for(int zdash = 0; zdash < numModelsM2; zdash++){
        long double tmpDenominator = 0.0;
        for(int z = 0; z < numModelsM2; z++){
            transP[zdash][z] = 0.0;
            for(int i = 0; i < numDyads; i++){
                for(int t = eBegin[i]+1; t < eEnd[i]+1; t++){
		// for(int t = 1; t < ntimes; t++){

                //for(int t = 0; t < ntimes; t++){
		  transP[zdash][z] += ksi[i][t-1][zdash][z] + OP_w_ij;
		  if(isnan(transP[zdash][z])){
		    printf("transP[%d][%d] = %Lf.\n",zdash,z,transP[zdash][z]);
		    printf("ksi[%d][%d][%d][%d] = %Lf.\n",i,t,zdash,z,ksi[i][t][zdash][z]);
		    printf("OP_w_ij = %Lf.\n",OP_w_ij);
		    exit(EXIT_FAILURE);
		  }
                }
            }
            tmpDenominator += transP[zdash][z];
        }

        for(int z = 0; z < numModelsM2; z++){
            transP[zdash][z] /= tmpDenominator;
            if(isnan(transP[zdash][z])){
                cout<<"zdash = "<<zdash<<" z = "<<z<<endl;
                cout<<"tmpDenominator = "<<tmpDenominator<<endl;
                cout<<"transP[zdash][z] = "<<transP[zdash][z]<<endl;
                exit(-1);
            }
        }
    }
    
/* HJ
#if 0
    //debug - set transP to debugTransP
    for(int zdash = 0; zdash < numModelsM2; zdash++){
        for(int z = 0; z < numModelsM2; z++){
            transP[zdash][z] = debugTransP[zdash][z];
        }
    }
#endif
*/

/* HJ
#if 0
transP[0][0] = 0.2;                                                     
transP[0][1] = 0.2;                                                            
transP[0][2] = 0.6;                                                            
transP[1][0] = 0.1;                                                            
transP[1][1] = 0.8;
transP[1][2] = 0.1;
transP[2][0] = 0.3;
transP[2][1] = 0.3;
transP[2][2] = 0.4;
transP[3][0] = 0.2;                                                     
transP[3][1] = 0.2;                                                            
transP[3][2] = 0.6;                                                            
transP[4][0] = 0.1;                                                            
transP[4][1] = 0.8;
transP[4][2] = 0.1;
transP[5][0] = 0.3;
transP[5][1] = 0.3;
transP[5][2] = 0.4;
#endif
*/

    //zero q matrix, muNew, sigmaSq, sigmaSqNew, kTmpSum
    for(int z = 0; z < numModelsM2; z++){
        for(int k = 0; k < numProdsK2; k++){
            q[z][k] = 0.0;
            muNew[z][k] = 0.0;
            sigmaSq[z][k] = 0.0;
            sigmaSqNew[z][k] =0.0;
            kTmpSum[z][k] = 0.0;
        }
        q[z][2*numProdsK - 1] = 0.0;
        q[z][2*numProdsK - 2] = 0.0;
    }

    //calculate q

#pragma omp parallel for
    for(int z = 0; z < numModelsM; z++){
        long double tmpDenominator = 0.0;
        for(int i = 0; i < numDyads; i++){
            for(int t = eBegin[i]; t < eEnd[i]+1; t++){
	    // for(int t = 0; t < ntimes; t++){
	      for(int k = 0; k < numProdsK; k++){
		q[z][k] += ((zeta[i][t][z] + OP_w_ij)*d[i][t][k] + (zeta[i][t][z+numModelsM]+OP_w_ij)*d[i][t][k+numProdsK]);
		q[z][k+numProdsK] += ((zeta[i][t][z] + OP_w_ij)*d[i][t][k+numProdsK] + (zeta[i][t][z+numModelsM]+OP_w_ij)*d[i][t][k]);
	      } //end of loop over k
	      tmpDenominator += ((zeta[i][t][z] + OP_w_ij) + (zeta[i][t][z+numModelsM] + OP_w_ij));
            } //end of loop over t
        } //end of loop over i

        for(int k = 0; k < 2*numProdsK; k++){
            if(tmpDenominator == 0.0){
                cout<<"Cluster "<<z<<" is empty."<<endl;
                logFile<<"Cluster "<<z<<" is empty."<<endl;
                abortFlag = 1;
            }
            else{
                q[z][k] = q[z][k]/tmpDenominator;
            }
            if(q[z][k] > 1.0000001){
                    cout<<"maximization: iter = "<<iter<<" Found q["<<z<<"]["<<k<<"] = "<<q[z][k]<<", greater than 1.0."<<endl;
                    logFile<<"maximization: iter = "<<iter<<" Found q["<<z<<"]["<<k<<"] = "<<q[z][k]<<", greater than 1.0."<<endl;
                    printf("q[z][k] = %10.10Lf\n",q[z][k]);
                    logFile<<"q[z][k] = "<<q[z][k]<<endl; 
                    cout<<"tmpDenominator = "<<tmpDenominator<<endl; 
                    cout<<"OP_w_ij = "<<OP_w_ij<<endl;
                    abortFlag = 9;
            }
        } //end of loop over k
    } //end of loop over z

    if(abortFlag == 1){
        cout<<"Found empty clusters."<<endl;
        logFile<<"Found empty clusters."<<endl;
        return(-1);
    }
    if(abortFlag == 9){
        cout<<"Found q[z][k] > 1. Simulation might be unstable."<<endl;
        logFile<<"Found q[z][k] > 1. Simulation might be unstable."<<endl;
        exit(-1);
    }

    //calculate mu
    //calcBigN(iter,logFile);

#pragma omp parallel for
    for(int z = 0; z < numModelsM; z++){
        //compute new values of mu[z][k]
        for(int i = 0; i < numDyads; i++){
            for(int t = eBegin[i]; t < eEnd[i]+1; t++){
	    // for(int t = 0; t < ntimes; t++){
	      //assert(accumulate(zeta[i][t].begin(),zeta[i][t].end(),0.0) > 0.999999999999 && accumulate(zeta[i][t].begin(),zeta[i][t].end(),0.0) < 1.00000000001);
	      for(int k = 0; k < numProdsK - 1; k++){
                    muNew[z][k] += (zeta[i][t][z]+OP_w_ij)*(long double)(1-dDash[i][t][k])*(long double)y[i][t][k];
                    muNew[z][k] += (zeta[i][t][z+numModelsM]+OP_w_ij)*(long double)(1-dDash[i][t][k+numProdsK-1])*(long double)y[i][t][k+numProdsK-1];

                    muNew[z][k+numProdsK-1] += (zeta[i][t][z]+OP_w_ij)*(long double)(1-dDash[i][t][k+numProdsK-1])*(long double)y[i][t][k+numProdsK-1];
                    muNew[z][k+numProdsK-1] += (zeta[i][t][z+numModelsM]+OP_w_ij)*(long double)(1-dDash[i][t][k])*(long double)y[i][t][k];

                    kTmpSum[z][k] += (zeta[i][t][z]+OP_w_ij)*(1-dDash[i][t][k]) + (zeta[i][t][z+numModelsM]+OP_w_ij)*(1-dDash[i][t][k+numProdsK-1]);
                    kTmpSum[z][k+numProdsK-1] += (zeta[i][t][z]+OP_w_ij)*(1-dDash[i][t][k+numProdsK-1]) + (zeta[i][t][z+numModelsM]+OP_w_ij)*(1-dDash[i][t][k]);
                } //end of k loop over numProdsK2
            } //end of t loop over ntimes
        } //end of i loop over numDyads

        for(int k = 0; k < numProdsK2; k++){
            //muNew[z][k] = (muNew[z][k] + BAYES_TAU*muAll[k]) / (kTmpSum[z][k] + BAYES_TAU);
            //mu[z][k] = (muNew[z][k] + BAYES_TAU*muAll[k]) / (kTmpSum[z][k] + BAYES_TAU);

            if(kTmpSum[z][k] == 0.0){
                //if(muNew[z][k] != 0.0){
                    cout<<"OP_w_ij = "<<OP_w_ij<<endl;
		    // HJ revised the following error message 2/28/2017
                    // cout<<"kTmpSum["<<k<<"] = "<<kTmpSum[z][k]<<" and muNew["<<z<<"]["<<k<<"] = "<<muNew[z][k]<<" not equal to 0."<<endl;
                    // logFile<<"kTmpSum["<<k<<"] = "<<kTmpSum[z][k]<<" and muNew["<<z<<"]["<<k<<"] = "<<muNew[z][k]<<" not equal to 0."<<endl;
                    cout << "kTmpSum[" << z << "]["<< k << "] = " << kTmpSum[z][k] << endl;
                    logFile << "kTmpSum[" << z << "][" << k << "] = "<< kTmpSum[z][k] << endl;
                    abortFlag = 9;
                //}
            }
        }

    }

    if(abortFlag == 9){
        // HJ revised the following error message 2/28/2017
        // cout<<"Found kTmpSum equal to 0 and corresponding muNew[z][k] not equal to 0."<<endl;
        // logFile<<"Found kTmpSum equal to 0 and corresponding muNew[z][k] not equal to 0."<<endl;
        cout<<"Found kTmpSum equal to 0. (hint: all 0s in some columns?)"<<endl;
        logFile<<"Found kTmpSum equal to 0. (hint: all 0s in some columns?"<<endl;
        exit(-1);
    }

#pragma omp parallel for
    for(int z = 0; z < numModelsM; z++){
        //compute new values of sigmaSq[z][k]
        for(int i = 0; i < numDyads; i++){
            for(int t = eBegin[i]; t < eEnd[i]+1; t++){
	    // for(int t = 0; t < ntimes; t++){
	      
/* HJ
#if 0
                if(!(accumulate(zeta[i][t].begin(),zeta[i][t].end(),0.0) > 0.9999999999 && accumulate(zeta[i][t].begin(),zeta[i][t].end(),0.0) < 1.00000000001)){
                    cout<<"zeta["<<i<<"] over all cluster z does not sum to 1. This simulation will stop."<<endl;
                    logFile<<"zeta["<<i<<"] over all cluster z does not sum to 1. This simulation will stop."<<endl;
                    abortFlag = 9;
                }
#endif
*/

                for(int k = 0; k < numProdsK - 1; k++){
                    //TODO: this can be made more efficient, same tmpDouble can be
                    //reused here
                    long double tmpDouble;

                    tmpDouble = ((long double)y[i][t][k] - mu[z][k]);
                    sigmaSqNew[z][k] += (zeta[i][t][z]+OP_w_ij)*(long double)(1-dDash[i][t][k])*tmpDouble*tmpDouble;
                    tmpDouble = ((long double)y[i][t][k+numProdsK-1] - mu[z][k+numProdsK-1]);
                    sigmaSqNew[z][k] += (zeta[i][t][z+numModelsM]+OP_w_ij)*(long double)(1-dDash[i][t][k+numProdsK-1])*tmpDouble*tmpDouble;

                    tmpDouble = ((long double)y[i][t][k+numProdsK-1] - mu[z][k+numProdsK-1]);
                    sigmaSqNew[z][k+numProdsK-1] += (zeta[i][t][z]+OP_w_ij)*(long double)(1-dDash[i][t][k+numProdsK-1])*tmpDouble*tmpDouble;
                    tmpDouble = ((long double)y[i][t][k] - mu[z][k]);
                    sigmaSqNew[z][k+numProdsK-1] += (zeta[i][t][z+numModelsM]+OP_w_ij)*(long double)(1-dDash[i][t][k])*tmpDouble*tmpDouble;
		    
/* HJ
#if 0
                    tmpDouble = ((long double)y[i][k] - muNew[z][k]);
                    sigmaSqNew[z][k] += ((long double)zeta[i][z])*(long double)(1-dDash[i][k])*tmpDouble*tmpDouble;
                    tmpDouble = ((long double)y[i][k+numProdsK-1] - muNew[z][k+numProdsK-1]);
                    sigmaSqNew[z][k] += ((long double)zeta[i][z+numModelsM])*(long double)(1-dDash[i][k+numProdsK-1])*tmpDouble*tmpDouble;

                    tmpDouble = ((long double)y[i][k+numProdsK-1] - muNew[z][k+numProdsK-1]);
                    sigmaSqNew[z][k+numProdsK-1] += ((long double)zeta[i][z])*(long double)(1-dDash[i][k+numProdsK-1])*tmpDouble*tmpDouble;
                    tmpDouble = ((long double)y[i][k] - muNew[z][k]);
                    sigmaSqNew[z][k+numProdsK-1] += ((long double)zeta[i][z+numModelsM])*(long double)(1-dDash[i][k])*tmpDouble*tmpDouble;
#endif
*/

                }
            }//end of loop over time
        }//end of loop over dyads

        for(int k = 0; k < numProdsK2; k++){
            sigmaSq[z][k] = sigmaSqNew[z][k]/(kTmpSum[z][k]); 

            assert(sigmaSq[z][k] != 0.0);
            if(isnan(sigmaSq[z][k])){
                    cout<<"sigmaSq[z][k] is nan: "<<sigmaSq[z][k]<<endl;
                    cout<<"sigmaSqNew[z][k] = "<<sigmaSqNew[z][k]<<endl;
                    //cout<<"tmpDouble = "<<tmpDouble<<endl;
                    cout<<"kTmpSum[z][k] = "<<kTmpSum[z][k]<<endl; 

                    logFile<<"sigmaSq[z][k] is nan: "<<sigmaSq[z][k]<<endl;
                    logFile<<"sigmaSqNew[z][k] = "<<sigmaSqNew[z][k]<<endl;
                    //logFile<<"tmpDouble = "<<tmpDouble<<endl;
                    logFile<<"kTmpSum[z][k] = "<<kTmpSum[z][k]<<endl; 
            }

            //mu[z][k] = (muNew[z][k] + BAYES_TAU*muAll[k]) / (kTmpSum[z][k] + BAYES_TAU);
            //mu[z][k] = (muNew[z][k] + BAYES_TAU*muClust[z][k]) / (kTmpSum[z][k] + BAYES_TAU);
            mu[z][k] = (muNew[z][k]) / (kTmpSum[z][k]);

        }

	}//end of loop over z

    if(abortFlag == 1){
        cout<<"Found empty clusters."<<endl;
        logFile<<"Found empty clusters."<<endl;
        return(-1);
    }
    if(abortFlag == 9){
        cout<<"Found q[z][k] > 1. Simulation might be unstable."<<endl;
        logFile<<"Found q[z][k] > 1. Simulation might be unstable."<<endl;
        exit(-1);
    }

	//unfold, the mu vector and sigma matrices so that f(z + M) = Af(z) for z = 0,....,M-1
	for(int z = 0; z < numModelsM; z++){

        for(int k = 0; k < numProdsK - 1; k++){
            mu[z+numModelsM][k] = mu[z][k+numProdsK-1];
            mu[z+numModelsM][k+numProdsK-1] = mu[z][k]; 

            sigmaSq[z+numModelsM][k] = sigmaSq[z][k+numProdsK-1];
            sigmaSq[z+numModelsM][k+numProdsK-1] = sigmaSq[z][k]; 
        }

        for(int k = 0; k < numProdsK; k++){
            q[z+numModelsM][k] = q[z][k+numProdsK];
            q[z+numModelsM][k+numProdsK]=q[z][k]; 
            assert(q[z][k+numProdsK] >= 0.0);
            assert(q[z][k+numProdsK] <= 1.0);
            assert(q[z][k] >= 0.0);
            assert(q[z][k] <= 1.0);
        }

	}
	
/* HJ
#if 0
    //debug - set true values of q, sigmasq, mu
    for(int z = 0; z < numModelsM2; z++){
        for(int k = 0; k < numProdsK2; k++){
           mu[z][k] = muTrue[z][k];
            sigmaSq[z][k] = sigmaSqTrue[z][k];
        }
    }
#endif
*/

/* HJ
#if 0
    for(int z = 0; z < numModelsM2; z++){
        for(int k = 0; k < 2*numProdsK; k++){
            q[z][k] = qTrue[z][k];
        }
    }
#endif
*/

    return 0;
}

void ZeroTradeModel::reorganizeProducts(){
    //vMap.insert(VMapType::value_type(mapKey,tmpVector));
    VMapType::iterator vit;
    vector<vector<int> > prodCount;
    vector<int>::iterator it;

    maxFreqProd = -999;

    //initialize count of products over all times to 0
    for(int t = 0; t < ntimes; t++){
        prodCount[t] = vector<int>(numProdsK, 0);
    }

    for(vit = vMap.begin(); vit != vMap.end(); vit++){
        for(int t = 0; t < ntimes; t++){
            for(int k = 0; k < numProdsK; k++){
                if(vit->second[t][k] > 0){
                    prodCount[t][k]++;
                }
            }
        }
    }

//TODO: is it OK to include products that are not traded at all

/* HJ
#if 0
    for(int k = 0; k < numProdsK; k++){
        if(prodCount[k] == 0){
            cout<<"Product "<<k<<" (0-index) is not traded at all."<<endl;
            cout<<"Product "<<prodNames[k]<<" (0-index) is not traded at all."<<endl;
            exit(-1);
        }
    }
#endif
*/

    int maxFreqProdAll = -9999;
    for(int t = 0; t < ntimes; t++){
        it = max_element(prodCount[t].begin(),prodCount[t].end());
        maxFreqProd = distance(prodCount[t].begin(),it);
        if(maxFreqProd > maxFreqProdAll){
            maxFreqProdAll = maxFreqProd;
        }
#if DEBUG
        cout<<"ZeroTradeModel::reorganizeProducts Index (0-index) of most frequently traded product is: "<<maxFreqProd<<endl;
        cout<<"ZeroTradeModel::reorganizeProducts Index (0-index) of most frequently traded product name is: "<<prodNames[maxFreqProd]<<endl;
#endif
    }
#if DEBUG
        cout<<"ZeroTradeModel::reorganizeProducts Index (0-index) of most frequently traded product is: "<<maxFreqProdAll<<endl;
#endif

    if(maxFreqProdAll != numProdsK-1){
        //swap the max_element with the last element in the product names vec
        swap(prodNames[maxFreqProdAll],prodNames[numProdsK-1]);

        //swap the max_element with the last element in each of the dyads
        for(vit = vMap.begin(); vit != vMap.end(); vit++){
            vector<vector<long double> > &tmpVector = vit->second;
            for(int t = 0; t < ntimes; t++){
                swap(tmpVector[t][maxFreqProd], tmpVector[t][numProdsK - 1]);
            }
        }
    }

    //after swapping print out product names again
#if DEBUG
    cout<<"reorganizeProducts: prodNames["<<maxFreqProd<<"] = "<<prodNames[maxFreqProd]<<endl;
    cout<<"reorganizePorducts: prodNames["<<numProdsK-1<<"] = "<<prodNames[numProdsK-1]<<endl;
#endif

    return;
}

/* HJ
#if 0
void ZeroTradeModel::analytics(int isim, ofstream &logFile){
   cout<<"Analyzing clusters for simulation "<<isim<<endl;
   logFile<<"Analyzing clusters for simulation "<<isim<<endl;

   //clusterOcc.resize(numModelsM);
   //maxZeta.resize(numDyads);
   for(int z = 0; z< numModelsM; z++){
       clusterOcc[z].clear();
   }

   int dyadCount = 0;
   long double  max_zeta, tmpSum;
   int maxZ;
   //Calculate Cluster Occupancy
   for(int i = 0;i < numDyads; i++){
            max_zeta = -9999.0;
            for(int z = 0; z < numModelsM; z++){
                  tmpSum = zeta[i][z] + zeta[i][z+numModelsM];
                  if(tmpSum > max_zeta){
                       max_zeta = tmpSum;
                       maxZ = z;
                   }
              }
              clusterOcc[maxZ].push_back(i);
              maxZeta[i] = maxZ;
     }

     logFile<<"The following clusters are not empty: "<<endl;
     int tmpDPCount=0;
     vector<int> nempty;
     for(int z = 0; z < numModelsM; z++){
             if(clusterOcc[z].size() != 0){
                  nempty.push_back(z);
                  tmpDPCount += clusterOcc[z].size();
                  //cout<<"Cluster "<<z<<" has "<<clusterOcc[z].size()<<" dyad-pairs."<<endl; 
                  logFile<<"Cluster "<<z<<" has "<<clusterOcc[z].size()<<" dyad-pairs."<<endl; 
                  for(int i = 0; i < clusterOcc[z].size(); i++){
                      int DPID = clusterOcc[z][i];

                      //!!!!!FLIP!!!!!!!
                      if(zeta[DPID][z] < zeta[DPID][z+numModelsM]){
                          //flip y
                          vector<long double> tmpVec(2*(numProdsK-1));
                          for(int k = 0; k < numProdsK - 1; k++){
                              tmpVec[k] = y[DPID][k+numProdsK-1];
                              tmpVec[k+numProdsK-1] = y[DPID][k];
                          }
                          y[DPID].assign(tmpVec.begin(),tmpVec.begin() + 2*(numProdsK - 1));

                          //flip V - basicaly you want to flip order in which
                          //V_ij and V_ji will be printed out
                          int impID = (revDyadMap.find(DPID)->second).first;
                          int expID = (revDyadMap.find(DPID)->second).second;
                          //VMapType::iterator vit = vMap.find(make_pair(expID,impID));
                          revDyadMap.find(DPID)->second.first = expID; 
                          revDyadMap.find(DPID)->second.second = impID; 

                          long double tmpVal = zeta[DPID][z+numModelsM];
                          zeta[DPID][z+numModelsM] = zeta[DPID][z];
                          zeta[DPID][z] = tmpVal; 

                      }

                  }

             }
     }

      assert(numDyads == tmpDPCount);
}
#endif
*/

/* HJ
#if 0
void ZeroTradeModel::chkpModelParams(int isim,ofstream &logFile){
    chkpObj.chkpModelParams(isim, logFile,
                                numProdsK,
                                numModelsM,
                                numDyads,
                                prodNames,
                                prodCodePrefix,
                                q,
                                mu,
                                sigmaSq,
                                zeta,
                                pi,
                                revDyadMap,
                                llVal,
                                bicVal);
}
#endif
*/

/* HJ
# if 0
void ZeroTradeModel::chkpModelParams(int isim,ofstream &logFile){

        {
            ofstream tmpFile;
            ostringstream tmpFileName;
            cout<<"Checkpointing simulation "<<isim<<"."<<endl;
            logFile<<"Checkpointing simulation "<<isim<<"."<<endl;
            tmpFileName.clear();
            tmpFileName.seekp(0);
            tmpFileName<<"MAXIMUM_CHECKPOINT_Q.txt";
            tmpFile.open(tmpFileName.str().c_str());
            for(int k = 0; k < numProdsK; k++){
                tmpFile<<prodNames[k]<<"\t";
            }
            for(int k = 0; k < numProdsK; k++){
                tmpFile<<prodNames[k]<<"\t";
            }
            tmpFile<<endl;
            for(int z = 0; z < numModelsM2; z++){
                std::copy(q[z].begin(),q[z].end(),ostream_iterator<long double>(tmpFile,"\t"));
                tmpFile<<endl;
            }
            tmpFile.close();
        }

#if SIM_CHKP
        {
            ofstream tmpFile2;
            ostringstream tmpFileName2; 
            tmpFileName2<<"simulation_"<<(isim)<<"_CHECKPOINT_Q.txt";
            tmpFile2.open(tmpFileName2.str().c_str());
            for(int k = 0; k < numProdsK; k++){
                tmpFile2<<prodNames[k]<<"\t";
            }
            for(int k = 0; k < numProdsK; k++){
                tmpFile2<<prodNames[k]<<"\t";
            }
            tmpFile2<<endl;
            for(int z = 0; z < numModelsM2; z++){
                std::copy(q[z].begin(),q[z].end(),ostream_iterator<long double>(tmpFile2,"\t"));
                tmpFile2<<endl;
            }
            tmpFile2.close();
        }
#endif

        {
            ostringstream tmpFileName;
            ofstream tmpFile; 
            tmpFileName.clear();
            tmpFileName.seekp(0);
            tmpFileName<<"MAXIMUM_CHECKPOINT_PI.txt";
            //tmpFileName<<"simulation_"<<(isim)<<"_CHECKPOINT_PI.txt";
            tmpFile.open(tmpFileName.str().c_str());
            std::copy(pi.begin(),pi.end(),ostream_iterator<long double>(tmpFile,"\n"));
            tmpFile.close();
        }

#if SIM_CHKP
        {
            ostringstream tmpFileName2;
            tmpFileName2.clear();
            tmpFileName2.seekp(0);
            tmpFileName2<<"simulation_"<<(isim)<<"_CHECKPOINT_PI.txt";
            ofstream tmpFile2; 
            tmpFile2.open(tmpFileName2.str().c_str());
            for(int k = 0; k < numProdsK - 1; k++){
                tmpFile2<<prodNames[k]<<"\t";
            }
            for(int k = 0; k < numProdsK - 1; k++){
                tmpFile2<<prodNames[k]<<"\t";
            }
            tmpFile2<<endl;
            std::copy(pi.begin(),pi.end(),ostream_iterator<long double>(tmpFile2,"\n"));
            tmpFile2.close();
        }
#endif

        {
            ostringstream tmpFileName;
            ofstream tmpFile; 
            tmpFileName.clear();
            tmpFileName.seekp(0);
            tmpFileName<<"MAXIMUM_CHECKPOINT_MU.txt";
            //tmpFileName<<"simulation_"<<(isim)<<"_CHECKPOINT_MU.txt";
            tmpFile.open(tmpFileName.str().c_str());
            for(int k = 0; k < numProdsK - 1; k++){
                tmpFile<<prodNames[k]<<"\t";
            }
            for(int k = 0; k < numProdsK - 1; k++){
                tmpFile<<prodNames[k]<<"\t";
            }
            tmpFile<<endl;
            for(int z = 0; z < numModelsM2; z++){
                std::copy(mu[z].begin(),mu[z].end(),ostream_iterator<long double>(tmpFile,"\t"));
                tmpFile<<endl;
            }
                tmpFile.close();
            }

#if SIM_CHKP
            {
                ostringstream tmpFileName2;
                ofstream tmpFile2; 
                tmpFileName2.clear();
                tmpFileName2.seekp(0);
                tmpFileName2<<"simulation_"<<(isim)<<"_CHECKPOINT_MU.txt";
                tmpFile2.open(tmpFileName2.str().c_str());
                for(int k = 0; k < numProdsK - 1; k++){
                    tmpFile2<<prodNames[k]<<"\t";
                }
                for(int k = 0; k < numProdsK - 1; k++){
                    tmpFile2<<prodNames[k]<<"\t";
                }
                tmpFile2<<endl;
                for(int z = 0; z < numModelsM2; z++){
                    std::copy(mu[z].begin(),mu[z].end(),ostream_iterator<long double>(tmpFile2,"\t"));
                    tmpFile2<<endl;
                }
                tmpFile2.close();
            }
#endif

            {
                ostringstream tmpFileName;
                ofstream tmpFile; 
                tmpFileName.clear();
                tmpFileName.seekp(0);
                tmpFileName<<"MAXIMUM_CHECKPOINT_SIGMASQ.txt";
                //tmpFileName<<"simulation_"<<(isim)<<"_CHECKPOINT_SIGMASQ.txt";
                tmpFile.open(tmpFileName.str().c_str());
                for(int k = 0; k < numProdsK - 1; k++){
                    tmpFile<<prodNames[k]<<"\t";
                }
                for(int k = 0; k < numProdsK - 1; k++){
                    tmpFile<<prodNames[k]<<"\t";
                }
                tmpFile<<endl;
                for(int z = 0; z < numModelsM2; z++){
                    std::copy(sigmaSq[z].begin(),sigmaSq[z].end(),ostream_iterator<long double>(tmpFile,"\t"));
                    tmpFile<<endl;
                }
                tmpFile.close();
            }

#if SIM_CHKP
            {
                ostringstream tmpFileName2;
                ofstream tmpFile2; 
                tmpFileName2.clear();
                tmpFileName2.seekp(0);
                tmpFileName2<<"MAXIMUM_CHECKPOINT_SIGMASQ.txt";
                //tmpFileName<<"simulation_"<<(isim)<<"_CHECKPOINT_SIGMASQ.txt";
                tmpFile2.open(tmpFileName2.str().c_str());
                for(int k = 0; k < numProdsK - 1; k++){
                    tmpFile2<<prodNames[k]<<"\t";
                }
                for(int k = 0; k < numProdsK - 1; k++){
                    tmpFile2<<prodNames[k]<<"\t";
                }
                tmpFile2<<endl;
                for(int z = 0; z < numModelsM2; z++){
                    std::copy(sigmaSq[z].begin(),sigmaSq[z].end(),ostream_iterator<long double>(tmpFile2,"\t"));
                    tmpFile2<<endl;
                }
                tmpFile2.close();
            }
#endif

            {
                ostringstream tmpFileName;
                ofstream tmpFile; 
                tmpFileName<<"MAXIMUM_CHECKPOINT_ZETA.txt";
                //tmpFileName<<"simulation_"<<(isim)<<"_ZETA.txt";
                tmpFile.open(tmpFileName.str().c_str());
                tmpFile<<"IMPID"<<"\t"<<"EXPID"<<"\t"; 
                for(int z = 0; z < numModelsM; z++){
                    tmpFile<<z<<"\t";
                }
                tmpFile<<endl;
                for(int i = 0; i < numDyads; i++){
                    RMapType::iterator vit = revDyadMap.find(i);
                    tmpFile<<vit->second.first<<"\t"<<vit->second.second<<"\t";
                    for(int z = 0; z < numModelsM; z++){
                            tmpFile<<zeta[i][z]+zeta[i][z+numModelsM]<<"\t";
                    }
#if 0
                    for(int z = 0; z < numModelsM2; z++){
                            tmpFile<<zeta[i][z]<<"\t";
                    }
#endif
                    tmpFile<<endl;
                }
                tmpFile.close();
            }

            {
                ostringstream tmpFileName;
                ofstream tmpFile; 
                //use the ostringstream object instantiated above to write out the
                //log-likelihood and BIC values for this simulation
                tmpFileName.clear();
                tmpFileName.seekp(0);
                tmpFileName<<"MAXIMUM_LOGLIKELIHOOD.txt";
                tmpFile.open(tmpFileName.str().c_str());
                tmpFile<<"Final logLikelihood = "<<llVal<<" Final BIC = "<<bicVal<<endl;
                tmpFile.close(); 
            }

}
#endif
*/

/* HJ
#if 0
void ZeroTradeModel::chkpClusters(int isim,ofstream &logFile){

            vector<int> maxZeta;
            vector<vector<int> > clusterOcc;
            vector<int> nemptyIdx;
            ofstream* clusterFiles = new ofstream[numModelsM];
            clusterOcc.resize(numModelsM);
            maxZeta.resize(numDyads);

            int dyadCount = 0;
            double  max_zeta, tmpSum;
            int maxZ;
#if 0
            for(int k = 0; k < numProdsK - 1; k++){
                tmpFile<<prodNames[k]<<"\t";
            }
            for(int k = 0; k < numProdsK - 1; k++){
                tmpFile<<prodNames[k]<<"\t";
            }
            tmpFile<<endl;
#endif
            for(int i = 0;i < numDyads; i++){
                    max_zeta = -1.0;
                    for(int z = 0; z < numModelsM; z++){
                        double tmpSum = zeta[i][z] + zeta[i][z+numModelsM];
                        if(tmpSum > max_zeta){
                            max_zeta = tmpSum;
                            maxZ = z;
                        }
                    }
                    clusterOcc[maxZ].push_back(i);
                    maxZeta[i] = maxZ;
            }


            //TODO: ssert that the y vectors are in the correct  order ie
            //zeta[m] > zeta[m+M]
            logFile<<"The following clusters are not empty: "<<endl;
            int tmpDPCount=0;
            for(int z = 0; z < numModelsM; z++){
                if(clusterOcc[z].size() != 0){
                    nemptyIdx.push_back(z);
                    tmpDPCount += clusterOcc[z].size();
                    //cout<<"Cluster "<<z<<" has "<<clusterOcc[z].size()<<" dyad-pairs."<<endl; 
                    logFile<<"Final: Cluster "<<z<<" has "<<clusterOcc[z].size()<<" dyad-pairs."<<endl; 
                }
            }
            logFile<<"Found "<<nemptyIdx.size()<<" non-empty clusters."<<endl;
            cout<<"Found "<<nemptyIdx.size()<<" non-empty clusters."<<endl;

            assert(numDyads == tmpDPCount);

#if SIM_CHKP
            ofstream* sim_clusterFiles = new ofstream[numModelsM];
#endif

            for(int z = 0; z < numModelsM; z++){
                if(clusterOcc[z].size() != 0){
                    ostringstream zFile, zFile2;
                    //ostringstream zEXPFile, zPROPFile;
                    //ostringstream flipZFile;
                    zFile<<"MAXIMUM_DYADS_cluster_"<<z<<".txt";
                    zFile2<<"simulation_"<<(isim)<<"_DYADS_cluster_"<<z<<".txt";
                    clusterFiles[z].open(zFile.str().c_str());
                    clusterFiles[z]<<"imp\texp\t";
#if SIM_CHKP
                    sim_clusterFiles[z].open(zFile2.str().c_str());
                    sim_clusterFiles[z]<<"imp\texp\t";
#endif

                    for(int k = 0; k < numProdsK; k++){
                        string currProd = prodNames[k];
                        string currProdDesc; 
                        SITCMapType::iterator sit = SITCMap.find(currProd);
                        if(sit == SITCMap.end()){
                            //cout<<"Product code "<<currProd<<" not found in SITC map."<<endl;
                            currProdDesc = currProd;
                        }
                        else{
                            currProdDesc = sit->second; 
                        }
                        //clusterFiles[z]<<prodNames[k]<<"\t";
                        clusterFiles[z]<<"\""<<currProdDesc<<"\""<<"\t";
#if SIM_CHKP
                        //sim_clusterFiles[z]<<prodNames[k]<<"\t";
                        sim_clusterFiles[z]<<"\""<<currProdDesc<<"\""<<"\t";
#endif
                    }
                    for(int k = 0; k < numProdsK; k++){
                        string currProd = prodNames[k];
                        string currProdDesc; 
                        SITCMapType::iterator sit = SITCMap.find(currProd);
                        if(sit == SITCMap.end()){
                            //cout<<"Product code "<<currProd<<" not found in SITC map."<<endl;
                            currProdDesc = currProd;
                        }
                        else{
                            currProdDesc = sit->second; 
                        }
                        //clusterFiles[z]<<prodNames[k]<<"\t";
                        clusterFiles[z]<<"\""<<currProdDesc<<"\""<<"\t";
#if SIM_CHKP
                        //sim_clusterFiles[z]<<prodNames[k]<<"\t";
                        sim_clusterFiles[z]<<"\""<<currProdDesc<<"\""<<"\t";
#endif
                    }
                    clusterFiles[z]<<endl;
#if SIM_CHKP
                    sim_clusterFiles[z]<<endl;
#endif

                    for(int i = 0; i < clusterOcc[z].size(); i++){
                        int tmpDyadPairID = clusterOcc[z][i]; 
                        RMapType::iterator rit = revDyadMap.find(tmpDyadPairID);
                        //int impID = (revDyadMap.find(tmpDyadPairID)->second).first;
                        //int expID = (revDyadMap.find(tmpDyadPairID)->second).second;
                        assert(rit != revDyadMap.end());
                        int impID = rit->second.first;
                        int expID = rit->second.second; 

                        VMapType::iterator vit = vMap.find(make_pair(impID,expID));
                        VMapType::iterator vitPair = vMap.find(make_pair(expID,impID)); 

                        assert(vit != vMap.end());
                        assert(vitPair != vMap.end());

                        vector<long double> &wijVec = vit->second;
                        vector<long double> &wijPairVec = vitPair->second; 
                        clusterFiles[z]<<impID<<"\t"<<expID<<"\t";
                        std::copy(wijVec.begin(),wijVec.end(),ostream_iterator<long double>(clusterFiles[z],"\t"));
                        std::copy(wijPairVec.begin(),wijPairVec.end(),ostream_iterator<long double>(clusterFiles[z],"\t"));
                        clusterFiles[z]<<endl; 
#if SIM_CHKP
                        sim_clusterFiles[z]<<impID<<"\t"<<expID<<"\t";
                        std::copy(wijVec.begin(),wijVec.end(),ostream_iterator<long double>(sim_clusterFiles[z],"\t"));
                        std::copy(wijPairVec.begin(),wijPairVec.end(),ostream_iterator<long double>(sim_clusterFiles[z],"\t"));
                        sim_clusterFiles[z]<<endl; 
#endif
                    }
                    clusterFiles[z].close();
#if SIM_CHKP
                    sim_clusterFiles[z].close();
#endif
                }
            }
#if 0
                zEXPFile.clear();
                zEXPFile.seekp(0);
                zEXPFile<<"MAXIMUM_EXPECTATION_cluster_"<<z<<".txt";
                //zEXPFile<<"simulation_"<<(isim)<<"_EXPECTATION_cluster_"<<z<<".txt";
                clusterFiles[z].open(zEXPFile.str().c_str());

                for(int k = 0; k < 2*(numProdsK-1); k++){
                    clusterFiles[z]<<expectedMu[z][k]<<endl;
                }
                clusterFiles[z].close();

                zPROPFile.clear();
                zPROPFile.seekp(0);
                zPROPFile<<"MAXIMUM_TFPROPORTION_cluster_"<<z<<".txt";
                //zPROPFile<<"simulation_"<<(isim)<<"_TFPROPORTION_cluster_"<<z<<".txt";
                clusterFiles[z].open(zPROPFile.str().c_str());

                for(int k = 0; k < 2*numProdsK; k++){
                    clusterFiles[z]<<prodProportion[z][k]<<endl;
                }
                clusterFiles[z].close();
                }
            }
#endif

            for(int z = 0; z < numModelsM; z++){
                if(clusterOcc[z].size() != 0){
                    ostringstream zFile;
                    zFile<<"MAXIMUM_PI_cluster_"<<z<<".txt";
                    //zFile<<"simulation_"<<(isim)<<"_PI_cluster_"<<z<<".txt";
                    clusterFiles[z].open(zFile.str().c_str());

                    clusterFiles[z]<<"pi["<<z<<"] = "<<pi[z]<<" pi["<<z+numModelsM<<"] = "<<pi[z+numModelsM]<<" "; 
                    clusterFiles[z]<<"aggregate pi["<<z<<"] = "<<pi[z] + pi[z+numModelsM]<<endl;

                    clusterFiles[z].close();
                }
            }

            for(int z = 0; z < numModelsM; z++){
                if(clusterOcc[z].size() != 0){
                    ostringstream zFile;
                    zFile<<"MAXIMUM_Q_cluster_"<<z<<".txt";
                    //zFile<<"simulation_"<<(isim)<<"_Q_cluster_"<<z<<".txt";
                    clusterFiles[z].open(zFile.str().c_str());

                    for(int k = 0; k < 2*numProdsK; k++){
                        clusterFiles[z]<<q[z][k]<<endl;
                    }

                    clusterFiles[z].close();
                }
            }

            for(int z = 0; z < numModelsM; z++){
                if(clusterOcc[z].size() != 0){
                    ostringstream zFile;
                    zFile<<"MAXIMUM_MU_cluster_"<<z<<".txt";
                    //zFile<<"simulation_"<<(isim)<<"_MU_cluster_"<<z<<".txt";
                    clusterFiles[z].open(zFile.str().c_str());

                    //clusterFiles[z]<<"pi["<<z<<"] = "<<pi[z]<<" pi["<<z+numModelsM<<"] = "<<pi[z+numModelsM]<<" "; 
                    //clusterFiles[z]<<"aggregate pi["<<z<<"] = "<<pi[z] + pi[z+numModelsM]<<endl;

                    for(int k = 0; k < numProdsK2; k++){
                        clusterFiles[z]<<mu[z][k]<<endl;
                    }

                    clusterFiles[z].close();
                }
            }
            
            for(int z = 0; z < numModelsM; z++){
                if(clusterOcc[z].size() != 0){
                    ostringstream zFile;
                    zFile<<"MAXIMUM_SIGMA_cluster_"<<z<<".txt";
                    //zFile<<"simulation_"<<(isim)<<"_SIGMA_cluster_"<<z<<".txt";
                    clusterFiles[z].open(zFile.str().c_str());

                    for(int k = 0; k < numProdsK2; k++){
                        clusterFiles[z]<<sigmaSq[z][k]<<endl;
                    }

                    clusterFiles[z].close();
                }
            }


}
#endif
*/

/* HJ
#if 0
void ZeroTradeModel::readCodeDescriptions(string filename){
    ifstream fin;
    string headerLine;
    string nextLine;
    fin.open(filename.c_str());

    getline(fin,headerLine);
    getline(fin,nextLine);

    while(!fin.eof()){
#if 0
        string tmpString1;
        string tmpString2;

        SITC_code<<"SITC_";

        fin>>tmpString1;
        fin>>tmpString2;
        SITC_code<<tmpString2;

        fin>>SITC_desc;
#endif
        char* col1, *col2, *col3;
        string tmpString;
        col2 = strtok(const_cast<char*>(nextLine.c_str()),"\t"); 
        //string SITC_code("SITC_");
        string prodCode(prodCodePrefix);
        prodCode.append("_");
        prodCode.append(col2);

        col3 = strtok(NULL,"\t");

        string prodDesc(col3); 
        ProdMap.insert(ProdMapType::value_type(prodCode,prodDesc));

        getline(fin,nextLine);

    }

    fin.close();
}
#endif
*/
