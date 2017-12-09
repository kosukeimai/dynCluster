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

Analytics::Analytics(int tmp_mc_trials,int tmp_numthreads,int tmp_numModelsM,
                        int tmp_ntimes, int tmp_numProdsK,string q_file,string mu_file,
                        string sigmasq_file,string pi_file,string zeta_file):
                            mc_trials(tmp_mc_trials),numThreads(tmp_numthreads),
                            numModelsM(tmp_numModelsM),ntimes(tmp_ntimes),
                            numProdsK(tmp_numProdsK),
                            q_filename(q_file),mu_filename(mu_file),
                            sigmasq_filename(sigmasq_file),pi_filename(pi_file),
                            zeta_filename(zeta_file){

        numModelsM2 = 2*numModelsM;
        numProdsK2 = 2*(numProdsK-1);

        srand(time(NULL));

        allocateMemory(pi,
                        sigmaSq,
                        mu,
                        q,
                        zeta,
                        clusterOcc,
                        &prodProportion);

        cout<<"Analytics::Analytics(): Before call to readInput()."<<endl;
        readInput();
        cout<<"Analytics::Analytics(): After call to readInput()."<<endl;

        //readCodeDescriptions("description.txt");
}

Analytics::~Analytics(){

}

void Analytics::allocateMemory(vector<long double>& piArg,
                        vector<vector<long double> >& sigmaSqArg,
                        vector<vector<long double> >& muArg,
                        vector<vector<long double> >& qArg,
                        vector<vector<vector<long double> > >& zetaArg,
                        vector<vector<vector<int > > >& clusterOccArg,
                        long double ****prodProportionPtrArg){

	piArg.resize(numModelsM2);
    sigmaSqArg.resize(numModelsM2);
	muArg.resize(numModelsM2);
	for(int i = 0; i < numModelsM2; i++){
        sigmaSqArg[i].resize(numProdsK2); 
		muArg[i].resize(numProdsK2);
	}

    qArg.resize(numModelsM2);
    for(int m = 0; m < numModelsM2; m++){
        qArg[m].resize(2*numProdsK);  
    }

#if 0
    //expectedMu has the same size as mu 
    expectedMuArg.resize(numModelsM); 
    for(int i = 0; i < numModelsM; i++){
        expectedMuArg[i].resize(2*(numProdsK-1));
    }
#endif

    //(*prodProportionPtrArg) = (long double**)malloc(numModelsM*sizeof(long double*));
    //long double **prodProportionArg = (*prodProportionPtrArg);
    //for(int z = 0; z < numModelsM; z++){
    //    prodProportionArg[z] = (long double*)malloc(2*numProdsK*sizeof(long double));
    //    for(int k = 0; k < 2*numProdsK; k++){
    //       prodProportion[z][k] = -9.0;
    //    }
    //}

    long double*** currPtr = (long double***)malloc(numThreads*sizeof(long double**));
    for(int ith = 0; ith < numThreads; ith++){
        currPtr[ith] = (long double**)malloc(numModelsM*sizeof(long double*));
        for(int z = 0; z < numModelsM; z++){
            currPtr[ith][z] = (long double*)malloc(2*numProdsK*sizeof(long double));
            for(int k = 0; k < 2*numProdsK; k++){
                currPtr[ith][z][k] = -9.0;
            }
        }
    }

    (*prodProportionPtrArg) = currPtr;

    zetaArg.resize(ntimes);

    clusterOccArg.resize(ntimes);
    for(int t = 0; t < ntimes; t++)
    {
        clusterOccArg[t].resize(numModelsM);
    }
}


void Analytics::readInput(){
    ifstream fin;
    fin.open(q_filename.c_str());
    string qHeader;
    char* colName;
    //regex regex_base(".*([[:digit:]]+).*");
    //regex regex_base("([:alpha:]+_[:digit:]+");
    regex regex_base("([[:alpha:]]+[[:digit:]]+)_([[:alnum:]]+)");
    match_results<string::const_iterator> sm;
    getline(fin,qHeader);

    //parse the qHeader to get the K product category names
    colName = strtok(const_cast<char*>(qHeader.c_str()),"\t");

    while(colName != NULL){
        string tmpStr(colName);
        if(regex_match(tmpStr,sm,regex_base)){
            prodNames.push_back(tmpStr);
            colName = strtok(NULL,"\t"); 
        //    cout<<"product "<<prodNames.size() + 1<<" "<<tmpStr<<endl;
        }
    }

    //set the string for product code prefix - HS/SITC
    prodCodePrefix = sm[0];
    #if 0
    if(regex_match(string(colName),sm,regex_base)){
        productCodePrefix = sm[0];
        cout<<"Setting productCodePrefix to "<<productCodePrefix<<endl;
    }
    else{
        cout<<"Input file contains bad format for product code."<<endl;
        cout<<"Product code should be [a-zA-Z]+_[0-9]+ or [a-zA-Z]+_[0-9]+. Exiting."<<endl;
        exit(-1);
    }
    #endif


    //only the first numProdsK values in prodNames will be used
    assert(prodNames.size() == 2*numProdsK);

    for(int z = 0; z < numModelsM2; z++){
        for(int k = 0; k < 2*numProdsK; k++){
            fin >> q[z][k];
        }
    }
    fin.close();
    cout<<"Analytics::readInput(): done reading q file."<<endl;

    fin.open(mu_filename.c_str());
    string muHeader;
    getline(fin,muHeader);
    for(int z = 0; z < numModelsM2; z++){
        for(int k = 0; k < 2*(numProdsK-1); k++){
            fin >> mu[z][k];
        }
    }
    fin.close();
    cout<<"Analytics::readInput(): done reading mu file."<<endl;

    fin.open(sigmasq_filename.c_str());
    string sigmasqHeader;
    getline(fin,sigmasqHeader);
    for(int z = 0; z < numModelsM2; z++){
        for(int k = 0; k < 2*(numProdsK-1); k++){
            fin >> sigmaSq[z][k];
        }
    }
    fin.close();
    cout<<"Analytics::readInput(): done reading sigmasq file."<<endl;

    fin.open(pi_filename.c_str());
    string piHeader;
    //getline(fin,piHeader); 
    for(int z = 0; z < numModelsM2; z++){
        fin >> pi[z];
    }
    fin.close();
    cout<<"Analytics::readInput(): done reading pi file."<<endl;

    for(int t = 0; t < ntimes; t++)
    {
        ostringstream tmpFileName;
        tmpFileName<<"CHECKPOINT_ZETA_T"<<t<<".txt";
        //fin.open(zeta_filename.c_str());
        fin.open(tmpFileName.str().c_str());
        string zetaHeader;
        getline(fin,zetaHeader);
        numDyads = 0;
        string tmpString;
        while(!fin.eof()){
            getline(fin,tmpString);
            if(fin.eof()) break;
            numDyads++;
        }
        fin.close();

        //fin.open(zeta_filename.c_str());
        fin.open(tmpFileName.str().c_str());
        getline(fin,zetaHeader);
        zeta[t].resize(numDyads);
        for(int i = 0; i < numDyads; i++){
            zeta[t][i].resize(numModelsM);
            int tmpImp, tmpExp;
            fin >> tmpImp;
            fin >> tmpExp;
            for(int z = 0; z < numModelsM; z++){
                fin >> zeta[t][i][z];
            }
        }
        fin.close(); 
    }
    cout<<"Analytics::readInput(): done reading zeta file."<<endl;

}

void Analytics::readInputQ(string filename,vector<vector<long double> >& qVecVec,vector<string>& namesVec){
    ifstream fin;
    fin.open(filename.c_str());
    string qHeader;
    char* colName;
    regex regex_base("(SITC_.*)");
    match_results<string::const_iterator> sm;
    getline(fin,qHeader); 

    //clear the namesVec
    namesVec.clear();

    //parse the qHeader to get the K product category names 
    colName = strtok(const_cast<char*>(qHeader.c_str()),"\t"); 
    while(colName != NULL){
        string tmpStr(colName);
        if(regex_match(tmpStr,sm,regex_base)){
            namesVec.push_back(tmpStr);
            colName = strtok(NULL,"\t"); 
        }
    }

    //only the first numProdsK values in prodNames will be used
    assert(namesVec.size() == 2*numProdsK);

    for(int z = 0; z < numModelsM2; z++){
        for(int k = 0; k < 2*numProdsK; k++){
            fin >> qVecVec[z][k];
        }
    }
    fin.close();
}

void Analytics::readInputMu(string filename,vector<vector<long double> >& muVecVec){
    ifstream fin; 
    fin.open(filename.c_str());
    string muHeader;
    getline(fin,muHeader);
    for(int z = 0; z < numModelsM2; z++){
        for(int k = 0; k < 2*(numProdsK-1); k++){
            fin >> muVecVec[z][k];
        }
    }
    fin.close();
}

void Analytics::readInputSigmaSq(string filename,vector<vector<long double> >& sigmaSqVecVec){
    ifstream fin; 
    fin.open(filename.c_str());
    string sigmaSqHeader;
    getline(fin,sigmaSqHeader);
    for(int z = 0; z < numModelsM2; z++){
        for(int k = 0; k < 2*(numProdsK-1); k++){
            fin >> sigmaSqVecVec[z][k];
        }
    }
    fin.close();
}

void Analytics::readInputPi(string filename,vector<long double>& piVec){
    ifstream fin;
    fin.open(filename.c_str());
    string piHeader;
    //getline(fin,piHeader); 
    for(int z = 0; z < numModelsM2; z++){
        fin >> piVec[z];
    }
    fin.close();
}

void Analytics::readInputZeta(string filename,vector<vector<long double> >& zetaVecVec){
    string zetaHeader;
    ifstream fin;
    fin.open(filename.c_str());
    getline(fin,zetaHeader); 
    zetaVecVec.resize(numDyads);
    for(int i = 0; i < numDyads; i++){
        zetaVecVec[i].resize(numModelsM); 
        int tmpImp, tmpExp;
        fin >> tmpImp;
        fin >> tmpExp;
        for(int z = 0; z < numModelsM; z++){
            fin >> zetaVecVec[i][z]; 
        }
    }
    fin.close(); 
    return; 
}

void Analytics::analytics(int rank){

   vector<vector<int> > maxZeta(ntimes);
   for(int t = 0; t < ntimes; t++)
   {
    maxZeta[t].resize(numDyads);
    //clusterOcc[t].resize(numModelsM);
   }

   for(int t = 0; t < ntimes; t++)
   {
    //find out which clusters are not "empty"
    int dyadCount = 0;
    long double  max_zeta, tmpSum;
    int maxZ;
    //Calculate Cluster Occupancy
    for(int i = 0;i < numDyads; i++){
                max_zeta = -9999.0;
                for(int z = 0; z < numModelsM; z++){
                    tmpSum = zeta[t][i][z];
                    if(tmpSum > max_zeta){
                        max_zeta = tmpSum;
                        maxZ = z;
                    }
                }
                clusterOcc[t][maxZ].push_back(i);
                maxZeta[t][i] = maxZ;
        }
   }

#if 0
    int tmpDPCount=0;
    //vector<int> nempty;
    for(int z = 0; z < numModelsM; z++){
        //cout<<"Occupancy of cluster "<<z<<" is "<<clusterOcc[z].size()<<"\n";
            if(clusterOcc[z].size() != 0){
                //nempty.push_back(z);
                tmpDPCount += clusterOcc[z].size();
                //cout<<"Cluster "<<z<<" has "<<clusterOcc[z].size()<<" dyad-pairs."<<endl; 
            }
    }
    assert(tmpDPCount == numDyads);
#endif
}

//doMC calculates the product proportion for cluster z 
void Analytics::doMC(int z,int rank){

    unsigned int iter;
    int numThreads,threadID;
    int ix; 
    double expectation = 0.0;
    double old_expectation = -999.0;
    double tol=0.0001;
    //vector<double> sampleMu;
    vector<double> sampleProp;
    base_generator_type generatorb2;
    threadID = omp_get_thread_num();

    generatorb2.seed(static_cast<unsigned int>(time(0) + rank*100));

    //an array of uniform distributions for each of the 2*(numProdsK-1)
    //variables 
    uniform_real<double> uni_dist(0.0,1.0);    //zeroProb is between 0 and 1.0 
    vector<normal_distribution<double> > norm_dist; //(mean,sqrt(varSq));
    for(int k = 0; k < 2*(numProdsK-1); k++){
        normal_distribution<double> tmpDist((double)mu[z][k],(double)sqrt(sigmaSq[z][k]));
        norm_dist.push_back(tmpDist);
    }

    //an array of 
    variate_generator<base_generator_type&,uniform_real<double> > uni(generatorb2,uni_dist); 
    vector<variate_generator<base_generator_type*,normal_distribution<double> > > norm_rnd; 
    for(int k = 0; k < 2*(numProdsK-1); k++){
        variate_generator<base_generator_type*,normal_distribution<double> > tmp_rnd(&generatorb2,norm_dist[k]); 
        norm_rnd.push_back(tmp_rnd); 
    }

    sampleProp.resize(2*numProdsK); 
    for(int k = 0; k < 2*numProdsK; k++){
        prodProportion[threadID][z][k] = 0.0;
    } 

    for(iter = 0; iter < mc_trials; iter++){
      double tmpDenominator = 0.0;
      for(int k = 0; k < numProdsK-1; k++){
        if(uni() < q[z][k]){
            sampleProp[k] = 0.0;
        }
        else{
            sampleProp[k] = exp(norm_rnd[k]());
            if(isnan(sampleProp[k]) || isinf(sampleProp[k])){
                cout<<"Rank: "<<rank<<"encountered nan."<<endl;
                cout<<"sampleProp["<<k<<"] = "<<sampleProp[k]<<endl;
                cout<<"mu["<<z<<"]["<<k<<"] = "<<mu[z][k]<<endl;
                cout<<"sigmaSq["<<z<<"]["<<k<<"] = "<<sigmaSq[z][k]<<endl;
                cout<<"q["<<z<<"]["<<k<<"] = "<<q[z][k]<<endl;
                cout<<"pi["<<z<<"] = "<<pi[z]<<endl;
                cout.flush(); 
                exit(-1);
            }
        }
        tmpDenominator += sampleProp[k];
      }
      sampleProp[numProdsK-1] = 1.0; 
      tmpDenominator += 1.0;

      for(int k = 0; k < numProdsK; k++){
            sampleProp[k] /= tmpDenominator; 
      }

      //cout<<"tmpDenominator is "<<tmpDenominator<<endl;
      
      tmpDenominator = 0.0; 
      for(int k = numProdsK-1; k < 2*(numProdsK-1); k++){
        if(uni() < q[z][k+1]){
      //      sampleMu[k] = 0.0;
            sampleProp[k+1] = 0.0;
        }
        else{
            //sampleMu[k] = norm_rnd[k]();
            sampleProp[k+1] = exp(norm_rnd[k]());
            if(isnan(sampleProp[k+1]) || isinf(sampleProp[k+1])){
                cout<<"Rank: "<<rank<<"encountered nan."<<endl;
                cout<<"sampleProp["<<k+1<<"] = "<<sampleProp[k+1]<<endl;
                cout<<"mu["<<z<<"]["<<k<<"] = "<<mu[z][k]<<endl;
                cout<<"sigmaSq["<<z<<"]["<<k<<"] = "<<sigmaSq[z][k]<<endl;
                cout<<"q["<<z<<"]["<<k+1<<"] = "<<q[z][k+1]<<endl;
                cout<<"pi["<<z<<"] = "<<pi[z]<<endl;
                cout.flush();
                exit(-1);
            }
        }
        tmpDenominator += sampleProp[k+1];
      }
      sampleProp[2*numProdsK-1] = 1.0;
      tmpDenominator += 1.0;

      for(int k = numProdsK; k < 2*numProdsK; k++){
        sampleProp[k] /= tmpDenominator;
            if(isnan(sampleProp[k])){
                cout<<"Rank: "<<rank<<"encountered nan."<<endl;
                cout<<"sampleProp["<<k<<"] = "<<sampleProp[k]<<endl;
                cout<<"mu["<<z<<"]["<<k<<"] = "<<mu[z][k]<<endl;
                cout<<"sigmaSq["<<z<<"]["<<k<<"] = "<<sigmaSq[z][k]<<endl;
                cout<<"q["<<z<<"]["<<k<<"] = "<<q[z][k]<<endl;
                cout<<"pi["<<z<<"] = "<<pi[z]<<endl;
                cout<<"tmpDenominator = "<<tmpDenominator<<endl;
                cout.flush();
                exit(-1);
            }
      }

      for(int k = 0; k < 2*numProdsK; k++){
        prodProportion[threadID][z][k] += sampleProp[k];
      }
   
    }

    for(int k = 0; k < 2*numProdsK; k++){
        prodProportion[threadID][z][k] /= mc_trials;
        if(isnan(prodProportion[threadID][z][k])){
                cout<<"Rank: "<<rank<<" encountered nan."<<endl;
                cout<<"prodProportion["<<threadID<<"]["<<z<<"]["<<k<<"] = "<<prodProportion[threadID][z][k]<<endl;
                cout<<"mc_trials = "<<mc_trials<<endl;
                exit(-1);
         }
    }

    return;
}




void Analytics::readCodeDescriptions(string filename){
    ifstream fin;
    string headerLine;
    string nextLine;
    fin.open(filename.c_str());

    getline(fin,headerLine);
    getline(fin,nextLine);

    int count = 0;
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
        string productCode(col2);
        //SITC_code.append(col2);

        col3 = strtok(NULL,"\t");

        string productDesc(col3); 
        CodeMap.insert(CodeMapType::value_type(productCode,productDesc));
        count++;

        getline(fin,nextLine);

    }

    fin.close();
    return;
}

