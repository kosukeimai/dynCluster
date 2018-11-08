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

#include "ChkpBase.h"

ChkpBase::ChkpBase(){
}

ChkpBase::~ChkpBase() {
}


void ChkpBase::chkpModelParams(int isim,ofstream &logFile,
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
			       vector<vector<long double> >& local_transP,
			       RMapType& local_revDyadMap,
			       long double& local_llVal,
			       long double& local_bicVal,
			       vector<vector<int> >& local_dExist,
			       vector<int>& local_eBegin,
			       vector<int>& local_eEnd){

        int local_numModelsM2 = local_numModelsM*2;

        {
            ofstream tmpFile; 
            ostringstream tmpFileName;
            logFile<<"Checkpointing simulation "<<isim<<"."<<endl;
            tmpFileName.clear();
            tmpFileName.seekp(0);
            tmpFileName<<"CHECKPOINT_Q.txt";
            tmpFile.open(tmpFileName.str().c_str());
            for(int k = 0; k < local_numProdsK; k++){
                tmpFile<<local_prodNames[k]<<"\t";
            }
            for(int k = 0; k < local_numProdsK; k++){
                tmpFile<<local_prodNames[k]<<"\t";
            }
            tmpFile<<endl;
            for(int z = 0; z < local_numModelsM2; z++){
                std::copy(local_q[z].begin(),local_q[z].end(),ostream_iterator<long double>(tmpFile,"\t"));
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
            for(int k = 0; k < local_numProdsK; k++){
                tmpFile2<<local_prodNames[k]<<"\t";
            }
            for(int k = 0; k < local_numProdsK; k++){
                tmpFile2<<local_prodNames[k]<<"\t";
            }
            tmpFile2<<endl;
            for(int z = 0; z < local_numModelsM2; z++){
                std::copy(local_q[z].begin(),local_q[z].end(),ostream_iterator<long double>(tmpFile2,"\t"));
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
            tmpFileName<<"CHECKPOINT_PI.txt";
            //tmpFileName<<"simulation_"<<(isim)<<"_CHECKPOINT_PI.txt";
            tmpFile.open(tmpFileName.str().c_str());
            std::copy(local_pi.begin(),local_pi.end(),ostream_iterator<long double>(tmpFile,"\n"));
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
            for(int k = 0; k < local_numProdsK - 1; k++){
                tmpFile2<<local_prodNames[k]<<"\t";
            }
            for(int k = 0; k < local_numProdsK - 1; k++){
                tmpFile2<<local_prodNames[k]<<"\t";
            }
            tmpFile2<<endl;
            std::copy(local_pi.begin(),local_pi.end(),ostream_iterator<long double>(tmpFile2,"\n"));
            tmpFile2.close();
        }
#endif

        {
            ostringstream tmpFileName;
            ofstream tmpFile; 
            tmpFileName.clear();
            tmpFileName.seekp(0);
            tmpFileName<<"CHECKPOINT_PMAT.txt";
            tmpFile.open(tmpFileName.str().c_str());
            for(int z = 0; z < local_numModelsM2; z++){
                for(int zz = 0; zz < local_numModelsM2; zz++){
                    tmpFile<<local_transP[z][zz]<<" ";
                }
                tmpFile<<endl;
            }
            tmpFile.close();
        }

        {
            ostringstream tmpFileName;
            ofstream tmpFile; 
            tmpFileName.clear();
            tmpFileName.seekp(0);
            tmpFileName<<"CHECKPOINT_MU.txt";
            //tmpFileName<<"simulation_"<<(isim)<<"_CHECKPOINT_MU.txt";
            tmpFile.open(tmpFileName.str().c_str());
            for(int k = 0; k < local_numProdsK - 1; k++){
                tmpFile<<local_prodNames[k]<<"\t";
            }
            for(int k = 0; k < local_numProdsK - 1; k++){
                tmpFile<<local_prodNames[k]<<"\t";
            }
            tmpFile<<endl;
            for(int z = 0; z < local_numModelsM2; z++){
                std::copy(local_mu[z].begin(),local_mu[z].end(),ostream_iterator<long double>(tmpFile,"\t"));
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
                for(int k = 0; k < local_numProdsK - 1; k++){
                    tmpFile2<<local_prodNames[k]<<"\t";
                }
                for(int k = 0; k < local_numProdsK - 1; k++){
                    tmpFile2<<local_prodNames[k]<<"\t";
                }
                tmpFile2<<endl;
                for(int z = 0; z < local_numModelsM2; z++){
                    std::copy(local_mu[z].begin(),local_mu[z].end(),ostream_iterator<long double>(tmpFile2,"\t"));
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
                tmpFileName<<"CHECKPOINT_SIGMASQ.txt";
                //tmpFileName<<"simulation_"<<(isim)<<"_CHECKPOINT_SIGMASQ.txt";
                tmpFile.open(tmpFileName.str().c_str());
                for(int k = 0; k < local_numProdsK - 1; k++){
                    tmpFile<<local_prodNames[k]<<"\t";
                }
                for(int k = 0; k < local_numProdsK - 1; k++){
                    tmpFile<<local_prodNames[k]<<"\t";
                }
                tmpFile<<endl;
                for(int z = 0; z < local_numModelsM2; z++){
                    std::copy(local_sigmaSq[z].begin(),local_sigmaSq[z].end(),ostream_iterator<long double>(tmpFile,"\t"));
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
                tmpFileName2<<"CHECKPOINT_SIGMASQ.txt";
                //tmpFileName<<"simulation_"<<(isim)<<"_CHECKPOINT_SIGMASQ.txt";
                tmpFile2.open(tmpFileName2.str().c_str());
                for(int k = 0; k < local_numProdsK - 1; k++){
                    tmpFile2<<local_prodNames[k]<<"\t";
                }
                for(int k = 0; k < local_numProdsK - 1; k++){
                    tmpFile2<<local_prodNames[k]<<"\t";
                }
                tmpFile2<<endl;
                for(int z = 0; z < local_numModelsM2; z++){
                    std::copy(local_sigmaSq[z].begin(),local_sigmaSq[z].end(),ostream_iterator<long double>(tmpFile2,"\t"));
                    tmpFile2<<endl;
                }
                tmpFile2.close();
            }
#endif

            {
                for(int t = 0; t < local_ntimes; t++)
                {
                    ostringstream tmpFileName;
                    ofstream tmpFile;
                    tmpFileName<<"CHECKPOINT_ZETA_T"<<t<<".txt";
                    tmpFile.open(tmpFileName.str().c_str());
                    tmpFile<<"IMPID"<<"\t"<<"EXPID"<<"\t";
                    for(int z = 0; z < local_numModelsM; z++){
                        tmpFile<<z<<"\t";
                    }
                    tmpFile<<endl;
                    for(int i = 0; i < local_numDyads; i++){
                        RMapType::iterator vit = local_revDyadMap.find(i);
                        tmpFile<<vit->second.first<<"\t"<<vit->second.second<<"\t";
                        for(int z = 0; z < local_numModelsM; z++){
                                tmpFile<<local_zeta[i][t][z]+local_zeta[i][t][z+local_numModelsM]<<"\t";
                        }
                        tmpFile<<endl;
                    }
                    tmpFile.close();
                }
            }

#if 0
            {
                ostringstream tmpFileName;
                ofstream tmpFile; 
                //use the ostringstream object instantiated above to write out the
                //log-likelihood and BIC values for this simulation
                tmpFileName.clear();
                tmpFileName.seekp(0);
                tmpFileName<<"LOGLIKELIHOOD.txt";
                tmpFile.open(tmpFileName.str().c_str());
                tmpFile<<"Final logLikelihood = "<<local_llVal<<" Final BIC = "<<local_bicVal<<endl;
                tmpFile.close(); 
            }
#endif

}

void ChkpBase::chkpClusters(int isim,ofstream &logFile,
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
			    vector<int>& local_eEnd){

            int local_numProdsK2 = (local_numProdsK - 1)*2;
            //int local_numModelsM2 = local_numModelsM*2; 

            vector<vector<int> > maxZeta;
            vector<vector<vector<int> > > clusterOcc;
            vector<vector<int> > nemptyIdx;

            maxZeta.resize(local_ntimes);
            for(int t = 0; t < local_ntimes; t++)
            {
                maxZeta[t].resize(local_numDyads);
            }

            clusterOcc.resize(local_ntimes);
            for(int t = 0; t < local_ntimes; t++)
            {
                clusterOcc[t].resize(local_numModelsM);
            }

            nemptyIdx.resize(local_ntimes);

#if 0
            vector<vector<ofstream> > clusterFiles;
            clusterFiles.resize(local_ntimes);
            for(int t = 0; t < local_ntimes; t++)
            {
                clusterFiles[t].resize(local_numModelsM);
            }
#endif

            int dyadCount = 0;
            double  max_zeta, tmpSum;
            int maxZ;
	    
	    // HJ: use this to track the number of non exist dyads
	    int nonExist = 0;

            for(int t = 0; t < local_ntimes; t++)
            {

	        // HJ: reset the  
	        nonExist = 0;

                for(int i = 0;i < local_numDyads; i++){
		  // HJ: Add the skip for the (dyad, t) combination that does not exist
		  if (local_dExist[i][t] == 1) {
                        max_zeta = -1.0;
                        for(int z = 0; z < local_numModelsM; z++){
                            double tmpSum = local_zeta[i][t][z] + local_zeta[i][t][z+local_numModelsM];
                            if(tmpSum > max_zeta){
                                max_zeta = tmpSum;
                                maxZ = z;
                            }
                        }
                        clusterOcc[t][maxZ].push_back(i);
			// HJ: Why do we need this? It is not used later!
                        // maxZeta[t][i] = maxZ;
		  } else {
		    nonExist += 1;
		  }
                }

                //DONE: ssert that the y vectors are in the correct  order ie
                //zeta[m] > zeta[m+M]
                logFile<<"The following clusters are not empty: "<<endl;
                int tmpDPCount=0;
                for(int z = 0; z < local_numModelsM; z++){
                    if(clusterOcc[t][z].size() != 0){
                        nemptyIdx[t].push_back(z);
                        tmpDPCount += clusterOcc[t][z].size();
                        logFile<<"Final: Cluster "<<z<<" has "<<clusterOcc[t][z].size()<<" dyad-pairs."<<endl;
                    }
                }
                logFile<<"Found "<<nemptyIdx[t].size()<<" non-empty clusters."<<endl;

                // HJ: assert(local_numDyads == tmpDPCount);
                assert(local_numDyads == tmpDPCount + nonExist);

            } //end of loop over t

#if 0
#if SIM_CHKP
            //ofstream* sim_clusterFiles = new ofstream[local_numModelsM*local_ntimes];
            vector<vector<ofstream> > sim_clusterFiles;
            sim_clusterFiles.resize(local_ntimes);
            for(int t = 0; t < local_ntimes; t++)
            {
                sim_clusterFiles[t].resize(local_numModelsM);
            }
#endif
#endif

            for(int t = 0; t < local_ntimes; t++){

                ofstream* clusterFiles = new ofstream[local_numModelsM*local_ntimes];
#if SIM_CHKP
                ofstream* sim_clusterFiles = new ofstream[local_numModelsM];
#endif
                for(int z = 0; z < local_numModelsM; z++){
                    if(clusterOcc[t][z].size() != 0){
                        ostringstream zFile, zFile2;
                        //ostringstream zEXPFile, zPROPFile;
                        //ostringstream flipZFile;
                        zFile<<"DYADS_cluster"<<"_T"<<t<<"_Z"<<z<<".txt";
                        zFile2<<"simulation_"<<(isim)<<"_DYADS_cluster"<<"_T"<<t<<"_Z"<<z<<".txt";
                        clusterFiles[z].open(zFile.str().c_str());
                        clusterFiles[z]<<"imp\texp\t";
#if SIM_CHKP
                        sim_clusterFiles[z].open(zFile2.str().c_str());
                        sim_clusterFiles[z]<<"imp\texp\t";
#endif

                        for(int k = 0; k < local_numProdsK; k++){
                            string currProd = local_prodNames[k];
                            string currProdDesc;
                            ProdMapType::iterator sit = local_ProdMap.find(currProd);
                            if(sit == local_ProdMap.end()){
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
                        for(int k = 0; k < local_numProdsK; k++){
                            string currProd = local_prodNames[k];
                            string currProdDesc;
                            ProdMapType::iterator sit = local_ProdMap.find(currProd);
                            if(sit == local_ProdMap.end()){
                                //cout<<"Product code "<<currProd<<" not found in SITC map."<<endl;
                                currProdDesc = currProd;
                            }
                            else{
                                currProdDesc = sit->second;
                            }
                            if(k == local_numProdsK - 1)
                            {
                                //clusterFiles[z]<<prodNames[k]<<"\t";
                                clusterFiles[z]<<"\""<<currProdDesc<<"\"";
#if SIM_CHKP
                                //sim_clusterFiles[z]<<prodNames[k]<<"\t";
                                sim_clusterFiles[z]<<"\""<<currProdDesc<<"\"";
#endif
                            }
                            else //if k < local_numProdsK - 1
                            {
                                //clusterFiles[z]<<prodNames[k]<<"\t";
                                clusterFiles[z]<<"\""<<currProdDesc<<"\""<<"\t";
#if SIM_CHKP
                                //sim_clusterFiles[z]<<prodNames[k]<<"\t";
                                sim_clusterFiles[z]<<"\""<<currProdDesc<<"\""<<"\t";
#endif
                            }
                        }
                        clusterFiles[z]<<endl;
#if SIM_CHKP
                        sim_clusterFiles[z]<<endl;
#endif

                        for(int i = 0; i < clusterOcc[t][z].size(); i++){
                            int tmpDyadPairID = clusterOcc[t][z][i];
                            RMapType::iterator rit = local_revDyadMap.find(tmpDyadPairID);
                            //int impID = (revDyadMap.find(tmpDyadPairID)->second).first;
                            //int expID = (revDyadMap.find(tmpDyadPairID)->second).second;
                            assert(rit != local_revDyadMap.end());
                            int impID = rit->second.first;
                            int expID = rit->second.second;

                            VMapType::iterator vit = local_vMap.find(make_pair(impID,expID));
                            VMapType::iterator vitPair = local_vMap.find(make_pair(expID,impID));

                            assert(vit != local_vMap.end());
                            assert(vitPair != local_vMap.end());

                            vector<vector<long double> > &mapVal = vit->second;
                            vector<long double> &wijVec = mapVal[t];
                            vector<vector<long double> > &mapValPair = vitPair->second;
                            vector<long double> &wijPairVec = mapValPair[t];

                            if(local_zeta[tmpDyadPairID][t][z] > local_zeta[tmpDyadPairID][t][z+local_numModelsM]){
                                clusterFiles[z]<<impID<<"\t"<<expID<<"\t";
                                //std::copy(wijVec.begin(),wijVec.end(),ostream_iterator<long double>(clusterFiles[z],"\t"));
                                //std::copy(wijPairVec.begin(),wijPairVec.end(),ostream_iterator<long double>(clusterFiles[z],"\t"));
                                for(int i = 0; i < wijVec.size(); i++)
                                {
                                    clusterFiles[z]<<wijVec[i]<<"\t";
                                }
                                for(int i = 0; i < wijPairVec.size() - 1; i++)
                                {
                                    clusterFiles[z]<<wijPairVec[i]<<"\t";
                                }
                                clusterFiles[z]<<wijPairVec[wijPairVec.size() - 1]<<endl;
                                //clusterFiles[z]<<endl;
                            }
                            else{
                                clusterFiles[z]<<expID<<"\t"<<impID<<"\t";
                                //std::copy(wijPairVec.begin(),wijPairVec.end(),ostream_iterator<long double>(clusterFiles[z],"\t"));
                                //std::copy(wijVec.begin(),wijVec.end(),ostream_iterator<long double>(clusterFiles[z],"\t"));
                                for(int i = 0; i < wijPairVec.size(); i++)
                                {
                                    clusterFiles[z]<<wijPairVec[i]<<"\t";
                                }
                                for(int i = 0; i < wijVec.size() - 1; i++)
                                {
                                    clusterFiles[z]<<wijVec[i]<<"\t";
                                }
                                clusterFiles[z]<<wijVec[wijVec.size() - 1]<<endl;
                                //clusterFiles[z]<<endl;
                            }
#if SIM_CHKP
                            if(local_zeta[tmpDyadPairID][t][z] > local_zeta[tmpDyadPairID][t][z+local_numModelsM]){
                                sim_clusterFiles[z]<<impID<<"\t"<<expID<<"\t";
                                std::copy(wijVec.begin(),wijVec.end(),ostream_iterator<long double>(sim_clusterFiles[z],"\t"));
                                std::copy(wijPairVec.begin(),wijPairVec.end(),ostream_iterator<long double>(sim_clusterFiles[z],"\t"));
                                sim_clusterFiles[z]<<endl;
                            }
                            else{
                                sim_clusterFiles[z]<<expID<<"\t"<<impID<<"\t";
                                std::copy(wijPairVec.begin(),wijPairVec.end(),ostream_iterator<long double>(sim_clusterFiles[z],"\t"));
                                std::copy(wijVec.begin(),wijVec.end(),ostream_iterator<long double>(sim_clusterFiles[z],"\t"));
                                sim_clusterFiles[z]<<endl;
                            }
#endif
                        }
                        clusterFiles[z].close();
#if SIM_CHKP
                        sim_clusterFiles[z].close();
#endif
                    }
                } //end of loop over z

                vector<long double> tmpPi(local_numModelsM,0.0);
                long double tmpPiDenominator = 0.0;
                for(int z = 0; z < local_numModelsM; z++)
                {
                        for(int i = 0; i < local_numDyads; i++)
                        {
                                tmpPi[z] += local_zeta[i][t][z];
                                tmpPi[z] += local_zeta[i][t][z+local_numModelsM];
                        }
                        tmpPiDenominator += tmpPi[z];
                }

                for(int z = 0; z < local_numModelsM; z++){
                    if(clusterOcc[t][z].size() != 0){
                        ostringstream zFile;
                        zFile<<"PI_cluster_T"<<t<<"_Z"<<z<<".txt";
                        clusterFiles[z].open(zFile.str().c_str());

                        clusterFiles[z]<<tmpPi[z]/tmpPiDenominator<<endl;

                        clusterFiles[z].close();
                    }
                }

                for(int z = 0; z < local_numModelsM; z++){
                    if(clusterOcc[t][z].size() != 0){
                        ostringstream zFile;
                        zFile<<"Q_cluster_T"<<t<<"_Z"<<z<<".txt";
                        //zFile<<"simulation_"<<(isim)<<"_Q_cluster_"<<z<<".txt";
                        clusterFiles[z].open(zFile.str().c_str());

#if 0
                        for(int k = 0; k < 2*local_numProdsK; k++){
                            clusterFiles[z]<<local_q[z][k]<<endl;
                        }
#endif
                        for(int k = 0; k < local_numProdsK; k++){
                            clusterFiles[z]<<local_prodNames[k]<<"\t"<<local_q[z][k]<<endl;
                        }
                        for(int k = local_numProdsK; k < 2*local_numProdsK; k++){
                            clusterFiles[z]<<local_prodNames[k-local_numProdsK]<<"\t"<<local_q[z][k]<<endl;
                        }

                        clusterFiles[z].close();
                    }
                }

                for(int z = 0; z < local_numModelsM; z++){
                    if(clusterOcc[t][z].size() != 0){
                        ostringstream zFile;
                        zFile<<"MU_cluster_T"<<t<<"_Z"<<z<<".txt";
                        //zFile<<"simulation_"<<(isim)<<"_MU_cluster_"<<z<<".txt";
                        clusterFiles[z].open(zFile.str().c_str());

                        //clusterFiles[z]<<"pi["<<z<<"] = "<<pi[z]<<" pi["<<z+numModelsM<<"] = "<<pi[z+numModelsM]<<" ";
                        //clusterFiles[z]<<"aggregate pi["<<z<<"] = "<<pi[z] + pi[z+numModelsM]<<endl;

#if 0
                        for(int k = 0; k < local_numProdsK2; k++){
                            clusterFiles[z]<<local_mu[z][k]<<endl;
                        }
#endif
                        for(int k = 0; k < local_numProdsK - 1; k++){
                            clusterFiles[z]<<local_prodNames[k]<<"\t"<<local_mu[z][k]<<endl;
                        }
                        for(int k = local_numProdsK - 1; k < 2*(local_numProdsK - 1); k++){
                            clusterFiles[z]<<local_prodNames[k-(local_numProdsK - 1)]<<"\t"<<local_mu[z][k]<<endl;
                        }

                        clusterFiles[z].close();
                    }
                }

                for(int z = 0; z < local_numModelsM; z++){
                    if(clusterOcc[t][z].size() != 0){
                        ostringstream zFile;
                        zFile<<"SIGMA_cluster_T"<<t<<"_Z"<<z<<".txt";
                        //zFile<<"simulation_"<<(isim)<<"_SIGMA_cluster_"<<z<<".txt";
                        clusterFiles[z].open(zFile.str().c_str());
#if 0
                        for(int k = 0; k < local_numProdsK2; k++){
                            clusterFiles[z]<<local_sigmaSq[z][k]<<endl;
                        }
#endif
                        for(int k = 0; k < local_numProdsK - 1; k++){
                            clusterFiles[z]<<local_prodNames[k]<<"\t"<<local_sigmaSq[z][k]<<endl;
                        }
                        for(int k = local_numProdsK - 1; k < 2*(local_numProdsK - 1); k++){
                            clusterFiles[z]<<local_prodNames[k-(local_numProdsK - 1)]<<"\t"<<local_sigmaSq[z][k]<<endl;
                        }

                        clusterFiles[z].close();
                    }
                }
        }// end of loop over t
}
