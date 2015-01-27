#include <memory>
#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
#include <map>
#include <TF1.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TNtuple.h>
#include <TChain.h>
#include <TFile.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TRandom.h>
#include <TGraphErrors.h>

#define BIGNUMBER 1000
#define MIDNUMBER 100

#include "DataFormats/HeavyIonEvent/interface/CentralityBins.h"

using namespace std;

double addQuad(double *x,int n){
   
   double sum = 0;
   for(int i = 0; i < n; ++i){
      sum+= x[i]*x[i];
   }
   return sqrt(sum);
 
}

void printError(double *var0,double uncSmear = 0.002,int uncEff = 1);

void analyze_systematics(){

   int nTable = 19;
   int nBins = 40;

   // Analysis bins defined by the PInG
   vector<vector<int> > bins;
   vector<int> b0, b1, b2;
   for(int i = 0; i < 4; ++i) b0.push_back(i);
   for(int i = 4; i < 12; ++i) b1.push_back(i);
   for(int i = 12; i < 40; ++i) b2.push_back(i);
   bins.push_back(b0);
   bins.push_back(b1);
   bins.push_back(b2);

   CentralityBins* table[20];

   string tags[200] = {

      "Eff100_sim00_g0",
      "Eff099_sim00_g0",
      "Eff098_sim00_g0",
      "Eff097_sim00_g0",
      "Eff096_sim00_g0",
      "Eff095_sim00_g0",

      "Eff100_sim01_g0",
      "Eff100_sim10_g0",
      "Eff100_sim11_g0",

      "Eff100_sim00_g1a",
      "Eff100_sim00_g1b",
      "Eff100_sim00_g2a",
      "Eff100_sim00_g2b",
      "Eff100_sim00_g3a",
      "Eff100_sim00_g3b",
      "Eff100_sim00_g4a",
      "Eff100_sim00_g4b",
      "Eff100_sim00_g5a",
      "Eff100_sim00_g5b"

   };


   string names[100] = {

      "Eff 100","Eff 99", "Eff 98", "Eff 97", "Eff 96", "Eff 95",
      "Changing SIM","Changing GEN","Changing GEN+SIM",
      "Nuclear Radius","Nuclear Radius",
      "Skin Depth %2","Skin Depth %2","Skin Depth %10","Skin Depth %10",
      "d Min","d Min",
      "Sigma NN","Sigma NN"
   };


   TFile* inf = new TFile("glaubertest_d0107_v0.root");

   CentralityBins* standardTable = (CentralityBins*)inf->Get("Eff100_sim00_g0/run1");
   double Ncoll[20] = {0,0,0,0,0,0,0,0,0,0,0,0,0};
   double Npart[20] = {0,0,0,0,0,0,0,0,0,0,0,0,0};
   double NpartCP2 = 0;
   double NcollCP2 = 0;


   double var0[100], var1[100],var2[100],varcp2[100];

   for(int t = 0; t < bins.size(); ++t){
      for(int j = 0; j < bins[t].size(); ++j){
	 //         Ncoll[t] += standardTable->NcollMeanOfBin((bins[t])[j]);
	 Ncoll[t] += standardTable->NpartMeanOfBin((bins[t])[j]);
      }
      Ncoll[t] /= bins[t].size();
   }

   NpartCP2 = (Npart[0])/(Npart[2]);
   NcollCP2 = (Ncoll[0])/(Ncoll[2]);

   cout<<"Ncoll 0-10%    :  "<<Ncoll[0]<<endl;
   cout<<"Ncoll 30-100%  :  "<<Ncoll[2]<<endl;
   cout<<"Ncoll 0-10% / Ncoll 30-100%: "<<NcollCP2<<endl;
   
   for(int i = 0; i < nTable; ++i){
      cout<<"Analyzing variation : "<<names[i].data()<<endl;
      table[i] = (CentralityBins*)inf->Get(Form("%s/run1",tags[i].data()));
      
      double ncoll[20] ={0,0,0,0,0,0,0,0,0,0,0,0};
      for(int t = 0; t < bins.size(); ++t){
	 for(int j = 0; j < bins[t].size(); ++j){
	    //            ncoll[t] += table[i]->NcollMeanOfBin((bins[t])[j]);
	    ncoll[t] += table[i]->NpartMeanOfBin((bins[t])[j]);
	 }
	 ncoll[t] /= bins[t].size();
      }

      cout<<"Results : "<<ncoll[0]<<endl;

      var0[i] = ((ncoll[0])-(Ncoll[0]))/Ncoll[0];
      var1[i] = ((ncoll[1])-(Ncoll[1]))/Ncoll[1];
      var2[i] = (ncoll[2]-Ncoll[2])/Ncoll[2];

      cout<<"Ncoll 0-10%    :  "<<ncoll[0]<<" Variation : "<<var0[i]*100<<"%"<<endl;
      cout<<"Ncoll 10-30%    :  "<<ncoll[1]<<" Variation : "<<var1[i]*100<<"%"<<endl;
      cout<<"Ncoll 30-100%  :  "<<ncoll[2]<<" Variation : "<<var2[i]*100<<"%"<<endl;
      double ncollcp1 = (ncoll[0])/(ncoll[1]);
      double ncollcp2 = ncoll[0]/ncoll[2];
      varcp2[i] = (ncollcp2 - NcollCP2)/NcollCP2;
      cout<<"Ncoll 0-10% / Ncoll 30-100%: "<<ncollcp2<<" Variation : "<<varcp2[i]*100<<"%"<<endl;

   }

   cout<<"0 - 10%"<<endl;
   printError(var0,0.0005);
   cout<<"10 - 30%"<<endl;

   printError(var1,0.005);
   cout<<"30 - 100%"<<endl;
   printError(var2,0.005);

   cout<<"0 - 10% / 30 -100%"<<endl;
   printError(varcp2);


   cout<<"Custom Calculations : "<<endl;
   if(1){
      double a[100] = {0.2,   7.4,   0.36,   2.7,   3 };
      cout<<addQuad(a,5)<<endl;
      
      double b[100] = {3,5.7,0.17,8.4,1.2};
      cout<<addQuad(b,5)<<endl;
      
      double c[100] = {3,1.8,0.5,6,2};
      cout<<addQuad(c,5)<<endl;
   }else{
      
      double a[100] = {0.1,   .65,   0.,   0.3,   0.2 };
      cout<<addQuad(a,5)<<endl;

      double b[100] = {2.4, 1.1, 0.1, 3.5, .2};
      cout<<addQuad(b,5)<<endl;

      double c[100] = { 2.3, .8, 0.1, 3.3 ,0.2};
      cout<<addQuad(c,5)<<endl;

   }




}

void printError(double *var0, double uncSmear, int uncEff){
   
   double glauberArray0[20] = {
      (fabs(var0[9])+fabs(var0[10]))/2,
      (fabs(var0[13])+fabs(var0[14]))/2,
      (fabs(var0[15])+fabs(var0[16]))/2,
      (fabs(var0[17])+fabs(var0[18]))/2,
   };
   
   double glauberError0 = addQuad(glauberArray0,4);
   double finalArray0[20] = {
      (fabs(var0[2-uncEff])+fabs(var0[2+uncEff]))/2,
      glauberError0,
      uncSmear
   };
   double error0 = addQuad(finalArray0,3);
   cout<<"Glauber Error : "<<glauberError0<<endl;
   cout<<"Final Error   : "<<error0<<endl;
  
}   


