///////////////////////////////////////////////////////////////////                                
// forest2yskim.C                                                //                                                 
// Creator : Yongsun Kim (MIT), jazzitup@mit.edu                 //                                                 
// Function : Transform hiForest files into yskim file           //
// yskims for MinBias1, Minbias2 and photon jet skims            //
///////////////////////////////////////////////////////////////////         

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TFile.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <TMath.h>
#include "../../HiForestAnalysis/hiForest.h"
//#include "../gammaJetAnalysis/CutAndBinCollection2012.h"
#include <time.h>


using namespace std;

static const long MAXTREESIZE = 10000000000;





void lumiCutter(int runN = 181611, int lumiFrom =1, int lumiTo =10) {
  //collisionType
  collisionType colli = cPbPb;
  TString inputFile_="/home/goyeonju/CMS/Files/centrality/HiForest_PbPb_minbias_DATA_20141011_53X_byKisoo.root";
  TString outname = Form("/home/goyeonju/CMS/Files/centrality/HiForest_PbPb_minbias_DATA_20141011_53X_byKisoo_run%d_lumi%d_%d.root",runN,lumiFrom,lumiTo);
  
  
  // start from here
  // path length histogram
  
  HiForest * t;
  t = new HiForest(inputFile_.Data(),"",colli);
   
  t->SetOutputFile(outname.Data());
  // LOOP!!
  t->InitTree();
    cout << "sdlf" <<endl;
  int nentries = t->GetEntries();
  //  nentries = 100000;
  cout << "number of entries = " << nentries << endl;
  for (Long64_t jentry = 0 ; jentry < nentries; jentry++) {
    if (jentry% 10000 == 0)  {
      cout <<jentry<<" / "<<nentries<<" "<<setprecision(2)<<(double)jentry/nentries*100<<endl;
    }
    t->GetEntry(jentry);
    
    if (t->hlt.HLT_HIMinBiasHfOrBSC_v1==0) 
        continue;
    if ( !( (t->evt.run = runN) && (t->evt.lumi>lumiFrom) && (t->evt.lumi<lumiTo) ) )
      continue;
    t->FillOutput();
  }
  delete t;
}

