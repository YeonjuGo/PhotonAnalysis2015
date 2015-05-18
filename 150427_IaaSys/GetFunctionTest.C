// Author : Yeonju Go

//basic c++ header, string ...
#include "../HiForestAnalysis/hiForest.h"
#include "../gammaJetAnalysis/CutAndBinCollection2012.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <TString.h>
#include <string>
#include <math.h>
//tree, hist, vector ...
#include <TROOT.h>
#include "TSystem.h"
#include <TObject.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TNtuple.h>
#include <TMath.h>
#include <math.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include "TClonesArray.h"
#include <TAxis.h>
#include <cmath>
#include <TLorentzRotation.h>
//canvas, legend, latex ... //cosmetic
#include <TCanvas.h>
#include <TStyle.h>
#include <TPaveStats.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLine.h>
#include <TAxis.h>
//random 
#include <TRandom.h>
#include <TStopwatch.h>
#include <ctime>

void GetFunctionTest(){
	TH1::SetDefaultSumw2();
	gStyle -> SetOptStat(0);
	gStyle -> SetTitleYSize(0.05);
    TFile* f1 = new TFile("/home/goyeonju/CMS/2015/gammaJetSystematics/resultHistograms/resultHistograms_ppSmeared10030.root");
    
    TCanvas c1
        ((TH1D *) f1 -> Get(“meanJetPt_pp;1”)) -> Draw()
        TCanvas c2
        ((TH1D *) _file0 -> Get(“meanJetPt_pp;4”)) -> Draw()
    TH1D* h_IaaBin_pp[4];
    TH1D* h_IaaBin_pp_ratio[4];
    TCanvas* c1=new TCanvas("c1","",1300,500);
    c1->Divide(4,1);
    for(int i=0; i<4; i++){
        h_IaaBin_pp[i] = (TH1D*) f1 -> Get(Form("dNdJetPt_pp_ptBin1;%d",i+1));
        h_IaaBin_pp[i]->
    }
 
    for(int i=0; i<4; i++){
        c1 ->cd(i+1);
        h_IaaBin_pp[i]->Draw();
    }

}
