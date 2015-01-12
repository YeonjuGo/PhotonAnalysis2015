//Author : Yeonju Go
//to derive pthat weight
//


#include "TROOT.h"
#include "TString.h"
#include <iostream>
#include "TFile.h"
#include "TSystem.h"
#include "TPad.h"
#include "TH1D.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include <stdlib.h>
using namespace std;

const int filepthat[]={30,50,80,120,170};
const int nFile = 5;

void xSectionCal(const char *ksp="pp"){

	//pp    
	/*    const char *fileName[nFile] = {
	      "/u/user/goyeonju/PRODUCTION/CMSSW_5_3_20/src/pp_JEC/merging-forest/ppAllQCD15/JEC_ppAllQCD15.root",        
	      "/u/user/goyeonju/PRODUCTION/CMSSW_5_3_20/src/pp_JEC/merging-forest/ppAllQCD30/JEC_ppAllQCD30.root",        
	      "/u/user/goyeonju/PRODUCTION/CMSSW_5_3_20/src/pp_JEC/merging-forest/ppAllQCD50/JEC_ppAllQCD50.root",        
	      "/u/user/goyeonju/PRODUCTION/CMSSW_5_3_20/src/pp_JEC/merging-forest/ppAllQCD80/JEC_ppAllQCD80.root",        
	      "/u/user/goyeonju/PRODUCTION/CMSSW_5_3_20/src/pp_JEC/merging-forest/ppAllQCD120/JEC_ppAllQCD120.root"        
	      }; 
	      */
	// pPb FINAL MC produced by Alex 
	const char *fileName[nFile] = {
		"/u/user/goyeonju/files/forest/pA/forJEC/HiForest_pPb_MIX_AllQCDPhoton30.root",
		"/u/user/goyeonju/files/forest/pA/forJEC/HiForest_pPb_MIX_AllQCDPhoton50.root",
		"/u/user/goyeonju/files/forest/pA/forJEC/HiForest_pPb_MIX_AllQCDPhoton80.root",
		"/u/user/goyeonju/files/forest/pA/forJEC/HiForest_pPb_MIX_AllQCDPhoton120.root",
		"/u/user/goyeonju/files/forest/pA/forJEC/HiForest_pPb_MIX_AllQCDPhoton170.root",
	};

	// pPb MC produced by yeonju
	/*    const char *fileName[nFile] = {
	      "/u/user/goyeonju/PRODUCTION/CMSSW_5_3_20/src/pPb_JEC/merging-forest/pPbAllQCD15/JEC_pPbAllQCD15.root",        
	      "/u/user/goyeonju/PRODUCTION/CMSSW_5_3_20/src/pPb_JEC/merging-forest/pPbAllQCD30/JEC_pPbAllQCD30.root",        
	      "/u/user/goyeonju/PRODUCTION/CMSSW_5_3_20/src/pPb_JEC/merging-forest/pPbAllQCD50/JEC_pPbAllQCD50.root",        
	      "/u/user/goyeonju/PRODUCTION/CMSSW_5_3_20/src/pPb_JEC/merging-forest/pPbAllQCD80/JEC_pPbAllQCD80.root",        
	      "/u/user/goyeonju/PRODUCTION/CMSSW_5_3_20/src/pPb_JEC/merging-forest/pPbAllQCD120/JEC_pPbAllQCD120.root"        
	      };
	      */
	// Pbp (reverse direction)
	/*    const char *fileName[nFile] = {
	      "/u/user/goyeonju/PRODUCTION/CMSSW_5_3_20/src/Pbp_JEC/merging-forest/PbpAllQCD15/JEC_PbpAllQCD15.root",        
	      "/u/user/goyeonju/PRODUCTION/CMSSW_5_3_20/src/Pbp_JEC/merging-forest/PbpAllQCD30/JEC_PbpAllQCD30.root",        
	      "/u/user/goyeonju/PRODUCTION/CMSSW_5_3_20/src/Pbp_JEC/merging-forest/PbpAllQCD50/JEC_PbpAllQCD50.root",        
	      "/u/user/goyeonju/PRODUCTION/CMSSW_5_3_20/src/Pbp_JEC/merging-forest/PbpAllQCD80/JEC_PbpAllQCD80.root",        
	      "/u/user/goyeonju/PRODUCTION/CMSSW_5_3_20/src/Pbp_JEC/merging-forest/PbpAllQCD120/JEC_PbpAllQCD120.root"        
	      };
	      */

	TFile *fin[nFile];
	TTree *jetTree[nFile];
	float pEnt[nFile][nFile];
	float totalEnt[nFile];
	for(int i=0; i<nFile ; i++){
		totalEnt[i]=0.0;
		for(int j=0; j<nFile ; j++){
			pEnt[i][j]=0.0;
		}
	}

	for(int ifile=0; ifile<nFile; ifile++){
		fin[ifile] = new TFile(fileName[ifile]);
		jetTree[ifile] = (TTree*) fin[ifile] -> Get("akPu3PFJetAnalyzer/t");
		Float_t pthat;
		TBranch *b_pthat; 
		jetTree[ifile]->SetBranchAddress("pthat",&pthat, &b_pthat);
		totalEnt[ifile] = jetTree[ifile] -> GetEntries();
		for(int ipthat=0; ipthat<nFile; ipthat++){
			if(ipthat==nFile-1) pEnt[ifile][ipthat] = jetTree[ifile] -> GetEntries(Form("pthat>=%d && pthat<9999", filepthat[ipthat]));
			else pEnt[ifile][ipthat] = jetTree[ifile] -> GetEntries(Form("pthat>=%d && pthat<%d", filepthat[ipthat], filepthat[ipthat+1]));
		}
		cout << jetTree[ifile] << ", : " << pEnt[ifile][0] << ", " << pEnt[ifile][1] <<  ", " << pEnt[ifile][2] <<  ", " << pEnt[ifile][3] <<  ", " << pEnt[ifile][4] <<  ", " << ", total Entries = " << totalEnt[ifile] << endl;
	}

	//
	// Calculate weighting factor
	//
	double weight[nFile];
	for(int ipart=0; ipart<nFile; ipart++){
		weight[ipart] = (pEnt[0][ipart])/(pEnt[0][ipart]+pEnt[1][ipart]+pEnt[2][ipart]+pEnt[3][ipart]+pEnt[4][ipart]+pEnt[5][ipart]);
		cout << ipart << "th partition weight : " << weight[ipart] << endl;
	}

	cout << "devide the range of pthat15 and each entries : " << endl;
	cout << "pthat "<<filepthat[0]<<" - "<<filepthat[1]<<" : " << pEnt[0][0]+pEnt[0][1]+pEnt[0][2]+pEnt[0][3]+pEnt[0][4] << endl;
	cout << "pthat "<<filepthat[1]<<" - "<<filepthat[2]<<" : " << pEnt[0][1]+pEnt[0][2]+pEnt[0][3]+pEnt[0][4]  << endl;
	cout << "pthat "<<filepthat[2]<<" - "<<filepthat[3]<<" : " << pEnt[0][2]+pEnt[0][3]+pEnt[0][4]  << endl;
	cout << "pthat "<<filepthat[3]<<" - "<<filepthat[4]<<" : " << pEnt[0][3]+pEnt[0][4] << endl;
	cout << "pthat "<<filepthat[4]<<" - 9999 : " << pEnt[0][4] << endl;

	TH1D *hpthat[nFile];
	for(int ifile=0; ifile<nFile; ifile++){
		hpthat[ifile] = new TH1D(Form("hpthat%d",ifile),"", 210, 10, 1000);
		jetTree[ifile]-> SetWeight(weight[ifile]);
		jetTree[ifile]-> Draw("pthat>>+hpthat[ifile]");
		hpthat[ifile] = (TH1D*)gDirectory->Get("hpthat[ifile]");
		// hpthat = (TH1F*)gPad->GetPrimitive("hpthat[ifile]");
	}
	/*
	   TCanvas *can = new TCanvas("can", "can", 300,300);

	   for(int ifile=0; ifile<nFile; ifile++){
	   if(ifile==0) hpthat[0] -> Draw();
	   else  hpthat[ifile] -> Draw("same");
	   cout << hpthat[ifile] -> GetBinContent(20) << endl; 
	   }
	   */
}
