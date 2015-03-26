// Author : Yeonju Go
// using forest file
// compare the isolation variables

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
#include <TCanvas.h>
#include <TStyle.h>
#include <TLine.h>
#include <TAxis.h>
//random 
#include <TRandom.h>
#include <TStopwatch.h>
#include <ctime>

void drawLongText(float x1, float y1, float x2, float y2, const char* head="", int col=kBlack);
void iso_data_mc(TString isoName="trkSumPtHollowConeDR04"){
	TH1::SetDefaultSumw2();
	gStyle -> SetOptStat(0);
	gStyle -> SetTitleYSize(0.04);

	//in this case, not mc but old file
	TString treeName = "multiPhotonAnalyzer/photon";
	TFile* fpho_mc = new TFile("/u/user/goyeonju/files/forest/pA/pA_photonSkimForest_v85_HLT_PAPhoton30_NoCaloIdVL_v1_highPtPhoton40.root");
	TTree* tpho_mc = (TTree*) fpho_mc->Get(treeName.Data());

	TFile* fpho_data= new TFile("/u/user/goyeonju/files/forest/pA/pPb_DATA_photon30trig_localJEC_v1.root");
	TTree* tpho_data= (TTree*) fpho_data->Get(treeName.Data());

//photon -> Draw("trkSumPtHollowConeDR04", "pt>40 && abs(eta)<1.44 && hadronicOverEm<0.1 && sigmaIphiIphi > 0.002 && sigmaIetaIeta > 0.002 && swissCrx < 0.9 && abs(seedTime) < 3 && rawEnergy/energy>0.5")
	const int nPhoCut = 3;
	const int nTotCut = 3;
	const TCut photonEtaCut = "abs(eta)<1.44";
	const TCut spikeCut = "sigmaIphiIphi > 0.002 && sigmaIetaIeta > 0.002 && swissCrx < 0.9 && abs(seedTime) < 3";
	const TCut hoeCut = "hadronicOverEm<0.1";
	const TCut resCut = "rawEnergy/energy>0.5";
	TCut phoEtCut[nPhoCut] = {"pt>40","pt>60","pt>80"} ;
 	const TCut candidateCut = "sigmaIetaIeta<0.01";
  	const TCut decayCut = "(sigmaIetaIeta>0.011) && (sigmaIetaIeta<0.017)";
	TString str[nPhoCut] = {"pt>40", "pt>60", "pt>80"};
	
	
	TCut totCut_data[nTotCut];
	totCut_data[0] = photonEtaCut;
	totCut_data[1] = photonEtaCut && spikeCut && hoeCut && resCut;
	totCut_data[2] = photonEtaCut && spikeCut && hoeCut && resCut && candidateCut;

	TCut totCut_mc[nTotCut];
	totCut_mc[0] = photonEtaCut;
	totCut_mc[1] = photonEtaCut && spikeCut && hoeCut && resCut;
	totCut_mc[2] = photonEtaCut && spikeCut && hoeCut && resCut && candidateCut;

	TH1D* hiso_data[nTotCut][nPhoCut];
	TH1D* hiso_mc[nTotCut][nPhoCut];
	int nbin, xmin, xmax;
	if(isoName=="trkSumPtHollowConeDR04") { nbin=100;xmin=0;xmax=13000; }
	else { nbin=100;xmin=0;xmax=300; }

	TCanvas* c_temp = new TCanvas("c_temp", "jet p_{T} distribution", 300, 300); 
	c_temp->cd();
	for(int i=0;i<nTotCut;i++){
		for(int j=0;j<nPhoCut;j++){
			hiso_data[i][j] = new TH1D(Form("hiso_data_tot%d_pho%d",i,j), Form("h%s_data_tot%d_pho%d;%s;",isoName.Data(),i,j,isoName.Data()), nbin,xmin,xmax);
			hiso_mc[i][j] = new TH1D(Form("hiso_mc_tot%d_pho%d",i,j),  Form("h%s_mc_tot%d_pho%d;%s;",isoName.Data(),i,j,isoName.Data()), nbin,xmin,xmax);
		}	
	}

	for(int i=0;i<nTotCut;i++){
		for(int j=0;j<nPhoCut;j++){
			tpho_data->Project(hiso_data[i][j]->GetName(),isoName.Data(), totCut_data[i] && phoEtCut[j]); 
			tpho_mc->Project(hiso_mc[i][j]->GetName(),isoName.Data(), totCut_data[i] && phoEtCut[j]); 
		}
	}	

	TCanvas* c_cut = new TCanvas("c_cut", Form("%s distribution",isoName.Data()), 1200, 700); 
	c_cut->Divide(3,3);

	TLegend *l1 = new TLegend(0.3365615,0.6445304,0.7577623,0.846736,NULL,"brNDC");
	easyLeg(l1);
	for(int i=0;i<nTotCut;i++){
		for(int j=0;j<nPhoCut;j++){
			c_cut->cd(3*i+j+1);
			hiso_data[i][j]->Scale(1./hiso_data[i][j]->Integral());
			hiso_mc[i][j]->Scale(1./hiso_mc[i][j]->Integral());
			//hist cosmetics
			double range = cleverRange(hiso_data[i][j],hiso_mc[i][j]);
                        hiso_data[i][j]->GetYaxis()->SetRangeUser(0.00001,range);
                        hiso_mc[i][j]->GetYaxis()->SetRangeUser(0.00001,range);
			gPad->SetLogy();
			hiso_data[i][j]->SetMarkerStyle(20);
			hiso_data[i][j]->SetMarkerSize(0.4);
			hiso_mc[i][j]->SetLineColor(4);
			hiso_mc[i][j]->Draw("hist e");
			hiso_data[i][j]->Draw("same p");
			if(i==0 && j==0) {
				l1->AddEntry(hiso_data[0][0],"2015 data","p");
				l1->AddEntry(hiso_mc[0][0],"2013 data","l");
				l1->Draw("same");
			}
			drawText(str[j],0.6,0.75);			
			/*	
				if(j==0){
				TLegend* leg = new TLegend(0.3,0.3,0.9,0.9);
				leg->AddEntry((TObject*)0,Form("%s",totCut_data[i].GetTitle()),"");
				leg->SetTextColor(kBlack);
				leg->Draw();
				}	
				*/
			//if(j==0)drawLongText(0.5,0.5,0.7,0.7,Form("%s",totCut_data[i].GetTitle()),kBlack);
			//if(j==1)drawLongText(0.5,0.5,0.7,0.7,Form("%s",totCut_mc[i].GetTitle()),kBlue);
		}
	}
//	c_cut -> Update();
	c_cut ->SaveAs(Form("pdf/%s.pdf",isoName.Data()));
}

void drawLongText(float x1, float y1, float x2, float y2, const char* head, int col){
	TLegend* leg = new TLegend(x1,y1,x2,y2,head);
	leg->SetTextColor(col);
	leg->SetTextSize(15);
	leg->Draw();
}
	
