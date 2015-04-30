// Author Yeonju Go
#include "TFile.h"
#include "TTree.h"
#include "TLine.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TMath.h"
#include "TDirectory.h"
#include "TProfile.h"
#include "TPad.h"
#include "stdio.h"
#include "../HiForestAnalysis/hiForest.h"

const int Ncut = 5;
Get2DEffPlots(TTree* t_evt=0, TStrinb* v1="hiHF", TString* v2="hiNpix", int xbin=200, double xmin=0, double xmax=4500, int ybin=200, double ymin=0, double ymax=10000, TCut* cut="", TCanvas* c_tot=0, TString* cap="");
void FilterEffi_data_2D()
{
	/*    const int onlineFilter = 0;//HLT_HIMinBiasHfOrBSC_v1
	      const int collCut = 0; //pcollisionEventSelection
	      const int vertexCut = 0; //pprimaryVertexFilter
	      const int pixelCut = 0; //phltPixelClusterShapeFilter
	      const int hfCoincCut = 0; //phfCoincFilter3
	      */
	const TCut runCut = "run==181611";
	const TCut lumiCut = "lumi>=1 && lumi<=895";
	const TCut eventCut = runCut && lumiCut;
	TH1::SetDefaultSumw2();
	gStyle -> SetOptStat(0);

	TString *fname = "/u/user/goyeonju/PRODUCTION/CMSSW_5_3_20/src/centralityMC/forest_officialMC/merge/Centrality_officialMC_Hydjet1p8_TuneDrum_Quenched_MinBias_2760GeV.root";

	TFile *fin = new TFile(fname.Data());
	TTree *t_evt = (TTree*) fin -> Get("hiEvtAnalyzer/HiTree");
	TTree *t_skim = (TTree*) fin -> Get("skimanalysis/HltTree");
	TTree *t_hlt = (TTree*) fin -> Get("hltanalysis/HltTree");
	t_evt -> AddFriend(t_hlt);
	t_evt -> AddFriend(t_skim);

	TCanvas *c_tot = new TCanvas("c_tot", "c_tot", 900,600);
	c_tot->Divide(3,2);

	TLine* t1 = new TLine(0,1,1000,1);
	t1->SetLineWidth(1);
	t1->SetLineStyle(7); // 7 is jumSun , 1 is onSun
	t1->SetLineColor(1); // 2 is red
	Get2DEffPlots(t_evt, "hiHF"
}
Get2DEffPlots(TTree* t_evt, TString* v1, TString* v2, int xbin, double xmin, double xmax, int ybin, double ymin, double ymax, TCut* cut, TCanvas* c_tot, TString* cap);
	TCut* totcut[Ncut];
	totcut[0] = cut;
	totcut[1] = cut && "pprimaryVertexFilter==1";
	totcut[2] = cut&& "phltPixelClusterShapeFilter==1";
	totcut[3] = cut&& "phfCoincFilter3==1";
	totcut[4] = cut&& "pcollisionEventSelection==1";
 /*    const int onlineFilter = 0;//HLT_HIMinBiasHfOrBSC_v1
 *                  const int collCut = 0; //pcollisionEventSelection
 *                                const int vertexCut = 0; //pprimaryVertexFilter
 *                                              const int pixelCut = 0; //phltPixelClusterShapeFilter
 *                                                            const int hfCoincCut = 0; //phfCoincFilter3
 *                                                                          */

	TCanvas *c_temp = new TCanvas("c_temp", "c_temp", 300,300);
	c_temp->cd();
	TH2D *h2D[5];
	for(int i=0; i<Ncut; i++){
		h2D[i] = new TH2D(Form("h2D_%d",i), Form(";%s;%s", v1.Data(), v2.Data()), xbin, xmin, xmax, ybin, ymin, ymax);
	}
	for(int i=0; i<Ncut; i++){
		t_evt->Draw(Form("%s:%s>>+%s",v2.Data(), v1.Data(), h2D[i]->GetName() ), totcut[i]);
		h2D[i]=(TH2D*)gDirectory->Get(h2D[i]->GetName());
	}
	for(int i=0;i<Ncut;i++){
                c_tot -> cd(i+1);
                h2D[i] -> Draw("colz");
        }
	c_tot->SaveAs(Form("pdf/h2D_%s_%s_%s",v1.Data(), v2.Data(), cap.Data()));
}
 
