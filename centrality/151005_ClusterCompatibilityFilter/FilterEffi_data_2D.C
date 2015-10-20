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
//#include "../HiForestAnalysis/hiForest.h"

const int Ncut = 5;
void Get2DEvtPlots(TTree* t_evt=0, TString v1="hiHF", TString v2="hiNpix", int xbin=200, double xmin=0, double xmax=4500, int ybin=200, double ymin=0, double ymax=10000, TCut cut="", TCanvas* c_tot=0, TString cap="", bool isPassed=0, bool isAOD=0);
void FilterEffi_data_2D(const char* fname="/afs/cern.ch/work/y/ygo/public/centrality/merged_centrality_MB_DATA750_RECO_150914.root", TString type = "RECO_750",bool isAOD=0)
{
    //pcollisionEventSelection
    //phfCoincFilter
    //phfCoincFilter3
    //pprimaryVertexFilter
    //phltPixelClusterShapeFilter
    //phiEcalRecHitSpikeFilter
    //HLT_HIMinBiasHfOrBSC_v1
	const TCut runCut = "run==181611";
	const TCut lumiCut = "lumi>=1 && lumi<=895";
	const TCut eventCut = runCut && lumiCut;
	TH1::SetDefaultSumw2();
	gStyle -> SetOptStat(0);

	TFile* fin = new TFile(fname);
	TTree* t_evt = (TTree*) fin -> Get("hiEvtAnalyzer/HiTree");
	TTree* t_skim = (TTree*) fin -> Get("skimanalysis/HltTree");
	TTree* t_hlt = (TTree*) fin -> Get("hltanalysis/HltTree");
	t_evt -> AddFriend(t_hlt);
	t_evt -> AddFriend(t_skim);

	TCanvas *c_tot = new TCanvas("c_tot", "c_tot", 900,600);
	c_tot->Divide(3,2);

	TLine* t1 = new TLine(0,1,1000,1);
	t1->SetLineWidth(1);
	t1->SetLineStyle(7); // 7 is jumSun , 1 is onSun
	t1->SetLineColor(1); // 2 is red
	Get2DEvtPlots(t_evt, "hiBin", "hiHF", 100,0,200, 100,0,5000,"",c_tot,type,0,(int)isAOD);
	Get2DEvtPlots(t_evt, "hiBin", "hiHFhit", 100,0,200, 100,0,50000,"",c_tot,type,0,(int)isAOD);
	Get2DEvtPlots(t_evt, "hiBin", "hiNpix", 100,0,200, 100,0,40000,"",c_tot,type,0,(int)isAOD);
	Get2DEvtPlots(t_evt, "hiBin", "hiZDC", 100,0,200, 100,0,80000,"",c_tot,type,0,(int)isAOD);
	Get2DEvtPlots(t_evt, "hiBin", "hiEE", 100,0,200, 100,0,3000,"",c_tot,type,0,(int)isAOD);
	Get2DEvtPlots(t_evt, "hiBin", "hiEB", 100,0,200, 100,0,4000,"",c_tot,type,0,(int)isAOD);
	Get2DEvtPlots(t_evt, "hiBin", "hiET", 100,0,200, 100,0,2000,"",c_tot,type,0,(int)isAOD);
	Get2DEvtPlots(t_evt, "hiHF", "hiNpix", 100,0,3000, 100,0,40000,"",c_tot,type,0,(int)isAOD);
	Get2DEvtPlots(t_evt, "hiEE", "hiEB", 100,0,4000, 100,0,4000,"",c_tot,type,0,(int)isAOD);

    Get2DEvtPlots(t_evt, "hiBin", "hiHF", 100,0,200, 100,0,5000,"",c_tot,type,1,(int)isAOD);
	Get2DEvtPlots(t_evt, "hiBin", "hiHFhit", 100,0,200, 100,0,50000,"",c_tot,type,1,(int)isAOD);
	Get2DEvtPlots(t_evt, "hiBin", "hiNpix", 100,0,200, 100,0,40000,"",c_tot,type,1,(int)isAOD);
	Get2DEvtPlots(t_evt, "hiBin", "hiZDC", 100,0,200, 100,0,80000,"",c_tot,type,1,(int)isAOD);
	Get2DEvtPlots(t_evt, "hiBin", "hiEE", 100,0,200, 100,0,3000,"",c_tot,type,1,(int)isAOD);
	Get2DEvtPlots(t_evt, "hiBin", "hiEB", 100,0,200, 100,0,4000,"",c_tot,type,1,(int)isAOD);
	Get2DEvtPlots(t_evt, "hiBin", "hiET", 100,0,200, 100,0,2000,"",c_tot,type,1,(int)isAOD);
	Get2DEvtPlots(t_evt, "hiHF", "hiNpix", 100,0,3000, 100,0,40000,"",c_tot,type,1,(int)isAOD);
	Get2DEvtPlots(t_evt, "hiEE", "hiEB", 100,0,4000, 100,0,4000,"",c_tot,type,1,(int)isAOD);
}
void Get2DEvtPlots(TTree* t_evt, TString v1, TString v2, int xbin, double xmin, double xmax, int ybin, double ymin, double ymax, TCut cut, TCanvas* c_tot, TString cap, bool isPassed, bool isAOD )
{
    TCut totcut[Ncut];
    totcut[0] = cut;
    totcut[1] = cut&& Form("pprimaryVertexFilter==%d",(int)isPassed);
    if(isAOD) totcut[2] = cut&& Form("pclusterCompatibilityFilter==%d",(int)isPassed);
    else totcut[2] = cut&& Form("phltPixelClusterShapeFilter==%d",(int)isPassed);
    totcut[3] = cut&& Form("phfCoincFilter3==%d",(int)isPassed);
    totcut[4] = cut&& Form("pcollisionEventSelection==%d",(int)isPassed);

    TCanvas *c_temp = new TCanvas("c_temp", "c_temp", 300,300);
    c_temp->cd();
    TH2D *h2D[10];
    for(int i=0; i<Ncut; i++){
        h2D[i] = new TH2D(Form("h2D_%d",i), Form(";%s;%s", v1.Data(), v2.Data()), xbin, xmin, xmax, ybin, ymin, ymax);
    }
    for(int i=0; i<Ncut; i++){
        t_evt->Draw(Form("%s:%s>>+%s",v2.Data(), v1.Data(), h2D[i]->GetName() ), totcut[i]);
        h2D[i]=(TH2D*)gDirectory->Get(h2D[i]->GetName());
    }

    TLegend* leg = new TLegend(0.4,0.2,0.95,0.6);
    //leg->
    for(int i=0;i<Ncut;i++){
        c_tot -> cd(i+1);
        h2D[i] -> Draw("colz");
    }
    c_tot->SaveAs(Form("png/h2D_%s_%s_isPassed%d_%s.png",v1.Data(), v2.Data(),(int)isPassed,cap.Data()));
    for(int i=0; i<Ncut; i++){
        delete h2D[i];
    }
    delete c_temp;
}

