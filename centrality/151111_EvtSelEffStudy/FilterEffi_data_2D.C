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
#include "../../yjUtility.h"

const int Ncut = 5;
void Get2DEffPlots(TTree* t_evt=0, TString v1="hiHF", TString v2="hiNpix", int xbin=200, double xmin=0, double xmax=4500, int ybin=200, double ymin=0, double ymax=10000, TCut cut="", TString cap="",bool isPassed=1);
void FilterEffi_data_2D(bool isMC=0)
{
    const TCut runCut = "run==181611";
    const TCut lumiCut = "lumi>=1 && lumi<=895";
    const TCut eventCut = runCut && lumiCut;
    TH1::SetDefaultSumw2();
    gStyle -> SetOptStat(0);

    TString fname;
    if(isMC) fname = "/u/user/goyeonju/files/centrality/Centrality_officialMC_Hydjet1p8_TuneDrum_Quenched_MinBias_2760GeV.root";
    else fname = "/u/user/goyeonju/files/centrality/PbPb_minbias_data_2760_HIRun2011-14Mar2014-v2_run181611_CMSSW5320_byYJ.root";

    TFile *fin = new TFile(fname.Data());
    TTree *t_evt = (TTree*) fin -> Get("hiEvtAnalyzer/HiTree");
    TTree *t_skim = (TTree*) fin -> Get("skimanalysis/HltTree");
    TTree *t_hlt = (TTree*) fin -> Get("hltanalysis/HltTree");
    t_evt -> AddFriend(t_hlt);
    t_evt -> AddFriend(t_skim);

    Get2DEffPlots(t_evt, "hiHF","hiNpix",100,0,5000,100,0,10000,"HLT_HIMinBiasHfOrBSC_v1==1","HLT_HIMinBiasHfOrBSC_v1",0);
    Get2DEffPlots(t_evt, "hiHF","hiZDC",100,0,5000,100,0,50000,"HLT_HIMinBiasHfOrBSC_v1==1","HLT_HIMinBiasHfOrBSC_v1",0);
    Get2DEffPlots(t_evt, "hiNpix","hiZDC",100,0,10000,100,0,50000,"HLT_HIMinBiasHfOrBSC_v1==1","HLT_HIMinBiasHfOrBSC_v1",0);
    Get2DEffPlots(t_evt, "hiHF","hiNpix",100,0,5000,100,0,10000,"HLT_HIMinBiasHfOrBSC_v1==1","HLT_HIMinBiasHfOrBSC_v1",1);
    Get2DEffPlots(t_evt, "hiHF","hiZDC",100,0,5000,100,0,50000,"HLT_HIMinBiasHfOrBSC_v1==1","HLT_HIMinBiasHfOrBSC_v1",1);
    Get2DEffPlots(t_evt, "hiNpix","hiZDC",100,0,10000,100,0,50000,"HLT_HIMinBiasHfOrBSC_v1==1","HLT_HIMinBiasHfOrBSC_v1",1);

    Get2DEffPlots(t_evt, "hiHF","hiNpix",100,0,5000,100,0,10000,"!HLT_HIMinBiasHfOrBSC_v1 && HLT_HIZeroBias_v1","noHfOrBSC_passZero",0);
    Get2DEffPlots(t_evt, "hiHF","hiZDC",100,0,5000,100,0,50000,"!HLT_HIMinBiasHfOrBSC_v1 && HLT_HIZeroBias_v1","noHfOrBSC_passZero",0);
    Get2DEffPlots(t_evt, "hiNpix","hiZDC",100,0,10000,100,0,50000,"!HLT_HIMinBiasHfOrBSC_v1 && HLT_HIZeroBias_v1","noHfOrBSC_passZero",0);
    Get2DEffPlots(t_evt, "hiHF","hiNpix",100,0,5000,100,0,10000,"!HLT_HIMinBiasHfOrBSC_v1 && HLT_HIZeroBias_v1","noHfOrBSC_passZero",1);
    Get2DEffPlots(t_evt, "hiHF","hiZDC",100,0,5000,100,0,50000,"!HLT_HIMinBiasHfOrBSC_v1 && HLT_HIZeroBias_v1","noHfOrBSC_passZero",1);
    Get2DEffPlots(t_evt, "hiNpix","hiZDC",100,0,10000,100,0,50000,"!HLT_HIMinBiasHfOrBSC_v1 && HLT_HIZeroBias_v1","noHfOrBSC_passZero",1);


    Get2DEffPlots(t_evt, "hiHF","hiNpix",100,0,5000,100,0,10000,"!HLT_HIMinBiasHfOrBSC_v1 && HLT_HIMinBiasHf_OR_v1","noHfOrBSC_passHfOr",0);
    Get2DEffPlots(t_evt, "hiHF","hiZDC",100,0,5000,100,0,50000,"!HLT_HIMinBiasHfOrBSC_v1 && HLT_HIMinBiasHf_OR_v1","noHfOrBSC_passHfOr",0);
    Get2DEffPlots(t_evt, "hiNpix","hiZDC",100,0,10000,100,0,50000,"!HLT_HIMinBiasHfOrBSC_v1 && HLT_HIMinBiasHf_OR_v1","noHfOrBSC_passHfOr",0);
    Get2DEffPlots(t_evt, "hiHF","hiNpix",100,0,5000,100,0,10000,"!HLT_HIMinBiasHfOrBSC_v1 && HLT_HIMinBiasHf_OR_v1","noHfOrBSC_passHfOr",1);
    Get2DEffPlots(t_evt, "hiHF","hiZDC",100,0,5000,100,0,50000,"!HLT_HIMinBiasHfOrBSC_v1 && HLT_HIMinBiasHf_OR_v1","noHfOrBSC_passHfOr",1);
    Get2DEffPlots(t_evt, "hiNpix","hiZDC",100,0,10000,100,0,50000,"!HLT_HIMinBiasHfOrBSC_v1 && HLT_HIMinBiasHf_OR_v1","noHfOrBSC_passHfOr",1);

}
void Get2DEffPlots(TTree* t_evt, TString v1, TString v2, int xbin, double xmin, double xmax, int ybin, double ymin, double ymax, TCut cut, TString cap,bool isPassed){
    TCut totcut[Ncut];
    totcut[0] = cut;
    totcut[1] = cut&& Form("pprimaryVertexFilter==%d",(int)isPassed);
    totcut[2] = cut&& Form("phltPixelClusterShapeFilter==%d",(int)isPassed);
    totcut[3] = cut&& Form("phfCoincFilter3==%d",(int)isPassed);
    totcut[4] = cut&& Form("pcollisionEventSelection==%d",(int)isPassed);

    TCanvas *c_temp = new TCanvas("c_temp", "c_temp", 300,300);
    c_temp->cd();
    TH2D *h2D[Ncut];
    for(int i=0; i<Ncut; i++){
        h2D[i] = new TH2D(Form("h2D_%s_%s_filter%d_%s_isPassed%d",v1.Data(),v2.Data(),i,cap.Data(),(int)isPassed), Form("%s;%s;%s",totcut[i].GetTitle(), v1.Data(), v2.Data()), xbin, xmin, xmax, ybin, ymin, ymax);
        t_evt->Draw(Form("%s:%s>>+%s",v2.Data(), v1.Data(), h2D[i]->GetName() ), totcut[i]);
        h2D[i]=(TH2D*)gDirectory->Get(h2D[i]->GetName());
    }
    TCanvas *c_tot = new TCanvas("c_temp", "c_temp", 900, 600);
    c_tot->Divide(3,2);
    for(int i=0;i<Ncut;i++){
        c_tot -> cd(i+1);
        h2D[i] -> Draw("colz");
    }
    c_tot->SaveAs(Form("pdf/h2D_%s_isPassed%d_%s_%s.pdf", cap.Data(),(int)isPassed, v1.Data(), v2.Data()));
}

