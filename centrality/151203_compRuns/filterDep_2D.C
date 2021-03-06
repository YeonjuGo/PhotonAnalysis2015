// to draw the ''''correlation plot'''' of centrality variables
// with the different event selection filters in one canvas.
// should modify evtfilter definition part
//
// Author Yeonju Go
// 23 Nov 2015
//
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
#include "yjUtility.h"

//////////// modify here as you want!! ////////////////
const int nfilter = 2;
//const int nfilter = 7;
const char* evtfilter[] = {"","phfCoincFilter3==1 && pprimaryVertexFilter"};
//const char* evtfilter[] = {"","phfCoincFilter3","pprimaryVertexFilter","phfCoincFilter3 && pprimaryVertexFilter"};
//const char* evtfilter[] = {"","pBeamScrapingFilter","pPAprimaryVertexFilter","pHBHENoiseFilterResultProducer","phfCoincFilter","phfCoincFilter3","PAcollisionEventSelection"};
const double dy= 0.5;

void Get2DEffPlots(TTree* t_evt=0, TString v1="hiHF", TString v2="hiNpix", int xbin=200, double xmin=0, double xmax=4500, int ybin=200, double ymin=0, double ymax=10000, TCut cut="", const char* cap="",bool isPassed=1);

void filterDep_2D(
        const char* fname="root://eoscms//eos/cms//store/group/phys_heavyions/velicanu/forest/Run2015E/HIExpressPhysics/Merged//HIForestExpress_run262836.root",
        //const char* fname="root://eoscms//eos/cms//store/group/phys_heavyions/velicanu/forest/Run2015E/HIExpressPhysics/Merged//HIForestExpress_run262836.root",
        //const char* fname="root://eoscms//eos/cms//store/group/phys_heavyions/velicanu/forest/Run2015E/HIExpressPhysics/Merged//HIForestExpress_run262548-v6.root",
        //const char* fname="root://eoscms//eos/cms//store/group/phys_heavyions/velicanu/forest/Run2015E/ExpressPhysics/Merged/ExpressHiForest_run262163-262172_1.4M.root",
        const char* type="Run262836",
        //TCut trig = "HLT_HIL1MinimumBiasHF1AND_v1"
    //TCut trig = "HLT_HIL1MinimumBiasHF1AND_v1"
        TCut trig = "HLT_HIL1HFplusANDminusTH0Express_v1"
        //TCut trig = "HLT_HIL1MinimumBiasHF1ANDExpress_v1"
        //TCut trig = "HLT_L1MinimumBiasHF1OR_part1_v1"
        )
{
    TH1::SetDefaultSumw2();
    gStyle -> SetOptStat(1);
    SetHistTitleStyle(0.06,0.04);
    SetyjPadStyle();

    TFile *fin = TFile::Open(fname);
    TTree *t_evt = (TTree*) fin -> Get("hiEvtAnalyzer/HiTree");
    TTree *t_skim = (TTree*) fin -> Get("skimanalysis/HltTree");
    TTree *t_hlt = (TTree*) fin -> Get("hltanalysis/HltTree");
    t_evt -> AddFriend(t_hlt);
    t_evt -> AddFriend(t_skim);

    double hiNpixMax = 45000;
    double hiNtracksMax= 4000;
    double hiNtracksPtCutMax= 3000;
    double hiNtracksEtaCutMax= 3000;
    double hiNtracksEtaPtCutMax= 2000;
    double hiZDCMax = 80000;
    double hiHFMax = 6000;
    double hiHFhitMax = 150000;
    double hiHFhitPlusMax = 80000;
    double hiBinMax = 200;
    int nbin = 100;
    const char* cap = type;
    //for(int i=0;i<2;i++){
    int i=1;
        Get2DEffPlots(t_evt, "hiHF","hiNpix",nbin,0,hiHFMax,nbin,0,hiNpixMax,trig,cap,i);
        Get2DEffPlots(t_evt, "hiNtracks","hiNpix",nbin,0,hiNtracksMax,nbin,0,hiNpixMax,trig,cap,i);
        Get2DEffPlots(t_evt, "hiNtracks","hiHF",nbin,0,hiNtracksMax,nbin,0,hiHFMax,trig,cap,i);
        //Get2DEffPlots(t_evt, "hiHF","hiBin",nbin,0,hiHFMax,nbin,0,hiBinMax,trig,cap,i);
        
        Get2DEffPlots(t_evt, "hiHFplus","hiHFminus",nbin,0,hiHFMax,nbin,0,hiHFMax,trig,cap,i);
        Get2DEffPlots(t_evt, "hiHFhitPlus","hiHFhitMinus",nbin,0,hiHFhitPlusMax,nbin,0,hiHFhitPlusMax,trig,cap,i);
        //Get2DEffPlots(t_evt, "hiHF","hiZDC",nbin,0,hiHFMax,nbin,0,hiZDCMax,trig,cap,trig,cap,i);
        //Get2DEffPlots(t_evt, "hiNpix","hiZDC",nbin,0,hiNpixMax,nbin,0,hiZDCMax,trig,cap,trig,cap,i);
    //}
}
void Get2DEffPlots(TTree* t_evt, TString v1, TString v2, int xbin, double xmin, double xmax, int ybin, double ymin, double ymax, TCut cut, const char* cap,bool isPassed){
    TCut totcut[nfilter];
    for(int i=0; i<nfilter;i++){
        if(i==0) totcut[i] = cut;
        else totcut[i] = cut && Form("%s==%d",evtfilter[i],(int)isPassed);
    } 

    TCanvas *c_temp = new TCanvas("c_temp", "c_temp", 300,300);
    c_temp->cd();
    TH2D *h2D[nfilter];
    TProfile *prof[nfilter];
    for(int i=0; i<nfilter; i++){
        h2D[i] = new TH2D(Form("h2D_%s_%s_filter%d_%s_isPassed%d",v1.Data(),v2.Data(),i,cap,(int)isPassed), Form("%s;%s;%s",totcut[i].GetTitle(), v1.Data(), v2.Data()), xbin, xmin, xmax, ybin, ymin, ymax);
        t_evt->Draw(Form("%s:%s>>+%s",v2.Data(), v1.Data(), h2D[i]->GetName() ), totcut[i]);
        h2D[i]=(TH2D*)gDirectory->Get(h2D[i]->GetName());
        prof[i] = h2D[i]->ProfileX();
    }
    TCanvas *c_tot = new TCanvas("c_temp", "c_temp", 400*nfilter, 400);
    c_tot->Divide(nfilter,1);
    for(int i=0;i<nfilter;i++){
        c_tot->cd(i+1);
        h2D[i]->Draw("colz");
        gPad->SetLogz();
        prof[i]->Draw("same");
    }
    c_tot->SaveAs(Form("pdf/filterDep_2D_%s_%s_%s_isPassed%d.pdf", v1.Data(),v2.Data(),cap,(int)isPassed));
}

