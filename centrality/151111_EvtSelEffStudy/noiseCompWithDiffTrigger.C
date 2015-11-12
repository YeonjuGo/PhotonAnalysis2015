// original macro by Yetkin?
// modified by yeonju
// to get number of the towers?

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

void noiseCompWithDiffTrigger(){

    const char* infile = "/u/user/goyeonju/files/centrality/PbPb_minbias_data_2760_HIRun2011-14Mar2014-v2_run181611_CMSSW5320_byYJ.root";
    const char* trig[]={"HLT_HIZeroBias_v1", "HLT_HIMinBiasZDC_Calo_v1","HLT_HIMinBiasHfOrBSC_v1", "!HLT_HIMinBiasHfOrBSC_v1 && HLT_HIZeroBias_v1"};
    //const char* trig[]={"HLT_HIZeroBias_v1==1", "HLT_HIMinBiasZDC_Calo_v1==1","HLT_HIMinBiasHfOrBSC_v1==1", "HLT_HIMinBiasHfOrBSC_v1!=1 && HLT_HIZeroBias_v1==1"};
    const int nTrig = 4;

    const char* outName = Form("histfiles/hist_noiseCompDiffTrigger.root");

    TFile * inf = new TFile(infile);
    TFile* outf = new TFile(outName,"recreate");

    TH1D* e[nTrig];
    TH1D* et[nTrig];
   
    TH1D* e_tower[nTrig];
    TH1D* et_tower[nTrig];
    for(int i = 0; i < nTrig; ++i){
        e[i] = new TH1D(Form("e%d",i),";HF tower Energy (GeV);",1600,0,800);
        et[i] = new TH1D(Form("et%d",i),";HF tower E_{T} (GeV);",1600,0,800);
        e_tower[i] = new TH1D(Form("e_tower%d",i),";HF tower Energy (GeV);",1600,0,800);
        et_tower[i] = new TH1D(Form("et_tower%d",i),";HF tower E_{T} (GeV);",1600,0,800);
 
    }

    TTree* t1 = (TTree*)inf->Get("hltanalysis/HltTree");
    TTree* t2 = (TTree*)inf->Get("rechitanalyzer/hf");
    TTree* t3 = (TTree*)inf->Get("rechitanalyzer/hbhe");
    TTree* t4 = (TTree*)inf->Get("rechitanalyzer/ee");
    TTree* t5 = (TTree*)inf->Get("rechitanalyzer/eb");
    TTree* t6 = (TTree*)inf->Get("rechitanalyzer/tower");
    TTree* t7 = (TTree*)inf->Get("hiEvtAnalyzer/HiTree");

    t1->AddFriend(t2);
    t1->AddFriend(t3);
    t1->AddFriend(t4);
    t1->AddFriend(t5);
    t1->AddFriend(t6);
    t1->AddFriend(t7);

    for(int i = 0; i < nTrig; ++i){
        t1->Draw(Form("hf.e>>e%d",i),trig[i]);
        t1->Draw(Form("hf.et>>et%d",i),trig[i]);
        t1->Draw(Form("tower.e>>e_tower%d",i),Form("abs(tower.eta)>2.87 && abs(tower.eta)<5.205 && %s", trig[i]));
        t1->Draw(Form("tower.et>>et_tower%d",i),Form("abs(tower.eta)>2.87 && abs(tower.eta)<5.205 && %s", trig[i]));
    }
    outf->Write();
}

