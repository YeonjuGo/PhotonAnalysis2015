// Author Yeonju Go
// last modification : 2015/01/27 
// 
#include "TFile.h"
#include "TTree.h"
#include "TLine.h"
#include "TH1F.h"
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
#include "TLatex.h"
#include "stdio.h"
#include "../../yjUtility.h"

void makeHist_nTower(float etThr=0.0, float eThr=1.0)
{
    const int Ncut = 5;
    int col[] = {1,2,4,6,8,28,46,41};
    int marker[] = {20,22,29,33,34};
    string FilterName[Ncut] = {"no evt filter","pprimaryVertexFilter","phltPixelClusterShapeFilter", "phfCoincFilter3", "pcollisionEventSelection"};

    TH1::SetDefaultSumw2();
    gStyle -> SetOptStat(0);

    // ===================================================================================
    // Get Trees from data & mc files.
    // ===================================================================================

    TFile *f[2]; // 0:data, 1:MC

    TTree *t_evt[2]; 
    TTree *t_skim[2];
    TTree *t_hlt[2];
    TTree *t_recTower[2];
    TTree *t_hf[2];
    int Entries[2];   
    for(int i=0;i<2;i++){ 
 //        if(i==0) f[0] = TFile::Open("root://cluster142.knu.ac.kr//store/user/ygo/PbPb_minbias_data_2760_HIRun2011-14Mar2014-v2_run181611.root");
 //       else f[1] = TFile::Open("root://cluster142.knu.ac.kr//store/user/ygo/officialMC_Hydjet1p8_TuneDrum_Quenched_MinBias_2760GeV.root");
//        if(i==0) f[0] = TFile::Open("/u/user/goyeonju/files/centrality/PbPb_minbias_data_2760_HIRun2011-14Mar2014-v2_run181611_CMSSW5320_byYJ.root");
//        else f[1] = TFile::Open("/u/user/goyeonju/files/centrality/Centrality_officialMC_Hydjet1p8_TuneDrum_Quenched_MinBias_2760GeV.root");
        if(i==0) f[0] = TFile::Open("/mnt/hadoop/cms/store/user/dgulhan/hiForest_HIMinBiasUPC_HIRun2011-v1_7_5_3_patch1/HiForest_HIMinBiasUPC_HIRun2011-v1_merged.root");
        else f[1] = TFile::Open("/mnt/hadoop/cms/store/user/dgulhan/hiForest_HydjetMB_2076GeV_FOREST_753p1_merged/HydjetMB_2076GeV_FOREST_753p1_v0_merged.root");
        t_evt[i] = (TTree*) f[i] -> Get("hiEvtAnalyzer/HiTree");
        t_skim[i] = (TTree*) f[i] -> Get("skimanalysis/HltTree");
        t_hlt[i] = (TTree*) f[i] -> Get("hltanalysis/HltTree");
        t_recTower[i] = (TTree*) f[i] -> Get("rechitanalyzer/tower");
        t_hf[i] = (TTree*) f[i] -> Get("rechitanalyzer/hf");
        t_recTower[i] -> AddFriend(t_hlt[i]);
        t_recTower[i] -> AddFriend(t_evt[i]);
        t_recTower[i] -> AddFriend(t_skim[i]);
        t_recTower[i] -> AddFriend(t_hf[i]);
        Entries[i] = t_recTower[i] -> GetEntries();
    }

    // ===============================================================================================
    // [et] Define et histograms (data/mc , pf tree/rechit tree) 
    // ===============================================================================================
    cout << "LET'S DEFINE HISTOGRAMS" << endl;

    TH1F* h1F_sumHF_data[10]; //[0:data, 1:mc][condition]
    TH1F* h1F_sumHF_mc; //[0:data, 1:mc][condition]

    for(int j=0; j<Ncut; j++){
        h1F_sumHF_data[j] = new TH1F(Form("sumHF_data_filter%d",j), ";# of HF towers;Events",200,0,1000);
        //h1F_sumHF_data[j] -> SetMarkerStyle(marker[i]);
        h1F_sumHF_data[j] -> SetMarkerSize(0.9);
        //h1F_sumHF_data[j] -> SetMarkerColor(col[i]);
        //h1F_sumHF_data[j] -> SetLineColor(col[i]);
        h1F_sumHF_data[j] -> SetLabelSize(0.03);
        //ratio_rec_n[i] -> SetAxisRange(0.5,1.5,"Y");
    }
    h1F_sumHF_mc = (TH1F*) h1F_sumHF_data[0]->Clone("sumHF_mc");

    // ===============================================================================================
    // [rechit] Fill histogram by using <Draw> function in TTree.
    // ===============================================================================================
    cout << "LET'S FILL HISTOGRAMS FROM TREE" << endl;
    //********************************************************
    //cut define
    TString towerCut = Form("hf.et > %.1f && hf.e > %.1f",etThr,eThr);
    //TString towerCut = Form("abs(tower.eta) > 2.87 && abs(tower.eta) < 5.02 && tower.et > %.1f && tower.e > %.1f",etThr,eThr);
    TCut cut="";
    TCut totcut[Ncut];
    totcut[0] = cut;
    totcut[1] = cut&& Form("pprimaryVertexFilter==1");
    //if(isAOD) totcut[2] = cut&& Form("pclusterCompatibilityFilter==1");
    totcut[2] = cut&& Form("phltPixelClusterShapeFilter==1");
    totcut[3] = cut&& Form("phfCoincFilter3==1");
    totcut[4] = cut&& Form("pcollisionEventSelection==1");
    //********************************************************
    for(int j=0; j<Ncut; j++){
        t_recTower[0]->Draw(Form("Sum$(%s)>>+%s",towerCut.Data(),h1F_sumHF_data[j]->GetName()), totcut[j] && "HLT_HIMinBiasHfOrBSC_v1==1");
        //cout << "i = " << i << ", j = " << j << " :::: finished " << endl;
        h1F_sumHF_data[j] = (TH1F*)gDirectory->Get(h1F_sumHF_data[j]->GetName());
    }
    t_recTower[1]->Draw(Form("Sum$(%s)>>+%s",towerCut.Data(),h1F_sumHF_mc->GetName()));
    h1F_sumHF_mc = (TH1F*)gDirectory->Get(h1F_sumHF_mc->GetName());

    TFile* outFile = new TFile(Form("histfiles/sumHF_etThr%.1f_eThr%.1f.root",etThr,eThr),"recreate");
    outFile->cd();
    for(int j=0; j<Ncut; j++){
        h1F_sumHF_data[j]->Write();
    }
    h1F_sumHF_mc->Write();
    outFile->Close();
} 
