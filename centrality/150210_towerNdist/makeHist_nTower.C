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

void makeHist_nTower(float etThr=1.5)
{
    const TCut threCut = Form("abs(eta)>2.87 && abs(eta)<5.2 && et>%.1f",etThr);
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
    TTree *t_pfTower[2];
    int Entries[2];   
    for(int i=0;i<2;i++){ 
        if(i==0) f[0] = TFile::Open("root://cluster142.knu.ac.kr//store/user/ygo/PbPb_minbias_data_2760_HIRun2011-14Mar2014-v2_run181611.root");
        else f[1] = TFile::Open("root://cluster142.knu.ac.kr//store/user/ygo/officialMC_Hydjet1p8_TuneDrum_Quenched_MinBias_2760GeV.root");
        t_evt[i] = (TTree*) f[i] -> Get("hiEvtAnalyzer/HiTree");
        t_skim[i] = (TTree*) f[i] -> Get("skimanalysis/HltTree");
        t_hlt[i] = (TTree*) f[i] -> Get("hltanalysis/HltTree");
        t_recTower[i] = (TTree*) f[i] -> Get("rechitanalyzer/tower");
        t_pfTower[i] = (TTree*) f[i] -> Get("pfTowers/tower");
        t_recTower[i] -> AddFriend(t_hlt[i]);
        t_recTower[i] -> AddFriend(t_evt[i]);
        t_recTower[i] -> AddFriend(t_skim[i]);
        t_pfTower[i] -> AddFriend(t_hlt[i]);
        t_pfTower[i] -> AddFriend(t_evt[i]);
        t_pfTower[i] -> AddFriend(t_skim[i]);
        Entries[i] = t_recTower[i] -> GetEntries();
    }

    // ===============================================================================================
    // [et] Define et histograms (data/mc , pf tree/rechit tree) 
    // ===============================================================================================
    cout << "LET'S DEFINE HISTOGRAMS" << endl;

    TH1F* h1F_rec_nTower[3][10]; //[0:data, 1:mc][condition]

    for(int i=0; i<2; i++){
        for(int j=0; j<Ncut; j++){
            h1F_rec_nTower[i][j] = new TH1F(Form("rec_nTower_sample%d_filter%d",i,j), ";# of HF towers;Events",200,0,1000);
            h1F_rec_nTower[i][j] -> SetMarkerStyle(marker[i]);
            h1F_rec_nTower[i][j] -> SetMarkerSize(0.9);
            h1F_rec_nTower[i][j] -> SetMarkerColor(col[i]);
            h1F_rec_nTower[i][j] -> SetLineColor(col[i]);
            h1F_rec_nTower[i][j] -> SetLabelSize(0.03);
            //ratio_rec_n[i] -> SetAxisRange(0.5,1.5,"Y");
        }
    }

    // ===============================================================================================
    // [rechit] Fill histogram by using <Draw> function in TTree.
    // ===============================================================================================
    cout << "LET'S FILL HISTOGRAMS FROM TREE" << endl;
    //********************************************************
    //cut define
    TString towerCut = Form("abs(tower.eta) > 2.87 && abs(tower.eta)<5.205 && tower.e > %.1f",etThr);
    TCut cut="";
    TCut totcut[Ncut];
    totcut[0] = cut;
    totcut[1] = cut&& Form("pprimaryVertexFilter==1");
    //if(isAOD) totcut[2] = cut&& Form("pclusterCompatibilityFilter==1");
    totcut[2] = cut&& Form("phltPixelClusterShapeFilter==1");
    totcut[3] = cut&& Form("phfCoincFilter3==1");
    totcut[4] = cut&& Form("pcollisionEventSelection==1");
    //********************************************************
    std::clock_t    start_loop, end_loop;
    start_loop = std::clock();


    for(int i=0; i<2; i++){
        for(int j=0; j<Ncut; j++){
            std::clock_t    startT_draw, endT_draw;
            startT_draw = std::clock();
            cout << "i = " << i << ", j = " << j << endl;
            //if(i==0) t_recTower[i]->Draw(Form("Sum$(abs(tower.eta) > 2.87 && tower.et > 0.5)>>+%s",h1F_rec_nTower[i][j]->GetName()));
            if(i==0) t_recTower[i]->Draw(Form("Sum$(%s)>>+%s",towerCut.Data(),h1F_rec_nTower[i][j]->GetName()), totcut[j] && "HLT_HIMinBiasHfOrBSC_v1==1");
            else if(i==1) t_recTower[i]->Draw(Form("Sum$(%s)>>+%s",towerCut.Data(),h1F_rec_nTower[i][j]->GetName()));
            h1F_rec_nTower[i][j] = (TH1F*)gDirectory->Get(h1F_rec_nTower[i][j]->GetName());
            endT_draw = std::clock();
            std::cout.precision(6);      // get back to default precision
            std::cout << "DRAW finished in             : " << (endT_draw - startT_draw) / (double)CLOCKS_PER_SEC << " seconds" << std::endl;
        }
    }

    end_loop = std::clock();
    std::cout.precision(6);      // get back to default precision
    std::cout << "LOOP finished in             : " << (end_loop - start_loop) / (double)CLOCKS_PER_SEC << " seconds" << std::endl;
    std::cout << "exited loop" << std::endl;

    TFile* outFile = new TFile(Form("recTower_n_etThr%.1f.root",etThr),"recreate");
    outFile->cd();
    for(int j=0; j<Ncut; j++){
        for(int i=0; i<2; i++){
            h1F_rec_nTower[i][j]->Write();
        }
    }
    outFile->Close();
} 
