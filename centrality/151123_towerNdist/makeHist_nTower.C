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

TH1* makeHist_nTower(float eThr=0, float eThr=0)
{

    TH1::SetDefaultSumw2();
    gStyle -> SetOptStat(0);

    TCut dataCut("HLT_HIL1MinimumBiasHF1AND_v1 && phfCoincFilter3");
    TCut mcCut("");

    // ===================================================================================
    // Get Trees from data & mc files.
    // ===================================================================================
    TFile *f[2]; // 0:data, 1:MC

    TTree *t_skim[2];
    TTree *t_hlt[2];
    TTree *t[2];
    int Entries[2];   
    const char* fname_data="root://eoscms//eos/cms/store/group/phys_heavyions/velicanu/forest/Run2015E/HIExpressPhysics/Merged/HiForest_Streamer_run262315.root";
    const char* fname_mc="/afs/cern.ch/work/y/ygo/public/centrality/HiForest_Centrality_Unpacker_Hydjet_Quenched_MinBias_5020GeV_750_RECODEBUG_v0_merged_forest_6.root";

    for(int i=0;i<2;i++){ 
        if(i==0) f[0] = TFile::Open(fname_data);
        else f[1] = TFile::Open(fname_mc);
        t_skim[i] = (TTree*) f[i] -> Get("skimanalysis/HltTree");
        t_hlt[i] = (TTree*) f[i] -> Get("hltanalysis/HltTree");
        t[i] = (TTree*) f[i] -> Get("rechitanalyzer/tower");
        t[i] -> AddFriend(t_hlt[i]);
        t[i] -> AddFriend(t_skim[i]);
        if(i==0) Entries[i] = t[i] -> GetEntries(dataCut);
        else Entries[i] = t[i] -> GetEntries();
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
