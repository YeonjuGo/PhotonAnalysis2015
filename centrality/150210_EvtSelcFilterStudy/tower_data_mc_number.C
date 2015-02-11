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
#include "../../HiForestAnalysis/hiForest.h"
#include "../../gammaJetAnalysis/commonUtility.h"

void tower_data_mc_number()
{
    const TCut runCut = "run==181611";
    const TCut lumiCut = "lumi>=1 && lumi<=895";
    // const TCut eventCut = "";
    const TCut eventCut = runCut && lumiCut;
    const TCut threCut = "abs(eta)>2.87 && abs(eta)<5.2 && et>1.4";
    const int Ncut = 5;
    double color[Ncut] = {12,8,9,41,46};
    double marker[Ncut] = {20,22,29,33,34};
    string FilterName[Ncut] = {"no evt filter","pprimaryVertexFilter","phltPixelClusterShapeFilter", "phfCoincFilter3", "pcollisionEventSelection"};

    TH1::SetDefaultSumw2();
    gStyle -> SetOptStat(0);

    TLatex* latex = new TLatex();
    latex->SetNDC();
    latex->SetTextAlign(12);
    latex->SetTextSize(0.04);
    // ===================================================================================
    // Get Trees from data & mc files.
    // ===================================================================================

    TFile *dataf = new TFile("/home/goyeonju/CMS/Files/centrality/HiForest_PbPb_minbias_DATA_byYJ.root");
    //TFile *dataf = new TFile("/home/goyeonju/CMS/Files/centrality/HiForest_PbPb_minbias_DATA_20141011_53X_byKisoo.root");
    TTree *datat_evt = (TTree*) dataf -> Get("hiEvtAnalyzer/HiTree");
    TTree *datat_skim = (TTree*) dataf -> Get("skimanalysis/HltTree");
    TTree *datat_hlt = (TTree*) dataf -> Get("hltanalysis/HltTree");
    TTree *datat_recTower = (TTree*) dataf -> Get("rechitanalyzer/tower");
    TTree *datat_pfTower = (TTree*) dataf -> Get("pfTowers/tower");
    datat_recTower -> AddFriend(datat_hlt);
    datat_recTower -> AddFriend(datat_evt);
    datat_recTower -> AddFriend(datat_skim);
    datat_pfTower -> AddFriend(datat_hlt);
    datat_pfTower -> AddFriend(datat_evt);
    datat_pfTower -> AddFriend(datat_skim);
    cout << datat_evt << ", " << datat_skim << ", " << datat_hlt << ", " << datat_recTower << ", " << datat_pfTower << ", " << endl;
    cout << "data pfTower tree Friends : " << datat_pfTower -> GetListOfFriends() << endl;
    cout << "data rechitTower tree Friends : " << datat_recTower -> GetListOfFriends() << endl;
    int Nrec_dataEntries = datat_recTower -> GetEntries();

    TFile *mcf = new TFile("/home/goyeonju/CMS/Files/centrality/centrality_PbPb_minbias_MC_53X_byYJ.root");
    //TFile *mcf = new TFile("/home/goyeonju/CMS/Files/centrality/HiForest_HydjetMB_730_53XBS_merged.root");
    TTree *mct_evt = (TTree*) mcf -> Get("hiEvtAnalyzer/HiTree");
    TTree *mct_skim = (TTree*) mcf -> Get("skimanalysis/HltTree");
    TTree *mct_hlt = (TTree*) mcf -> Get("hltanalysis/HltTree");
    TTree *mct_recTower = (TTree*) mcf -> Get("rechitanalyzer/tower");
    TTree *mct_pfTower = (TTree*) mcf -> Get("pfTowers/tower");
    mct_recTower -> AddFriend(mct_hlt);
    mct_recTower -> AddFriend(mct_evt);
    mct_recTower -> AddFriend(mct_skim);
    mct_pfTower -> AddFriend(mct_hlt);
    mct_pfTower -> AddFriend(mct_evt);
    mct_pfTower -> AddFriend(mct_skim);
    int Nrec_mcEntries = mct_recTower -> GetEntries();

    // int Nevt_mct = mct_evt -> GetEntries();
    // cout << "# of MC events = " << Nevt_mct << endl;

    // ***********************************************************************************************************    

    // ===============================================================================================
    // [et] Define et histograms (data/mc , pf tree/rechit tree) 
    // ===============================================================================================
    cout << "LET'S DEFINE HISTOGRAMS" << endl;

    TH1F* data_rec_n[5];
    TH1F* mc_rec_n[5];
    TH1F* ratio_rec_n[5];// ratio = data/mc

    for(int i=0; i<Ncut; i++)
    {
        data_rec_n[i] = new TH1F(Form("data_rec_n%d",i), ";# of HF towers;Events",100,0,1000);
        data_rec_n[i] -> SetMarkerStyle(20);
        data_rec_n[i] -> SetMarkerSize(0.9);
        data_rec_n[i] -> SetMarkerColor(color[i]);
        data_rec_n[i] -> SetLineColor(color[i]);
        data_rec_n[i] -> SetLabelSize(0.03);

        mc_rec_n[i] = (TH1F*)data_rec_n[i]->Clone(Form("mc_rec_n%d",i));
        //mc_rec_n[i] = new TH1F(Form("mc_rec_n%d",i), ";# of HF towers;Events",200,0,1000);
        mc_rec_n[i] -> SetMarkerStyle(marker[i]);
        mc_rec_n[i] -> SetMarkerSize(0.9);
        mc_rec_n[i] -> SetMarkerColor(color[i]); //marker color
        mc_rec_n[i] -> SetLineColor(color[i]); //line color
        mc_rec_n[i] -> SetLabelSize(0.03);

        ratio_rec_n[i] = (TH1F*)data_rec_n[i]->Clone(Form("mc_rec_n%d",i));
        //ratio_rec_n[i] = new TH1F(Form("ratio_rec_n%d",i), ";# of HF towers;DATA/MC",200,0,1000);
        ratio_rec_n[i] -> SetMarkerStyle(20);
        ratio_rec_n[i] -> SetMarkerSize(0.5);
        ratio_rec_n[i] -> SetMarkerColor(color[i]);
        ratio_rec_n[i] -> SetLineColor(color[i]);
        ratio_rec_n[i] -> SetLabelSize(0.03);
        //ratio_rec_n[i] -> SetAxisRange(0.5,1.5,"Y");

    }

    // ===============================================================================================
    // [rechit] Fill histogram by using <Draw> function in TTree.
    // ===============================================================================================
    cout << "LET'S FILL HISTOGRAMS FROM TREE" << endl;

    datat_recTower -> Draw("Sum$(abs(tower.eta) > 2.87 && tower.et > 0.5)>>+ data_rec_n0",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1");
    data_rec_n[0] = (TH1F*)gDirectory->Get("data_rec_n0");
    datat_recTower -> Draw("Sum$(abs(tower.eta) > 2.87 && tower.et > 0.5)>>+ data_rec_n1",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && pprimaryVertexFilter==1");
    data_rec_n[1] = (TH1F*)gDirectory->Get("data_rec_n1");
    datat_recTower -> Draw("Sum$(abs(tower.eta) > 2.87 && tower.et > 0.5)>>+ data_rec_n2",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && phltPixelClusterShapeFilter==1");
    data_rec_n[2] = (TH1F*)gDirectory->Get("data_rec_n2");
    datat_recTower -> Draw("Sum$(abs(tower.eta) > 2.87 && tower.et > 0.5)>>+ data_rec_n3",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && phfCoincFilter3==1");
    data_rec_n[3] = (TH1F*)gDirectory->Get("data_rec_n3");
    datat_recTower -> Draw("Sum$(abs(tower.eta) > 2.87 && tower.et > 0.5)>>+ data_rec_n4",eventCut && "HLT_HIMinBiasHfOrBSC_v1==1 && pcollisionEventSelection==1");
    data_rec_n[4] = (TH1F*)gDirectory->Get("data_rec_n4");

    mct_recTower -> Draw("Sum$(abs(tower.eta) > 2.87 && tower.et > 0.5)>>+ mc_rec_n0");
    mc_rec_n[0] = (TH1F*)gDirectory->Get("mc_rec_n0");
    mct_recTower -> Draw("Sum$(abs(tower.eta) > 2.87 && tower.et > 0.5)>>+ mc_rec_n1","pprimaryVertexFilter==1");
    mc_rec_n[1] = (TH1F*)gDirectory->Get("mc_rec_n1");
    mct_recTower -> Draw("Sum$(abs(tower.eta) > 2.87 && tower.et > 0.5)>>+ mc_rec_n2","phltPixelClusterShapeFilter==1");
    mc_rec_n[2] = (TH1F*)gDirectory->Get("mc_rec_n2");
    mct_recTower -> Draw("Sum$(abs(tower.eta) > 2.87 && tower.et > 0.5)>>+ mc_rec_n3","phfCoincFilter3==1");
    mc_rec_n[3] = (TH1F*)gDirectory->Get("mc_rec_n3");
    mct_recTower -> Draw("Sum$(abs(tower.eta) > 2.87 && tower.et > 0.5)>>+ mc_rec_n4", "pcollisionEventSelection==1");
    mc_rec_n[4] = (TH1F*)gDirectory->Get("mc_rec_n4");

    for(int i=0; i<Ncut; i++){
        int cutBinFrom = data_rec_n[i]->FindBin(200); 
        int cutBinTo = data_rec_n[i]->FindBin(700); 

        data_rec_n[i] -> Scale(mc_rec_n[i]->Integral(cutBinFrom,cutBinTo)/data_rec_n[i]->Integral(cutBinFrom,cutBinTo));
        ratio_rec_n[i] -> Divide(data_rec_n[i],mc_rec_n[i]);
        cout << "total integral of MC in the whole range: " << mc_rec_n[i]->Integral()<< endl;
//        cout << "total integral of MC in the cutted range : " << mc_rec_n[i]->Integral(cutBinFrom,cutBinTo)<< endl;
    }


// ===============================================================================================
// [# of HF tower] Draw histograms in Canvas. 
// ===============================================================================================
    cout << "LET'S DRAW HISTOGRAMS IN CANVAS" << endl;
    TCanvas *c_tot = new TCanvas("c_tot","c_tot", 1500,350);
    makeMultiPanelCanvas(c_tot,5,1,0.0,0.0,0.2,0.15,0.02);
    for(int i=0; i<Ncut; i++){
         c_tot->cd(i+1);
         mc_rec_n[i]->DrawCopy("hist");
         data_rec_n[i]->Draw("same");
         if(i==0)latex->DrawLatex(0.46,0.88,"rechitTowers");
         latex->DrawLatex( 0.46, 0.78, Form("%s", FilterName[i].c_str()) );
    }

    TCanvas *c_mc = new TCanvas("c_mc", "c_mc", 400,800);
//    makeMultiPanelCanvas(c_mc,1,2,0.0,0.0,0.2,0.15,0.02);
    c_mc -> Divide(1,2);

    c_mc->cd(1);
    gPad->SetLogy();
    mc_rec_n[0]->GetXaxis()->SetRangeUser(0.,250.);
    mc_rec_n[0] -> Draw("ep");
    mc_rec_n[1] -> Draw("same ep");
    mc_rec_n[2] -> Draw("same ep");
    mc_rec_n[3] -> Draw("same ep");
    mc_rec_n[4] -> Draw("same ep");
  //  latex->DrawLatex(0.46,0.68,"rechitTowers MC");
     
    c_mc->cd(2);
    gPad->SetLogy();
    TH1F* hFilterRatio = (TH1F*)mc_rec_n[0]->Clone("hFilterRatio");
    hFilterRatio->GetYaxis()->SetTitle("Filter Efficiency");
    for(int i=1;i<Ncut;i++){
        hFilterRatio->Reset();
        hFilterRatio->GetXaxis()->SetRangeUser(0.,250.);
        hFilterRatio->Divide(mc_rec_n[i],mc_rec_n[0]);
        hFilterRatio->SetMarkerColor(color[i]);
        hFilterRatio->SetMarkerStyle(marker[i]);
        if(i==1) hFilterRatio->DrawCopy();
        hFilterRatio->DrawCopy("same");
    }

    TLegend *leg = new TLegend(0.6,0.6,0.9,0.9);
    leg->SetFillColor(0);
    leg->SetTextFont(43);
    leg->SetTextSize(15);
    leg->SetHeader("rechitTowers MC");
    leg->AddEntry(mc_rec_n[1],Form("%s",FilterName[1].c_str()),"p");
    leg->AddEntry(mc_rec_n[2],Form("%s",FilterName[2].c_str()),"p");
    leg->AddEntry(mc_rec_n[3],Form("%s",FilterName[3].c_str()),"p");
    leg->AddEntry(mc_rec_n[4],Form("%s",FilterName[4].c_str()),"p");
    leg->Draw();

    c_mc -> SaveAs("pdf/rectower_mc_n.pdf");
    c_tot-> SaveAs("pdf/rectower_data_mc_n_tot.pdf");
    TFile* outFile = new TFile("outfile.root","recreate");
    outFile->cd();
    for(int i=0;i<Ncut;i++){
           mc_rec_n[i]->Write();
           data_rec_n[i]->Write();
           ratio_rec_n[i]->Write();
    }
    outFile->Close();
 } 
