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
//#include "TPad.h"
#include "stdio.h"
#include "../HiForestAnalysis/hiForest.h"

void FilterEffi()
{
    /*    const int onlineFilter = 0;//HLT_HIMinBiasHfOrBSC_v1
          const int collCut = 0; //pcollisionEventSelection
          const int vertexCut = 0; //pprimaryVertexFilter
          const int pixelCut = 0; //phltPixelClusterShapeFilter
          const int hfCoincCut = 0; //phfCoincFilter3
          */
    const int Ncut = 5;
    TH1::SetDefaultSumw2();
    gStyle -> SetOptStat(0);

    TFile *dataf = new TFile("PbPb_Centrality_data.root");
    TTree *datat_evt = (TTree*) dataf -> Get("hiEvtAnalyzer/HiTree");
    TTree *datat_skim = (TTree*) dataf -> Get("skimanalysis/HltTree");
    TTree *datat_hlt = (TTree*) dataf -> Get("hltanalysis/HltTree");
    datat_evt -> AddFriend(datat_hlt);
    datat_evt -> AddFriend(datat_skim);

    TFile *mcf = new TFile("PbPb_Centrality_mc.root");
    TTree *mct_evt = (TTree*) mcf -> Get("hiEvtAnalyzer/HiTree");
    TTree *mct_skim = (TTree*) mcf -> Get("skimanalysis/HltTree");
    TTree *mct_hlt = (TTree*) mcf -> Get("hltanalysis/HltTree");
    mct_evt -> AddFriend(mct_hlt);
    mct_evt -> AddFriend(mct_skim);
    cout << mcf << endl;

    //
    //HF sum!!
    //
    TLine* t1 = new TLine(0,1,30,1);
    t1->SetLineWidth(1);
    t1->SetLineStyle(7); // 7 is jumSun , 1 is onSun
    t1->SetLineColor(1); // 2 is red
    /*
       TH1D *HFsum_mc[5];
       for(int i=0; i<Ncut; i++)
       {
       HFsum_mc[i] = new TH1D(Form("HFsum_mc%d",i), ";#Sigma E_{T}^{HF};Events",5000,0,5000);
       HFsum_mc[i] -> SetMarkerStyle(20);
       HFsum_mc[i] -> SetMarkerSize(1.0);
       HFsum_mc[i] -> SetMarkerColor(i+1);
       }

       mct_evt -> Draw("hiHF >>+ HFsum_mc0", "HLT_HIMinBiasHfOrBSC_v1!=1");
       HFsum_mc[0] = (TH1D*)gDirectory->Get("HFsum_mc0");
       mct_evt -> Draw("hiHF >>+ HFsum_mc1", "pprimaryVertexFilter==1");
       HFsum_mc[1] = (TH1D*)gDirectory->Get("HFsum_mc1");
       mct_evt -> Draw("hiHF >>+ HFsum_mc2", "phltPixelClusterShapeFilter==1");
       HFsum_mc[2] = (TH1D*)gDirectory->Get("HFsum_mc2");
       mct_evt -> Draw("hiHF >>+ HFsum_mc3", "phfCoincFilter3==1");
       HFsum_mc[3] = (TH1D*)gDirectory->Get("HFsum_mc3");
       mct_evt -> Draw("hiHF >>+ HFsum_mc4", "pcollisionEventSelection==1");
       HFsum_mc[4] = (TH1D*)gDirectory->Get("HFsum_mc4");

       TLegend* l1 = new TLegend(0.4, 0.65, 0.85, 0.80);
    //l1 -> AddEntry((TObject*)0, "PbPb Minbias MC"); 
    l1 -> AddEntry(HFsum_mc[0], "No filters");
    l1 -> AddEntry(HFsum_mc[1], "primay vertex filter");
    l1 -> AddEntry(HFsum_mc[2], "pixel cluster shape filter");
    l1 -> AddEntry(HFsum_mc[3], "HF coinc. 3 filter");
    l1 -> AddEntry(HFsum_mc[4], "collision event filter");

    TLegend* l2 = new TLegend(0.4, 0.65, 0.85, 0.80);
    l2 -> AddEntry(HFsum_mc[1], "primay vertex filter");
    l2 -> AddEntry(HFsum_mc[2], "pixel cluster shape filter");
    l2 -> AddEntry(HFsum_mc[3], "HF coinc. 3 filter");
    l2 -> AddEntry(HFsum_mc[4], "collision event filter");

    TCanvas *c_HFsum = new TCanvas("c_HFsum", "c_HFsum", 400,400);
    c_HFsum -> SetLogy();
    l1 -> Draw();
    HFsum_mc[0] -> Draw("e1p");
    HFsum_mc[1] -> Draw("same&&e1p");
    HFsum_mc[2] -> Draw("same&&e1p");
    HFsum_mc[3] -> Draw("same&&e1p");
    HFsum_mc[4] -> Draw("same&&e1p");

    TH1D *HFsum_mc_effi[5];
    for(int i=0; i<Ncut; i++){
    //   HFsum_mc_effi[i] = new TH1D(Form("HFsum_mc_effi%d",i), "", 500,0, 5000);
    HFsum_mc_effi[i] = (TH1D*)HFsum_mc[i]->Clone(Form("HFsum_mc_effi%d",i));
    HFsum_mc_effi[i] -> SetTitle(";#Sigma E_{T}^{HF};Filter Efficiency");
    if(i!=0)
    HFsum_mc_effi[i] -> Divide(HFsum_mc[i],HFsum_mc[0]);
    HFsum_mc_effi[i] -> SetMarkerStyle(20);
    HFsum_mc_effi[i] -> SetMarkerSize(1.0);
    HFsum_mc_effi[i] -> SetMarkerColor(i+1);
    HFsum_mc_effi[i] -> SetAxisRange(0.0,1.1,"Y");
    HFsum_mc_effi[i] -> SetAxisRange(0.0,30.0,"X");
    }

    TCanvas *c_HFsum_effi = new TCanvas("c_HFsum_effi", "c_HFsum_effi", 400, 400);
    t1 -> Draw();
    l2 -> Draw();
    HFsum_mc_effi[1] -> Draw("elp"); 
    HFsum_mc_effi[2] -> Draw("same elp"); 
    HFsum_mc_effi[3] -> Draw("same elp"); 
    HFsum_mc_effi[4] -> Draw("same elp"); 
    */
    //
    // hiNpix
    //

    TH1D *HFNpix_mc[5];
    for(int i=0; i<Ncut; i++)
    {
        HFNpix_mc[i] = new TH1D(Form("HFNpix_mc%d",i), ";hiNpix;Events",500,0,10000);
        HFNpix_mc[i] -> SetMarkerStyle(20);
        HFNpix_mc[i] -> SetMarkerSize(1.0);
        HFNpix_mc[i] -> SetMarkerColor(i+1);
    }

    mct_evt -> Draw("hiNpix >>+ HFNpix_mc0", "HLT_HIMinBiasHfOrBSC_v1!=1");
    HFNpix_mc[0] = (TH1D*)gDirectory->Get("HFNpix_mc0");
    mct_evt -> Draw("hiNpix >>+ HFNpix_mc1", "pprimaryVertexFilter==1");
    HFNpix_mc[1] = (TH1D*)gDirectory->Get("HFNpix_mc1");
    mct_evt -> Draw("hiNpix >>+ HFNpix_mc2", "phltPixelClusterShapeFilter==1");
    HFNpix_mc[2] = (TH1D*)gDirectory->Get("HFNpix_mc2");
    mct_evt -> Draw("hiNpix >>+ HFNpix_mc3", "phfCoincFilter3==1");
    HFNpix_mc[3] = (TH1D*)gDirectory->Get("HFNpix_mc3");
    mct_evt -> Draw("hiNpix >>+ HFNpix_mc4", "pcollisionEventSelection==1");
    HFNpix_mc[4] = (TH1D*)gDirectory->Get("HFNpix_mc4");

    TCanvas *c_HFNpix = new TCanvas("c_HFNpix", "c_HFNpix", 400,400);
    c_HFNpix -> SetLogy();
    HFNpix_mc[0] -> Draw("e1p");
    HFNpix_mc[1] -> Draw("same&&e1p");
    HFNpix_mc[2] -> Draw("same&&e1p");
    HFNpix_mc[3] -> Draw("same&&e1p");
    HFNpix_mc[4] -> Draw("same&&e1p");
    //  l1 -> Draw();

    TH1D *HFNpix_mc_effi[5];
    for(int i=0; i<Ncut; i++){
        //   HFNpix_mc_effi[i] = new TH1D(Form("HFNpix_mc_effi%d",i), "", 500,0, 5000);
        HFNpix_mc_effi[i] = (TH1D*)HFNpix_mc[i]->Clone(Form("HFNpix_mc_effi%d",i));
        HFNpix_mc_effi[i] -> SetTitle("hiNpix;Filter Efficiency");
        if(i!=0)
            HFNpix_mc_effi[i] -> Divide(HFNpix_mc[i],HFNpix_mc[0]);
        HFNpix_mc_effi[i] -> SetMarkerStyle(20);
        HFNpix_mc_effi[i] -> SetMarkerSize(1.0);
        HFNpix_mc_effi[i] -> SetMarkerColor(i+1);
        HFNpix_mc_effi[i] -> SetAxisRange(0.0,1.1,"Y");
        HFNpix_mc_effi[i] -> SetAxisRange(0.0,100.0,"X");
    }

    TCanvas *c_HFNpix_effi = new TCanvas("c_HFNpix_effi", "c_HFNpix_effi", 400, 400);
    //  l2 -> Draw();
    t1 -> Draw();
    HFNpix_mc_effi[1] -> Draw("elp"); 
    HFNpix_mc_effi[2] -> Draw("same elp"); 
    HFNpix_mc_effi[3] -> Draw("same elp"); 
    HFNpix_mc_effi[4] -> Draw("same elp"); 



    //
    // hiBin
    //
    /*
       TH1D *HFBin_mc[5];
       for(int i=0; i<Ncut; i++)
       {
       HFBin_mc[i] = new TH1D(Form("HFBin_mc%d",i), ";hiBin;Events",105,0,210);
       HFBin_mc[i] -> SetMarkerStyle(20);
       HFBin_mc[i] -> SetMarkerSize(1.0);
       HFBin_mc[i] -> SetMarkerColor(i+1);
       }

       mct_evt -> Draw("hiBin >>+ HFBin_mc0", "HLT_HIMinBiasHfOrBSC_v1!=1");
       HFBin_mc[0] = (TH1D*)gDirectory->Get("HFBin_mc0");
       mct_evt -> Draw("hiBin >>+ HFBin_mc1", "pprimaryVertexFilter==1");
       HFBin_mc[1] = (TH1D*)gDirectory->Get("HFBin_mc1");
       mct_evt -> Draw("hiBin >>+ HFBin_mc2", "phltPixelClusterShapeFilter==1");
       HFBin_mc[2] = (TH1D*)gDirectory->Get("HFBin_mc2");
       mct_evt -> Draw("hiBin >>+ HFBin_mc3", "phfCoincFilter3==1");
       HFBin_mc[3] = (TH1D*)gDirectory->Get("HFBin_mc3");
       mct_evt -> Draw("hiBin >>+ HFBin_mc4", "pcollisionEventSelection==1");
       HFBin_mc[4] = (TH1D*)gDirectory->Get("HFBin_mc4");

       TCanvas *c_HFBin = new TCanvas("c_HFBin", "c_HFBin", 400,400);
    //c_HFBin -> SetLogy();
    HFBin_mc[0] -> Draw("e1p");
    HFBin_mc[1] -> Draw("same&&e1p");
    HFBin_mc[2] -> Draw("same&&e1p");
    HFBin_mc[3] -> Draw("same&&e1p");
    HFBin_mc[4] -> Draw("same&&e1p");
    //  l1 -> Draw();

    TH1D *HFBin_mc_effi[5];
    for(int i=0; i<Ncut; i++){
    //   HFBin_mc_effi[i] = new TH1D(Form("HFBin_mc_effi%d",i), "", 500,0, 5000);
    HFBin_mc_effi[i] = (TH1D*)HFBin_mc[i]->Clone(Form("HFBin_mc_effi%d",i));
    HFBin_mc_effi[i] -> SetTitle("hiBin;Filter Efficiency");
    if(i!=0)
    HFBin_mc_effi[i] -> Divide(HFBin_mc[i],HFBin_mc[0]);
    HFBin_mc_effi[i] -> SetMarkerStyle(20);
    HFBin_mc_effi[i] -> SetMarkerSize(1.0);
    HFBin_mc_effi[i] -> SetMarkerColor(i+1);
    HFBin_mc_effi[i] -> SetAxisRange(0.0,1.1,"Y");
    HFBin_mc_effi[i] -> SetAxisRange(140.0,210.0,"X");
    }

    TCanvas *c_HFBin_effi = new TCanvas("c_HFBin_effi", "c_HFBin_effi", 400, 400);
    //    l2 -> Draw();
    t1 -> Draw();
    HFBin_mc_effi[1] -> Draw("elp"); 
    HFBin_mc_effi[2] -> Draw("same elp"); 
    HFBin_mc_effi[3] -> Draw("same elp"); 
    HFBin_mc_effi[4] -> Draw("same elp"); 

*/
    /*
    //
    // HF pixel hits
    // 

    TH1D *HFpixelhit_mc[5];
    for(int i=0; i<Ncut; i++){
    HFpixelhit_mc[i] = new TH1D(Form("HFpixelhit_mc%d",i), ";#Sigma E_{T}^{HF pixel hit};Events",5000,0,5000);
    HFpixelhit_mc[i] -> SetMarkerStyle(20);
    HFpixelhit_mc[i] -> SetMarkerSize(1.0);
    HFpixelhit_mc[i] -> SetMarkerColor(i+1);
    }

    mct_evt -> Draw("hiHF >>+ HFpixelhit_mc0", "HLT_HIMinBiasHfOrBSC_v1!=1");
    HFpixelhit_mc[0] = (TH1D*)gDirectory->Get("HFpixelhit_mc0");
    mct_evt -> Draw("hiHF >>+ HFpixelhit_mc1", "pprimaryVertexFilter==1");
    HFpixelhit_mc[1] = (TH1D*)gDirectory->Get("HFpixelhit_mc1");
    mct_evt -> Draw("hiHF >>+ HFpixelhit_mc2", "phltPixelClusterShapeFilter==1");
    HFpixelhit_mc[2] = (TH1D*)gDirectory->Get("HFpixelhit_mc2");
    mct_evt -> Draw("hiHF >>+ HFpixelhit_mc3", "phfCoincFilter3==1");
    HFpixelhit_mc[3] = (TH1D*)gDirectory->Get("HFpixelhit_mc3");
    mct_evt -> Draw("hiHF >>+ HFpixelhit_mc4", "pcollisionEventSelection==1");
    HFpixelhit_mc[4] = (TH1D*)gDirectory->Get("HFpixelhit_mc4");

    TCanvas *c_HFpixelhit = new TCanvas("c_HFpixelhit", "c_HFpixelhit", 400,400);
    c_HFpixelhit -> SetLogy();
    HFpixelhit_mc[0] -> Draw("e1p");
    HFpixelhit_mc[1] -> Draw("same&&e1p");
    HFpixelhit_mc[2] -> Draw("same&&e1p");
    HFpixelhit_mc[3] -> Draw("same&&e1p");
    HFpixelhit_mc[4] -> Draw("same&&e1p");
    l1 -> Draw();

    TH1D *HFpixelhit_mc_effi[5];
    for(int i=0; i<Ncut; i++){
    HFpixelhit_mc_effi[i] = (TH1D*)HFpixelhit_mc[i]->Clone(Form("HFpixelhit_mc_effi%d",i));
    HFpixelhit_mc_effi[i] -> SetTitle(";#Sigma E_{T}^{HF pixel hit};Filter Efficiency");
    if(i!=0)
    HFpixelhit_mc_effi[i] -> Divide(HFpixelhit_mc[i],HFpixelhit_mc[0]);
    HFpixelhit_mc_effi[i] -> SetMarkerStyle(20);
    HFpixelhit_mc_effi[i] -> SetMarkerSize(1.0);
    HFpixelhit_mc_effi[i] -> SetMarkerColor(i+1);
    HFpixelhit_mc_effi[i] -> SetAxisRange(0.0,1.1,"Y");
    HFpixelhit_mc_effi[i] -> SetAxisRange(0.0,30.0,"X");
    }

    TCanvas *c_HFpixelhit_effi = new TCanvas("c_HFpixelhit_effi", "c_HFpixelhit_effi", 400, 400);
    HFpixelhit_mc_effi[1] -> Draw("elp"); 
    HFpixelhit_mc_effi[2] -> Draw("same elp"); 
    HFpixelhit_mc_effi[3] -> Draw("same elp"); 
    HFpixelhit_mc_effi[4] -> Draw("same elp"); 


    //
    // HF pixel hits vs. HFsum 2D plot
    //

    //gPad-> SetLogy(1);
    TH2D *HFPixelSum_mc[5]; 
    for(int i=0; i<Ncut; i++){
    HFPixelSum_mc[i] = new TH2D(Form("HFPixelSum_mc%d",i), ";#Sigma E_{T}^{HF}; #Sigma E_{HF pixel hits}",30,0,30, 30,0,1400);
    HFpixelhit_mc[i] -> SetMarkerStyle(20);
    HFpixelhit_mc[i] -> SetMarkerSize(1.0);
    HFpixelhit_mc[i] -> SetMarkerColor(i+1);
    }

    mct_evt -> Draw("hiHFhit:hiHF >>+ HFPixelSum_mc0");
    HFPixelSum_mc[0] = (TH2D*)gDirectory->Get("HFPixelSum_mc0");
    mct_evt -> Draw("hiHFhit:hiHF >>+ HFPixelSum_mc1", "pprimaryVertexFilter==1");
    HFPixelSum_mc[1] = (TH2D*)gDirectory->Get("HFPixelSum_mc1");
    mct_evt -> Draw("hiHFhit:hiHF >>+ HFPixelSum_mc2", "phltPixelClusterShapeFilter==1");
    HFPixelSum_mc[2] = (TH2D*)gDirectory->Get("HFPixelSum_mc2");
    mct_evt -> Draw("hiHFhit:hiHF >>+ HFPixelSum_mc3", "phfCoincFilter3==1");
    HFPixelSum_mc[3] = (TH2D*)gDirectory->Get("HFPixelSum_mc3");
    mct_evt -> Draw("hiHFhit:hiHF >>+ HFPixelSum_mc4", "pcollisionEventSelection==1");
    HFPixelSum_mc[4] = (TH2D*)gDirectory->Get("HFPixelSum_mc4");

    TCanvas *c_HFPixelSum = new TCanvas("c_HFPixelSum", "c_HFPixelSum", 1500,300);
    c_HFPixelSum -> Divide(5,1);
    for(int i=0;i<Ncut;i++){
        c_HFPixelSum -> cd(i+1);
        HFPixelSum_mc[i] -> Draw("colz");
    }
    */
        //
        // HF E plus vs. HF E minus 2D plot
        //
        /*
           TH2D *HFPlMa_data = new TH2D("HFPlMa_data", "; #Sigma E_{HF minus  hits}; #Sigma E_{HF plus  hits}",200,0,60000,200,0,60000);
           TH2D *HFPlMa_mc = new TH2D("HFPlMa_mc", "; #Sigma E_{HF minus  hits}; #Sigma E_{HF plus  hits}",200,0,60000,200,0,60000);

           if(onlineFilter==1){
           datat_evt -> Draw("hiHFhitPlus:hiHFhitMinus >>+ HFPlMa_data", "HLT_HIMinBiasHfOrBSC_v1==1");
           HFPlMa_data = (TH2D*)gDirectory->Get("HFPlMa_data");
           mct_evt -> Draw("hiHFhitPlus:hiHFhitMinus >>+ HFPlMa_mc", "HLT_HIMinBiasHfOrBSC_v1==1");
           HFPlMa_mc = (TH2D*)gDirectory->Get("HFPlMa_mc");
           } else {
           datat_evt -> Draw("hiHFhitPlus:hiHFhitMinus >>+ HFPlMa_data");
           HFPlMa_data = (TH2D*)gDirectory->Get("HFPlMa_data");
           mct_evt -> Draw("hiHFhitPlus:hiHFhitMinus >>+ HFPlMa_mc");
           HFPlMa_mc = (TH2D*)gDirectory->Get("HFPlMa_mc");
           }

           TCanvas *c_HFPlMa = new TCanvas("c_HFPlMa", "c_HFPlMa", 800,400);
           c_HFPlMa -> SetLogy();
           c_HFPlMa -> Divide(2,1);
           c_HFPlMa -> cd(1);
           HFPlMa_data -> Draw("colz");
           c_HFPlMa -> cd(2);
           HFPlMa_mc -> Draw("colz");
           */

        //
        // HF pixel hits plus vs. HF pixel hits minus 2D plot
        //
        //
        /*
           TH2D *HFPixelPlMa_data = new TH2D("HFPixelPlMa_data", "; #Sigma E_{HF minus pixel hits}; #Sigma E_{HF plus pixel hits}",200,0,60000,600,0,60000);
           TH2D *HFPixelPlMa_mc = new TH2D("HFPixelPlMa_mc", "; #Sigma E_{HF minus pixel hits}; #Sigma E_{HF plus pixel hits}",200,0,60000,600,0,60000);
           TH2D *ghost = (TH2D*) HFPixelPlMa_data->Clone("ghost"); 

           if(onlineFilter==1){
           datat_evt -> Draw("hiHFhitPlus:hiHFhitMinus >>+ HFPixelPlMa_data", "HLT_HIMinBiasHfOrBSC_v1==1");
           HFPixelPlMa_data = (TH2D*)gDirectory->Get("HFPixelPlMa_data");
           mct_evt -> Draw("hiHFhitPlus:hiHFhitMinus >>+ HFPixelPlMa_mc", "HLT_HIMinBiasHfOrBSC_v1==1");
           HFPixelPlMa_mc = (TH2D*)gDirectory->Get("HFPixelPlMa_mc");
           } else {
           datat_evt -> Draw("hiHFhitPlus:hiHFhitMinus >>+ HFPixelPlMa_data");
           HFPixelPlMa_data = (TH2D*)gDirectory->Get("HFPixelPlMa_data");
           mct_evt -> Draw("hiHFhitPlus:hiHFhitMinus >>+ HFPixelPlMa_mc");
           HFPixelPlMa_mc = (TH2D*)gDirectory->Get("HFPixelPlMa_mc");
           }

           TCanvas *c_HFPixelPlMa = new TCanvas("c_HFPixelPlMa", "c_HFPixelPlMa", 800,400);
           c_HFPixelPlMa -> Divide(2,1);
           c_HFPixelPlMa -> SetLogy(); 
           c_HFPixelPlMa -> cd(1);
           HFPixelPlMa_data -> Draw("colz");
           c_HFPixelPlMa -> cd(2);
           ghost -> Draw();
           HFPixelPlMa_mc -> Draw("same && colz");
           */
} 
