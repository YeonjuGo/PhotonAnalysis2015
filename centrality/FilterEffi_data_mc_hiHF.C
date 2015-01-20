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

void FilterEffi_data_mc_1D()
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
    const int Ncut = 5;
    TH1::SetDefaultSumw2();
    gStyle -> SetOptStat(0);
    
    TFile *dataf = new TFile("/u/user/goyeonju/PRODUCTION/CMSSW_5_3_20/src/centralityDATA/merging-forest/HiForest_100_1_pKq.root");
    //TFile *dataf = new TFile("files/HiForest_20141011.root");
    TTree *datat_evt = (TTree*) dataf -> Get("hiEvtAnalyzer/HiTree");
    TTree *datat_skim = (TTree*) dataf -> Get("skimanalysis/HltTree");
    TTree *datat_hlt = (TTree*) dataf -> Get("hltanalysis/HltTree");
    datat_evt -> AddFriend(datat_hlt);
    datat_evt -> AddFriend(datat_skim);
    double Nevt_datat = datat_evt -> GetEntries();
    cout << "# of DATA events = " << Nevt_datat << endl;

    TFile *mcf = new TFile("/u/user/goyeonju/files/centrality/HiForest_HydjetMB_730_53XBS_merged.root");//in KNU server
    TTree *mct_evt = (TTree*) mcf -> Get("hiEvtAnalyzer/HiTree");
    TTree *mct_skim = (TTree*) mcf -> Get("skimanalysis/HltTree");
    TTree *mct_hlt = (TTree*) mcf -> Get("hltanalysis/HltTree");
    mct_evt -> AddFriend(mct_hlt);
    mct_evt -> AddFriend(mct_skim);
    int Nevt_mct = mct_evt -> GetEntries();
    cout << "# of MC events = " << Nevt_mct << endl;


//======================================
//HF sum!!
//======================================

    const double HFsum_bins[] = {0,2,5,10,15,20,30,40,50,60,70,80,90,100,200,300,400,500,1000,2000,3000,4000,5000};
    const int n_HFsum_bins = sizeof(HFsum_bins)/sizeof(double) - 1;
    
    TLine* t1 = new TLine(0,1,1000,1);
    t1->SetLineWidth(1);
    t1->SetLineStyle(7); // 7 is jumSun , 1 is onSun
    t1->SetLineColor(1); // 2 is red

    TH1D *HFsum_data[5];
    TH1D *HFsum_mc[5];
    for(int i=0; i<Ncut; i++)
    {
       // HFsum_data[i] = new TH1D(Form("HFsum_data%d",i), ";hiHF;Normalized Events",n_HFsum_bins, HFsum_bins);
        HFsum_data[i] = new TH1D(Form("HFsum_data%d",i), ";hiHF;Normalized Events",100,0,5000);
        HFsum_data[i] -> SetMarkerStyle(20+i);
        HFsum_data[i] -> SetMarkerSize(0.7);
        HFsum_data[i] -> SetMarkerColor(kRed+i);
        HFsum_data[i] -> SetLabelSize(0.03);
   	
	HFsum_mc[i] = new TH1D(Form("HFsum_mc%d",i), ";hiHF;Normalized Events",n_HFsum_bins, HFsum_bins);
	HFsum_mc[i] = new TH1D(Form("HFsum_mc%d",i), ";hiHF;Normalized Events",100,0,5000);
        //HFsum_mc[i] -> SetMarkerStyle(20);
        //HFsum_mc[i] -> SetMarkerSize(1.0);
        HFsum_mc[i] -> SetLineColor(2+i); //line color
        HFsum_mc[i] -> SetLabelSize(0.03);
    }
    TCanvas *c_temp = new TCanvas("c_temp", "c_temp", 300,300);
	
    c_temp -> cd();
    datat_evt -> Draw("hiHF >>+ HFsum_data0",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1");
    HFsum_data[0] = (TH1D*)gDirectory->Get("HFsum_data0");
    datat_evt -> Draw("hiHF >>+ HFsum_data1",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && pprimaryVertexFilter==1");
    HFsum_data[1] = (TH1D*)gDirectory->Get("HFsum_data1");
    datat_evt -> Draw("hiHF >>+ HFsum_data2",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && phltPixelClusterShapeFilter==1");
    HFsum_data[2] = (TH1D*)gDirectory->Get("HFsum_data2");
    datat_evt -> Draw("hiHF >>+ HFsum_data3",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && phfCoincFilter3==1");
    HFsum_data[3] = (TH1D*)gDirectory->Get("HFsum_data3");
    datat_evt -> Draw("hiHF >>+ HFsum_data4",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && pcollisionEventSelection==1");
    HFsum_data[4] = (TH1D*)gDirectory->Get("HFsum_data4");

    mct_evt -> Draw("hiHF >>+ HFsum_mc0");
    HFsum_mc[0] = (TH1D*)gDirectory->Get("HFsum_mc0");
    mct_evt -> Draw("hiHF >>+ HFsum_mc1","pprimaryVertexFilter==1");
    HFsum_mc[1] = (TH1D*)gDirectory->Get("HFsum_mc1");
    mct_evt -> Draw("hiHF >>+ HFsum_mc2","phltPixelClusterShapeFilter==1");
    HFsum_mc[2] = (TH1D*)gDirectory->Get("HFsum_mc2");
    mct_evt -> Draw("hiHF >>+ HFsum_mc3","phfCoincFilter3==1");
    HFsum_mc[3] = (TH1D*)gDirectory->Get("HFsum_mc3");
    mct_evt -> Draw("hiHF >>+ HFsum_mc4","pcollisionEventSelection==1");
    HFsum_mc[4] = (TH1D*)gDirectory->Get("HFsum_mc4");

    TLegend* l1 = new TLegend(0.3, 0.65, 0.6, 0.80, "PbPb Minbias rereco DATA");
   // l1 -> AddEntry((TObject*)0, "PbPb Minbias rereco DATA"); 
    l1 -> AddEntry(HFsum_data[0], "No filters");
    l1 -> AddEntry(HFsum_data[1], "primay vertex filter");
    l1 -> AddEntry(HFsum_data[2], "pixel cluster shape filter");
    l1 -> AddEntry(HFsum_data[3], "HF coinc. 3 filter");
    l1 -> AddEntry(HFsum_data[4], "collision event filter");
  
    TLegend* l1_mc = new TLegend(0.6, 0.65, 0.9, 0.80, "PbPb Minbias MC");
    //l1 -> AddEntry((TObject*)0, "PbPb Minbias MC"); 
    l1_mc -> AddEntry(HFsum_mc[0], "No filters");
    l1_mc -> AddEntry(HFsum_mc[1], "primay vertex filter");
    l1_mc -> AddEntry(HFsum_mc[2], "pixel cluster shape filter");
    l1_mc -> AddEntry(HFsum_mc[3], "HF coinc. 3 filter");
    l1_mc -> AddEntry(HFsum_mc[4], "collision event filter");
   
    TLegend* l2 = new TLegend(0.4, 0.65, 0.85, 0.80,"PbPb Minbias rereco DATA");
    //l2 -> AddEntry((TObject*)0, "PbPb Minbias rereco DATA"); 
    l2 -> AddEntry(HFsum_data[1], "primay vertex filter");
    l2 -> AddEntry(HFsum_data[2], "pixel cluster shape filter");
    l2 -> AddEntry(HFsum_data[3], "HF coinc. 3 filter");
    l2 -> AddEntry(HFsum_data[4], "collision event filter");

    TLegend* l2_mc = new TLegend(0.4, 0.65, 0.85, 0.80, "PbPb Minbias MC");
    //l2_mc -> AddEntry((TObject*)0, "PbPb Minbias MC"); 
    l2_mc -> AddEntry(HFsum_mc[1], "primay vertex filter");
    l2_mc -> AddEntry(HFsum_mc[2], "pixel cluster shape filter");
    l2_mc -> AddEntry(HFsum_mc[3], "HF coinc. 3 filter");
    l2_mc -> AddEntry(HFsum_mc[4], "collision event filter");

    TCanvas *c_HFsum = new TCanvas("c_HFsum", "c_HFsum", 400,400);
    c_HFsum -> SetLogy();

    double norm_data = HFsum_data[0]->Integral("width");
    double norm_mc = HFsum_mc[0]->Integral("width");
    for(int i=0; i<Ncut; i++){
        if(i==0) HFsum_data[i] -> Draw("ep");
        else HFsum_data[i] -> Draw("same&&ep");
        HFsum_mc[i] -> Draw("same&&ehist");

        //HFsum_data[i] -> Scale(1./Nevt_datat);
        //HFsum_mc[i] -> Scale(1./Nevt_mct);
        HFsum_data[i] -> Scale(1./norm_data);
        HFsum_mc[i] -> Scale(1./norm_mc);

        cout << "The integral of HFsum data " << i << " : " << HFsum_data[i]->Integral("width") << endl;
        cout << "The integral of HFsum mc " << i << " : " << HFsum_mc[i]->Integral("width") << endl;
    }
   
    l1 -> Draw();
    l1_mc -> Draw();
    
    c_HFsum -> SaveAs("pdf/HFsum.pdf");


    TH1D *HFsum_data_effi[5];
    for(int i=0; i<Ncut; i++){
     //   HFsum_data_effi[i] = new TH1D(Form("HFsum_data_effi%d",i), "", 500,0, 5000);
     //  HFsum_data_effi[i] = new TH1D(Form("HFsum_data_effi%d",i), ";hiHF;Normalized Events",n_HFsum_bins,HFsum_bins);
        HFsum_data_effi[i] = (TH1D*)HFsum_data[i]->Clone(Form("HFsum_data_effi%d",i));
	   HFsum_data_effi[i] -> SetTitle(";hiHF;Filter Efficiency");
        if(i!=0)
            HFsum_data_effi[i] -> Divide(HFsum_data[i],HFsum_data[0]);

        HFsum_data_effi[i] -> SetMarkerStyle(20+i);
        HFsum_data_effi[i] -> SetMarkerSize(0.7);
        HFsum_data_effi[i] -> SetMarkerColor(kRed+i);
        HFsum_data_effi[i] -> SetAxisRange(0.0,1.1,"Y");
        //HFsum_data_effi[i] -> SetAxisRange(0.0,300.0,"X");
    }
    TCanvas *c_HFsum_data_effi = new TCanvas("c_HFsum_data_effi", "c_HFsum_data_effi", 400, 400);
    HFsum_data_effi[1] -> Draw("elp"); 
    HFsum_data_effi[2] -> Draw("same ep"); 
    HFsum_data_effi[3] -> Draw("same ep"); 
    HFsum_data_effi[4] -> Draw("same ep"); 
    t1 -> Draw();
    l2 -> Draw();
    
    c_HFsum_data_effi -> SetLogx();
    c_HFsum_data_effi -> SaveAs("pdf/HFsum_data_effi.pdf");

   TH1D *HFsum_mc_effi[5];
    for(int i=0; i<Ncut; i++){
     //   HFsum_mc_effi[i] = new TH1D(Form("HFsum_mc_effi%d",i), "", 500,0, 5000);
        HFsum_mc_effi[i] = (TH1D*)HFsum_mc[i]->Clone(Form("HFsum_mc_effi%d",i));
        HFsum_mc_effi[i] -> SetTitle(";hiHF;Filter Efficiency");
        if(i!=0)
            HFsum_mc_effi[i] -> Divide(HFsum_mc[i],HFsum_mc[0]);
        HFsum_mc_effi[i] -> SetMarkerStyle(20+i);
        HFsum_mc_effi[i] -> SetMarkerSize(0.7);
        HFsum_mc_effi[i] -> SetMarkerColor(2+i);
        HFsum_mc_effi[i] -> SetAxisRange(0.0,1.1,"Y");
        //HFsum_mc_effi[i] -> SetAxisRange(0.0,300.0,"X");
    }

    TCanvas *c_HFsum_mc_effi = new TCanvas("c_HFsum_mc_effi", "c_HFsum_mc_effi", 400, 400);
    HFsum_mc_effi[1] -> Draw("elp"); 
    HFsum_mc_effi[2] -> Draw("same ep"); 
    HFsum_mc_effi[3] -> Draw("same ep"); 
    HFsum_mc_effi[4] -> Draw("same ep"); 
    t1 -> Draw();
    l2_mc -> Draw();
    
    c_HFsum_mc_effi -> SetLogx();
    c_HFsum_mc_effi -> SaveAs("pdf/HFsum_mc_effi.pdf");

    TCanvas *c_HFsum_effi = new TCanvas("c_HFsum_effi","c_HFsum_effi", 400,400);
    c_HFsum_effi -> Divide(1,2,0);
    c_HFsum_effi -> SetLogx(); 
    c_HFsum_effi -> cd(1);
    HFsum_data_effi[0] -> Draw("ep");
    HFsum_data_effi[1] -> Draw("same ep");
    HFsum_data_effi[2] -> Draw("same ep");
    HFsum_data_effi[3] -> Draw("same ep");
    HFsum_data_effi[4] -> Draw("same ep");
    
    c_HFsum_effi -> cd(2);
    HFsum_mc_effi[0] -> Draw("ep");
    HFsum_mc_effi[1] -> Draw("same ep");
    HFsum_mc_effi[2] -> Draw("same ep");
    HFsum_mc_effi[3] -> Draw("same ep");
    HFsum_mc_effi[4] -> Draw("same ep");

} 
