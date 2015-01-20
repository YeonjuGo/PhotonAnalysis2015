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

void FilterEffi_data_mc_hiNpix()
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
    
    TFile *dataf = new TFile("/u/user/goyeonju/PRODUCTION/CMSSW_5_3_20/src/centralityDATA/merging-forest/HiForest_100_1_pKq.root");//in KNU server
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

    const double hiNpix_bins[] = {0,2,5,10,15,20,30,40,50,60,70,80,90,100,200,300,400,500,1000,2000,3000,4000,5000};
    const int n_hiNpix_bins = sizeof(hiNpix_bins)/sizeof(double) - 1;
    
    TLine* t1 = new TLine(0,1,1000,1);
    t1->SetLineWidth(1);
    t1->SetLineStyle(7); // 7 is jumSun , 1 is onSun
    t1->SetLineColor(1); // 2 is red

    TH1D *hiNpix_data[5];
    TH1D *hiNpix_mc[5];
    for(int i=0; i<Ncut; i++)
    {
       // hiNpix_data[i] = new TH1D(Form("hiNpix_data%d",i), ";hiNpix;Normalized Events",n_hiNpix_bins, hiNpix_bins);
        hiNpix_data[i] = new TH1D(Form("hiNpix_data%d",i), ";hiNpix;Normalized Events",100,0,5000);
        hiNpix_data[i] -> SetMarkerStyle(20+i);
        hiNpix_data[i] -> SetMarkerSize(0.7);
        hiNpix_data[i] -> SetMarkerColor(kRed+i);
        hiNpix_data[i] -> SetLabelSize(0.03);
   	
	hiNpix_mc[i] = new TH1D(Form("hiNpix_mc%d",i), ";hiNpix;Normalized Events",n_hiNpix_bins, hiNpix_bins);
	hiNpix_mc[i] = new TH1D(Form("hiNpix_mc%d",i), ";hiNpix;Normalized Events",100,0,5000);
        //hiNpix_mc[i] -> SetMarkerStyle(20);
        //hiNpix_mc[i] -> SetMarkerSize(1.0);
        hiNpix_mc[i] -> SetLineColor(2+i); //line color
        hiNpix_mc[i] -> SetLabelSize(0.03);
    }
    TCanvas *c_temp = new TCanvas("c_temp", "c_temp", 300,300);
	
    c_temp -> cd();
    datat_evt -> Draw("hiNpix >>+ hiNpix_data0",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1");
    hiNpix_data[0] = (TH1D*)gDirectory->Get("hiNpix_data0");
    datat_evt -> Draw("hiNpix >>+ hiNpix_data1",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && pprimaryVertexFilter==1");
    hiNpix_data[1] = (TH1D*)gDirectory->Get("hiNpix_data1");
    datat_evt -> Draw("hiNpix >>+ hiNpix_data2",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && phltPixelClusterShapeFilter==1");
    hiNpix_data[2] = (TH1D*)gDirectory->Get("hiNpix_data2");
    datat_evt -> Draw("hiNpix >>+ hiNpix_data3",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && phfCoincFilter3==1");
    hiNpix_data[3] = (TH1D*)gDirectory->Get("hiNpix_data3");
    datat_evt -> Draw("hiNpix >>+ hiNpix_data4",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && pcollisionEventSelection==1");
    hiNpix_data[4] = (TH1D*)gDirectory->Get("hiNpix_data4");

    mct_evt -> Draw("hiNpix >>+ hiNpix_mc0");
    hiNpix_mc[0] = (TH1D*)gDirectory->Get("hiNpix_mc0");
    mct_evt -> Draw("hiNpix >>+ hiNpix_mc1","pprimaryVertexFilter==1");
    hiNpix_mc[1] = (TH1D*)gDirectory->Get("hiNpix_mc1");
    mct_evt -> Draw("hiNpix >>+ hiNpix_mc2","phltPixelClusterShapeFilter==1");
    hiNpix_mc[2] = (TH1D*)gDirectory->Get("hiNpix_mc2");
    mct_evt -> Draw("hiNpix >>+ hiNpix_mc3","phfCoincFilter3==1");
    hiNpix_mc[3] = (TH1D*)gDirectory->Get("hiNpix_mc3");
    mct_evt -> Draw("hiNpix >>+ hiNpix_mc4","pcollisionEventSelection==1");
    hiNpix_mc[4] = (TH1D*)gDirectory->Get("hiNpix_mc4");

    TLegend* l1 = new TLegend(0.3, 0.65, 0.6, 0.80, "PbPb Minbias rereco DATA");
   // l1 -> AddEntry((TObject*)0, "PbPb Minbias rereco DATA"); 
    l1 -> AddEntry(hiNpix_data[0], "No filters");
    l1 -> AddEntry(hiNpix_data[1], "primay vertex filter");
    l1 -> AddEntry(hiNpix_data[2], "pixel cluster shape filter");
    l1 -> AddEntry(hiNpix_data[3], "HF coinc. 3 filter");
    l1 -> AddEntry(hiNpix_data[4], "collision event filter");
  
    TLegend* l1_mc = new TLegend(0.6, 0.65, 0.9, 0.80, "PbPb Minbias MC");
    //l1 -> AddEntry((TObject*)0, "PbPb Minbias MC"); 
    l1_mc -> AddEntry(hiNpix_mc[0], "No filters");
    l1_mc -> AddEntry(hiNpix_mc[1], "primay vertex filter");
    l1_mc -> AddEntry(hiNpix_mc[2], "pixel cluster shape filter");
    l1_mc -> AddEntry(hiNpix_mc[3], "HF coinc. 3 filter");
    l1_mc -> AddEntry(hiNpix_mc[4], "collision event filter");
   
    TLegend* l2 = new TLegend(0.4, 0.65, 0.85, 0.80,"PbPb Minbias rereco DATA");
    //l2 -> AddEntry((TObject*)0, "PbPb Minbias rereco DATA"); 
    l2 -> AddEntry(hiNpix_data[1], "primay vertex filter");
    l2 -> AddEntry(hiNpix_data[2], "pixel cluster shape filter");
    l2 -> AddEntry(hiNpix_data[3], "HF coinc. 3 filter");
    l2 -> AddEntry(hiNpix_data[4], "collision event filter");

    TLegend* l2_mc = new TLegend(0.4, 0.65, 0.85, 0.80, "PbPb Minbias MC");
    //l2_mc -> AddEntry((TObject*)0, "PbPb Minbias MC"); 
    l2_mc -> AddEntry(hiNpix_mc[1], "primay vertex filter");
    l2_mc -> AddEntry(hiNpix_mc[2], "pixel cluster shape filter");
    l2_mc -> AddEntry(hiNpix_mc[3], "HF coinc. 3 filter");
    l2_mc -> AddEntry(hiNpix_mc[4], "collision event filter");

    TCanvas *c_hiNpix = new TCanvas("c_hiNpix", "c_hiNpix", 400,400);
    c_hiNpix -> SetLogy();

    double norm_data = hiNpix_data[0]->Integral("width");
    double norm_mc = hiNpix_mc[0]->Integral("width");
    for(int i=0; i<Ncut; i++){
        if(i==0) hiNpix_data[i] -> Draw("ep");
        else hiNpix_data[i] -> Draw("same&&ep");
        hiNpix_mc[i] -> Draw("same&&ehist");

        //hiNpix_data[i] -> Scale(1./Nevt_datat);
        //hiNpix_mc[i] -> Scale(1./Nevt_mct);
        hiNpix_data[i] -> Scale(1./norm_data);
        hiNpix_mc[i] -> Scale(1./norm_mc);

        cout << "The integral of hiNpix data " << i << " : " << hiNpix_data[i]->Integral("width") << endl;
        cout << "The integral of hiNpix mc " << i << " : " << hiNpix_mc[i]->Integral("width") << endl;
    }
   
    l1 -> Draw();
    l1_mc -> Draw();
    
    c_hiNpix -> SaveAs("pdf/hiNpix.pdf");


    TH1D *hiNpix_data_effi[5];
    for(int i=0; i<Ncut; i++){
     //   hiNpix_data_effi[i] = new TH1D(Form("hiNpix_data_effi%d",i), "", 500,0, 5000);
     //  hiNpix_data_effi[i] = new TH1D(Form("hiNpix_data_effi%d",i), ";hiNpix;Normalized Events",n_hiNpix_bins,hiNpix_bins);
        hiNpix_data_effi[i] = (TH1D*)hiNpix_data[i]->Clone(Form("hiNpix_data_effi%d",i));
	   hiNpix_data_effi[i] -> SetTitle(";hiNpix;Filter Efficiency");
        if(i!=0)
            hiNpix_data_effi[i] -> Divide(hiNpix_data[i],hiNpix_data[0]);

        hiNpix_data_effi[i] -> SetMarkerStyle(20+i);
        hiNpix_data_effi[i] -> SetMarkerSize(0.7);
        hiNpix_data_effi[i] -> SetMarkerColor(kRed+i);
        hiNpix_data_effi[i] -> SetAxisRange(0.0,1.1,"Y");
        //hiNpix_data_effi[i] -> SetAxisRange(0.0,300.0,"X");
    }
    TCanvas *c_hiNpix_data_effi = new TCanvas("c_hiNpix_data_effi", "c_hiNpix_data_effi", 400, 400);
    hiNpix_data_effi[1] -> Draw("elp"); 
    hiNpix_data_effi[2] -> Draw("same ep"); 
    hiNpix_data_effi[3] -> Draw("same ep"); 
    hiNpix_data_effi[4] -> Draw("same ep"); 
    t1 -> Draw();
    l2 -> Draw();
    
    c_hiNpix_data_effi -> SetLogx();
    c_hiNpix_data_effi -> SaveAs("pdf/hiNpix_data_effi.pdf");

   TH1D *hiNpix_mc_effi[5];
    for(int i=0; i<Ncut; i++){
     //   hiNpix_mc_effi[i] = new TH1D(Form("hiNpix_mc_effi%d",i), "", 500,0, 5000);
        hiNpix_mc_effi[i] = (TH1D*)hiNpix_mc[i]->Clone(Form("hiNpix_mc_effi%d",i));
        hiNpix_mc_effi[i] -> SetTitle(";hiNpix;Filter Efficiency");
        if(i!=0)
            hiNpix_mc_effi[i] -> Divide(hiNpix_mc[i],hiNpix_mc[0]);
        hiNpix_mc_effi[i] -> SetMarkerStyle(20+i);
        hiNpix_mc_effi[i] -> SetMarkerSize(0.7);
        hiNpix_mc_effi[i] -> SetMarkerColor(2+i);
        hiNpix_mc_effi[i] -> SetAxisRange(0.0,1.1,"Y");
        //hiNpix_mc_effi[i] -> SetAxisRange(0.0,300.0,"X");
    }

    TCanvas *c_hiNpix_mc_effi = new TCanvas("c_hiNpix_mc_effi", "c_hiNpix_mc_effi", 400, 400);
    hiNpix_mc_effi[1] -> Draw("elp"); 
    hiNpix_mc_effi[2] -> Draw("same ep"); 
    hiNpix_mc_effi[3] -> Draw("same ep"); 
    hiNpix_mc_effi[4] -> Draw("same ep"); 
    t1 -> Draw();
    l2_mc -> Draw();
    
    c_hiNpix_mc_effi -> SetLogx();
    c_hiNpix_mc_effi -> SaveAs("pdf/hiNpix_mc_effi.pdf");

    TCanvas *c_hiNpix_effi = new TCanvas("c_hiNpix_effi","c_hiNpix_effi", 400,400);
    c_hiNpix_effi -> Divide(1,2,0);
    c_hiNpix_effi -> SetLogx(); 
    c_hiNpix_effi -> cd(1);
    hiNpix_data_effi[0] -> Draw("ep");
    hiNpix_data_effi[1] -> Draw("same ep");
    hiNpix_data_effi[2] -> Draw("same ep");
    hiNpix_data_effi[3] -> Draw("same ep");
    hiNpix_data_effi[4] -> Draw("same ep");
    
    c_hiNpix_effi -> cd(2);
    hiNpix_mc_effi[0] -> Draw("ep");
    hiNpix_mc_effi[1] -> Draw("same ep");
    hiNpix_mc_effi[2] -> Draw("same ep");
    hiNpix_mc_effi[3] -> Draw("same ep");
    hiNpix_mc_effi[4] -> Draw("same ep");

} 
