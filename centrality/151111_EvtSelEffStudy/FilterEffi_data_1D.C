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
#include <iostream>
#include "../../yjUtility.h"

const double dy= 0.5;
const int Ncut = 5;
void Get1DEffPlots(TTree* t_evt=0, TString v1="hiHF",int xbin=200, double xmin=0, double xmax=4500, TCut cut="", TCanvas* c_tot=0, TString cap="", bool isPassed=0, double eff_ymin=0.50, bool isAOD=0);

void FilterEffi_data_1D(const char* fname="/afs/cern.ch/work/y/ygo/public/centrality/merged_centrality_MB_DATA750_RECO_150914.root", TString type="RECO_750", bool isAOD=0)
{
   const TCut runCut = "run==181611";
    const TCut lumiCut = "lumi>=1 && lumi<=895";
    const TCut eventCut = runCut && lumiCut;
    TH1::SetDefaultSumw2();
    gStyle -> SetOptStat(0);

    TFile *fin = new TFile(fname);
    TTree *t_evt = (TTree*) fin -> Get("hiEvtAnalyzer/HiTree");
    TTree *t_skim = (TTree*) fin -> Get("skimanalysis/HltTree");
    TTree *t_hlt = (TTree*) fin -> Get("hltanalysis/HltTree");
    t_evt -> AddFriend(t_hlt);
    t_evt -> AddFriend(t_skim);
    double Nevt_t = t_evt -> GetEntries();
    cout << "# of DATA events = " << Nevt_t << endl;


    TCanvas *c_tot = new TCanvas("c_tot", "c_tot", 300,600);
    c_tot->Divide(1,2);
    //    const double HFsum_bins[] = {0,2,5,10,15,20,30,40,50,60,70,80,90,100,200,300,400,500,1000,2000,3000,4000,5000};
    //    const int n_HFsum_bins = sizeof(HFsum_bins)/sizeof(double) - 1;

    Get1DEffPlots(t_evt, "hiHF",100,0,5000,"HLT_HIMinBiasHfOrBSC_v1",c_tot,type,1,0.70,(int)isAOD);
    Get1DEffPlots(t_evt, "hiNpix",100,0,10000,"HLT_HIMinBiasHfOrBSC_v1",c_tot,type,1,0.95,(int)isAOD);
    Get1DEffPlots(t_evt, "hiBin",105,0,210,"HLT_HIMinBiasHfOrBSC_v1",c_tot,type,1,0.0,(int)isAOD);
    Get1DEffPlots(t_evt, "hiZDC",100,0,50000,"HLT_HIMinBiasHfOrBSC_v1",c_tot,type,1,0.0,(int)isAOD);
    Get1DEffPlots(t_evt, "hiET",100,0,2000,"HLT_HIMinBiasHfOrBSC_v1",c_tot,type,1,0.994,(int)isAOD);
    Get1DEffPlots(t_evt, "hiEE",100,0,4000,"HLT_HIMinBiasHfOrBSC_v1",c_tot,type,1,0.992,(int)isAOD);
    Get1DEffPlots(t_evt, "hiEB",100,0,5000,"HLT_HIMinBiasHfOrBSC_v1",c_tot,type,1,0.2,(int)isAOD);
    Get1DEffPlots(t_evt, "vz",100,-50,50,"HLT_HIMinBiasHfOrBSC_v1",c_tot,type,1,0.2,(int)isAOD);

    Get1DEffPlots(t_evt, "hiHF",100,0,5000,"HLT_HIMinBiasHfOrBSC_v1",c_tot,type,0,0.70,(int)isAOD);
    Get1DEffPlots(t_evt, "hiNpix",100,0,10000,"HLT_HIMinBiasHfOrBSC_v1",c_tot,type,0,0.95,(int)isAOD);
    Get1DEffPlots(t_evt, "hiBin",105,0,210,"HLT_HIMinBiasHfOrBSC_v1",c_tot,type,0,0.0,(int)isAOD);
    Get1DEffPlots(t_evt, "hiZDC",100,0,50000,"HLT_HIMinBiasHfOrBSC_v1",c_tot,type,0,0.0,(int)isAOD);
    Get1DEffPlots(t_evt, "hiET",100,0,2000,"HLT_HIMinBiasHfOrBSC_v1",c_tot,type,0,0.994,(int)isAOD);
    Get1DEffPlots(t_evt, "hiEE",100,0,4000,"HLT_HIMinBiasHfOrBSC_v1",c_tot,type,0,0.992,(int)isAOD);
    Get1DEffPlots(t_evt, "hiEB",100,0,5000,"HLT_HIMinBiasHfOrBSC_v1",c_tot,type,0,0.2,(int)isAOD);
    Get1DEffPlots(t_evt, "vz",100,-50,50,"HLT_HIMinBiasHfOrBSC_v1",c_tot,type,0,0.2,(int)isAOD);
}

void Get1DEffPlots(TTree* t_evt, TString v1, int xbin, double xmin, double xmax, TCut cut, TCanvas* c_tot, TString cap, bool isPassed, double eff_ymin, bool isAOD)
{

    TCut totcut[Ncut];
    totcut[0] = cut;
    totcut[1] = cut&& Form("pprimaryVertexFilter==%d",(int)isPassed);
    if(isAOD) totcut[2] = cut&& Form("pclusterCompatibilityFilter==%d",(int)isPassed);
    else totcut[2] = cut&& Form("phltPixelClusterShapeFilter==%d",(int)isPassed);
    totcut[3] = cut&& Form("phfCoincFilter3==%d",(int)isPassed);
    totcut[4] = cut&& Form("pcollisionEventSelection==%d",(int)isPassed);

    TCanvas *c_temp = new TCanvas("c_temp", "c_temp", 300,300);
    c_temp->cd();
    TH1D *h1D[10];
    TH1D *h1D_eff[10];
    for(int i=0; i<Ncut; i++){
        h1D[i] = new TH1D(Form("h1D_%d",i), Form(";%s;Events", v1.Data()), xbin, xmin,xmax );
        h1D_eff[i] = (TH1D*)h1D[i]->Clone(Form("h1D_eff_%d",i));
        h1D_eff[i] -> SetTitle(Form(";%s;Filter Rate",v1.Data()));
        h1D[i] -> SetMarkerStyle(20);
        //h1D[i] -> SetMarkerSize();
        h1D[i] -> SetMarkerColor(col[i]);
        //h1D[i] -> SetLabelSize(0.03);
        h1D_eff[i] -> SetMarkerStyle(20);
        //h1D_eff[i] -> SetMarkerSize(0.7);
        h1D_eff[i] -> SetMarkerColor(col[i]);
        h1D_eff[i] -> SetAxisRange(eff_ymin,1.0005,"Y");

    }
    for(int i=0; i<Ncut; i++){
        t_evt->Draw(Form("%s>>+%s",v1.Data(), h1D[i]->GetName() ), totcut[i]);
        h1D[i]=(TH1D*)gDirectory->Get(h1D[i]->GetName());

    }

    TLegend* l1 = new TLegend(0.45, 0.7, 0.9, 0.95, "PbPb Minbias rereco DATA");
    // TLegend* l1 = new TLegend(0.1, 0.15, 0.6, 0.50, "PbPb Minbias rereco DATA");
    // l1 -> AddEntry((TObject*)0, "PbPb Minbias rereco DATA"); 
    legStyle(l1);
    l1 -> AddEntry(h1D[0], "No filters");
    l1 -> AddEntry(h1D[1], "primay vertex filter");
    if(isAOD) l1 -> AddEntry(h1D[2], "cluster compatibility filter");
    else l1 -> AddEntry(h1D[2], "pixel cluster shape filter");
    l1 -> AddEntry(h1D[3], "HF coinc. 3 filter");
    l1 -> AddEntry(h1D[4], "collision event filter");

    TLine* t1 = new TLine(0,1,100000,1);
    t1->SetLineWidth(1);
    t1->SetLineStyle(7); // 7 is jumSun , 1 is onSun
    t1->SetLineColor(1); // 2 is red

    for(int i=0;i<Ncut;i++){
        c_tot -> cd(1);
        //gPad -> SetLogx();
        //gPad -> SetLogy();
        if(i==0) h1D[i] -> Draw("ep");
        else h1D[i] -> Draw("ep same");
        l1 -> Draw();
        
        if(i!=0){
            c_tot->cd(2);
            gPad -> SetLogx();
            h1D_eff[i] -> Divide(h1D[i],h1D[0],1,1,"B");
            h1D_eff[i] -> GetYaxis()-> SetRangeUser(eff_ymin,1.005);
            if(i==1) {
                h1D_eff[i] -> Draw("ep"); 
                t1 -> Draw("same"); 
                drawText(cap, 0.5,0.6);
                drawText("2011 PbPb minbias", 0.5,0.68);
            }
            else h1D_eff[i] -> Draw("ep same"); 
        }
    }
    for(int i=0;i<Ncut;i++){
       //h1D[i]-> Scale(1./t_evt->GetEntries());
    }

    c_tot->SaveAs(Form("png/h1D_%s_isPassed%d_%s.png",v1.Data(),(int)isPassed,cap.Data()));

    for(int i=0; i<Ncut; i++){
        delete h1D[i];
        delete h1D_eff[i];
    }
    delete c_temp;

}
