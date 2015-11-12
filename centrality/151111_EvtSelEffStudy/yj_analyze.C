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

void yj_analyze(bool isMC=0){
    const char* infile;
    if(isMC) infile = "/u/user/goyeonju/files/centrality/Centrality_officialMC_Hydjet1p8_TuneDrum_Quenched_MinBias_2760GeV.root";
    else infile = "/u/user/goyeonju/files/centrality/PbPb_minbias_data_2760_HIRun2011-14Mar2014-v2_run181611_CMSSW5320_byYJ.root";
    const char* trig[]={""};
    //const char* trig = Form("HltTree.pcollisionEventSelection==1");//"L1Tech_BSC_minBias_threshold1.v0";
    double towerCut = 3.0;
    double ebCut = 3;

    const char* outName = Form("histfiles/allHist_isMC%d_towerCut%d_%s.root",(int)isMC,(int)towerCut,trig);

    TFile * inf = new TFile(infile);
    TFile* outf = new TFile(outName,"recreate");

    TH1D* ha[100];
    TH1D* hb[100];
    TH1D* hc[100];
    TH1D* he[100];
    TH1D* hn[100];

    TH1D* hneb[100];
    TH1D* hnee[100];

    TH2D* h2o[100];
    TH2D* h2a[100];
    TH2D* h2b[100];
    TH2D* h2c[100];
    TH2D* h2d[100];
    TH2D* h2pm[100];
    TH2D* h2ebhf[100];

    for(int i = 0; i < 20; ++i){
        ha[i] = new TH1D(Form("ha%d",i),";N HF towers above throshold;Events",20,0,1000);
        hb[i] = new TH1D(Form("hb%d",i),";Sum Tower E_{T} [GeV];Events",20,0,4000);
        hc[i] = new TH1D(Form("hc%d",i),"",1000,0,1000);
        he[i] = new TH1D(Form("he%d",i),"",1000,0,4000);
        hn[i] = new TH1D(Form("hn%d",i),"",1000,0,10000);

        hneb[i] = new TH1D(Form("hneb%d",i),";N_{hits} EE > 0.3 GeV; Events",20,0,10000);
        hnee[i] = new TH1D(Form("hnee%d",i),";N_{hits} EE > 0.3 GeV; Events",20,0,10000);

        h2o[i] = new TH2D(Form("h2o%d",i),"",420,0,420,600,0,6000);
        h2a[i] = new TH2D(Form("h2a%d",i),"",420,0,420,900,0,9000);
        h2b[i] = new TH2D(Form("h2b%d",i),"",900,0,9000,600,0,6000);
        h2c[i] = new TH2D(Form("h2c%d",i),"",420,0,420,370,0,3700);
        h2d[i] = new TH2D(Form("h2d%d",i),"",370,0,3700,600,0,6000);

        h2pm[i] = new TH2D(Form("h2pm%d",i),"; Energy Sum HF- [GeV]; Energy Sum HF+ [GeV]",20,0,60000,20,0,60000);
        h2ebhf[i] = new TH2D(Form("h2ebhf%d",i),"; Energy Sum HF [GeV]; E_{T} Sum EB [GeV]",20,0,120000,20,0,80);

    }

    ha[1]->SetLineColor(2);
    ha[1]->SetMarkerColor(2);

    TTree* t1 = (TTree*)inf->Get("hltanalysis/HltTree");
    TTree* t2 = (TTree*)inf->Get("rechitanalyzer/hf");
    TTree* t3 = (TTree*)inf->Get("rechitanalyzer/hbhe");
    TTree* t4 = (TTree*)inf->Get("rechitanalyzer/ee");
    TTree* t5 = (TTree*)inf->Get("rechitanalyzer/eb");
    //TTree* t6 = (TTree*)inf->Get("rechitanalyzer/bkg");
    TTree* t7 = (TTree*)inf->Get("rechitanalyzer/tower");
    TTree* t8 = (TTree*)inf->Get("hiEvtAnalyzer/HiTree");

    if(isMC){
//        TTree* t9 = (TTree*)inf->Get("HiGenParticleAna/hi");
//        t1->AddFriend(t9);
    }

    t1->AddFriend(t2);
    t1->AddFriend(t3);
    t1->AddFriend(t4);
    t1->AddFriend(t5);
    //t1->AddFriend(t6);
    t1->AddFriend(t7);
    t1->AddFriend(t8);

    TCanvas* c1 = new TCanvas("c1","",1200,800);
    c1->Divide(5,3);

    c1->cd(1);
    t1->Draw(Form("Sum$(hf.e > %f)",towerCut),trig);
    c1->cd(2);
    t1->Draw(Form("Sum$(tower.et)>>hb0"),trig);

    c1->cd(3);

    t1->Draw(Form("Sum$(eb.et > 0.1)>>hneb0"),trig);

    c1->cd(4);

    t1->Draw(Form("Sum$(ee.et > 0.3)>>hnee0"),trig);

    c1->cd(5);
    t1->Draw(Form("hiHFhitPlus:hiHFhitMinus>>h2pm0"),trig,"colz");

    c1->cd(6);
    t1->Draw(Form("hiHFhit:hiEB>>hebhf0"),trig,"colz");

    c1->cd(7);
    t1->Draw(Form("hiHFhit:hiNpix>>h2ebhf0"),trig,"colz");

    c1->cd(8);
    c1->cd(9);
    c1->cd(10);
    //  t1->Draw(Form("Sum$(abs(tower.eta) > 2.87 && tower.e > %f)>>ha0",towerCut));
    t1->Draw(Form("Sum$(abs(tower.eta) > 2.87 && tower.e > %f)>>ha1",towerCut),trig,"same");

    c1->cd(11);
    t1->Draw("Sum$(tower.et * (abs(tower.eta) > 2.87))",trig);

    c1->cd(12);
    t1->Draw(Form("Sum$(abs(tower.eta) > 2.87 && tower.e > %f)>>ha1",towerCut),trig);

    if(isMC){
        TCanvas * c2 = new TCanvas("c2","",900,600);
        c2->Divide(3,2);
        c2->cd(1);
        t1->Draw("hiHF:npart>>h2o0",trig,"colz");

        c2->cd(2);

        t1->Draw("Sum$(hi.eta > 2.87 && hi.eta < 5.205):npart>>h2a0",trig,"colz");
        c2->cd(3);
        t1->Draw("hiHF:Sum$(hi.eta > 2.87 && hi.eta < 5.205)>>h2b0",trig,"colz");

        c2->cd(4);
        t1->Draw("Sum$(hi.pt* (hi.eta > 2.87 && hi.eta < 5.205)):npart>>h2c0",trig,"colz");
        c2->cd(5);
        t1->Draw("hiHF:Sum$(hi.pt* (hi.eta > 2.87 && hi.eta < 5.205))>>h2d0",trig,"colz");

        c2->cd(6);
        t1->Draw("Sum$(hi.eta > 2.87 && hi.eta < 5.205)>>hn0",trig);
        t1->Draw("Sum$(hi.pt* (hi.eta > 2.87 && hi.eta < 5.205))>>he0",trig);

        c2 -> SaveAs("pdf/c2.pdf");
    }
    c1 -> SaveAs("pdf/c1.pdf");
    outf->Write();

}

