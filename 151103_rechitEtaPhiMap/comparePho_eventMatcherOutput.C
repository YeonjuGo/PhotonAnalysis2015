#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TLatex.h"
#include "stdio.h"
#include "../yjUtility.h"
//#include "../HIUtils/histoUtil.h"
//last forward run is 211256


//const Double_t hfBins[] = {0, 20, 30, 1000}; //last entry is upper bound on last bin
//const Int_t nhfBins = 3;

//int returnHFBin(double hf);

void compareTwo(TTree* t1=0 ,TTree* t2=0,TString var="pt", int nBins=10, double xMin=0, double xMax=10, const char* theCut="");
void compareTwo(TH1* h1=0, TH1* h2=0,double xmin=0.0, double xmax=200.0);
void comparePho_eventMatcherOutput()
{
    TH1::SetDefaultSumw2();

    TFile *inf1 = new TFile("/afs/cern.ch/work/y/ygo/private/PhotonAnalysis2015/151103_rechitEtaPhiMap/skimFiles/eventMatched_53X_75X_biransSample.root");
    //TFile *inf1 = new TFile("/afs/cern.ch/work/y/ygo/private/PhotonAnalysis2015/151103_rechitEtaPhiMap/skimFiles/eventMatched_53X_75X.root");
    TTree* tevt = (TTree*)inf1->Get("HiTree");
    TTree* tpho[2];
    const char* ver = "53";
    for(int i=0;i<2;i++){
        if(i==1) ver = "75";
        tpho[i] = (TTree*)inf1->Get(Form("pho%sx",ver));
    }

    const double xMin(0.0), xMax(200.0);
    const int xBin=100;

    const char* ptCut = "";
    const char* eCut = "";
    const char* totCut[2][2];//[sample][eta]
//    totCut[0][0] = "abs(eta)>1.44";
//    totCut[1][0] = "abs(phoEta)>1.44";

        totCut[0][0] = "abs(eta)<1.44";
    totCut[1][0] = "abs(phoEta)<1.44";
//    totCut[0][1] = "abs(eta)>1.44";
//    totCut[1][1] = "abs(phoEta)>1.44";
//    totCut[0][1] = "abs(eta)>1.44";
//    totCut[1][1] = "abs(phoEta)>1.44";
    //const char* totCut = Form("%s && %s && %s",ptCut eCut);
    //cout << "totCut :" << totCut->GetName() << endl; 
    TH1D* phopt[2];
    TH1D* phoE[2];
    TH1D* phoRawE[2];
    TH1D* phoEoverRawE[2];
    TH1D* phosigmaIetaIeta[2];
    TH1D* phor9[2];
    TH1D* phophiwidth[2];
    TH1D* phoetawidth[2];
    TH1D* phobrem[2];

    for(int i=0; i<2;i++){
        phopt[i] = new TH1D(Form("phopt_%d",i),";Photon p_{T} (GeV);",xBin,xMin,xMax);
        phoE[i] = new TH1D(Form("phoE_%d",i),";Photon Energy (GeV);",xBin,xMin,xMax);
        phoRawE[i] = new TH1D(Form("phoRawE_%d",i),";Photon Raw Energy (GeV);",xBin,xMin,xMax);
        phoEoverRawE[i] = new TH1D(Form("phoEOverRawE_%d",i),";Photon Energy/RawEnergy;",xBin,0.95,1.2);
    }
    tpho[0]->Draw(Form("pt>>%s",phopt[0]->GetName()), Form("%s",totCut[0][0]));
    tpho[1]->Draw(Form("phoEt>>%s",phopt[1]->GetName()),Form("%s",totCut[1][0]));
    tpho[0]->Draw(Form("energy>>%s",phoE[0]->GetName()), Form("%s",totCut[0][0]));
    tpho[1]->Draw(Form("phoE>>%s",phoE[1]->GetName()), Form("%s",totCut[1][0]));
    tpho[0]->Draw(Form("rawEnergy>>%s",phoRawE[0]->GetName()), Form("%s",totCut[0][0]));
    tpho[1]->Draw(Form("phoSCRawE>>%s",phoRawE[1]->GetName()), Form("%s",totCut[1][0]));
    tpho[0]->Draw(Form("energy/rawEnergy>>%s",phoEoverRawE[0]->GetName()), Form("%s",totCut[0][0]));
    tpho[1]->Draw(Form("phoE/phoSCRawE>>%s",phoEoverRawE[1]->GetName()), Form("%s",totCut[1][0]));

    compareTwo(phopt[0], phopt[1],xMin,xMax);
    compareTwo(phoE[0], phoE[1],xMin,xMax);
    compareTwo(phoRawE[0], phoRawE[1],xMin,xMax);
    compareTwo(phoEoverRawE[0], phoEoverRawE[1],0.95,1.2);

/*
    TH1D* phoeta[2];
    TH1D* phophi[2];
    TH1D* phohoe[2];
    TH1D* phoecalIso[2];
    TH1D* phohcalIso[2];
    TH1D* photrackerIso[2];
    TH1D* phosigmaIetaIeta[2];
    TH1D* phor9[2];
    TH1D* phophiwidth[2];
    TH1D* phoetawidth[2];
    TH1D* phobrem[2];
*/



//    compareTwo(t1, t2, "run",100, 210400, 212000, "abs(photonEta)<1.44 && photonEt>40");
//    compareTwo(t1, t2, "hf4Sum", 20, 0, 150, "abs(photonEta)<1.44 && photonEt>40");
//    compareTwo(t1, t2, "hovere", 20, 0, 0.2, "abs(photonEta)<1.44 && photonEt>40");
 //   compareTwo(t1, t2, "ecalIso", 20,-5,20, "abs(photonEta)<1.44 && photonEt>40");
  //  compareTwo(t1, t2, "hcalIso", 20,-5,20, "abs(photonEta)<1.44 && photonEt>40");
   // compareTwo(t1, t2, "trackIso", 20, -5, 20, "abs(photonEta)<1.44 && photonEt>40");
    //compareTwo(t1, t2, "sigmaIetaIeta", 20, 0, 0.025, "abs(photonEta)<1.44 && photonEt>40 && ecalIso < 4.2  &&  hcalIso < 2.2  &&  trackIso < 2");
//    compareTwo(t1, t2, "photonEt", 20, 0, 200,  "abs(photonEta)<1.44 && photonEt>40 && ecalIso < 4.2  &&  hcalIso < 2.2  &&  trackIso < 2");



}

void compareTwo(TH1* h1, TH1* h2, double xmin, double xmax)  {
    TCanvas* c=  new TCanvas(Form("c_%s",h1->GetName()),"", 400,800);
    c->Divide(1,2);
    c->cd(1);
    h1->Sumw2();
    h2->Sumw2();
    //  h1->Scale( 1. / h1->Integral());
    //  h2->Scale( 1. / h2->Integral());

    // default plotting options
    h1->SetMarkerStyle(kFullSquare);
    h2->SetMarkerStyle(kFullCircle);

    h1->SetMarkerColor(kBlack);
    h2->SetMarkerColor(kRed);

    h2->SetMarkerSize(h1->GetMarkerSize()*0.8);     // to distinguish between points when they overlap

    h1->SetLineColor(4);
    h1->SetMarkerColor(4);
    h2->SetLineColor(2);
    h2->SetMarkerColor(2);

    // set Y-axis ranges
    // default Y-axis range is that of h1
    // make sure that both plots will not run out of y-axis
    double max1 = h1->GetMaximum();
    double max2 = h2->GetMaximum();
    double min1 = h1->GetMinimum();
    double min2 = h2->GetMinimum();

    if (max2 > max1)
        h1->SetMaximum(max2+TMath::Abs(max2)*0.2);

    if (min2 < min1)
        h1->SetMinimum(min2-TMath::Abs(min2)*0.2);

    h1->DrawCopy("hist");
    h2->DrawCopy("hist same");
    c->cd(2);
    h2->SetAxisRange(0,3,"Y");
    h2->Divide(h1);
    h2->SetYTitle("75X/53X Ratio");
    h2->DrawCopy();
    jumSun(xmin,1,xmax,1);
    //c -> SaveAs(Form("pdf/%s_biransSample_endcap.pdf",c->GetName())); 
    c -> SaveAs(Form("pdf/%s_biransSample.pdf",c->GetName())); 
    //c -> SaveAs(Form("pdf/%s_biransSample_barrel.pdf",c->GetName())); 
}

void compareTwo(TTree* t1, TTree* t2, TString var, int nBins, double xMin, double xMax, const char* theCut)  {
//TCanvas* compareTwo(TTree* t1, TTree* t2, TString var, int nBins, double xMin, double xMax, const char* theCut)  {
    TCanvas* c=  new TCanvas(Form("c_%s_%s",var.Data(),t1->GetName()),"", 400,800);
    c->Divide(1,2);
    c->cd(1);
    TH1D* h1 = new TH1D(Form("h1_%s_%s",var.Data(),t1->GetName()), Form(";%s;",var.Data()), nBins,xMin,xMax);
    TH1D* h2 = (TH1D*)h1->Clone(Form("h2_%s_%s",var.Data(),t1->GetName()));
    h1->Sumw2();
    h2->Sumw2();
    t1->Draw(Form("%s>>%s",var.Data(),h1->GetName()), theCut);
    t2->Draw(Form("%s>>%s",var.Data(),h2->GetName()), theCut);
    //  h1->Scale( 1. / h1->Integral());
    //  h2->Scale( 1. / h2->Integral());

    // default plotting options
    h1->SetMarkerStyle(kFullSquare);
    h2->SetMarkerStyle(kFullCircle);

    h1->SetMarkerColor(kBlack);
    h2->SetMarkerColor(kRed);

    h2->SetMarkerSize(h1->GetMarkerSize()*0.8);     // to distinguish between points when they overlap

    h1->SetLineColor(4);
    h1->SetMarkerColor(4);
    h2->SetLineColor(2);
    h2->SetMarkerColor(2);

    // set Y-axis ranges
    // default Y-axis range is that of h1
    // make sure that both plots will not run out of y-axis
    double max1 = h1->GetMaximum();
    double max2 = h2->GetMaximum();
    double min1 = h1->GetMinimum();
    double min2 = h2->GetMinimum();

    if (max2 > max1)
        h1->SetMaximum(max2+TMath::Abs(max2)*0.2);

    if (min2 < min1)
        h1->SetMinimum(min2-TMath::Abs(min2)*0.2);

    h2->DrawCopy();
    h1->DrawCopy("hist same");
    c->cd(2);
    h2->SetAxisRange(0,3,"Y");
    h2->Divide(h1);
    h2->SetYTitle("New/Old Ratio");
    h2->DrawCopy();
    jumSun(xMin,1,xMax,1);
    c -> SaveAs(Form("pdf/%s.pdf",c->GetName())); 
//    return c; 
}
