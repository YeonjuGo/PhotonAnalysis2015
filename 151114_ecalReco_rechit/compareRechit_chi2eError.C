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
#include "TStyle.h"
#include "TPaveStats.h"
#include "TROOT.h"
#include "TLatex.h"
#include "stdio.h"
#include "../yjUtility.h"
#include "../HIUtils/histoUtil.h"
#include "../yjCMSstyle.C"
//last forward run is 211256

void compareTwo_diff(TTree* t1=0 ,TTree* t2=0,TString var="pt", int nBins=10, double xMin=0, double xMax=10, const char* theCut="", const char* cap="");
void compareTwo_diff_eta(TTree* t1_ee=0 ,TTree* t2_ee=0, TTree* t1_eb=0 ,TTree* t2_eb=0, TString var="eta",int nBins=10, double xMin=0, double xMax=10, const char* theCut="", const char* cap="");
void compareTwo_ratio_eta(TTree* t1_ee=0 ,TTree* t2_ee=0, TTree* t1_eb=0 ,TTree* t2_eb=0, TString var="eta",int nBins=10, double xMin=0, double xMax=10, const char* theCut="", const char* cap="");

void compareTwo_ratio(TTree* t1=0 ,TTree* t2=0,TString var="pt", int nBins=10, double xMin=0, double xMax=10, const char* theCut="", TString type="photonMultifit", TString sam1="DATA", TString sam2="MC");
void draw_2D(TTree* t1=0, const char* xvar="chi2", const char* yvar="eError", int xbins=70, double xMin=0, double xMax=70, int ybins=50, double yMin=0, double yMax=1.5, const char* theCut="", TString cap="");

void compareRechit_chi2eError(const char* f1="/afs/cern.ch/work/y/ygo/private/PRODUCTION/CMSSW_7_5_5_patch2/src/EcalReco/test/HiForestAOD_photonData_multifit_chisquare.root",
        const char* f2="/afs/cern.ch/user/y/ygo/workspace/public/ecalLocalReco/forest_allqcdphoton30_multifit_ch2eError_755p1.root",
        TString type="photonMultifit", TString sam1="DATA",TString sam2="MC")
{
   // gROOT->Reset();
    //CMSstyle();
    SetHistTitleStyle();
    SetyjPadStyle();
    TH1::SetDefaultSumw2();
    gStyle->SetOptStat(0);  //0: don’t show statistic
    const char* fname[2];
    fname[0] = f1; 
    fname[1] = f2; 
    TFile* inf[2];
    TTree* tee[2];
    TTree* teb[2];
    TTree* tevt[2];
    TTree* tjet[2];
    TTree* tpho[2];
    TTree* tskim[2];
    for(int i=0;i<2;i++){
        inf[i] = TFile::Open(fname[i]);
        //inf[i] = new TFile(fname[i], "READ");
        tee[i] = (TTree*)inf[i]->Get(Form("rechitanalyzer/ee"));
        teb[i] = (TTree*)inf[i]->Get(Form("rechitanalyzer/eb"));
        tskim[i] = (TTree*)inf[i]->Get(Form("skimanalysis/HltTree"));
        //thbhe[i] = (TTree*)inf[i]->Get(Form("rechitAnalyzer/hbhe"));
        tevt[i] = (TTree*)inf[i]->Get("hiEvtAnalyzer/HiTree");
        tjet[i] = (TTree*)inf[i]->Get(Form("akPu4CaloJetAnalyzer/t"));
        tpho[i] = (TTree*)inf[i]->Get(Form("ggHiNtuplizer/EventTree"));
        tee[i]->AddFriend(tevt[i]);
        teb[i]->AddFriend(tevt[i]);
        tee[i]->AddFriend(tjet[i]);
        teb[i]->AddFriend(tjet[i]);
        tee[i]->AddFriend(tpho[i]);
        teb[i]->AddFriend(tpho[i]);
        tee[i]->AddFriend(tskim[i]);
        teb[i]->AddFriend(tskim[i]);
        tevt[i]->AddFriend(tjet[i]);
        tevt[i]->AddFriend(tpho[i]);
    }
    cout << "dd" << endl;
    int chi2Bin = 35;
    double chi2Max = 35; 
    int eErrBin = 20;
    double eErrMax_barrel = 0.1;
    double eErrMax_endcap = 0.3;
    int eBin = 100;
    double eMax = 300;
 //   int eBin = 100;
 //   double eMax = 50;

    const char* cut0 = "pcollisionEventSelection>0 && Max$(phoEt)>30";
    TString addcap0 = "_evtSel_leadingPho30";
    compareTwo_ratio(tee[0],tee[1],"e",eBin,0,eMax,Form("%s",cut0),type+addcap0, sam1, sam2);
    compareTwo_ratio(teb[0],teb[1],"e",eBin,0,eMax,Form("%s",cut0),type+addcap0, sam1, sam2);
 
/*    const char* c2ut = "pcollisionEventSelection>0 && Max$(phoEt)>30 && chi2==64";
    TString add2cap = "_evtSel_leadingPho30_chisq64";
    compareTwo_ratio(tee[0],tee[1],"e",eBin,0,eMax,Form("%s",c2ut),type+add2cap, sam1, sam2);
    compareTwo_ratio(teb[0],teb[1],"e",eBin,0,eMax,Form("%s",cut1),type+addcap1, sam1, sam2);
    compareTwo_ratio(tee[0],tee[1],"eError",eErrBin,0,eErrMax_endcap,Form("%s",cut1),type+addcap1, sam1, sam2);
    compareTwo_ratio(teb[0],teb[1],"eError",eErrBin,0,eErrMax_barrel,Form("%s",cut1),type+addcap1, sam1, sam2);
*/
    eBin = 100;
    eMax = 50;
     
    const char* cut2bar = "pcollisionEventSelection>0 && Max$(phoEt)>30 && eError>0.01 && eError<0.03";
    const char* cut2end = "pcollisionEventSelection>0 && Max$(phoEt)>30 && eError>0.01 && eError<0.3";
    TString addcap2 = "_evtSel_leadingPho30_eErrAroundBump";
    compareTwo_ratio(tee[0],tee[1],"e",eBin,0,eMax,Form("%s",cut2end),type+addcap2, sam1, sam2);
    compareTwo_ratio(teb[0],teb[1],"e",eBin,0,eMax,Form("%s",cut2bar),type+addcap2, sam1, sam2);
    compareTwo_ratio(tee[0],tee[1],"chi2",chi2Bin,0,chi2Max,Form("%s",cut2end),type+addcap2, sam1, sam2);
    compareTwo_ratio(teb[0],teb[1],"chi2",chi2Bin,0,chi2Max,Form("%s",cut2bar),type+addcap2, sam1, sam2);

    const char* cut3bar = "pcollisionEventSelection>0 && Max$(phoEt)>30 && (eError==0 || chi2==0 || chi2==64)";
    const char* cut3end = "pcollisionEventSelection>0 && Max$(phoEt)>30 && (eError==0 || chi2==0 || chi2==64)";
    TString addcap3 = "_evtSel_leadingPho30_eError0_OR_chisq0_OR_chisq64";
    compareTwo_ratio(tee[0],tee[1],"chi2",chi2Bin,0,chi2Max,Form("%s",cut3end),type+addcap3, sam1, sam2);
    compareTwo_ratio(teb[0],teb[1],"chi2",chi2Bin,0,chi2Max,Form("%s",cut3bar),type+addcap3, sam1, sam2);
#if 0
    compareTwo_ratio(tee[0],tee[1],"chi2",chi2Bin,0,chi2Max,Form("%s",cut0),type+addcap1, sam1, sam2);
    compareTwo_ratio(teb[0],teb[1],"chi2",chi2Bin,0,chi2Max,Form("%s",cut0),type+addcap1, sam1, sam2);
    compareTwo_ratio(tee[0],tee[1],"eError",eErrBin,0,eErrMax_endcap,Form("%s",cut0),type+addcap1, sam1, sam2);
    compareTwo_ratio(teb[0],teb[1],"eError",eErrBin,0,eErrMax_barrel,Form("%s",cut0),type+addcap1, sam1, sam2);
#endif
#if 0
//    addcap = "evtSel_leadingPho30_hiNtrack50";
//    const char* cut1 = "hiNtracks<50 && hiNtracks>0 && pcollisionEventSelection>0 && Max$(phoEt)>30";
    TString addcap1 = "evtSel_leadingPho30_hiHF1000";
    const char* cut1 = "hiHF<1000 && pcollisionEventSelection>0 && Max$(phoEt)>30";
    compareTwo_ratio(tee[0],tee[1],"chi2",chi2Bin,0,chi2Max,Form("%s",cut1),Form("%s_%s",cap, Form("%s",addcap1.Data()) ));
    compareTwo_ratio(teb[0],teb[1],"chi2",chi2Bin,0,chi2Max,Form("%s",cut1),Form("%s_%s",cap, Form("%s",addcap1.Data()) ));
    compareTwo_ratio(tee[0],tee[1],"eError",eErrBin,0,eErrMax_endcap,Form("%s",cut1),Form("%s_%s",cap, Form("%s",addcap1.Data()) ));
    compareTwo_ratio(teb[0],teb[1],"eError",eErrBin,0,eErrMax_barrel,Form("%s",cut1),Form("%s_%s",cap, Form("%s",addcap1.Data()) ));
  
//    addcap = "evtSel_leadingPho30_hiNtrack1000";
//    const char* cut2 = "hiNtracks>1000 && pcollisionEventSelection>0 && Max$(phoEt)>30";
    TString addcap2 = "evtSel_leadingPho30_hiHF2500";
    const char* cut2 = "hiHF>2500 && pcollisionEventSelection>0 && Max$(phoEt)>30";
    compareTwo_ratio(tee[0],tee[1],"chi2",chi2Bin,0,chi2Max,Form("%s",cut2),Form("%s_%s",cap, Form("%s",addcap2.Data()) ));
    compareTwo_ratio(teb[0],teb[1],"chi2",chi2Bin,0,chi2Max,Form("%s",cut2),Form("%s_%s",cap, Form("%s",addcap2.Data()) ));
    compareTwo_ratio(tee[0],tee[1],"eError",eErrBin,0,eErrMax_endcap,Form("%s",cut2),Form("%s_%s",cap, Form("%s",addcap2.Data()) ));
    compareTwo_ratio(teb[0],teb[1],"eError",eErrBin,0,eErrMax_barrel,Form("%s",cut2),Form("%s_%s",cap, Form("%s",addcap2.Data()) ));
#endif
    eErrBin = 100;
    eErrMax_barrel = 0.5;
    eErrMax_endcap = 2.0;
    eBin = 100;
    eMax = 150;
#if 0
    for(int i=0;i<2;i++){
        TString sample = "_DATA";
        if(i==1) sample = "_MC";
        TString tmpSt = addcap0 + sample;
        draw_2D(teb[i],"chi2","eError", 70,0,70, eErrBin,0,eErrMax_barrel,cut0,tmpSt);
        draw_2D(tee[i],"chi2","eError", 70,0,70, eErrBin,0,eErrMax_endcap,cut0,tmpSt);
        draw_2D(teb[i],"e","eError", eBin,0,eMax, eErrBin,0,eErrMax_barrel,cut0,tmpSt);
        draw_2D(tee[i],"e","eError", eBin,0,eMax, eErrBin,0,eErrMax_endcap,cut0,tmpSt);
        tmpSt = addcap1 + sample;
        draw_2D(teb[i],"chi2","eError", 70,0,70, eErrBin,0,eErrMax_barrel,cut1,tmpSt);
        draw_2D(tee[i],"chi2","eError", 70,0,70, eErrBin,0,eErrMax_endcap,cut1,tmpSt);
        draw_2D(teb[i],"e","eError", eBin,0,eMax, eErrBin,0,eErrMax_barrel,cut1,tmpSt);
        draw_2D(tee[i],"e","eError", eBin,0,eMax, eErrBin,0,eErrMax_endcap,cut1,tmpSt);
        tmpSt = addcap2 + sample;
        draw_2D(teb[i],"chi2","eError", 70,0,70, eErrBin,0,eErrMax_barrel,cut2,tmpSt);
        draw_2D(tee[i],"chi2","eError", 70,0,70, eErrBin,0,eErrMax_endcap,cut2,tmpSt);
        draw_2D(teb[i],"e","eError", eBin,0,eMax, eErrBin,0,eErrMax_barrel,cut2,tmpSt);
        draw_2D(tee[i],"e","eError", eBin,0,eMax, eErrBin,0,eErrMax_endcap,cut2,tmpSt);
    }
#endif
}
void draw_2D(TTree* t1, const char* xvar, const char* yvar, int xbins, double xMin, double xMax, int ybins, double yMin, double yMax, const char* theCut, TString cap){
    TCanvas* c1=  new TCanvas(Form("c2D_%s_%s_%s_%s",xvar,yvar,t1->GetName(),cap.Data()),"", 400,400);
    c1->cd();
    TH2D* h2D = new TH2D(Form("h2D_%s_%s_%s_%s",xvar,yvar,t1->GetName(),cap.Data()),Form("%s_%s;%s;%s",t1->GetName(),cap.Data(),xvar,yvar),xbins,xMin,xMax,ybins,yMin,yMax);
    t1->Draw(Form("%s:%s>>%s",yvar,xvar,h2D->GetName()),Form("%s",theCut),"colz");
    gPad->SetLogz();
    c1-> SaveAs(Form("pdf/%s.pdf",c1->GetName())); 
}

void compareTwo_ratio(TTree* t1, TTree* t2, TString var, int nBins, double xMin, double xMax, const char* theCut, TString type, TString sam1, TString sam2)  {
    TString tempcap = type + "_" + sam1 + "_" + sam2;
    const char* cap = tempcap.Data();
    TCanvas* cc=  new TCanvas(Form("c_%s_%s_%s",var.Data(),t1->GetName(),cap),"", 400,800);
    cc->Divide(1,2);
    cc->cd(1);
    TH1D* h1 = new TH1D(Form("h1_%s_%s_%s",var.Data(),t1->GetName(),cap), Form(";%s;",var.Data()), nBins,xMin,xMax);
    TH1D* h2 = (TH1D*)h1->Clone(Form("h2_%s_%s_%s",var.Data(),t1->GetName(),cap));
    h1->Sumw2();
    h2->Sumw2();
    t1->Draw(Form("%s>>%s",var.Data(),h1->GetName()), theCut);
    t2->Draw(Form("%s>>%s",var.Data(),h2->GetName()), theCut);
    h1->Scale( 1. / h1->Integral());
    h2->Scale( 1. / h2->Integral());

    // default plotting options
    h1->SetMarkerStyle(kFullSquare);
    h2->SetMarkerStyle(kFullCircle);

    h1->SetMarkerColor(kRed);
    h2->SetMarkerColor(kBlack);

    h2->SetMarkerSize(h1->GetMarkerSize()*0.8);     // to distinguish between points when they overlap

    h1->SetLineColor(2);
    h1->SetMarkerColor(2);
    h2->SetLineColor(1);
    h2->SetMarkerColor(1);

    // set Y-axis ranges
    // default Y-axis range is that of h1
    // make sure that both plots will not run out of y-axis
    double max1 = h1->GetMaximum();
    double max2 = h2->GetMaximum();
    double min1 = h1->GetMinimum();
    double min2 = h2->GetMinimum();

    if (max2 > max1)
        h1->SetMaximum(max2+TMath::Abs(max2)*0.2);

    //if (min2 < min1)
        //h1->SetMinimum(min2-TMath::Abs(min2)*0.2);
    gPad->SetLogy();
    //if(!(min1<=0 || min2<=0)) gPad->SetLogy();
    //cleverRange(h1,h2);
    h1->DrawCopy("hist");
    h2->DrawCopy("same");
    TLegend* leg = new TLegend(0.6,0.75,0.9,0.9);
    legStyle(leg);
  //  leg->AddEntry(h1,Form("Global"), "l");
  //  leg->AddEntry(h2,Form("Multifit"));
//   leg->AddEntry(h1,Form("Jet80Data"), "l");
//    leg->AddEntry(h2,Form("Photon20Data"));
    leg->AddEntry(h1,Form("%s",sam1.Data()), "l");
    leg->AddEntry(h2,Form("%s",sam2.Data()));
    leg->Draw("same");
    drawText(Form("%s",t1->GetName()),0.15,0.95);
    drawText(Form("%s",cap),0.25,0.95);
    cc->cd(2);
    h2->SetAxisRange(0.0,2.0,"Y");
    h2->Divide(h1);
    h2->SetYTitle(Form("%s/%s Ratio",sam1.Data(),sam2.Data()));
    //h2->SetYTitle("Jet80Data/Photon20Data Ratio");
    //h2->SetYTitle("GlobalReco/MultiFit Ratio");
    h2->DrawCopy();
    jumSun(xMin,1,xMax,1);
    cc-> SaveAs(Form("pdf/%s.pdf",cc->GetName())); 
    //return c; 
}


void compareTwo_diff_eta(TTree* t1_ee, TTree* t2_ee, TTree* t1_eb, TTree* t2_eb, TString var, int nBins, double xMin, double xMax, const char* theCut, const char* cap)  {
    cout << "s" << endl;
    TCanvas* cc=  new TCanvas(Form("c_%s_%s_%s",var.Data(),t1_ee->GetName(),cap),"", 400,800);
    cc->Divide(1,2);
    cout << "s" << endl;
    cc->cd(1);
    TH1D* h1 = new TH1D(Form("h1_%s_%s_%s",var.Data(),t1_ee->GetName(),cap), Form(";%s;",var.Data()), nBins,xMin,xMax);
    TH1D* h2 = (TH1D*)h1->Clone(Form("h2_%s_%s_%s",var.Data(),t1_ee->GetName(),cap));
    h1->Sumw2();
    h2->Sumw2();
    cout << "s" << endl;
    t1_ee->Draw(Form("%s>>%s",var.Data(),h1->GetName()), theCut);
    t2_ee->Draw(Form("%s>>%s",var.Data(),h2->GetName()), theCut);
    t1_eb->Draw(Form("%s>>+%s",var.Data(),h1->GetName()), theCut);
    cout << "s" << endl;
    t2_eb->Draw(Form("%s>>+%s",var.Data(),h2->GetName()), theCut);
    cout << "s" << endl;

    cout << "s" << endl;
    Int_t t1eeEntries = t1_ee->GetEntries(Form("%s",theCut));
    cout << "s" << endl;
    Int_t t2eeEntries = t2_ee->GetEntries(Form("%s",theCut));
    cout << "s" << endl;
    Int_t t1ebEntries = t1_eb->GetEntries(Form("%s",theCut));
    cout << "s" << endl;
    Int_t t2ebEntries = t2_eb->GetEntries(Form("%s",theCut));
    cout << "s" << endl;
    h1->Scale( 1. / (t1eeEntries+t1ebEntries));
    cout << "s" << endl;
    h2->Scale( 1. / (t2eeEntries+t2ebEntries));
    cout << "s" << endl;
    //h1->Scale( 1. / h1->Integral());
    //h2->Scale( 1. / h2->Integral());

    // default plotting options
    h1->SetMarkerStyle(kFullSquare);
    h2->SetMarkerStyle(kFullCircle);

    h1->SetMarkerColor(kRed);
    h2->SetMarkerColor(kBlack);

    h2->SetMarkerSize(h1->GetMarkerSize()*0.8);     // to distinguish between points when they overlap

    h1->SetLineColor(2);
    h1->SetMarkerColor(2);
    h2->SetLineColor(1);
    h2->SetMarkerColor(1);

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

    cleverRange(h1,h2);
    h2->DrawCopy();
    h1->DrawCopy("hist same");
/*
    gPad->Update();
    TPaveStats *statsbox1 = (TPaveStats*)gPad->GetPrimitive("stats");
    double sth1 = statsbox1->GetY1NDC();
    double sth2 = statsbox1->GetY2NDC();
    double stnewh1 = 2 * sth1 - sth2;
    double stnewh2 = sth1;
    statsbox->SetY1NDC(stnewh1);
    statsbox->SetY2NDC(stnewh2);
   // hZmassglbEm20->Draw("sames");
*/

    drawText(Form("%s",cap),0.45,0.95);
    TLegend* leg = new TLegend(0.6,0.75,0.9,0.9);
    //legStyle(leg);
    leg->AddEntry(h1,Form("Global"), "l");
    leg->AddEntry(h2,Form("Multifit"));
    leg->Draw();
    cc->cd(2);
    //h2->SetAxisRange(0.5,1.5,"Y");
    //h1->SetYTitle("GlobalReco - MultiFit");
    for(int i=0;i< h2->GetNbinsX();i++){
        h2->SetBinError(i+1,0.00000001);
        //cout << "h1 : " << h1->GetBinContent(i+1)<< ", h2 : " << h2->GetBinContent(i+1) << endl;
    }
    h2->Add(h1,-1);
    h2->SetYTitle("MultiFit - GlobalReco");
    h2->DrawCopy("p");
    jumSun(xMin,0,xMax,0);
    cc-> SaveAs(Form("pdf/%s.pdf",cc->GetName())); 
//    delete h1;
//    delete h2;
    //return c; 
}


void compareTwo_diff(TTree* t1, TTree* t2, TString var, int nBins, double xMin, double xMax, const char* theCut, const char* cap)  {
    TCanvas* cc=  new TCanvas(Form("c_%s_%s_%s",var.Data(),t1->GetName(),cap),"", 400,800);
    cc->Divide(1,2);
    cc->cd(1);
    TH1D* h1 = new TH1D(Form("h1_%s_%s_%s",var.Data(),t1->GetName(),cap), Form(";%s;",var.Data()), nBins,xMin,xMax);
    TH1D* h2 = (TH1D*)h1->Clone(Form("h2_%s_%s_%s",var.Data(),t1->GetName(),cap));
    h1->Sumw2();
    h2->Sumw2();
    t1->Draw(Form("%s>>%s",var.Data(),h1->GetName()), theCut);
    t2->Draw(Form("%s>>%s",var.Data(),h2->GetName()), theCut);
    Int_t t1Entries = t1->GetEntries(Form("%s",theCut));
    Int_t t2Entries = t2->GetEntries(Form("%s",theCut));
    h1->Scale( 1. / (t1Entries));
    h2->Scale( 1. / (t2Entries));
    

    h1->Scale( 1. / h1->Integral());
    h2->Scale( 1. / h2->Integral());

    // default plotting options
    h1->SetMarkerStyle(kFullSquare);
    h2->SetMarkerStyle(kFullCircle);

    h1->SetMarkerColor(kRed);
    h2->SetMarkerColor(kBlack);

    h2->SetMarkerSize(h1->GetMarkerSize()*0.8);     // to distinguish between points when they overlap

    h1->SetLineColor(2);
    h1->SetMarkerColor(2);
    h2->SetLineColor(1);
    h2->SetMarkerColor(1);

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

    cleverRange(h1,h2);
    h2->DrawCopy();
    h1->DrawCopy("hist same");
    drawText(Form("%s",cap),0.45,0.95);
    drawText(Form("%s",t1->GetName()),0.25,0.95);
    TLegend* leg = new TLegend(0.6,0.75,0.9,0.9);
    //legStyle(leg);
//    leg->AddEntry(h1,Form("Global"), "l");
//    leg->AddEntry(h2,Form("Multifit"));
//    leg->AddEntry(h1,Form("Jet80Data"), "l");
//    leg->AddEntry(h2,Form("Photon20Data"));
    leg->AddEntry(h1,Form("PhotonData"), "l");
    leg->AddEntry(h2,Form("PhotonMC"));
    leg->Draw();
    cc->cd(2);
    //h2->SetAxisRange(0.5,1.5,"Y");
    h2->Add(h1,-1);
    //h1->SetYTitle("GlobalReco - MultiFit");
    for(int i=0;i< h2->GetNbinsX();i++){
        h2->SetBinError(i+1,0.000000000000000000001);
    }
    //h2->SetYTitle("MultiFit - GlobalReco");
    //h2->SetYTitle("Photon20Data - Jet80Data");
    h2->SetYTitle("PhotonData - PhotonMC");
    h2->DrawCopy("p");
    jumSun(xMin,0,xMax,0);
    cc-> SaveAs(Form("pdf/%s.pdf",cc->GetName())); 
    //delete h1;
    //delete h2;
    //return c; 
}

void compareTwo_ratio_eta(TTree* t1_ee, TTree* t2_ee, TTree* t1_eb, TTree* t2_eb, TString var, int nBins, double xMin, double xMax, const char* theCut, const char* cap)  {
    TCanvas* cc=  new TCanvas(Form("c_%s_%s",var.Data(),t1_ee->GetName()),"", 400,800);
    cc->Divide(1,2);
    cc->cd(1);
    TH1D* h1 = new TH1D(Form("h1_%s_%s",var.Data(),t1_ee->GetName()), Form("%s;%s;",cap,var.Data()), nBins,xMin,xMax);
    TH1D* h2 = (TH1D*)h1->Clone(Form("h2_%s_%s",var.Data(),t1_ee->GetName()));
    h1->Sumw2();
    h2->Sumw2();
    t1_ee->Draw(Form("%s>>%s",var.Data(),h1->GetName()), theCut);
    t2_ee->Draw(Form("%s>>%s",var.Data(),h2->GetName()), theCut);
    t1_eb->Draw(Form("%s>>+%s",var.Data(),h1->GetName()), theCut);
    t2_eb->Draw(Form("%s>>+%s",var.Data(),h2->GetName()), theCut);
    h1->Scale( 1. / h1->Integral());
    h2->Scale( 1. / h2->Integral());

    // default plotting options
    h1->SetMarkerStyle(kFullSquare);
    h2->SetMarkerStyle(kFullCircle);

    h1->SetMarkerColor(kRed);
    h2->SetMarkerColor(kBlack);

    h2->SetMarkerSize(h1->GetMarkerSize()*0.8);     // to distinguish between points when they overlap

    h1->SetLineColor(2);
    h1->SetMarkerColor(2);
    h2->SetLineColor(1);
    h2->SetMarkerColor(1);

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

    cleverRange(h1,h2);
    h2->DrawCopy();
    h1->DrawCopy("hist same");
/*
    gPad->Update();
    TPaveStats *statsbox1 = (TPaveStats*)gPad->GetPrimitive("stats");
    double sth1 = statsbox1->GetY1NDC();
    double sth2 = statsbox1->GetY2NDC();
    double stnewh1 = 2 * sth1 - sth2;
    double stnewh2 = sth1;
    statsbox->SetY1NDC(stnewh1);
    statsbox->SetY2NDC(stnewh2);
    // hZmassglbEm20->Draw("sames");
*/
    TLegend* leg = new TLegend(0.6,0.75,0.9,0.9);
    //legStyle(leg);
    leg->AddEntry(h1,Form("Global"), "l");
    leg->AddEntry(h2,Form("Multifit"));
    leg->Draw();
    cc->cd(2);
    h2->SetAxisRange(0.5,1.5,"Y");
    h2->Divide(h1);
    h2->SetYTitle("GlobalReco/MultiFit Ratio");
    h2->DrawCopy();
    jumSun(xMin,1,xMax,1);
    cc-> SaveAs(Form("pdf/%s_%s.pdf",cc->GetName(),cap));
}


