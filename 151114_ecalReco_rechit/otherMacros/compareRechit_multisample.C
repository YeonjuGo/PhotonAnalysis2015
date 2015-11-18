#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TChain.h"
#include "TString.h"
#include "TROOT.h"
#include "TLatex.h"
#include "stdio.h"
#include <vector>
#include "TBranch.h"
#include "TList.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "THnSparse.h"
#include "Math/QuantFuncMathCore.h"
#include "Math/SpecFuncMathCore.h"
#include "Math/ProbFunc.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include <vector>
#include <iostream>
#include <fstream>
#include "../yjUtility.h"
#include "../HIUtils/histoUtil.h"
#include "../CMSstyle.C"

//last forward run is 211256
using namespace std;

void compareTwo_diff(TTree* t1=0 ,TTree* t2=0,TString var="pt", int nBins=10, double xMin=0, double xMax=10, const char* cap="", const char* theCut="");
void compareTwo_ratio(TTree* t1=0 ,TTree* t2=0,TString var="pt", int nBins=10, double xMin=0, double xMax=10, const char* cap="", const char* theCut="");
void compareTwo_diff_eta(TTree* t1_ee=0 ,TTree* t2_ee=0, TTree* t1_eb=0 ,TTree* t2_eb=0, TString var="eta",int nBins=10, double xMin=0, double xMax=10, const char* cap="", const char* theCut="");

void compareRechit_multisample(std::vector<std::string> &infnames1,
                                std::vector<std::string> &infnames2, const char* cap="")
{

    CMSstyle();
    TChain* tevt[2];
    TChain* tee[2];
    TChain* teb[2];
    for(int i=0;i<2;i++){
        tevt[i]= new TChain("hiEvtAnalyzer/HiTree");
        tee[i]= new TChain("rechitanalyzer/ee");
        teb[i]= new TChain("rechitanalyzer/eb");
    }

    size_t startFile = 0;
    size_t nfiles1 = infnames1.size();
    size_t nfiles2 = infnames2.size();

    Printf("start file: %zu  nfiles: %zu",startFile,nfiles1);
    for(size_t i=startFile; i<nfiles1; i++) {
        TString tmp(infnames1[i].c_str());
        if(tmp.IsNull()) continue;
        tee[0]->Add(infnames1[i].c_str());
        teb[0]->Add(infnames1[i].c_str());
        tevt[0]->Add(infnames1[i].c_str());
        Printf("i: %zu  %s",i,infnames1[i].c_str());
    }
    
    Printf("start file: %zu  nfiles: %zu",startFile,nfiles2);
    for(size_t i=startFile; i<nfiles2; i++) {
        TString tmp(infnames2[i].c_str());
        if(tmp.IsNull()) continue;
        tee[1]->Add(infnames2[i].c_str());
        teb[1]->Add(infnames2[i].c_str());
        tevt[1]->Add(infnames2[i].c_str());
        Printf("i: %zu  %s",i,infnames2[i].c_str());
    }

    //compareTwo_diff(tee[0],tee[1],"e",50,0,3, cap);
    //compareTwo_diff(tee[0],tee[1],"et",50,0,3, cap);
    compareTwo_diff(tee[0],tee[1],"phi",50,-3.14,3.14,cap);
    compareTwo_diff_eta(tee[0],tee[1],teb[0],teb[1],"eta",50,-3,3,cap);
    //compareTwo_diff(teb[0],teb[1],"e",50,0,1.5,"",cap);
    //compareTwo_diff(teb[0],teb[1],"et",50,0,1.5,cap);
    //compareTwo_diff(teb[0],teb[1],"phi",50,-3.14,3.14,cap);

}

void compareTwo_diff_eta(TTree* t1_ee, TTree* t2_ee, TTree* t1_eb, TTree* t2_eb, TString var, int nBins, double xMin, double xMax, const char* cap, const char* theCut)  {
    TCanvas* cc=  new TCanvas(Form("c_%s_%s",var.Data(),t1_ee->GetName()),"", 400,800);
    cc->Divide(1,2);
    cc->cd(1);
    TH1D* h1 = new TH1D(Form("h1_%s_%s",var.Data(),t1_ee->GetName()), Form(";%s;",var.Data()), nBins,xMin,xMax);
    TH1D* h2 = (TH1D*)h1->Clone(Form("h2_%s_%s",var.Data(),t1_ee->GetName()));
    h1->Sumw2();
    h2->Sumw2();
    t1_ee->Draw(Form("%s>>%s",var.Data(),h1->GetName()), theCut);
    t2_ee->Draw(Form("%s>>%s",var.Data(),h2->GetName()), theCut);
    t1_eb->Draw(Form("%s>>+%s",var.Data(),h1->GetName()), theCut);
    t2_eb->Draw(Form("%s>>+%s",var.Data(),h2->GetName()), theCut);
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
        cout << "h1 : " << h1->GetBinContent(i+1)<< ", h2 : " << h2->GetBinContent(i+1) << endl;
    }
    h2->Add(h1,-1);
    h2->SetYTitle("MultiFit - GlobalReco");
    h2->DrawCopy("p");
    jumSun(xMin,0,xMax,0);
    cc-> SaveAs(Form("pdf/%s_%s.pdf",cap,cc->GetName())); 
    //return c; 
}


void compareTwo_diff(TTree* t1, TTree* t2, TString var, int nBins, double xMin, double xMax, const char* cap, const char* theCut)  {
    TCanvas* cc=  new TCanvas(Form("c_%s_%s",var.Data(),t1->GetName()),"", 400,800);
    cc->Divide(1,2);
    cc->cd(1);
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
    TLegend* leg = new TLegend(0.6,0.75,0.9,0.9);
    //legStyle(leg);
    leg->AddEntry(h1,Form("Global"), "l");
    leg->AddEntry(h2,Form("Multifit"));
    leg->Draw();
    cc->cd(2);
    //h2->SetAxisRange(0.5,1.5,"Y");
    h2->Add(h1,-1);
    //h1->SetYTitle("GlobalReco - MultiFit");
    for(int i=0;i< h2->GetNbinsX();i++){
        h2->SetBinError(i+1,0.000000000000000000001);
    }
    h2->SetYTitle("MultiFit - GlobalReco");
    h2->DrawCopy("p");
    jumSun(xMin,0,xMax,0);
    cc-> SaveAs(Form("pdf/%s_%s.pdf",cap,cc->GetName())); 
    //return c; 
}


void compareTwo_ratio(TTree* t1, TTree* t2, TString var, int nBins, double xMin, double xMax, const char* cap, const char* theCut)  {
    TCanvas* cc=  new TCanvas(Form("c_%s_%s",var.Data(),t1->GetName()),"", 400,800);
    cc->Divide(1,2);
    cc->cd(1);
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
    cc-> SaveAs(Form("pdf/%s_%s.pdf",cap,cc->GetName())); 
    //return c; 
}

