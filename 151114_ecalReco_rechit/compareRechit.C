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
void compareTwo_ratio(TTree* t1=0 ,TTree* t2=0,TString var="pt", int nBins=10, double xMin=0, double xMax=10, const char* theCut="", const char* cap="");
void compareTwo_diff_eta(TTree* t1_ee=0 ,TTree* t2_ee=0, TTree* t1_eb=0 ,TTree* t2_eb=0, TString var="eta",int nBins=10, double xMin=0, double xMax=10, const char* theCut="", const char* cap="");
void compareTwo_ratio_eta(TTree* t1_ee=0 ,TTree* t2_ee=0, TTree* t1_eb=0 ,TTree* t2_eb=0, TString var="eta",int nBins=10, double xMin=0, double xMax=10, const char* theCut="", const char* cap="");


void compareRechit(const char* f1="/afs/cern.ch/work/y/ygo/public/ecalLocalReco/forest_jet80_2011data_755p1_global.root", 
        const char* f2="/afs/cern.ch/work/y/ygo/public/ecalLocalReco/forest_jet80_2011data_755p1_multifit.root",
        const char* cap="jet80Data")
{
    gROOT->Reset();
    CMSstyle();
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
    const char* ver = "GlobalReco";
    for(int i=0;i<2;i++){
        if(i==1) ver = "MultiFit";
        inf[i] = TFile::Open(fname[i]);
        //inf[i] = new TFile(fname[i], "READ");
        tee[i] = (TTree*)inf[i]->Get(Form("rechitanalyzer/ee"));
        teb[i] = (TTree*)inf[i]->Get(Form("rechitanalyzer/eb"));
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
        tevt[i]->AddFriend(tjet[i]);
        tevt[i]->AddFriend(tpho[i]);
    }

    compareTwo_diff(tee[0],tee[1],"chi2",50,0,70,"",Form("%s",cap));
    compareTwo_diff(teb[0],teb[1],"chi2",50,0,70,"",Form("%s",cap));
    compareTwo_diff(tee[0],tee[1],"eError",50,0,0.2,"",Form("%s",cap));
    compareTwo_diff(teb[0],teb[1],"eError",50,0,0.2,"",Form("%s",cap));
   /* 
    const char* jetptcut = "Max$(jtpt)>80";
    compareTwo_diff(tee[0],tee[1],"chi2",50,0,1,Form("%s",cap));
    compareTwo_diff(teb[0],teb[1],"chi2",50,0,70,Form("%s",cap));
    compareTwo_diff(tee[0],tee[1],"eError",50,0,1,Form("%s",cap));
    compareTwo_diff(teb[0],teb[1],"eError",50,0,0.6,Form("%s",cap));
 */

#if 0
    double eMax=1.0;
    double etMax=0.5;
    const char* jetptcut = "Max$(jtpt)>90";
    compareTwo_diff(tpho[0],tpho[1],"phoEt",40,0,200,Form("%s",jetptcut),Form("%s_ptcut90",cap));
    compareTwo_diff(tpho[0],tpho[1],"phoEta",50,-3,3,Form("%s",jetptcut),Form("%s_ptcut90",cap));
    compareTwo_diff(tpho[0],tpho[1],"phoPhi",20,-3.14,3.14,Form("%s",jetptcut),Form("%s_ptcut90",cap));
    compareTwo_diff(tpho[0],tjet[1],"jtpt",40,0,200,Form("%s",jetptcut),Form("%s_ptcut90",cap));
    compareTwo_diff(tjet[0],tjet[1],"jteta",50,-3,3,Form("%s",jetptcut),Form("%s_ptcut90",cap));
    compareTwo_diff(tjet[0],tjet[1],"jtphi",20,-3.14,3.14,Form("%s",jetptcut),Form("%s_ptcut90",cap));
    compareTwo_diff(tevt[0],tevt[1],"hiNtracks",50,0,2000,Form("%s",jetptcut),Form("%s_ptcut90",cap));
    compareTwo_diff(tevt[0],tevt[1],"hiBin",50,0,200,Form("%s",jetptcut),Form("%s_ptcut90",cap));
    compareTwo_diff(tevt[0],tevt[1],"hiHF",40,0,4000,Form("%s",jetptcut),Form("%s_ptcut90",cap));
    compareTwo_diff(tevt[0],tevt[1],"hiHFhit",40,0,50000,Form("%s",jetptcut),Form("%s_ptcut90",cap));
// no cut 
    compareTwo_diff(tee[0],tee[1],"e",50,0,2,Form("%s",jetptcut),Form("%s_ptcut90",cap));
    compareTwo_diff(tee[0],tee[1],"et",50,0,1,Form("%s",jetptcut),Form("%s_ptcut90",cap));
    compareTwo_diff(tee[0],tee[1],"phi",20,-3.14,3.14,Form("%s",jetptcut),Form("%s_ptcut90",cap));
    compareTwo_diff_eta(tee[0],tee[1],teb[0],teb[1],"eta",50,-3,3,Form("%s",jetptcut),Form("%s_ptcut90",cap));
    compareTwo_diff(teb[0],teb[1],"e",50,0,0.6,Form("%s",jetptcut),Form("%s_ptcut90",cap));
    compareTwo_diff(teb[0],teb[1],"et",50,0,0.6,Form("%s",jetptcut),Form("%s_ptcut90",cap));
    compareTwo_diff(teb[0],teb[1],"phi",20,-3.14,3.14,Form("%s",jetptcut),Form("%s_ptcut90",cap));
// hiNtracks < 50 
    compareTwo_diff(tee[0],tee[1],"e",50,0,1.5,Form("hiNtracks<50 && %s",jetptcut),Form("%s_hiNtracks50_ptcut90",cap));
    compareTwo_diff(tee[0],tee[1],"et",50,0,1.5,Form("hiNtracks<50 && %s",jetptcut),Form("%s_hiNtracks50_ptcut90",cap));
    compareTwo_diff(tee[0],tee[1],"phi",20,-3.14,3.14,Form("hiNtracks<50 && %s",jetptcut),Form("%s_hiNtracks50_ptcut90",cap));
    compareTwo_diff_eta(tee[0],tee[1],teb[0],teb[1],"eta",50,-3,3,Form("hiNtracks<50 && %s",jetptcut),Form("%s_hiNtracks50_ptcut90",cap));
    compareTwo_diff(teb[0],teb[1],"e",50,0,0.5,Form("hiNtracks<50 && %s",jetptcut),Form("%s_hiNtracks50_ptcut90",cap));
    compareTwo_diff(teb[0],teb[1],"et",50,0,0.5,Form("hiNtracks<50 && %s",jetptcut),Form("%s_hiNtracks50_ptcut90",cap));
    compareTwo_diff(teb[0],teb[1],"phi",20,-3.14,3.14,Form("hiNtracks<50 && %s",jetptcut),Form("%s_hiNtracks50_ptcut90",cap));
// hiNtracks >1000

    compareTwo_diff(tee[0],tee[1],"e",50,0,2,Form("hiNtracks>1000 && %s",jetptcut),Form("%s_hiNtracks1000_ptcut90",cap));
    compareTwo_diff(tee[0],tee[1],"et",50,0,2,Form("hiNtracks>1000 && %s",jetptcut),Form("%s_hiNtracks1000_ptcut90",cap));
    compareTwo_diff(tee[0],tee[1],"phi",20,-3.14,3.14,Form("hiNtracks>1000 && %s",jetptcut),Form("%s_hiNtracks1000_ptcut90",cap));
    compareTwo_diff_eta(tee[0],tee[1],teb[0],teb[1],"eta",50,-3,3,Form("hiNtracks>1000 && %s",jetptcut),Form("%s_hiNtracks1000_ptcut90",cap));
    compareTwo_diff(teb[0],teb[1],"e",50,0,0.6,Form("hiNtracks>1000 && %s",jetptcut),Form("%s_hiNtracks1000_ptcut90",cap));
    compareTwo_diff(teb[0],teb[1],"et",50,0,0.6,Form("hiNtracks>1000 && %s",jetptcut),Form("%s_hiNtracks1000_ptcut90",cap));
    compareTwo_diff(teb[0],teb[1],"phi",20,-3.14,3.14,Form("hiNtracks>1000 && %s",jetptcut),Form("%s_hiNtracks1000_ptcut90",cap));
#endif
#if 0
    const char* jetptcut = "Max$(jtpt)>80";
    compareTwo_diff(tpho[0],tpho[1],"phoEt",40,0,200,Form("%s",jetptcut),Form("%s_ptcut80",cap));
    compareTwo_diff(tpho[0],tpho[1],"phoEta",50,-3,3,Form("%s",jetptcut),Form("%s_ptcut80",cap));
    compareTwo_diff(tpho[0],tpho[1],"phoPhi",20,-3.14,3.14,Form("%s",jetptcut),Form("%s_ptcut80",cap));
    compareTwo_diff(tjet[0],tjet[1],"jtpt",40,0,200,Form("%s",jetptcut),Form("%s_ptcut80",cap));
    compareTwo_diff(tjet[0],tjet[1],"jteta",50,-3,3,Form("%s",jetptcut),Form("%s_ptcut80",cap));
    compareTwo_diff(tjet[0],tjet[1],"jtphi",20,-3.14,3.14,Form("%s",jetptcut),Form("%s_ptcut80",cap));
    compareTwo_diff(tevt[0],tevt[1],"hiNtracks",50,0,2000,Form("%s",jetptcut),Form("%s_ptcut80",cap));
    compareTwo_diff(tevt[0],tevt[1],"hiBin",50,0,200,Form("%s",jetptcut),Form("%s_ptcut80",cap));
    compareTwo_diff(tevt[0],tevt[1],"hiHF",40,0,4000,Form("%s",jetptcut),Form("%s_ptcut80",cap));
    compareTwo_diff(tevt[0],tevt[1],"hiHFhit",40,0,50000,Form("%s",jetptcut),Form("%s_ptcut80",cap));
// no cut 
    compareTwo_diff(tee[0],tee[1],"e",50,0,2,Form("%s",jetptcut),Form("%s_ptcut80",cap));
    compareTwo_diff(tee[0],tee[1],"et",50,0,1,Form("%s",jetptcut),Form("%s_ptcut80",cap));
    compareTwo_diff(tee[0],tee[1],"phi",20,-3.14,3.14,Form("%s",jetptcut),Form("%s_ptcut80",cap));
    compareTwo_diff_eta(tee[0],tee[1],teb[0],teb[1],"eta",50,-3,3,Form("%s",jetptcut),Form("%s_ptcut80",cap));
    compareTwo_diff(teb[0],teb[1],"e",50,0,0.6,Form("%s",jetptcut),Form("%s_ptcut80",cap));
    compareTwo_diff(teb[0],teb[1],"et",50,0,0.6,Form("%s",jetptcut),Form("%s_ptcut80",cap));
    compareTwo_diff(teb[0],teb[1],"phi",20,-3.14,3.14,Form("%s",jetptcut),Form("%s_ptcut80",cap));
// hiNtracks < 50 
    compareTwo_diff(tee[0],tee[1],"e",50,0,1.5,Form("hiNtracks<50 && %s",jetptcut),Form("%s_hiNtracks50_ptcut80",cap));
    compareTwo_diff(tee[0],tee[1],"et",50,0,1.5,Form("hiNtracks<50 && %s",jetptcut),Form("%s_hiNtracks50_ptcut80",cap));
    compareTwo_diff(tee[0],tee[1],"phi",20,-3.14,3.14,Form("hiNtracks<50 && %s",jetptcut),Form("%s_hiNtracks50_ptcut80",cap));
    compareTwo_diff_eta(tee[0],tee[1],teb[0],teb[1],"eta",50,-3,3,Form("hiNtracks<50 && %s",jetptcut),Form("%s_hiNtracks50_ptcut80",cap));
    compareTwo_diff(teb[0],teb[1],"e",50,0,0.5,Form("hiNtracks<50 && %s",jetptcut),Form("%s_hiNtracks50_ptcut80",cap));
    compareTwo_diff(teb[0],teb[1],"et",50,0,0.5,Form("hiNtracks<50 && %s",jetptcut),Form("%s_hiNtracks50_ptcut80",cap));
    compareTwo_diff(teb[0],teb[1],"phi",20,-3.14,3.14,Form("hiNtracks<50 && %s",jetptcut),Form("%s_hiNtracks50_ptcut80",cap));
// hiNtracks >1000

    compareTwo_diff(tee[0],tee[1],"e",50,0,2,Form("hiNtracks>1000 && %s",jetptcut),Form("%s_hiNtracks1000_ptcut80",cap));
    compareTwo_diff(tee[0],tee[1],"et",50,0,2,Form("hiNtracks>1000 && %s",jetptcut),Form("%s_hiNtracks1000_ptcut80",cap));
    compareTwo_diff(tee[0],tee[1],"phi",20,-3.14,3.14,Form("hiNtracks>1000 && %s",jetptcut),Form("%s_hiNtracks1000_ptcut80",cap));
    compareTwo_diff_eta(tee[0],tee[1],teb[0],teb[1],"eta",50,-3,3,Form("hiNtracks>1000 && %s",jetptcut),Form("%s_hiNtracks1000_ptcut80",cap));
    compareTwo_diff(teb[0],teb[1],"e",50,0,0.6,Form("hiNtracks>1000 && %s",jetptcut),Form("%s_hiNtracks1000_ptcut80",cap));
    compareTwo_diff(teb[0],teb[1],"et",50,0,0.6,Form("hiNtracks>1000 && %s",jetptcut),Form("%s_hiNtracks1000_ptcut80",cap));
    compareTwo_diff(teb[0],teb[1],"phi",20,-3.14,3.14,Form("hiNtracks>1000 && %s",jetptcut),Form("%s_hiNtracks1000_ptcut80",cap));
#endif

#if 0
// basic spectra (jet, centrality, photon)
    compareTwo_diff(tpho[0],tpho[1],"phoEt",40,0,200,Form("%s",cap));
    compareTwo_diff(tpho[0],tpho[1],"phoEta",50,-3,3,Form("%s",cap));
    compareTwo_diff(tpho[0],tpho[1],"phoPhi",20,-3.14,3.14,Form("%s",cap));
    compareTwo_diff(tjet[0],tjet[1],"jtpt",40,0,200,"",Form("%s",cap));
    compareTwo_diff(tjet[0],tjet[1],"jteta",50,-3,3,"",Form("%s",cap));
    compareTwo_diff(tjet[0],tjet[1],"jtphi",20,-3.14,3.14,"",Form("%s",cap));
    compareTwo_diff(tevt[0],tevt[1],"hiNtracks",50,0,2000,"",Form("%s",cap));
    compareTwo_diff(tevt[0],tevt[1],"hiBin",50,0,200,"",Form("%s",cap));
    compareTwo_diff(tevt[0],tevt[1],"hiHF",40,0,4000,"",Form("%s",cap));
// no cut 
    compareTwo_diff(tee[0],tee[1],"e",50,0,2,"",Form("%s",cap));
    compareTwo_diff(tee[0],tee[1],"et",50,0,1,"",Form("%s",cap));
    compareTwo_diff(tee[0],tee[1],"phi",20,-3.14,3.14,"",Form("%s",cap));
    compareTwo_diff_eta(tee[0],tee[1],teb[0],teb[1],"eta",50,-3,3,"",Form("%s",cap));
    compareTwo_diff(teb[0],teb[1],"e",50,0,0.6,"",Form("%s",cap));
    compareTwo_diff(teb[0],teb[1],"et",50,0,0.6,"",Form("%s",cap));
    compareTwo_diff(teb[0],teb[1],"phi",20,-3.14,3.14,"",Form("%s",cap));
// hiNtracks < 50 
    compareTwo_diff(tee[0],tee[1],"e",50,0,1.5,"hiNtracks<50",Form("%s_hiNtracks50",cap));
    compareTwo_diff(tee[0],tee[1],"et",50,0,1.5,"hiNtracks<50",Form("%s_hiNtracks50",cap));
    compareTwo_diff(tee[0],tee[1],"phi",20,-3.14,3.14,"hiNtracks<50",Form("%s_hiNtracks50",cap));
    compareTwo_diff_eta(tee[0],tee[1],teb[0],teb[1],"eta",50,-3,3,"hiNtracks<50",Form("%s_hiNtracks50",cap));
    compareTwo_diff(teb[0],teb[1],"e",50,0,0.5,"hiNtracks<50",Form("%s_hiNtracks50",cap));
    compareTwo_diff(teb[0],teb[1],"et",50,0,0.5,"hiNtracks<50",Form("%s_hiNtracks50",cap));
    compareTwo_diff(teb[0],teb[1],"phi",20,-3.14,3.14,"hiNtracks<50",Form("%s_hiNtracks50",cap));
// hiNtracks >1000
    compareTwo_diff(tee[0],tee[1],"e",50,0,2,"hiNtracks>1000",Form("%s_hiNtracks1000",cap));
    compareTwo_diff(tee[0],tee[1],"et",50,0,2,"hiNtracks>1000",Form("%s_hiNtracks1000",cap));
    compareTwo_diff(tee[0],tee[1],"phi",20,-3.14,3.14,"hiNtracks>1000",Form("%s_hiNtracks1000",cap));
    compareTwo_diff_eta(tee[0],tee[1],teb[0],teb[1],"eta",50,-3,3,"hiNtracks>1000",Form("%s_hiNtracks1000",cap));
    compareTwo_diff(teb[0],teb[1],"e",50,0,0.6,"hiNtracks>1000",Form("%s_hiNtracks1000",cap));
    compareTwo_diff(teb[0],teb[1],"et",50,0,0.6,"hiNtracks>1000",Form("%s_hiNtracks1000",cap));
    compareTwo_diff(teb[0],teb[1],"phi",20,-3.14,3.14,"hiNtracks>1000",Form("%s_hiNtracks1000",cap));
#endif

}

void compareTwo_diff_eta(TTree* t1_ee, TTree* t2_ee, TTree* t1_eb, TTree* t2_eb, TString var, int nBins, double xMin, double xMax, const char* theCut, const char* cap)  {
    TCanvas* cc=  new TCanvas(Form("c_%s_%s_%s",var.Data(),t1_ee->GetName(),cap),"", 400,800);
    cc->Divide(1,2);
    cc->cd(1);
    TH1D* h1 = new TH1D(Form("h1_%s_%s_%s",var.Data(),t1_ee->GetName(),cap), Form(";%s;",var.Data()), nBins,xMin,xMax);
    TH1D* h2 = (TH1D*)h1->Clone(Form("h2_%s_%s_%s",var.Data(),t1_ee->GetName(),cap));
    h1->Sumw2();
    h2->Sumw2();
    t1_ee->Draw(Form("%s>>%s",var.Data(),h1->GetName()), theCut);
    t2_ee->Draw(Form("%s>>%s",var.Data(),h2->GetName()), theCut);
    t1_eb->Draw(Form("%s>>+%s",var.Data(),h1->GetName()), theCut);
    t2_eb->Draw(Form("%s>>+%s",var.Data(),h2->GetName()), theCut);

    Int_t t1eeEntries = t1_ee->GetEntries(Form("%s",theCut));
    Int_t t2eeEntries = t2_ee->GetEntries(Form("%s",theCut));
    Int_t t1ebEntries = t1_eb->GetEntries(Form("%s",theCut));
    Int_t t2ebEntries = t2_eb->GetEntries(Form("%s",theCut));
    h1->Scale( 1. / (t1eeEntries+t1ebEntries));
    h2->Scale( 1. / (t2eeEntries+t2ebEntries));
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


void compareTwo_ratio(TTree* t1, TTree* t2, TString var, int nBins, double xMin, double xMax, const char* theCut, const char* cap)  {
    TCanvas* cc=  new TCanvas(Form("c_%s_%s",var.Data(),t1->GetName()),"", 400,800);
    cc->Divide(1,2);
    cc->cd(1);
    TH1D* h1 = new TH1D(Form("h1_%s_%s",var.Data(),t1->GetName()), Form("%s;%s;",cap,var.Data()), nBins,xMin,xMax);
    TH1D* h2 = (TH1D*)h1->Clone(Form("h2_%s_%s",var.Data(),t1->GetName()));
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

    if (min2 < min1)
        h1->SetMinimum(min2-TMath::Abs(min2)*0.2);

    cleverRange(h1,h2);
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
    cc-> SaveAs(Form("pdf/%s_%s.pdf",cc->GetName(),cap)); 
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


int main(){
    compareRechit();
    return 0;
}
