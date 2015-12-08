#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TGraphAsymmErrors.h>
#include <TLegend.h>
#include <TChain.h>
#include <iostream>
#include "../../yjUtility.h"
void plotCentRatio()
{
    SetHistTitleStyle();
    SetyjPadStyle();
    gStyle->SetOptStat(0);
    TH1::SetDefaultSumw2();
    
    int Nrun = 2;
    int run[] = {262548,262695, 262694, 262703, 262735};
    //int run[] = {262620,262656,262694,262695,262703,262735,262811,262816,263035};
    //int run[] = {262695,262548,262620,262656,262694,262811,262816};
    //string var = "hiBin";
    //int nomL(0), nomH(5), denL(0), denH(60);
    string var = "hiHF";
    int nomL(3840), nomH(6000), denL(240), denH(6000);


    const char* infname[Nrun];

    int nBin(50);
    double binL(0), binH(100);
    const char* cut1 = Form("%s>=%d && %s<=%d",var.data(),nomL,var.data(),nomH);
    const char* cut2 = Form("%s>=%d && %s<=%d",var.data(),denL,var.data(),denH);
    const char* presel = "pcollisionEventSelection && HLT_HIL1MinimumBiasHF1AND_v1";

    TGraphAsymmErrors *g = new TGraphAsymmErrors;
    for(int i=0;i<Nrun;i++){
        cout << "processing run " << run[i] << " ..." << endl;
        if(run[i]==262620 || run[i]==262656) infname[i]=Form("root://eoscms//eos/cms//store/group/phys_heavyions/velicanu/forest/HIRun2015/HIMinimumBias2/Merged/HIForestExpress_run%d.root",run[i]);
        else if(run[i]==263233) infname[i]=(Form("root://eoscms//eos/cms/store/group/cmst3/group/hintt/mverweij/PbPbReco/Forest/000/263/233/merge/HiForestRun263233PrivateRecoMB.root"));
        else if(run[i]==262694 || run[i]==262695 || run[i]==262703 || run[i]==262735) infname[i]=Form("root://eoscms//eos/cms//store/group/phys_heavyions/velicanu/forest/HIRun2015/HIMinimumBias2/Merged/HiForestPromptReco_%d.root",run[i]);
        else if(run[i]==262548) infname[i]=Form("root://eoscms//eos/cms//store/group/phys_heavyions/velicanu/forest/HIRun2015/HIExpressPhysics/Merged/HIForestMinbiasUPC_run262548.root");
        else infname[i]=Form("root://eoscms//eos/cms//store/group/phys_heavyions/velicanu/forest/HIRun2015/HIMinimumBias2/Merged/HIMinimumBias2_run%d.root",run[i]);
        TFile *inf = TFile::Open(infname[i]);
        TTree *tEvt = (TTree*)inf->Get("hiEvtAnalyzer/HiTree");
        TTree *tHLT = (TTree*)inf->Get("hltanalysis/HltTree");
        TTree *tSkim = (TTree*)inf->Get("skimanalysis/HltTree");
        tEvt->AddFriend(tHLT);
        tEvt->AddFriend(tSkim);
        int nom = tEvt->GetEntries(Form("(%s)&&(%s)",presel,cut1));
        int den = tEvt->GetEntries(Form("(%s)&&(%s)",presel,cut2));
        double ratio = (double)nom/den;
        //double ratio = getRatio(tEvt,var,cut1,cut2,presel,nBin,binL,binH);
        double err = (1./den)*(sqrt(nom*(1.-nom/(double)den)));
        //double err = (1./den)*(sqrt(nom*(1.-nom/(double)den)));
        cout << "ratio : " << ratio << ", err : " << err << endl;
        g->SetPoint(i,run[i],ratio);
        g->SetPointError(i,0.0001,0.0001,err,err);
    }
    double runL(262545), runH(263040);
    int runBin=runH-runL;
    TH1D *h= new TH1D("h","",runBin,runL,runH);
    h->SetXTitle("Run");
    h->SetYTitle(Form("%s #(%d-%d)/#(%d-%d)",var.data(),nomL,nomH,denL,denH));
    h->GetYaxis()->SetRangeUser(0,0.12);
    TCanvas *c = new TCanvas("c","",500,500);
    h->Draw();
    gPad->SetGridy();
    g->Draw("p same");
    g->SetMarkerStyle(20);
    g->SetMarkerColor(2);
    /*
       gHF2And->Draw("p same");
       gHF2And->SetLineColor(2);
       gHF2And->SetMarkerColor(2);
       gHF2And->SetMarkerStyle(24);
    //gHFHFplusANDminusTH0->Draw("p same");
    gHFHFplusANDminusTH0->SetLineColor(4);
    gHFHFplusANDminusTH0->SetMarkerColor(4);
    gHFHFplusANDminusTH0->SetMarkerStyle(24);

    TLegend *leg = new TLegend (0.44,0.26,0.9,0.53);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(g,"HIL1MinimumBiasHF1ANDExpress_v1","pl");  
    leg->AddEntry(gHF2And,"HIL1MinimumBiasHF2ANDExpress_v1","pl");
    //leg->AddEntry(gHFHFplusANDminusTH0,"L1_HFplusANDminusTH0","pl");
    leg->Draw();
    */
}
