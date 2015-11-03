/*
 * test_GammaJetAnalyzer.C
 *
 * macro to test GammaJetAnalyzer class
 *  1. histograms by GammayJetAnalyzer
 */

#include <TFile.h>
#include <TTree.h>
#include <TCut.h>

#include <iostream>

const int MAXPHOTONS = 500;
const int MAXJETS = 500;

void ySkimAnalyzerMacros()
{
    const char* inputfileName_20131113_photon35 = "/export/d00/scratch/luck/yskimmedFiles/yskim_HiForestPhoton-v7-noDuplicate_20131113_photon35.root";
    const char* inputfileName_PAS               = "/export/d00/scratch/tatar/yskimmedFiles2013/yskim_HiForestPhoton-v7-noDuplicate.root";

    TFile *file_20131113_photon35 = new TFile(inputfileName_20131113_photon35, "READ");
    TFile *file_PAS = new TFile(inputfileName_PAS, "READ");

    TTree* tgj_20131113 = (TTree*)file_20131113_photon35->Get("tgj");
    TTree* tgj_PAS      = (TTree*)file_PAS->Get("tgj");

    TTree* yJet_20131113 = (TTree*)file_20131113_photon35->Get("yJet");
    TTree* yJet_PAS      = (TTree*)file_PAS->Get("yJet");

    TTree* mJet_20131113 = (TTree*)file_20131113_photon35->Get("mJet");
    TTree* mJet_PAS      = (TTree*)file_PAS->Get("mJet");

    mJet_20131113->AddFriend(tgj_20131113);
    mJet_PAS->AddFriend(tgj_PAS);

    int centBin_index = 0;
    int pTBin_index   = 0;

    int cent_gt[2] = {0,  12};
    int cent_lt[2] = {11, 38};

    int photonEt_gt[4] = {40, 50, 60, 80};
    int photonEt_lt[4] = {50, 60, 80, 9999};

    TCut centCut  = Form("cBin >= %d && cBin<= %d",        cent_gt[centBin_index],   cent_lt[centBin_index]);
    TCut ptPhoCut = Form("photonEt > %d && photonEt < %d", photonEt_gt[pTBin_index], photonEt_lt[pTBin_index]);

//  https://github.com/CmsHI/gammaJetAnalysis/blob/527d9a4c5366679ce8afdf8404fe0ba22cfc6790/histogramProducer/gammaJetHistProducer.C#L95
    TCut caloIso;
    int colType = 2;    // 1=pp, 2=HI, 3=pA
    if ( colType == 1 )
      caloIso = "(ecalIso < 4.2  &&  hcalIso < 2.2  &&  trackIso < 2) && hovere<0.1";
    else if ( colType == 2 )
      caloIso = "(sumIso<1) && hovere<0.1";
    else {
      caloIso = "ecalIso < 4.2  &&  hcalIso < 2.2  &&  trackIso < 2 && hovere<0.1";
    }
    TCut sbIso = "(sumIso>10) && (sumIso<20) && hovere<0.1";

    TCut basicPhoCut  = centCut && ptPhoCut && caloIso;
    TCut sbPhoCut     = centCut && ptPhoCut && sbIso;
    TCut evtSeltCut   = basicPhoCut;
    TCut sbSeltCut    = sbPhoCut ;
    TCut phoCandCut   = "sigmaIetaIeta<0.010";
    TCut phoDecayCut  = "(sigmaIetaIeta>0.011) && (sigmaIetaIeta<0.017)";
    TCut jetCut       =  Form("abs(eta)<%f && pt>%f", (float)1.6, (float)30 );

    std::cout<<"Scan started"<<std::endl;

//    tgj_20131113->SetScanField(0);
//    tgj_20131113->Scan("run:evt:lumi:hiBin", (evtSeltCut && phoCandCut).GetTitle());

    std::cout<<"Scan finished"<<std::endl;


    std::cout<<"Draw started"<<std::endl;

    TFile* out=new TFile("ySkimAnalyzerMacros.root","recreate");

    mJet_20131113->Draw("nJet >> h_nJet_20131113", (evtSeltCut && phoCandCut && jetCut).GetTitle());
    mJet_PAS     ->Draw("nJet >> h_nJet_PAS",      (evtSeltCut && phoCandCut && jetCut).GetTitle());

    mJet_20131113->Draw("nJet >> h_nJet_20131113_decay", (evtSeltCut && phoDecayCut && jetCut).GetTitle());
    mJet_PAS     ->Draw("nJet >> h_nJet_PAS_decay",      (evtSeltCut && phoDecayCut && jetCut).GetTitle());

    mJet_20131113->Draw("dphi >> h_dphi_20131113", (evtSeltCut && phoCandCut && jetCut).GetTitle());
    mJet_PAS     ->Draw("dphi >> h_dphi_PAS",      (evtSeltCut && phoCandCut && jetCut).GetTitle());

    mJet_20131113->Draw("dphi >> h_dphi_20131113_decay", (evtSeltCut && phoDecayCut && jetCut).GetTitle());
    mJet_PAS     ->Draw("dphi >> h_dphi_PAS_decay",      (evtSeltCut && phoDecayCut && jetCut).GetTitle());

    std::cout<<"Draw finished"<<std::endl;

    out->Write();
    out->Close();
}

int main(int argc, char** argv)
{
    ySkimAnalyzerMacros();
}
