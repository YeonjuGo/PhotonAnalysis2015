/*
 * test_GammaJetAnalyzer.C
 *
 * macro to test GammaJetAnalyzer class
 *  1. histograms by GammayJetAnalyzer
 */

#include "GammaJetAnalyzer.h"
#include "GammaJetAnalyzer.cc"   // need to use this include if this macro and "GammaJetAnalyzer.h" are not in the same directory.
#include "histoUtil.h"
#include "smallPhotonUtil.h"

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TCanvas.h>

#include <iostream>
#include <ctime>
#include <iomanip>

const int MAXPHOTONS = 500;
const int MAXJETS = 500;

void GammaJetAnalyzerMacros_yskim()
{
    const char* inputfileName="/mnt/hadoop/cms/store/user/luck/PbPb2011_photons_Data/HiForestPhoton-v7-noDuplicate.root";
    collisionType collision = HI;

    TFile *inputFile = new TFile(inputfileName, "READ");

    GammaJetAnalyzer* gja = new GammaJetAnalyzer(inputFile, collision);

    TString cond_eta_spike     = GammaJetAnalyzer::mergeSelections(gja->cond_pt_eta,  gja->cond_spike);
    TString cond_eta_spike_iso = GammaJetAnalyzer::mergeSelections(cond_eta_spike, gja->cond_iso);

//  forest2yskim_jetSkim_forestV3.C
    TCut forest2yskimCut_event = "pcollisionEventSelection !=0 && abs(vz)<15";
    TCut forest2yskimCut_photon = "pt > 35 && abs(eta)<1.44";
    TCut forest2yskimCut_photon_spike = gja->cond_spike.Data();
    forest2yskimCut_photon = forest2yskimCut_photon && forest2yskimCut_photon_spike && TCut("hadronicOverEm<=0.2 && isEle==0 && seedTime!=0");
    forest2yskimCut_photon = forest2yskimCut_photon && TCut("rawEnergy/energy >= 0.5");
    TCut forest2yskimCut = forest2yskimCut_event && forest2yskimCut_photon;

//    gammaJetHistProducer.c
//    CUTS :
//    centCut      =  "cBin >= 0 && cBin<= 11"
//    ptPhoCut     =  "photonEt>40 && photonEt<50"
//    caloIso      =  "(sumIso<1) && hovere<0.1"
//    sbIso        =  "(sumIso>10) && (sumIso<20) && hovere<0.1"
//    basicPhoCut  =  centCut && ptPhoCut && caloIso
//    sbPhoCut     =  centCut && ptPhoCut && sbIso
    TCut centCut      =  "hiBin >= 12 && hiBin<= 38";
    TCut ptPhoCut     =  "pt>60 && pt<80";
    TCut caloIso      =  "((cr4+cc4+ct4PtCut20)/0.9<1) && hadronicOverEm<0.1";
    TCut sbIso        =  "((cr4+cc4+ct4PtCut20)/0.9>10) && ((cr4+cc4+ct4PtCut20)/0.9<20) && hadronicOverEm<0.1";
    TCut basicPhoCut  =  centCut && ptPhoCut && caloIso;
    TCut sbPhoCut     =  centCut && ptPhoCut && sbIso   ;

//    evtSeltCut   =  basicPhoCut
//    sbSeltCut    =  sbPhoCut
//    phoCandCut   =  "sigmaIetaIeta<0.010";
//    phoDecayCut  =  "(sigmaIetaIeta>0.011) && (sigmaIetaIeta<0.017)";
    TCut evtSeltCut = basicPhoCut;
    TCut sbSeltCut  = sbPhoCut ;
    TString phoCandCut   =  "sigmaIetaIeta<0.010";
    TString phoDecayCut  =  "(sigmaIetaIeta>0.011) && (sigmaIetaIeta<0.017)";

    std::cout<<"GetEntries() started"<<std::endl;

    // evtSeltCut
    Long64_t entries_evtSeltCut = gja->tree->GetEntries((forest2yskimCut && evtSeltCut).GetTitle());

    // sbPhoCut
    Long64_t entries_sbPhoCut = gja->tree->GetEntries((forest2yskimCut && sbPhoCut).GetTitle());

    // evtSeltCut && phoCandCut  ==> # of candidate photons
    Long64_t entries_phoCandCut = gja->tree->GetEntries((forest2yskimCut && evtSeltCut && phoCandCut).GetTitle());

    // evtSeltCut && phoDecayCut ==> # of sideband photons
    Long64_t entries_phoDecayCut = gja->tree->GetEntries((forest2yskimCut && evtSeltCut && phoDecayCut).GetTitle());

    // jetCut : evtSeltCut && phoCandCut && "abs(eta)<1.6 && pt>30" ==> # of raw jets
    // jetCut : ((cBin >= 12 && cBin<= 38 && photonEt>60 && photonEt<80 && (sumIso<1)  && hovere<0.1)&&(sigmaIetaIeta<0.010))&&(abs(eta)<1.6 && pt>30)
    Long64_t entries_jetCut = gja->tree->GetEntries((forest2yskimCut && evtSeltCut && phoCandCut
                                                                     && TCut("abs(jteta)<1.6 && jtpt>30")).GetTitle());



    // jetCutDphi : evtSeltCut && phoCandCut && "abs(eta)<1.6 && pt>30" && Form("abs(dphi)>%f", PI * 7./8.) ==> # of signal jets
    // jetCutDphi : (((cBin >= 12 && cBin<= 38 && photonEt>60 && photonEt<80 && (sumIso<1)  && hovere<0.1)&&(sigmaIetaIeta<0.010))&&(abs(eta)<1.6 && pt>30))&&(abs(dphi)>2.748894)
//    Long64_t entries_phoCandCut = gja->tree->GetEntries((forest2yskimCut && evtSeltCut && phoCandCut
//                                                                         && TCut("abs(jteta)<1.6 && jtpt>30")
//                                                                         && Form("abs(dphi)>%f", PI * 7./8.)).GetTitle());

    std::cout<<"GetEntries() finished"<<std::endl;

    std::cout<<"entries_evtSeltCut  = "<<entries_evtSeltCut<<std::endl;
    std::cout<<"entries_sbPhoCut    = "<<entries_sbPhoCut<<std::endl;
    std::cout<<"entries_phoCandCut  = "<<entries_phoCandCut<<std::endl;
    std::cout<<"entries_phoDecayCut = "<<entries_phoDecayCut<<std::endl;
    std::cout<<"entries_jetCut      = "<<entries_jetCut<<std::endl;

    std::cout<<"drawing started"<<std::endl;

    /////

    gja->tree->SetScanField(0);
    gja->tree->Scan("run:evt:lumi:hiBin",(forest2yskimCut && evtSeltCut && phoCandCut && TCut("abs(jteta)<1.6 && jtpt>30")).GetTitle());
    /////

//    gja->drawMax("pt", "pt", GammaJetAnalyzer::mergeSelections(gja->cond_photon, (ptPhoCut && TCut("rawEnergy/energy >= 0.5") && TCut("hadronicOverEm<=0.2 && isEle==0 && seedTime!=0")).GetTitle()), GammaJetAnalyzer::mergeSelections(gja->cond_event, (forest2yskimCut_event && centCut).GetTitle()));
    // gja->drawMaxJet("jtpt", "jtpt",
    //         GammaJetAnalyzer::mergeSelections("abs(jteta)<1.6 && jtpt>30",
    //                                            Form("%s >= 0.5",gja->constructFormula_dR(alias_maxphoton_eta ,alias_maxphoton_phi,"jteta","jtphi").Data())),
    //         GammaJetAnalyzer::mergeSelections(gja->cond_photon,
    //                                           (ptPhoCut && TCut("rawEnergy/energy >= 0.5") && TCut("hadronicOverEm<=0.2 && isEle==0 && seedTime!=0")).GetTitle()),
    //         GammaJetAnalyzer::mergeSelections(gja->cond_event,
    //                                           (forest2yskimCut_event && centCut).GetTitle()));

//    gja->drawMax("pt/jtpt", "pt", gja->cond_photon, gja->cond_event);
//    gja->tree->Draw(Form("jtpt/Max$(pt*(%s))", gja->cond_photon.Data()), Form("Sum$(%s)>0",GammaJetAnalyzer::mergeSelections(gja->cond_photon,  gja->cond_event).Data()));
//    TCanvas* c2=new TCanvas();
//    c2->cd();
//    gja->tree->Draw("nref", Form("Sum$(%s)>0",GammaJetAnalyzer::mergeSelections(gja->cond_photon,  gja->cond_event).Data()));
//
//    TCanvas* c3=new TCanvas();
//    c3->cd();
//    gja->tree->Draw(Form("Max$(jtpt)/Max$(pt*(%s))", gja->cond_photon.Data()), Form("Sum$(%s)>0",GammaJetAnalyzer::mergeSelections(gja->cond_photon,  gja->cond_event).Data()));
//
//    TCanvas* c4=new TCanvas();
//    c4->cd();
//    gja->tree->Draw(Form("Sum$(phi*(pt==Max$(pt*(%s)))):Sum$(jtphi*(jtpt==Max$(jtpt)))>>hcorr", gja->cond_photon.Data()), Form("Sum$(%s)>0",GammaJetAnalyzer::mergeSelections(gja->cond_photon,  gja->cond_event).Data()));
//
//    std::cout<<gja->tree->GetEntries(GammaJetAnalyzer::mergeSelections(gja->cond_photon,  gja->cond_event))<<std::endl;

    std::cout<<"drawing finished"<<std::endl;
}

int main(int argc, char** argv)
{
    GammaJetAnalyzerMacros_yskim();
}
