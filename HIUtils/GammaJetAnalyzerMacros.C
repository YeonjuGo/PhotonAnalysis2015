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
#include <TString.h>

#include <iostream>
#include <string>

#define PI 3.141592653589

//const int MAXPHOTONS = 500;
//const int MAXJETS = 500;

void GammaJetAnalyzerMacros();
void GammaJetAnalyzerMacrosv1();
void GammaJetAnalyzerMacrosv2();

void GammaJetAnalyzerMacrosv2()
{
    const char* inputfileName="/mnt/hadoop/cms/store/user/luck/PbPb2011_photons_Data/HiForestPhoton-v7-noDuplicate.root";
    collisionType collision = HI;

    TFile *inputFile  = new TFile(inputfileName, "READ");
    TFile *outputFile = new TFile("GammaJetAnalyzerMacrosv2.root", "RECREATE");

    GammaJetAnalyzer* gja = new GammaJetAnalyzer(inputFile);
    if(collision == PP || collision == PA)   {
        gja->setJetTree(ak3PFJets);
    }
    else {
        gja->setJetTree(akPu3PFJets);
    }

    TString cond_event = gja->cond_event;
    if(collision == HI)  {
        gja->cut_hiBin_gt = 11;
        gja->cut_hiBin_lt = 39;
        gja->updateSelections();
        cond_event = GammaJetAnalyzer::mergeSelections(gja->cond_event, gja->cond_hiBin);
    }

    TString cond_eta_spike     = GammaJetAnalyzer::mergeSelections(gja->cond_pt_eta,  gja->cond_spike);
    TString cond_eta_spike_iso = GammaJetAnalyzer::mergeSelections(cond_eta_spike, gja->cond_iso);

    TString cond_jet_pt_eta    = GammaJetAnalyzer::mergeSelections(gja->cond_jet_pt, gja->cond_jet_eta);

    std::string h_formula;
    std::string h_cond_photon;
    std::string h_cond_jet;
    std::string h_selection;
    std::string h_cond_event;
    std::string h_histName;
    TString dR, dphi, deta;
    std::cout<<"drawing started"<<std::endl;

    // xjg, dphi > 7/8 * PI
    gja->resetCuts();
    h_formula     = "MAXJETPT/MAXPHOTONPT";
    h_cond_photon = gja->cond_photon;
    h_cond_jet    = GammaJetAnalyzer::mergeSelections(cond_jet_pt_eta, gja->cond_jet_dphi);
    h_selection   = "1";
    h_cond_event  = cond_event;
    TH1D* h_Xjg_MaxGMaxJ_dphi78 = new TH1D("h_Xjg_MaxGMaxJ_dphi78",h_formula.c_str() ,16 ,0 , 2);
    h_histName    = h_Xjg_MaxGMaxJ_dphi78->GetName();
    std::cout << "--------------"   << std::endl;
    std::cout << "h_formula     = " << h_formula     <<std::endl;
    std::cout << "h_cond_photon = " << h_cond_photon <<std::endl;
    std::cout << "h_cond_jet    = " << h_cond_jet    <<std::endl;
    std::cout << "h_cond_event  = " << h_cond_event  <<std::endl;
    gja->drawPhotonJet(h_formula.c_str(),h_cond_photon.c_str(),h_cond_jet.c_str(),h_selection.c_str(),h_cond_event.c_str(),h_histName.c_str());

    // xjg, dphi > 7/8 * PI
    gja->resetCuts();
    h_formula     = "JETPT/MAXPHOTONPT";
    h_cond_photon = gja->cond_photon;
    h_cond_jet    = GammaJetAnalyzer::mergeSelections(cond_jet_pt_eta, gja->cond_jet_dphi);
    h_selection   = h_cond_jet.c_str();
    h_cond_event  = cond_event;
    TH1D* h_Xjg_MaxG_dphi78 = new TH1D("h_Xjg_MaxG_dphi78",h_formula.c_str() ,16 ,0 , 2);
    h_histName    = h_Xjg_MaxG_dphi78->GetName();
    std::cout << "--------------"   << std::endl;
    std::cout << "h_formula     = " << h_formula     <<std::endl;
    std::cout << "h_cond_photon = " << h_cond_photon <<std::endl;
    std::cout << "h_cond_jet    = " << h_cond_jet    <<std::endl;
    std::cout << "h_cond_event  = " << h_cond_event  <<std::endl;
    gja->drawPhotonJet(h_formula.c_str(),h_cond_photon.c_str(),h_cond_jet.c_str(),h_selection.c_str(),h_cond_event.c_str(),h_histName.c_str());

    // xjg, dphi > 5/8 * PI
    gja->resetCuts();
    gja->cut_jet_photon_deltaPhi = PI * 5./8.;
    gja->updateSelections();
    h_formula     = "MAXJETPT/MAXPHOTONPT";
    h_cond_photon = gja->cond_photon;
    h_cond_jet    = GammaJetAnalyzer::mergeSelections(cond_jet_pt_eta, gja->cond_jet_dphi);
    h_selection   = "1";
    h_cond_event  = cond_event;
    TH1D* h_Xjg_MaxGMaxJ_dphi58 = new TH1D("h_Xjg_MaxGMaxJ_dphi58",h_formula.c_str() ,16 ,0 , 2);
    h_histName    = h_Xjg_MaxGMaxJ_dphi58->GetName();
    std::cout << "--------------"   << std::endl;
    std::cout << "h_formula     = " << h_formula     <<std::endl;
    std::cout << "h_cond_photon = " << h_cond_photon <<std::endl;
    std::cout << "h_cond_jet    = " << h_cond_jet    <<std::endl;
    std::cout << "h_cond_event  = " << h_cond_event  <<std::endl;
    gja->drawPhotonJet(h_formula.c_str(),h_cond_photon.c_str(),h_cond_jet.c_str(),h_selection.c_str(),h_cond_event.c_str(),h_histName.c_str());

    // xjg, dphi > 5/8 * PI
    gja->resetCuts();
    gja->cut_jet_photon_deltaPhi = PI * 5./8.;
    gja->updateSelections();
    h_formula     = "JETPT/MAXPHOTONPT";
    h_cond_photon = gja->cond_photon;
    h_cond_jet    = GammaJetAnalyzer::mergeSelections(cond_jet_pt_eta, gja->cond_jet_dphi);
    h_selection   = h_cond_jet.c_str();
    h_cond_event  = cond_event;
    TH1D* h_Xjg_MaxG_dphi58 = new TH1D("h_Xjg_MaxG_dphi58",h_formula.c_str() ,16 ,0 , 2);
    h_histName    = h_Xjg_MaxG_dphi58->GetName();
    std::cout << "--------------"   << std::endl;
    std::cout << "h_formula     = " << h_formula     <<std::endl;
    std::cout << "h_cond_photon = " << h_cond_photon <<std::endl;
    std::cout << "h_cond_jet    = " << h_cond_jet    <<std::endl;
    std::cout << "h_cond_event  = " << h_cond_event  <<std::endl;
    gja->drawPhotonJet(h_formula.c_str(),h_cond_photon.c_str(),h_cond_jet.c_str(),h_selection.c_str(),h_cond_event.c_str(),h_histName.c_str());

    // draw deltaR (MAXPHOTON, MAXJET) : MAXPHOTONPT
    dR = gja->constructFormula_dR("MAXPHOTONETA", "MAXPHOTONPHI", "MAXJETETA", "MAXJETPHI");
    gja->resetCuts();
    gja->updateSelections();
    h_formula     = Form("%s:MAXPHOTONPT", dR.Data());
    h_cond_photon = gja->cond_photon;
    h_cond_jet    = cond_jet_pt_eta;
    h_selection   = "1";
    h_cond_event  = cond_event;
    TH2D* h2D_dR_MaxGMaxJ = new TH2D("h2D_dR_MaxGMaxJ",h_formula.c_str(),44 ,0 ,220 ,54 ,0 ,5.4);
    h_histName    = h2D_dR_MaxGMaxJ->GetName();
    std::cout << "--------------"   << std::endl;
    std::cout << "h_formula     = " << h_formula     <<std::endl;
    std::cout << "h_cond_photon = " << h_cond_photon <<std::endl;
    std::cout << "h_cond_jet    = " << h_cond_jet    <<std::endl;
    std::cout << "h_cond_event  = " << h_cond_event  <<std::endl;
    gja->drawPhotonJet(h_formula.c_str(),h_cond_photon.c_str(),h_cond_jet.c_str(),h_selection.c_str(),h_cond_event.c_str(),h_histName.c_str());

    // draw deltaR (MAXPHOTON, JET) : MAXPHOTONPT
    dR = gja->constructFormula_dR("MAXPHOTONETA", "MAXPHOTONPHI", "JETETA", "JETPHI");
    gja->resetCuts();
    gja->updateSelections();
    h_formula     = Form("%s:MAXPHOTONPT", dR.Data());
    h_cond_photon = gja->cond_photon;
    h_cond_jet    = cond_jet_pt_eta;
    h_selection   = cond_jet_pt_eta;
    h_cond_event  = cond_event;
    TH2D* h2D_dR_MaxG = new TH2D("h2D_dR_MaxG",h_formula.c_str(),44 ,0 ,220 ,54 ,0 ,5.4);
    h_histName    = h2D_dR_MaxG->GetName();
    std::cout << "--------------"   << std::endl;
    std::cout << "h_formula     = " << h_formula     <<std::endl;
    std::cout << "h_cond_photon = " << h_cond_photon <<std::endl;
    std::cout << "h_cond_jet    = " << h_cond_jet    <<std::endl;
    std::cout << "h_cond_event  = " << h_cond_event  <<std::endl;
    gja->drawPhotonJet(h_formula.c_str(),h_cond_photon.c_str(),h_cond_jet.c_str(),h_selection.c_str(),h_cond_event.c_str(),h_histName.c_str());

/*  "trackN" does not exist in /mnt/hadoop/cms/store/user/luck/PbPb2011_photons_Data/HiForestPhoton-v7-noDuplicate.root
    // draw trackN (MAXJET) : MAXPHOTONPT
    TString trackN_MAXJET = Form("Sum$(trackN*(jtpt == Max$(jtpt*(%s))))", cond_jet_pt_eta.Data());
    gja->resetCuts();
    gja->updateSelections();
    h_formula     = Form("%s:MAXPHOTONPT", trackN_MAXJET.Data());
    h_cond_photon = gja->cond_photon;
    h_cond_jet    = cond_jet_pt_eta;
    h_selection   = "1";
    h_cond_event  = cond_event;
    h_histName    = "h2D_trackN_MaxGMaxJ";
    std::cout << "--------------"   << std::endl;
    std::cout << "h_formula     = " << h_formula     <<std::endl;
    std::cout << "h_cond_photon = " << h_cond_photon <<std::endl;
    std::cout << "h_cond_jet    = " << h_cond_jet    <<std::endl;
    std::cout << "h_cond_event  = " << h_cond_event  <<std::endl;
    gja->drawPhotonJet(h_formula.c_str(),h_cond_photon.c_str(),h_cond_jet.c_str(),h_selection.c_str(),h_cond_event.c_str(),h_histName.c_str());
*/

/*  "trackN" does not exist in /mnt/hadoop/cms/store/user/luck/PbPb2011_photons_Data/HiForestPhoton-v7-noDuplicate.root
    // draw trackN (JET) : MAXPHOTONPT
    TString trackN_JET = Form("trackN*(%s)", cond_jet_pt_eta.Data());
    gja->resetCuts();
    gja->updateSelections();
    h_formula     = Form("%s:MAXPHOTONPT", trackN_JET.Data());
    h_cond_photon = gja->cond_photon;
    h_cond_jet    = cond_jet_pt_eta;
    h_selection   = cond_jet_pt_eta;
    h_cond_event  = cond_event;
    h_histName    = "h2D_trackN_MaxG";
    std::cout << "--------------"   << std::endl;
    std::cout << "h_formula     = " << h_formula     <<std::endl;
    std::cout << "h_cond_photon = " << h_cond_photon <<std::endl;
    std::cout << "h_cond_jet    = " << h_cond_jet    <<std::endl;
    std::cout << "h_cond_event  = " << h_cond_event  <<std::endl;
    gja->drawPhotonJet(h_formula.c_str(),h_cond_photon.c_str(),h_cond_jet.c_str(),h_selection.c_str(),h_cond_event.c_str(),h_histName.c_str());
*/

    // draw deltaPhi (MAXPHOTON, MAXJET) : MAXPHOTONPT
    dphi = gja->constructFormula_dphi("MAXPHOTONPHI","MAXJETPHI");
    gja->resetCuts();
    gja->updateSelections();
    h_formula     = Form("%s:MAXPHOTONPT", dphi.Data());
    h_cond_photon = gja->cond_photon;
    h_cond_jet    = cond_jet_pt_eta;
    h_selection   = "1";
    h_cond_event  = cond_event;
    TH2D* h2D_dphi_MaxGMaxJ = new TH2D("h2D_dphi_MaxGMaxJ",h_formula.c_str(),44 ,0 ,220 ,35 ,0 ,3.5);
    h_histName    = h2D_dphi_MaxGMaxJ->GetName();
    std::cout << "--------------"   << std::endl;
    std::cout << "h_formula     = " << h_formula     <<std::endl;
    std::cout << "h_cond_photon = " << h_cond_photon <<std::endl;
    std::cout << "h_cond_jet    = " << h_cond_jet    <<std::endl;
    std::cout << "h_cond_event  = " << h_cond_event  <<std::endl;
    gja->drawPhotonJet(h_formula.c_str(),h_cond_photon.c_str(),h_cond_jet.c_str(),h_selection.c_str(),h_cond_event.c_str(),h_histName.c_str());

    // draw deltaEta (MAXPHOTON, MAXJET) : MAXPHOTONPT
    deta = gja->constructFormula_deta("MAXPHOTONETA","MAXJETETA");
    gja->resetCuts();
    gja->updateSelections();
    h_formula     = Form("%s:MAXPHOTONPT", deta.Data());
    h_cond_photon = gja->cond_photon;
    h_cond_jet    = cond_jet_pt_eta;
    h_selection   = "1";
    h_cond_event  = cond_event;
    TH2D* h2D_deta_MaxGMaxJ = new TH2D("h2D_deta_MaxGMaxJ",h_formula.c_str(),44 ,0 ,220 ,35 ,0 ,3.5);
    h_histName    = h2D_deta_MaxGMaxJ->GetName();
    std::cout << "--------------"   << std::endl;
    std::cout << "h_formula     = " << h_formula     <<std::endl;
    std::cout << "h_cond_photon = " << h_cond_photon <<std::endl;
    std::cout << "h_cond_jet    = " << h_cond_jet    <<std::endl;
    std::cout << "h_cond_event  = " << h_cond_event  <<std::endl;
    gja->drawPhotonJet(h_formula.c_str(),h_cond_photon.c_str(),h_cond_jet.c_str(),h_selection.c_str(),h_cond_event.c_str(),h_histName.c_str());

    // draw deltaPhi (MAXPHOTON, MAXJET) : deltaEta (MAXPHOTON, MAXJET)
    deta = gja->constructFormula_deta("MAXPHOTONETA","MAXJETETA");
    dphi = gja->constructFormula_dphi("MAXPHOTONPHI","MAXJETPHI");
    gja->resetCuts();
    gja->updateSelections();
    h_formula     = Form("%s:%s", deta.Data(), dphi.Data());
    h_cond_photon = gja->cond_photon;
    h_cond_jet    = cond_jet_pt_eta;
    h_selection   = "1";
    h_cond_event  = cond_event;
    TH2D* h2D_dphi_vs_deta_MaxGMaxJ = new TH2D("h2D_dphi_vs_deta_MaxGMaxJ",h_formula.c_str(),35 ,0 ,3.5 ,35 ,0 ,3.5);
    h_histName    = h2D_dphi_vs_deta_MaxGMaxJ->GetName();
    std::cout << "--------------"   << std::endl;
    std::cout << "h_formula     = " << h_formula     <<std::endl;
    std::cout << "h_cond_photon = " << h_cond_photon <<std::endl;
    std::cout << "h_cond_jet    = " << h_cond_jet    <<std::endl;
    std::cout << "h_cond_event  = " << h_cond_event  <<std::endl;
    gja->drawPhotonJet(h_formula.c_str(),h_cond_photon.c_str(),h_cond_jet.c_str(),h_selection.c_str(),h_cond_event.c_str(),h_histName.c_str());

    // draw nJets (JET) : MAXPHOTONPT
    TString nJets_JET = Form("Sum$((%s))", cond_jet_pt_eta.Data());
    gja->resetCuts();
    gja->updateSelections();
    h_formula     = Form("%s:MAXPHOTONPT", nJets_JET.Data());
    h_cond_photon = gja->cond_photon;
    h_cond_jet    = cond_jet_pt_eta;
    h_selection   = cond_jet_pt_eta;
    h_cond_event  = cond_event;
    TH2D* h2D_nJets_MaxG = new TH2D("h2D_nJets_MaxG",h_formula.c_str(),44 ,0 ,220 ,7 ,0 ,7);
    h_histName    = h2D_nJets_MaxG->GetName();
    std::cout << "--------------"   << std::endl;
    std::cout << "h_formula     = " << h_formula     <<std::endl;
    std::cout << "h_cond_photon = " << h_cond_photon <<std::endl;
    std::cout << "h_cond_jet    = " << h_cond_jet    <<std::endl;
    std::cout << "h_cond_event  = " << h_cond_event  <<std::endl;
    gja->drawPhotonJet(h_formula.c_str(),h_cond_photon.c_str(),h_cond_jet.c_str(),h_selection.c_str(),h_cond_event.c_str(),h_histName.c_str());

    // draw leading photon
    // MAXPHOTONETA:MAXPHOTONPT
    gja->resetCuts();
    gja->updateSelections();
    h_formula     = "MAXPHOTONETA:MAXPHOTONPT";
    h_cond_photon = gja->cond_photon;
    h_cond_jet    = cond_jet_pt_eta;
    h_selection   = "1";
    h_cond_event  = cond_event;
    TH2D* h2D_phoEta_MaxG = new TH2D("h2D_phoEta_MaxG" ,h_formula.c_str() ,44 ,0 ,220 ,16 ,-gja->cut_eta ,gja->cut_eta);
    h_histName    = h2D_phoEta_MaxG->GetName();
    std::cout << "--------------"   << std::endl;
    std::cout << "h_formula     = " << h_formula     <<std::endl;
    std::cout << "h_cond_photon = " << h_cond_photon <<std::endl;
    std::cout << "h_cond_jet    = " << h_cond_jet    <<std::endl;
    std::cout << "h_cond_event  = " << h_cond_event  <<std::endl;
    gja->drawPhotonJet(h_formula.c_str(),h_cond_photon.c_str(),h_cond_jet.c_str(),h_selection.c_str(),h_cond_event.c_str(),h_histName.c_str());

    // draw leading photon
    // MAXPHOTONPHI:MAXPHOTONPT
    gja->resetCuts();
    gja->updateSelections();
    h_formula     = "MAXPHOTONPHI:MAXPHOTONPT";
    h_cond_photon = gja->cond_photon;
    h_cond_jet    = cond_jet_pt_eta;
    h_selection   = "1";
    h_cond_event  = cond_event;
    TH2D* h2D_phoPhi_MaxG = new TH2D("h2D_phoPhi_MaxG" ,h_formula.c_str() ,44 ,0 ,220 ,16 ,-PI ,PI);
    h_histName    = h2D_phoPhi_MaxG->GetName();
    std::cout << "--------------"   << std::endl;
    std::cout << "h_formula     = " << h_formula     <<std::endl;
    std::cout << "h_cond_photon = " << h_cond_photon <<std::endl;
    std::cout << "h_cond_jet    = " << h_cond_jet    <<std::endl;
    std::cout << "h_cond_event  = " << h_cond_event  <<std::endl;
    gja->drawPhotonJet(h_formula.c_str(),h_cond_photon.c_str(),h_cond_jet.c_str(),h_selection.c_str(),h_cond_event.c_str(),h_histName.c_str());

    // draw leading photon
    // MAXPHOTONETA:hiEvtPlanes[0]
    gja->resetCuts();
    gja->updateSelections();
    h_formula     = "MAXPHOTONETA:hiEvtPlanes[0]";
    h_cond_photon = gja->cond_photon;
    h_cond_jet    = cond_jet_pt_eta;
    h_selection   = "1";
    h_cond_event  = cond_event;
//    TH2D* h2D_phoEta_vs_hiEvtPlanes_MaxG = new TH2D("h2D_phoEta_vs_hiEvtPlanes_MaxG" ,h_formula.c_str() ,44 ,0 ,220 ,16 ,-4 ,4);
    h_histName    = "h2D_phoEta_vs_hiEvtPlanes_MaxG";
    std::cout << "--------------"   << std::endl;
    std::cout << "h_formula     = " << h_formula     <<std::endl;
    std::cout << "h_cond_photon = " << h_cond_photon <<std::endl;
    std::cout << "h_cond_jet    = " << h_cond_jet    <<std::endl;
    std::cout << "h_cond_event  = " << h_cond_event  <<std::endl;
    gja->drawPhotonJet(h_formula.c_str(),h_cond_photon.c_str(),h_cond_jet.c_str(),h_selection.c_str(),h_cond_event.c_str(),h_histName.c_str());

    // draw leading photon
    // MAXPHOTONPHI:hiEvtPlanes[0]
    gja->resetCuts();
    gja->updateSelections();
    h_formula     = "MAXPHOTONPHI:hiEvtPlanes[0]";
    h_cond_photon = gja->cond_photon;
    h_cond_jet    = cond_jet_pt_eta;
    h_selection   = "1";
    h_cond_event  = cond_event;
//    TH2D* h2D_phoPhi_vs_hiEvtPlanes_MaxG = new TH2D("h2D_phoPhi_vs_hiEvtPlanes_MaxG" ,h_formula.c_str() ,44 ,0 ,220 ,16 ,-PI ,PI);
    h_histName    = "h2D_phoPhi_vs_hiEvtPlanes_MaxG";
    std::cout << "--------------"   << std::endl;
    std::cout << "h_formula     = " << h_formula     <<std::endl;
    std::cout << "h_cond_photon = " << h_cond_photon <<std::endl;
    std::cout << "h_cond_jet    = " << h_cond_jet    <<std::endl;
    std::cout << "h_cond_event  = " << h_cond_event  <<std::endl;
    gja->drawPhotonJet(h_formula.c_str(),h_cond_photon.c_str(),h_cond_jet.c_str(),h_selection.c_str(),h_cond_event.c_str(),h_histName.c_str());

    std::cout<<"drawing finished"<<std::endl;

    outputFile->Write();

    inputFile->Close();
    outputFile->Close();
}

void GammaJetAnalyzerMacrosv1()
{
    const char* inputfileName="/mnt/hadoop/cms/store/user/luck/PbPb2011_photons_Data/HiForestPhoton-v7-noDuplicate.root";
    collisionType collision = HI;

    TFile *inputFile = new TFile(inputfileName, "READ");

    GammaJetAnalyzer* gja = new GammaJetAnalyzer(inputFile);
    if(collision == PP || collision == PA)   {
        gja->setJetTree(ak3PFJets);
    }
    else {
        gja->setJetTree(akPu3PFJets);
    }

    TString cond_eta_spike     = GammaJetAnalyzer::mergeSelections(gja->cond_pt_eta,  gja->cond_spike);
    TString cond_eta_spike_iso = GammaJetAnalyzer::mergeSelections(cond_eta_spike, gja->cond_iso);

    std::cout<<"drawing started"<<std::endl;

    //gja->drawMax("pt/jtpt", "pt", gja->cond_photon, gja->cond_event);
    gja->tree->Draw(Form("jtpt/Max$(pt*(%s))", gja->cond_photon.Data()), Form("Sum$(%s)>0",GammaJetAnalyzer::mergeSelections(gja->cond_photon,  gja->cond_event).Data()));
    TCanvas* c2=new TCanvas();
    c2->cd();
    gja->tree->Draw("nref", Form("Sum$(%s)>0",GammaJetAnalyzer::mergeSelections(gja->cond_photon,  gja->cond_event).Data()));

    TCanvas* c3=new TCanvas();
    c3->cd();
    gja->tree->Draw(Form("Max$(jtpt)/Max$(pt*(%s))", gja->cond_photon.Data()), Form("Sum$(%s)>0",GammaJetAnalyzer::mergeSelections(gja->cond_photon,  gja->cond_event).Data()));

    TCanvas* c4=new TCanvas();
    c4->cd();
    gja->tree->Draw(Form("Sum$(phi*(pt==Max$(pt*(%s)))):Sum$(jtphi*(jtpt==Max$(jtpt)))>>hcorr", gja->cond_photon.Data()), Form("Sum$(%s)>0",GammaJetAnalyzer::mergeSelections(gja->cond_photon,  gja->cond_event).Data()));

    std::cout<<gja->tree->GetEntries(GammaJetAnalyzer::mergeSelections(gja->cond_photon,  gja->cond_event))<<std::endl;

    std::cout<<"drawing finished"<<std::endl;

    inputFile->Close();
}

void GammaJetAnalyzerMacros()
{
//    std::cout << "running GammaJetAnalyzerMacrosv1()" << std::endl;
//    GammaJetAnalyzerMacrosv1();
    std::cout << "running GammaJetAnalyzerMacrosv2()" << std::endl;
    GammaJetAnalyzerMacrosv2();
}

int main(int argc, char** argv)
{
    GammaJetAnalyzerMacros();
}
