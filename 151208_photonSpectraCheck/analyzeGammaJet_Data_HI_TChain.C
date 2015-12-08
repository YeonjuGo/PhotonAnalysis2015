/*
 * macro to analyze and save di-Photon, di-Electron and di-Muon spectrum
 */

#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TCut.h>
#include <TH1.h>
#include <TH1D.h>
#include <TList.h>
#include <TLorentzVector.h>

#include <vector>
#include <string>
#include <iostream>
#include <iomanip>

#include "../treeUtil.h"
#include "../histoUtil.h"
#include "../EventMatcher.h"

const int MAXJETS = 500;
const float cutPt = 10;
const float cutEta = 1.4791;
const float cutDeltaR = 0.15;
const long MAXTREESIZE = 200000000000; // set maximum tree size from 10 GB to 100 GB, so that the code does not switch to a new file after 10 GB7

const double eleMass = 0.000511;
const double muMass  = 0.105658;

const bool looseIsolation  = true;
#define PI 3.141592653589
const float awayRange = PI * 7./8.;

void analyzeGammaJet_Data_HI_TChain(const char* inputFile, const char* outputFile = "out_analyzeDiPhoEleMu_Data_HI.root");
Double_t getDR  ( Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2);
Double_t getDPHI( Double_t phi1, Double_t phi2);
Double_t getDETA(Double_t eta1, Double_t eta2);

void analyzeGammaJet_Data_HI_TChain(const char* inputFile, const char* outputFile)
{
       std::cout<<"running analyzeGammaJet_Data_HI_TChain()"<<std::endl;
       std::cout<<"inputFile   = "<< inputFile <<std::endl;
       std::cout<<"outputFile  = "<< outputFile <<std::endl;

       TChain* treeHLT   = new TChain("hltanalysis/HltTree");
       TChain* tree      = new TChain("ggHiNtuplizer/EventTree");
       TChain* treePho   = new TChain("ggHiNtuplizer/EventTree");
       TChain* treeEvent = new TChain("ggHiNtuplizer/EventTree");
       TChain* treeJet   = new TChain("akPu4PFJetAnalyzer/t");
       TChain* treeHiEvt = new TChain("hiEvtAnalyzer/HiTree");
       TChain* treeSkim  = new TChain("skimanalysis/HltTree");

       int numFiles = 30;
       const char* fileNames[numFiles] =
       {
               // "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/velicanu/forest/Run2015E/HIExpressPhysics/Merged/HiForestStreamer_Run262314-262315.root"

//               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/velicanu/forest/HIRun2015/HIExpressPhysics/Merged/HIForestExpress_run262314-262315-v3.root"
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/velicanu/forest/HIRun2015/HIExpressPhysics/Merged/HIForestExpress_run262548-v6.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/velicanu/forest/HIRun2015/HIExpressPhysics/Merged/HIForestExpress_run262620-v6.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/velicanu/forest/HIRun2015/HIExpressPhysics/Merged/HIForestExpress_run262656.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/velicanu/forest/HIRun2015/HIExpressPhysics/Merged/HIForestExpress_run262694.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/velicanu/forest/HIRun2015/HIExpressPhysics/Merged/HIForestExpress_run262695.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/velicanu/forest/HIRun2015/HIExpressPhysics/Merged/HIForestExpress_run262697.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/velicanu/forest/HIRun2015/HIExpressPhysics/Merged/HIForestExpress_run262703.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/velicanu/forest/HIRun2015/HIExpressPhysics/Merged/HIForestExpress_run262726.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/velicanu/forest/HIRun2015/HIExpressPhysics/Merged/HIForestExpress_run262768.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/velicanu/forest/HIRun2015/HIExpressPhysics/Merged/HIForestExpress_run262777.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/velicanu/forest/HIRun2015/HIExpressPhysics/Merged/HIForestExpress_run262784.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/velicanu/forest/HIRun2015/HIExpressPhysics/Merged/HIForestExpress_run262811.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/velicanu/forest/HIRun2015/HIExpressPhysics/Merged/HIForestExpress_run262816.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/velicanu/forest/HIRun2015/HIExpressPhysics/Merged/HIForestExpress_run262817.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/velicanu/forest/HIRun2015/HIExpressPhysics/Merged/HIForestExpress_run262818.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/velicanu/forest/HIRun2015/HIExpressPhysics/Merged/HIForestExpress_run262819.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/velicanu/forest/HIRun2015/HIExpressPhysics/Merged/HIForestExpress_run262834.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/velicanu/forest/HIRun2015/HIExpressPhysics/Merged/HIForestExpress_run262836.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/velicanu/forest/HIRun2015/HIExpressPhysics/Merged/HIForestExpress_run262837_preHBHE.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/velicanu/forest/HIRun2015/HIExpressPhysics/Merged/HIForestExpress_run263035.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/velicanu/forest/HIRun2015/HIExpressPhysics/Merged/HIForestExpress_run263212.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/velicanu/forest/HIRun2015/HIExpressPhysics/Merged/HIForestExpress_run263213.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/velicanu/forest/HIRun2015/HIExpressPhysics/Merged/HIForestExpress_run262921.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/velicanu/forest/HIRun2015/HIExpressPhysics/Merged/HIForestExpress_run263233.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/velicanu/forest/HIRun2015/HIExpressPhysics/Merged/HIForestExpress_run263234.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/velicanu/forest/HIRun2015/HIExpressPhysics/Merged/HIForestExpress_run263261.root",

               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/velicanu/forest/HIRun2015/HIExpressPhysics/Merged/HIForestExpress_run263284.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/velicanu/forest/HIRun2015/HIExpressPhysics/Merged/HIForestExpress_run263286.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/velicanu/forest/HIRun2015/HIExpressPhysics/Merged/HIForestExpress_run263293.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/velicanu/forest/HIRun2015/HIExpressPhysics/Merged/HIForestExpress_run263322.root"
       };

       for (int i=0; i<numFiles; ++i)
       {
           treeHLT->Add(fileNames[i]);
           treePho->Add(fileNames[i]);
           treeEvent->Add(fileNames[i]);
           treeJet->Add(fileNames[i]);
           treeHiEvt->Add(fileNames[i]);
           treeSkim->Add(fileNames[i]);
           tree->Add(fileNames[i]);
       }

       treeHLT->SetBranchStatus("*",0);     // disable all branches
       treeHLT->SetBranchStatus("HLT_HI*SinglePhoton*Eta*v1*",1);     // enable photon branches
       treeHLT->SetBranchStatus("HLT_HI*DoublePhoton*Eta*v1*",1);     // enable photon branches
       treeHLT->SetBranchStatus("*DoubleMu*",1);     // enable muon branches

       treeEvent->SetBranchStatus("*",0);
       treeEvent->SetBranchStatus("run",1);
       treeEvent->SetBranchStatus("event",1);
       treeEvent->SetBranchStatus("lumis",1);

       treePho->SetBranchStatus("*",0);
       treePho->SetBranchStatus("run",1);
       treePho->SetBranchStatus("event",1);
       treePho->SetBranchStatus("lumis",1);
       treePho->SetBranchStatus("nPho",1);
       treePho->SetBranchStatus("pho*",1);

       treeJet->SetBranchStatus("*",0);        // disable all branches
       treeJet->SetBranchStatus("nref",1);     // enable jet branches
       treeJet->SetBranchStatus("jtpt",1);     // enable jet branches
       treeJet->SetBranchStatus("jteta",1);     // enable jet branches
       treeJet->SetBranchStatus("jtphi",1);     // enable jet branches
       treeJet->SetBranchStatus("track*",1);
       treeJet->SetBranchStatus("charged*",1);
       treeJet->SetBranchStatus("photon*",1);
       treeJet->SetBranchStatus("neutral*",1);
       treeJet->SetBranchStatus("eMax*",1);
       treeJet->SetBranchStatus("eSum*",1);
       treeJet->SetBranchStatus("eN*",1);

       Int_t nref;
       Float_t jtpt[MAXJETS];
       Float_t jteta[MAXJETS];
       Float_t jtphi[MAXJETS];

       treeJet->SetBranchAddress("nref",&nref);
       treeJet->SetBranchAddress("jtpt",jtpt);
       treeJet->SetBranchAddress("jteta",jteta);
       treeJet->SetBranchAddress("jtphi",jtphi);

       treeHiEvt->SetBranchStatus("*",0);
       treeHiEvt->SetBranchStatus("hiBin",1);
       treeHiEvt->SetBranchStatus("vz",1);
       
       treeSkim->SetBranchStatus("*",0);
       treeSkim->SetBranchStatus("ana_step",1);
       treeSkim->SetBranchStatus("pcollisionEventSelection",1);
       treeSkim->SetBranchStatus("pHBHENoiseFilterResultProducer",1);
       treeSkim->SetBranchStatus("HBHE*",1);

       // event information
       Int_t run, lumis;
       Long64_t event;
       treeEvent->SetBranchAddress("run", &run);
       treeEvent->SetBranchAddress("event", &event);
       treeEvent->SetBranchAddress("lumis", &lumis);

       // RECO photons
       Int_t nPho;
       std::vector<float>* phoE=0;
       std::vector<float>* phoEt=0;
       std::vector<float>* phoEta=0;
       std::vector<float>* phoPhi=0;
       std::vector<float>* pho_ecalClusterIsoR4=0;
       std::vector<float>* pho_hcalRechitIsoR4=0;
       std::vector<float>* pho_trackIsoR4PtCut20=0;
       std::vector<float>* phoR9=0;
       std::vector<float>* phoHoverE=0;
       std::vector<float>* phoSigmaIEtaIEta=0;
       std::vector<float>* phoSigmaIEtaIEta_2012=0;
       std::vector<float>* phoSCRawE=0;
       std::vector<float>* phoE5x5=0;
       std::vector<float>* phoSCEtaWidth=0;
       std::vector<float>* phoSCPhiWidth=0;
       std::vector<float>* pho_swissCrx=0;
       std::vector<float>* pho_seedTime=0;
       std::vector<float>* pfcIso4=0;
       std::vector<float>* pfnIso4=0;

       // RECO electrons
       Int_t nEle;
       std::vector<int>*   eleCharge=0;
       std::vector<float>* eleEn=0;
       std::vector<float>* elePt=0;
       std::vector<float>* eleEta=0;
       std::vector<float>* elePhi=0;
       std::vector<float>* eleHoverE=0;
       std::vector<float>* eleSigmaIEtaIEta=0;
       std::vector<float>* eleSigmaIEtaIEta_2012=0;
       std::vector<float>* eleSigmaIPhiIPhi=0;
       std::vector<float>* eleEoverPInv=0;
       std::vector<float>* eledEtaAtVtx=0;
       std::vector<float>* eledPhiAtVtx=0;
       std::vector<float>* eleD0=0;
       std::vector<float>* eleDz=0;
       std::vector<float>* eleTrkPt=0;

       // RECO muons
       Int_t nMu;
       std::vector<int>*   muCharge=0;
       std::vector<float>* muPt=0;
       std::vector<float>* muEta=0;
       std::vector<float>* muPhi=0;
       std::vector<int>*   muIsGood=0;
       std::vector<float>* muD0=0;
       std::vector<float>* muDz=0;
       std::vector<float>* muChi2NDF=0;
       std::vector<float>* muInnerD0=0;
       std::vector<float>* muInnerDz=0;
       std::vector<float>* muIsoTrk=0;
       std::vector<float>* muPFChIso=0;
       std::vector<float>* muPFPhoIso=0;
       std::vector<float>* muPFNeuIso=0;
       std::vector<float>* muPFPUIso=0;

       tree->SetBranchStatus("*",0);     // disable all branches
       tree->SetBranchStatus("nPho",1);     // enable photon branches
       tree->SetBranchStatus("pho*",1);     // enable photon branches
       tree->SetBranchAddress("nPho",&nPho);
       tree->SetBranchAddress("phoE",&phoE);
       tree->SetBranchAddress("phoEt",&phoEt);
       tree->SetBranchAddress("phoEta",&phoEta);
       tree->SetBranchAddress("phoPhi",&phoPhi);
       tree->SetBranchAddress("pho_ecalClusterIsoR4",&pho_ecalClusterIsoR4);
       tree->SetBranchAddress("pho_hcalRechitIsoR4",&pho_hcalRechitIsoR4);
       tree->SetBranchAddress("pho_trackIsoR4PtCut20",&pho_trackIsoR4PtCut20);
       tree->SetBranchAddress("phoR9",&phoR9);
       tree->SetBranchAddress("phoHoverE",&phoHoverE);
       tree->SetBranchAddress("phoSigmaIEtaIEta",&phoSigmaIEtaIEta);
       tree->SetBranchAddress("phoSigmaIEtaIEta_2012",&phoSigmaIEtaIEta_2012);
       tree->SetBranchAddress("phoSCRawE",&phoSCRawE);
       tree->SetBranchAddress("phoE5x5",&phoE5x5);
       tree->SetBranchAddress("phoSCEtaWidth",&phoSCEtaWidth);
       tree->SetBranchAddress("phoSCPhiWidth",&phoSCPhiWidth);
       tree->SetBranchAddress("pho_swissCrx",&pho_swissCrx);
       tree->SetBranchAddress("pho_seedTime",&pho_seedTime);
       tree->SetBranchAddress("pfcIso4",&pfcIso4);
       tree->SetBranchAddress("pfnIso4",&pfnIso4);

       tree->SetBranchStatus("nEle",1);     // enable electron branches
       tree->SetBranchStatus("ele*",1);     // enable electron branches
       tree->SetBranchAddress("nEle",&nEle);
       tree->SetBranchAddress("eleEn",&eleEn);
       tree->SetBranchAddress("eleCharge",&eleCharge);
       tree->SetBranchAddress("elePt",&elePt);
       tree->SetBranchAddress("eleEta",&eleEta);
       tree->SetBranchAddress("elePhi",&elePhi);
       tree->SetBranchAddress("eleHoverE",&eleHoverE);
       tree->SetBranchAddress("eleSigmaIEtaIEta",&eleSigmaIEtaIEta);
       tree->SetBranchAddress("eleSigmaIEtaIEta_2012",&eleSigmaIEtaIEta_2012);
       tree->SetBranchAddress("eleSigmaIPhiIPhi",&eleSigmaIPhiIPhi);
       tree->SetBranchAddress("eleEoverPInv",&eleEoverPInv);
       tree->SetBranchAddress("eledEtaAtVtx",&eledEtaAtVtx);
       tree->SetBranchAddress("eledPhiAtVtx",&eledPhiAtVtx);
       tree->SetBranchAddress("eleD0",&eleD0);
       tree->SetBranchAddress("eleDz",&eleDz);
       tree->SetBranchAddress("eleTrkPt",&eleTrkPt);

       tree->SetBranchStatus("nMu",1);     // enable muon branches
       tree->SetBranchStatus("mu*",1);     // enable muon branches
       tree->SetBranchAddress("nMu",&nMu);
       tree->SetBranchAddress("muCharge",&muCharge);
       tree->SetBranchAddress("muPt",&muPt);
       tree->SetBranchAddress("muEta",&muEta);
       tree->SetBranchAddress("muPhi",&muPhi);
       tree->SetBranchAddress("muIsGood",&muIsGood);
       tree->SetBranchAddress("muD0",&muD0);
       tree->SetBranchAddress("muDz",&muDz);
       tree->SetBranchAddress("muChi2NDF",&muChi2NDF);
       tree->SetBranchAddress("muInnerD0",&muInnerD0);
       tree->SetBranchAddress("muInnerDz",&muInnerDz);
       tree->SetBranchAddress("muIsoTrk",&muIsoTrk);
       tree->SetBranchAddress("muPFChIso",&muPFChIso);
       tree->SetBranchAddress("muPFPhoIso",&muPFPhoIso);
       tree->SetBranchAddress("muPFNeuIso",&muPFNeuIso);
       tree->SetBranchAddress("muPFPUIso",&muPFPUIso);

       TFile* output = new TFile(outputFile,"RECREATE");
       // output tree variables

       TTree *outputTreeHLT    = treeHLT->CloneTree(0);
       outputTreeHLT->SetName("hltTree");
       TTree *outputTreePho    = treePho->CloneTree(0);
       outputTreePho->SetName("photons");
       outputTreePho->SetTitle("Event data + photons");
       TTree *outputTreeJet    = treeJet->CloneTree(0);
       outputTreeJet->SetName("jets");
       TTree *outputTreeHiEvt  = treeHiEvt->CloneTree(0);
       outputTreeHiEvt->SetName("HiEvt");
       outputTreeHiEvt->SetTitle("subbranches of hltanalysis/HltTree");
       TTree *outputTreeSkim   = treeSkim->CloneTree(0);
       outputTreeSkim->SetName("skim");
       outputTreeSkim->SetTitle("subbranches of skimanalysis/HltTree");

       TTree *outputTreeGammaJet = new TTree("gammaJet","leading photon-jet correlations");
       TTree *outputTreediPho  = new TTree("diphoton","");
       TTree *outputTreediEle  = new TTree("dielectron","");
       TTree *outputTreediMu   = new TTree("dimuon","");

       outputTreeHLT->SetMaxTreeSize(MAXTREESIZE);
       outputTreediPho->SetMaxTreeSize(MAXTREESIZE);
       outputTreediEle->SetMaxTreeSize(MAXTREESIZE);
       outputTreediMu->SetMaxTreeSize(MAXTREESIZE);

       Int_t phoIdx;
       std::vector<int>   jetIdx;
       std::vector<float> xjg;
       std::vector<float> deta;
       std::vector<float> dphi;
       std::vector<float> dR;
       std::vector<int>   insideJet;
       Int_t nJetin7over8;

       outputTreeGammaJet->Branch("phoIdx",&phoIdx);
       outputTreeGammaJet->Branch("jetIdx",&jetIdx);
       outputTreeGammaJet->Branch("xjg",&xjg);
       outputTreeGammaJet->Branch("deta",&deta);
       outputTreeGammaJet->Branch("dphi",&dphi);
       outputTreeGammaJet->Branch("dR",&dR);
       outputTreeGammaJet->Branch("insideJet",&insideJet);
       outputTreeGammaJet->Branch("nJetin7over8",&nJetin7over8);

       // output tree variables
       std::vector<float> phoEt_1;
       std::vector<float> phoEta_1;
       std::vector<float> phoPhi_1;
       std::vector<float> pho_ecalClusterIsoR4_1;
       std::vector<float> pho_hcalRechitIsoR4_1;
       std::vector<float> pho_trackIsoR4PtCut20_1;
       std::vector<float> phoR9_1;
       std::vector<float> phoHoverE_1;
       std::vector<float> phoSigmaIEtaIEta_1;
       std::vector<float> phoSigmaIEtaIEta_2012_1;
       std::vector<float> phoSCRawE_1;
       std::vector<float> phoE5x5_1;
       std::vector<float> phoSCEtaWidth_1;
       std::vector<float> phoSCPhiWidth_1;
       std::vector<float> pho_swissCrx_1;
       std::vector<float> pho_seedTime_1;
       std::vector<int>   matched_eleCharge_1;
       std::vector<int>   matched_eleIndex_1;
       std::vector<float> matched_elePt_1;
       std::vector<float> matched_eleEta_1;
       std::vector<float> matched_elePhi_1;
       std::vector<float> matched_eleTrkPt_1;

       std::vector<float> phoEt_2;
       std::vector<float> phoEta_2;
       std::vector<float> phoPhi_2;
       std::vector<float> pho_ecalClusterIsoR4_2;
       std::vector<float> pho_hcalRechitIsoR4_2;
       std::vector<float> pho_trackIsoR4PtCut20_2;
       std::vector<float> phoR9_2;
       std::vector<float> phoHoverE_2;
       std::vector<float> phoSigmaIEtaIEta_2;
       std::vector<float> phoSigmaIEtaIEta_2012_2;
       std::vector<float> phoSCRawE_2;
       std::vector<float> phoE5x5_2;
       std::vector<float> phoSCEtaWidth_2;
       std::vector<float> phoSCPhiWidth_2;
       std::vector<float> pho_swissCrx_2;
       std::vector<float> pho_seedTime_2;
       std::vector<int>   matched_eleCharge_2;
       std::vector<int>   matched_eleIndex_2;
       std::vector<float> matched_elePt_2;
       std::vector<float> matched_eleEta_2;
       std::vector<float> matched_elePhi_2;
       std::vector<float> matched_eleTrkPt_2;

       std::vector<float> vSum_gg_M;
       std::vector<float> vSum_gg_Energy;
       std::vector<float> vSum_gg_Pt;
       std::vector<float> vSum_gg_Eta;
       std::vector<float> vSum_gg_Phi;

       std::vector<float> vSum_gg_M_elePt;
       std::vector<float> vSum_gg_M_eleTrkPt;

       outputTreediPho->Branch("nPho",&nPho);
       outputTreediPho->Branch("phoEt_1",&phoEt_1);
       outputTreediPho->Branch("phoEta_1",&phoEta_1);
       outputTreediPho->Branch("phoPhi_1",&phoPhi_1);
       outputTreediPho->Branch("pho_ecalClusterIsoR4_1",&pho_ecalClusterIsoR4_1);
       outputTreediPho->Branch("pho_hcalRechitIsoR4_1",&pho_hcalRechitIsoR4_1);
       outputTreediPho->Branch("pho_trackIsoR4PtCut20_1",&pho_trackIsoR4PtCut20_1);
       outputTreediPho->Branch("phoR9_1",&phoR9_1);
       outputTreediPho->Branch("phoHoverE_1",&phoHoverE_1);
       outputTreediPho->Branch("phoSigmaIEtaIEta_1",&phoSigmaIEtaIEta_1);
       outputTreediPho->Branch("phoSigmaIEtaIEta_2012_1",&phoSigmaIEtaIEta_2012_1);
       outputTreediPho->Branch("phoSCRawE_1",&phoSCRawE_1);
       outputTreediPho->Branch("phoE5x5_1",&phoE5x5_1);
       outputTreediPho->Branch("phoSCEtaWidth_1",&phoSCEtaWidth_1);
       outputTreediPho->Branch("phoSCPhiWidth_1",&phoSCPhiWidth_1);
       outputTreediPho->Branch("pho_swissCrx_1",&pho_swissCrx_1);
       outputTreediPho->Branch("pho_seedTime_1",&pho_seedTime_1);
       outputTreediPho->Branch("matched_eleCharge_1",&matched_eleCharge_1);
       outputTreediPho->Branch("matched_eleIndex_1",&matched_eleIndex_1);
       outputTreediPho->Branch("matched_elePt_1",&matched_elePt_1);
       outputTreediPho->Branch("matched_eleEta_1",&matched_eleEta_1);
       outputTreediPho->Branch("matched_elePhi_1",&matched_elePhi_1);
       outputTreediPho->Branch("matched_eleTrkPt_1",&matched_eleTrkPt_1);

       outputTreediPho->Branch("phoEt_2",&phoEt_2);
       outputTreediPho->Branch("phoEta_2",&phoEta_2);
       outputTreediPho->Branch("phoPhi_2",&phoPhi_2);
       outputTreediPho->Branch("pho_ecalClusterIsoR4_2",&pho_ecalClusterIsoR4_2);
       outputTreediPho->Branch("pho_hcalRechitIsoR4_2",&pho_hcalRechitIsoR4_2);
       outputTreediPho->Branch("pho_trackIsoR4PtCut20_2",&pho_trackIsoR4PtCut20_2);
       outputTreediPho->Branch("phoR9_2",&phoR9_2);
       outputTreediPho->Branch("phoHoverE_2",&phoHoverE_2);
       outputTreediPho->Branch("phoSigmaIEtaIEta_2",&phoSigmaIEtaIEta_2);
       outputTreediPho->Branch("phoSigmaIEtaIEta_2012_2",&phoSigmaIEtaIEta_2012_2);
       outputTreediPho->Branch("phoSCRawE_2",&phoSCRawE_2);
       outputTreediPho->Branch("phoE5x5_2",&phoE5x5_2);
       outputTreediPho->Branch("phoSCEtaWidth_2",&phoSCEtaWidth_2);
       outputTreediPho->Branch("phoSCPhiWidth_2",&phoSCPhiWidth_2);
       outputTreediPho->Branch("pho_swissCrx_2",&pho_swissCrx_2);
       outputTreediPho->Branch("pho_seedTime_2",&pho_seedTime_2);
       outputTreediPho->Branch("matched_eleCharge_2",&matched_eleCharge_2);
       outputTreediPho->Branch("matched_eleIndex_2",&matched_eleIndex_2);
       outputTreediPho->Branch("matched_elePt_2",&matched_elePt_2);
       outputTreediPho->Branch("matched_eleEta_2",&matched_eleEta_2);
       outputTreediPho->Branch("matched_elePhi_2",&matched_elePhi_2);
       outputTreediPho->Branch("matched_eleTrkPt_2",&matched_eleTrkPt_2);

       outputTreediPho->Branch("vSum_M",&vSum_gg_M);
       outputTreediPho->Branch("vSum_Energy",&vSum_gg_Energy);
       outputTreediPho->Branch("vSum_Pt",&vSum_gg_Pt);
       outputTreediPho->Branch("vSum_Eta",&vSum_gg_Eta);
       outputTreediPho->Branch("vSum_Phi",&vSum_gg_Phi);

       outputTreediPho->Branch("vSum_M_elePt",&vSum_gg_M_elePt);
       outputTreediPho->Branch("vSum_M_eleTrkPt",&vSum_gg_M_eleTrkPt);

       std::vector<int>   eleCharge_1;
       std::vector<float> elePt_1;
       std::vector<float> eleEta_1;
       std::vector<float> elePhi_1;
       std::vector<float> eleHoverE_1;
       std::vector<float> eleSigmaIEtaIEta_1;
       std::vector<float> eleSigmaIEtaIEta_2012_1;
       std::vector<float> eleSigmaIPhiIPhi_1;
       std::vector<float> eleEoverPInv_1;
       std::vector<float> eledEtaAtVtx_1;
       std::vector<float> eledPhiAtVtx_1;
       std::vector<float> eleD0_1;
       std::vector<float> eleDz_1;

       std::vector<int>   eleCharge_2;
       std::vector<float> elePt_2;
       std::vector<float> eleEta_2;
       std::vector<float> elePhi_2;
       std::vector<float> eleHoverE_2;
       std::vector<float> eleSigmaIEtaIEta_2;
       std::vector<float> eleSigmaIEtaIEta_2012_2;
       std::vector<float> eleSigmaIPhiIPhi_2;
       std::vector<float> eleEoverPInv_2;
       std::vector<float> eledEtaAtVtx_2;
       std::vector<float> eledPhiAtVtx_2;
       std::vector<float> eleD0_2;
       std::vector<float> eleDz_2;

       std::vector<float> vSum_ee_M;
       std::vector<float> vSum_ee_Energy;
       std::vector<float> vSum_ee_Pt;
       std::vector<float> vSum_ee_Eta;
       std::vector<float> vSum_ee_Phi;

       outputTreediEle->Branch("nEle",&nEle);
       outputTreediEle->Branch("eleCharge_1",&eleCharge_1);
       outputTreediEle->Branch("elePt_1",&elePt_1);
       outputTreediEle->Branch("eleEta_1",&eleEta_1);
       outputTreediEle->Branch("elePhi_1",&elePhi_1);
       outputTreediEle->Branch("eleHoverE_1",&eleHoverE_1);
       outputTreediEle->Branch("eleSigmaIEtaIEta_1",&eleSigmaIEtaIEta_1);
       outputTreediEle->Branch("eleSigmaIEtaIEta_2012_1",&eleSigmaIEtaIEta_2012_1);
       outputTreediEle->Branch("eleSigmaIPhiIPhi_1",&eleSigmaIPhiIPhi_1);
       outputTreediEle->Branch("eleEoverPInv_1",&eleEoverPInv_1);
       outputTreediEle->Branch("eledEtaAtVtx_1",&eledEtaAtVtx_1);
       outputTreediEle->Branch("eledPhiAtVtx_1",&eledPhiAtVtx_1);
       outputTreediEle->Branch("eleD0_1",&eleD0_1);
       outputTreediEle->Branch("eleDz_1",&eleDz_1);

       outputTreediEle->Branch("eleCharge_2",&eleCharge_2);
       outputTreediEle->Branch("elePt_2",&elePt_2);
       outputTreediEle->Branch("eleEta_2",&eleEta_2);
       outputTreediEle->Branch("elePhi_2",&elePhi_2);
       outputTreediEle->Branch("eleHoverE_2",&eleHoverE_2);
       outputTreediEle->Branch("eleSigmaIEtaIEta_2",&eleSigmaIEtaIEta_2);
       outputTreediEle->Branch("eleSigmaIEtaIEta_2012_2",&eleSigmaIEtaIEta_2012_2);
       outputTreediEle->Branch("eleSigmaIPhiIPhi_2",&eleSigmaIPhiIPhi_2);
       outputTreediEle->Branch("eleEoverPInv_2",&eleEoverPInv_2);
       outputTreediEle->Branch("eledEtaAtVtx_2",&eledEtaAtVtx_2);
       outputTreediEle->Branch("eledPhiAtVtx_2",&eledPhiAtVtx_2);
       outputTreediEle->Branch("eleD0_2",&eleD0_2);
       outputTreediEle->Branch("eleDz_2",&eleDz_2);

       outputTreediEle->Branch("vSum_M",&vSum_ee_M);
       outputTreediEle->Branch("vSum_Energy",&vSum_ee_Energy);
       outputTreediEle->Branch("vSum_Pt",&vSum_ee_Pt);
       outputTreediEle->Branch("vSum_Eta",&vSum_ee_Eta);
       outputTreediEle->Branch("vSum_Phi",&vSum_ee_Phi);

       // RECO muons
       std::vector<int>   muCharge_1;
       std::vector<float> muPt_1;
       std::vector<float> muEta_1;
       std::vector<float> muPhi_1;
       std::vector<int>   muIsGood_1;
       std::vector<float> muD0_1;
       std::vector<float> muDz_1;
       std::vector<float> muChi2NDF_1;
       std::vector<float> muInnerD0_1;
       std::vector<float> muInnerDz_1;
       std::vector<float> muIsoTrk_1;
       std::vector<float> muPFChIso_1;
       std::vector<float> muPFPhoIso_1;
       std::vector<float> muPFNeuIso_1;
       std::vector<float> muPFPUIso_1;

       std::vector<int>   muCharge_2;
       std::vector<float> muPt_2;
       std::vector<float> muEta_2;
       std::vector<float> muPhi_2;
       std::vector<int>   muIsGood_2;
       std::vector<float> muD0_2;
       std::vector<float> muDz_2;
       std::vector<float> muChi2NDF_2;
       std::vector<float> muInnerD0_2;
       std::vector<float> muInnerDz_2;
       std::vector<float> muIsoTrk_2;
       std::vector<float> muPFChIso_2;
       std::vector<float> muPFPhoIso_2;
       std::vector<float> muPFNeuIso_2;
       std::vector<float> muPFPUIso_2;

       std::vector<float> vSum_mumu_M;
       std::vector<float> vSum_mumu_Energy;
       std::vector<float> vSum_mumu_Pt;
       std::vector<float> vSum_mumu_Eta;
       std::vector<float> vSum_mumu_Phi;

       outputTreediMu->Branch("nMu",&nMu);
       outputTreediMu->Branch("muCharge_1",&muCharge_1);
       outputTreediMu->Branch("muPt_1",&muPt_1);
       outputTreediMu->Branch("muEta_1",&muEta_1);
       outputTreediMu->Branch("muPhi_1",&muPhi_1);
       outputTreediMu->Branch("muIsGood_1",&muIsGood_1);
       outputTreediMu->Branch("muD0_1",&muD0_1);
       outputTreediMu->Branch("muDz_1",&muDz_1);
       outputTreediMu->Branch("muChi2NDF_1",&muChi2NDF_1);
       outputTreediMu->Branch("muInnerD0_1",&muInnerD0_1);
       outputTreediMu->Branch("muInnerDz_1",&muInnerDz_1);
       outputTreediMu->Branch("muIsoTrk_1",&muIsoTrk_1);
       outputTreediMu->Branch("muPFChIso_1",&muPFChIso_1);
       outputTreediMu->Branch("muPFPhoIso_1",&muPFPhoIso_1);
       outputTreediMu->Branch("muPFNeuIso_1",&muPFNeuIso_1);
       outputTreediMu->Branch("muPFPUIso_1",&muPFPUIso_1);

       outputTreediMu->Branch("muCharge_2",&muCharge_2);
       outputTreediMu->Branch("muPt_2",&muPt_2);
       outputTreediMu->Branch("muEta_2",&muEta_2);
       outputTreediMu->Branch("muPhi_2",&muPhi_2);
       outputTreediMu->Branch("muIsGood_2",&muIsGood_2);
       outputTreediMu->Branch("muD0_2",&muD0_2);
       outputTreediMu->Branch("muDz_2",&muDz_2);
       outputTreediMu->Branch("muChi2NDF_2",&muChi2NDF_2);
       outputTreediMu->Branch("muInnerD0_2",&muInnerD0_2);
       outputTreediMu->Branch("muInnerDz_2",&muInnerDz_2);
       outputTreediMu->Branch("muIsoTrk_2",&muIsoTrk_2);
       outputTreediMu->Branch("muPFChIso_2",&muPFChIso_2);
       outputTreediMu->Branch("muPFPhoIso_2",&muPFPhoIso_2);
       outputTreediMu->Branch("muPFNeuIso_2",&muPFNeuIso_2);
       outputTreediMu->Branch("muPFPUIso_2",&muPFPUIso_2);

       outputTreediMu->Branch("vSum_M",&vSum_mumu_M);
       outputTreediMu->Branch("vSum_Energy",&vSum_mumu_Energy);
       outputTreediMu->Branch("vSum_Pt",&vSum_mumu_Pt);
       outputTreediMu->Branch("vSum_Eta",&vSum_mumu_Eta);
       outputTreediMu->Branch("vSum_Phi",&vSum_mumu_Phi);

       EventMatcher* em = new EventMatcher();
       Long64_t duplicateEntries = 0;

       Long64_t entries = treeEvent->GetEntries();
       Long64_t entriesAnalyzed = 0;
       std::cout << "entries         = " << entries << std::endl;
       std::cout<< "Loop : ggHiNtuplizer/EventTree" <<std::endl;
       for (Long64_t j_entry=0; j_entry<entries; ++j_entry)
       {
           if (j_entry % 2000 == 0)  {
             std::cout << "current entry = " <<j_entry<<" out of "<<entries<<" : "<<std::setprecision(2)<<(double)j_entry/entries*100<<" %"<<std::endl;
           }

           treeHLT->GetEntry(j_entry);
           treeEvent->GetEntry(j_entry);
           treePho->GetEntry(j_entry);
           treeJet->GetEntry(j_entry);
           treeHiEvt->GetEntry(j_entry);
           treeSkim->GetEntry(j_entry);
           tree->GetEntry(j_entry);

           bool eventAdded = em->addEvent(run,lumis,event,j_entry);
           if(!eventAdded) // this event is duplicate, skip this one.
           {
               duplicateEntries++;
               continue;
           }

           // photon-jet block
           phoIdx = -1;
           double maxPhoEt = -1;
           for(int i=0; i<nPho; ++i)
           {
               bool passedEtaCut = TMath::Abs(phoEta->at(i)) < cutEta ;
               bool passedSpikeRejection = (phoSigmaIEtaIEta->at(i)  > 0.002  &&
                                            pho_swissCrx->at(i)      < 0.9    &&
                                            TMath::Abs(pho_seedTime->at(i) ) < 3);

               bool passedIsolation;
               if (looseIsolation){
                   passedIsolation = (pho_ecalClusterIsoR4->at(i)  < 4.2 &&
                                      pho_hcalRechitIsoR4->at(i)   < 2.2 &&
                                      pho_trackIsoR4PtCut20->at(i) < 2   &&
                                      phoHoverE->at(i) < 0.1);
               }
               else
               {
                   passedIsolation = ((pho_ecalClusterIsoR4->at(i) + pho_hcalRechitIsoR4->at(i) + pho_trackIsoR4PtCut20->at(i)) < 1 &&
                           phoHoverE->at(i) < 0.1) ;
               }

               bool passedPurity = (phoSigmaIEtaIEta->at(i) < 0.01);

               if (!(passedEtaCut && passedSpikeRejection && passedIsolation && passedPurity)) continue;

               if (phoEt->at(i) > maxPhoEt)
               {
                   maxPhoEt = phoEt->at(i);
                   phoIdx = i;
               }
           }
           if (phoIdx == -1) continue;
           entriesAnalyzed++;
           // photon-jet correlation)
           nJetin7over8 = 0;

           jetIdx.clear();
           xjg.clear();
           deta.clear();
           dphi.clear();
           dR.clear();
           insideJet.clear();
           for (int i=0; i<nref; ++i)
           {
               // cuts on jets will be applied during plotting

               float tmp_deta = getDETA(phoEta->at(phoIdx), jteta[i]);
               float tmp_dphi = getDPHI(phoPhi->at(phoIdx), jtphi[i]);
               if (tmp_dphi > awayRange)
                   nJetin7over8++;
               float tmp_dR   = getDR(phoEta->at(phoIdx), phoPhi->at(phoIdx), jteta[i], jtphi[i]);

               int tmp_insideJet;
               if(tmp_dR < 0.4) tmp_insideJet = 1;
               else             tmp_insideJet = 0;

               jetIdx.push_back(i);
               xjg.push_back((float)jtpt[i]/phoEt->at(phoIdx));
               deta.push_back(tmp_deta);
               dphi.push_back(tmp_dphi);
               dR.push_back(tmp_dR);
               insideJet.push_back(tmp_insideJet);
           }

           // diphoton block
           phoEt_1.clear();
           phoEta_1.clear();
           phoPhi_1.clear();
           pho_ecalClusterIsoR4_1.clear();
           pho_hcalRechitIsoR4_1.clear();
           pho_trackIsoR4PtCut20_1.clear();
           phoR9_1.clear();
           phoHoverE_1.clear();
           phoSigmaIEtaIEta_1.clear();
           phoSigmaIEtaIEta_2012_1.clear();
           phoSCRawE_1.clear();
           phoE5x5_1.clear();
           phoSCEtaWidth_1.clear();
           phoSCPhiWidth_1.clear();
           pho_swissCrx_1.clear();
           pho_seedTime_1.clear();
           matched_eleCharge_1.clear();
           matched_eleIndex_1.clear();
           matched_elePt_1.clear();
           matched_eleEta_1.clear();
           matched_elePhi_1.clear();
           matched_eleTrkPt_1.clear();

           phoEt_2.clear();
           phoEta_2.clear();
           phoPhi_2.clear();
           pho_ecalClusterIsoR4_2.clear();
           pho_hcalRechitIsoR4_2.clear();
           pho_trackIsoR4PtCut20_2.clear();
           phoR9_2.clear();
           phoHoverE_2.clear();
           phoSigmaIEtaIEta_2.clear();
           phoSigmaIEtaIEta_2012_2.clear();
           phoSCRawE_2.clear();
           phoE5x5_2.clear();
           phoSCEtaWidth_2.clear();
           phoSCPhiWidth_2.clear();
           pho_swissCrx_2.clear();
           pho_seedTime_2.clear();
           matched_eleCharge_2.clear();
           matched_eleIndex_2.clear();
           matched_elePt_2.clear();
           matched_eleEta_2.clear();
           matched_elePhi_2.clear();
           matched_eleTrkPt_2.clear();

           vSum_gg_M.clear();
           vSum_gg_Energy.clear();
           vSum_gg_Pt.clear();
           vSum_gg_Eta.clear();
           vSum_gg_Phi.clear();

           vSum_gg_M_elePt.clear();
           vSum_gg_M_eleTrkPt.clear();

           std::vector<int>   matched_eleCharge_tmp;
           std::vector<int>   matched_eleIndex_tmp;
           std::vector<float> matched_elePt_tmp;
           std::vector<float> matched_eleEta_tmp;
           std::vector<float> matched_elePhi_tmp;
           std::vector<float> matched_eleTrkPt_tmp;

           matched_eleCharge_tmp.clear();
           matched_eleIndex_tmp.clear();
           matched_elePt_tmp.clear();
           matched_eleEta_tmp.clear();
           matched_elePhi_tmp.clear();
           matched_eleTrkPt_tmp.clear();

           // electron-photon matching
           // based on  https://github.com/CmsHI/HiForestAnalysis/blob/master/PhotonUtilities.C#L154-L186
           for(int i=0; i<nPho; ++i) {
               float currentMaxPt=-1;
               int matched_Index =-1;
               for (int j=0; j<nEle; ++j) {
                   if (elePt->at(j) < currentMaxPt) continue;
                   if (getDR(eleEta->at(j), elePhi->at(j), phoEta->at(i), phoPhi->at(i))>cutDeltaR) continue;

                   matched_Index=j;
                   currentMaxPt=elePt->at(j);
               }
               matched_eleIndex_tmp.push_back(matched_Index);

               if (matched_Index > -1) {
                   matched_eleCharge_tmp.push_back(eleCharge->at(matched_Index));
                   matched_elePt_tmp.push_back(elePt->at(matched_Index));
                   matched_eleEta_tmp.push_back(eleEta->at(matched_Index));
                   matched_elePhi_tmp.push_back(elePhi->at(matched_Index));
                   matched_eleTrkPt_tmp.push_back(eleTrkPt->at(matched_Index));
               } else {
                   matched_eleCharge_tmp.push_back(0);
                   matched_elePt_tmp.push_back(0);
                   matched_eleEta_tmp.push_back(0);
                   matched_elePhi_tmp.push_back(0);
                   matched_eleTrkPt_tmp.push_back(0);
               }
           }

           for(int i=0; i<nPho; ++i)
           {
               for (int j=i+1; j<nPho; ++j)
               {
                   phoEt_1.push_back(phoEt->at(i));
                   phoEta_1.push_back(phoEta->at(i));
                   phoPhi_1.push_back(phoPhi->at(i));
                   pho_ecalClusterIsoR4_1.push_back(pho_ecalClusterIsoR4->at(i));
                   pho_hcalRechitIsoR4_1.push_back(pho_hcalRechitIsoR4->at(i));
                   pho_trackIsoR4PtCut20_1.push_back(pho_trackIsoR4PtCut20->at(i));
                   phoR9_1.push_back(phoR9->at(i));
                   phoHoverE_1.push_back(phoHoverE->at(i));
                   phoSigmaIEtaIEta_1.push_back(phoSigmaIEtaIEta->at(i));
                   phoSigmaIEtaIEta_2012_1.push_back(phoSigmaIEtaIEta_2012->at(i));
                   phoSCRawE_1.push_back(phoSCRawE->at(i));
                   phoE5x5_1.push_back(phoE5x5->at(i));
                   phoSCEtaWidth_1.push_back(phoSCEtaWidth->at(i));
                   phoSCPhiWidth_1.push_back(phoSCPhiWidth->at(i));
                   pho_swissCrx_1.push_back(pho_swissCrx->at(i));
                   pho_seedTime_1.push_back(pho_seedTime->at(i));

                   matched_eleIndex_1.push_back(matched_eleIndex_tmp.at(i));
                   matched_eleCharge_1.push_back(matched_eleCharge_tmp.at(i));
                   matched_elePt_1.push_back(matched_elePt_tmp.at(i));
                   matched_eleEta_1.push_back(matched_eleEta_tmp.at(i));
                   matched_elePhi_1.push_back(matched_elePhi_tmp.at(i));
                   matched_eleTrkPt_1.push_back(matched_eleTrkPt_tmp.at(i));

                   phoEt_2.push_back(phoEt->at(j));
                   phoEta_2.push_back(phoEta->at(j));
                   phoPhi_2.push_back(phoPhi->at(j));
                   pho_ecalClusterIsoR4_2.push_back(pho_ecalClusterIsoR4->at(j));
                   pho_hcalRechitIsoR4_2.push_back(pho_hcalRechitIsoR4->at(j));
                   pho_trackIsoR4PtCut20_2.push_back(pho_trackIsoR4PtCut20->at(j));
                   phoR9_2.push_back(phoR9->at(j));
                   phoHoverE_2.push_back(phoHoverE->at(j));
                   phoSigmaIEtaIEta_2.push_back(phoSigmaIEtaIEta->at(j));
                   phoSigmaIEtaIEta_2012_2.push_back(phoSigmaIEtaIEta_2012->at(j));
                   phoSCRawE_2.push_back(phoSCRawE->at(j));
                   phoE5x5_2.push_back(phoE5x5->at(j));
                   phoSCEtaWidth_2.push_back(phoSCEtaWidth->at(j));
                   phoSCPhiWidth_2.push_back(phoSCPhiWidth->at(j));
                   pho_swissCrx_2.push_back(pho_swissCrx->at(j));
                   pho_seedTime_2.push_back(pho_seedTime->at(j));

                   matched_eleIndex_2.push_back(matched_eleIndex_tmp.at(j));
                   matched_eleCharge_2.push_back(matched_eleCharge_tmp.at(j));
                   matched_elePt_2.push_back(matched_elePt_tmp.at(j));
                   matched_eleEta_2.push_back(matched_eleEta_tmp.at(j));
                   matched_elePhi_2.push_back(matched_elePhi_tmp.at(j));
                   matched_eleTrkPt_2.push_back(matched_eleTrkPt_tmp.at(j));

                   // diphoton inv mass
                   TLorentzVector v1, v2, vSum;
                   v1.SetPtEtaPhiE( phoEt->at(i), phoEta->at(i),
                           phoPhi->at(i), phoE->at(i));
                   v2.SetPtEtaPhiE( phoEt->at(j), phoEta->at(j),
                           phoPhi->at(j), phoE->at(j));
                   vSum = v1+v2;

                   vSum_gg_M.push_back(vSum.M());
                   vSum_gg_Energy.push_back(vSum.Energy());
                   vSum_gg_Pt.push_back(vSum.Pt());
                   vSum_gg_Eta.push_back(vSum.Eta());
                   vSum_gg_Phi.push_back(vSum.Phi());

                   if (matched_eleIndex_tmp.at(i) > -1 && matched_eleIndex_tmp.at(j) > -1) {
                       // diphoton inv mass with elePt
                       v1.SetPtEtaPhiM( elePt->at(matched_eleIndex_tmp.at(i)), eleEta->at(matched_eleIndex_tmp.at(i)),
                               elePhi->at(matched_eleIndex_tmp.at(i)), eleMass);
                       v2.SetPtEtaPhiM( elePt->at(matched_eleIndex_tmp.at(j)), eleEta->at(matched_eleIndex_tmp.at(j)),
                               elePhi->at(matched_eleIndex_tmp.at(j)), eleMass);
                       vSum = v1+v2;

                       vSum_gg_M_elePt.push_back(vSum.M());

                       // diphoton inv mass with eleTrkPt
                       v1.SetPtEtaPhiM( eleTrkPt->at(matched_eleIndex_tmp.at(i)), eleEta->at(matched_eleIndex_tmp.at(i)),
                               elePhi->at(matched_eleIndex_tmp.at(i)), eleMass);
                       v2.SetPtEtaPhiM( eleTrkPt->at(matched_eleIndex_tmp.at(j)), eleEta->at(matched_eleIndex_tmp.at(j)),
                               elePhi->at(matched_eleIndex_tmp.at(j)), eleMass);
                       vSum = v1+v2;

                       vSum_gg_M_eleTrkPt.push_back(vSum.M());
                   }
                   else {
                       vSum_gg_M_elePt.push_back(0);
                       vSum_gg_M_eleTrkPt.push_back(0);
                   }
               }
           }

           // dielectron block
           eleCharge_1.clear();
           elePt_1.clear();
           eleEta_1.clear();
           elePhi_1.clear();
           eleHoverE_1.clear();
           eleSigmaIEtaIEta_1.clear();
           eleSigmaIEtaIEta_2012_1.clear();
           eleSigmaIPhiIPhi_1.clear();
           eleEoverPInv_1.clear();
           eledEtaAtVtx_1.clear();
           eledPhiAtVtx_1.clear();
           eleD0_1.clear();
           eleDz_1.clear();

           eleCharge_2.clear();
           elePt_2.clear();
           eleEta_2.clear();
           elePhi_2.clear();
           eleHoverE_2.clear();
           eleSigmaIEtaIEta_2.clear();
           eleSigmaIEtaIEta_2012_2.clear();
           eleSigmaIPhiIPhi_2.clear();
           eleEoverPInv_2.clear();
           eledEtaAtVtx_2.clear();
           eledPhiAtVtx_2.clear();
           eleD0_2.clear();
           eleDz_2.clear();

           vSum_ee_M.clear();
           vSum_ee_Energy.clear();
           vSum_ee_Pt.clear();
           vSum_ee_Eta.clear();
           vSum_ee_Phi.clear();

           for(int i=0; i<nEle; ++i)
           {
               for (int j=i+1; j<nEle; ++j)
               {
                   eleCharge_1.push_back(eleCharge->at(i));
                   elePt_1.push_back(elePt->at(i));
                   eleEta_1.push_back(eleEta->at(i));
                   elePhi_1.push_back(elePhi->at(i));
                   eleHoverE_1.push_back(eleHoverE->at(i));
                   eleSigmaIEtaIEta_1.push_back(eleSigmaIEtaIEta->at(i));
                   eleSigmaIEtaIEta_2012_1.push_back(eleSigmaIEtaIEta_2012->at(j));
                   eleSigmaIPhiIPhi_1.push_back(eleSigmaIPhiIPhi->at(i));
                   eleEoverPInv_1.push_back(eleEoverPInv->at(i));
                   eledEtaAtVtx_1.push_back(eledEtaAtVtx->at(i));
                   eledPhiAtVtx_1.push_back(eledPhiAtVtx->at(i));
                   eleD0_1.push_back(eleD0->at(i));
                   eleDz_1.push_back(eleDz->at(i));

                   eleCharge_2.push_back(eleCharge->at(j));
                   elePt_2.push_back(elePt->at(j));
                   eleEta_2.push_back(eleEta->at(j));
                   elePhi_2.push_back(elePhi->at(j));
                   eleHoverE_2.push_back(eleHoverE->at(j));
                   eleSigmaIEtaIEta_2.push_back(eleSigmaIEtaIEta->at(j));
                   eleSigmaIEtaIEta_2012_2.push_back(eleSigmaIEtaIEta_2012->at(j));
                   eleSigmaIPhiIPhi_2.push_back(eleSigmaIPhiIPhi->at(j));
                   eleEoverPInv_2.push_back(eleEoverPInv->at(j));
                   eledEtaAtVtx_2.push_back(eledEtaAtVtx->at(j));
                   eledPhiAtVtx_2.push_back(eledPhiAtVtx->at(j));
                   eleD0_2.push_back(eleD0->at(j));
                   eleDz_2.push_back(eleDz->at(j));

                   TLorentzVector v1, v2, vSum;
                   v1.SetPtEtaPhiM( elePt->at(i), eleEta->at(i),
                           elePhi->at(i), eleMass);
                   v2.SetPtEtaPhiM( elePt->at(j), eleEta->at(j),
                           elePhi->at(j), eleMass);
                   vSum = v1+v2;

                   vSum_ee_M.push_back(vSum.M());
                   vSum_ee_Energy.push_back(vSum.Energy());
                   vSum_ee_Pt.push_back(vSum.Pt());
                   vSum_ee_Eta.push_back(vSum.Eta());
                   vSum_ee_Phi.push_back(vSum.Phi());
               }
           }

           // dimuon block
           muCharge_1.clear();
           muPt_1.clear();
           muEta_1.clear();
           muPhi_1.clear();
           muIsGood_1.clear();
           muD0_1.clear();
           muDz_1.clear();
           muChi2NDF_1.clear();
           muInnerD0_1.clear();
           muInnerDz_1.clear();
           muIsoTrk_1.clear();
           muPFChIso_1.clear();
           muPFPhoIso_1.clear();
           muPFNeuIso_1.clear();
           muPFPUIso_1.clear();

           muCharge_2.clear();
           muPt_2.clear();
           muEta_2.clear();
           muPhi_2.clear();
           muIsGood_2.clear();
           muD0_2.clear();
           muDz_2.clear();
           muChi2NDF_2.clear();
           muInnerD0_2.clear();
           muInnerDz_2.clear();
           muIsoTrk_2.clear();
           muPFChIso_2.clear();
           muPFPhoIso_2.clear();
           muPFNeuIso_2.clear();
           muPFPUIso_2.clear();

           vSum_mumu_M.clear();
           vSum_mumu_Energy.clear();
           vSum_mumu_Pt.clear();
           vSum_mumu_Eta.clear();
           vSum_mumu_Phi.clear();

           for(int i=0; i<nMu; ++i)
             {
                 for (int j=i+1; j<nMu; ++j)
                 {
                     muCharge_1.push_back(muCharge->at(i));
                     muPt_1.push_back(muPt->at(i));
                     muEta_1.push_back(muEta->at(i));
                     muPhi_1.push_back(muPhi->at(i));
                     muIsGood_1.push_back(muIsGood->at(i));
                     muD0_1.push_back(muD0->at(i));
                     muDz_1.push_back(muDz->at(j));
                     muChi2NDF_1.push_back(muChi2NDF->at(i));
                     muInnerD0_1.push_back(muInnerD0->at(i));
                     muInnerDz_1.push_back(muInnerDz->at(i));
                     muIsoTrk_1.push_back(muIsoTrk->at(i));
                     muPFChIso_1.push_back(muPFChIso->at(i));
                     muPFPhoIso_1.push_back(muPFPhoIso->at(i));
                     muPFNeuIso_1.push_back(muPFNeuIso->at(i));
                     muPFPUIso_1.push_back(muPFPUIso->at(i));

                     muCharge_2.push_back(muCharge->at(j));
                     muPt_2.push_back(muPt->at(j));
                     muEta_2.push_back(muEta->at(j));
                     muPhi_2.push_back(muPhi->at(j));
                     muIsGood_2.push_back(muIsGood->at(j));
                     muD0_2.push_back(muD0->at(j));
                     muDz_2.push_back(muDz->at(j));
                     muChi2NDF_2.push_back(muChi2NDF->at(j));
                     muInnerD0_2.push_back(muInnerD0->at(j));
                     muInnerDz_2.push_back(muInnerDz->at(j));
                     muIsoTrk_2.push_back(muIsoTrk->at(j));
                     muPFChIso_2.push_back(muPFChIso->at(j));
                     muPFPhoIso_2.push_back(muPFPhoIso->at(j));
                     muPFNeuIso_2.push_back(muPFNeuIso->at(j));
                     muPFPUIso_2.push_back(muPFPUIso->at(j));

                     TLorentzVector v1, v2, vSum;
                     v1.SetPtEtaPhiM( muPt->at(i), muEta->at(i),
                             muPhi->at(i), muMass);
                     v2.SetPtEtaPhiM( muPt->at(j), muEta->at(j),
                             muPhi->at(j), muMass);
                     vSum = v1+v2;

                     vSum_mumu_M.push_back(vSum.M());
                     vSum_mumu_Energy.push_back(vSum.Energy());
                     vSum_mumu_Pt.push_back(vSum.Pt());
                     vSum_mumu_Eta.push_back(vSum.Eta());
                     vSum_mumu_Phi.push_back(vSum.Phi());
                 }
             }

           outputTreeHLT->Fill();
           outputTreePho->Fill();
           outputTreeJet->Fill();
           outputTreeHiEvt->Fill();
           outputTreeSkim->Fill();
           
           outputTreeGammaJet->Fill();
           outputTreediPho->Fill();
           outputTreediEle->Fill();
           outputTreediMu->Fill();
       }
       std::cout<<  "Loop ENDED : ggHiNtuplizer/EventTree" <<std::endl;
       std::cout << "entries            = " << entries << std::endl;
       std::cout << "duplicateEntries   = " << duplicateEntries << std::endl;
       std::cout << "entriesAnalyzed    = " << entriesAnalyzed << std::endl;
       std::cout << "outputTreeHLT->GetEntries()   = " << outputTreeHLT->GetEntries() << std::endl;
       std::cout << "outputTreePho->GetEntries()      = " << outputTreePho->GetEntries() << std::endl;
       std::cout << "outputTreeJet->GetEntries()   = " << outputTreeJet->GetEntries() << std::endl;
       std::cout << "outputTreeHiEvt->GetEntries() = " << outputTreeHiEvt->GetEntries() << std::endl;
       std::cout << "outputTreeSkim->GetEntries()  = " << outputTreeSkim->GetEntries() << std::endl;
    
       std::cout << "outputTreeGammaJet->GetEntries() = " << outputTreeGammaJet->GetEntries() << std::endl;
       std::cout << "outputTreediPho->GetEntries()    = " << outputTreediPho->GetEntries() << std::endl;
       std::cout << "outputTreediEle->GetEntries()    = " << outputTreediEle->GetEntries() << std::endl;
       std::cout << "outputTreediMu->GetEntries()     = " << outputTreediMu->GetEntries() << std::endl;
  
       output->Write();
       output->Close();
}

int main(int argc, char** argv)
{
    if (argc == 3)    analyzeGammaJet_Data_HI_TChain(argv[1], argv[2]);
    if (argc == 2)    analyzeGammaJet_Data_HI_TChain(argv[1]);
    if (argc == 1)    return 1;
    return 0;
}


Double_t getDR( Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2){
  Double_t theDphi = getDPHI( phi1, phi2);
  Double_t theDeta = eta1 - eta2;
  return TMath::Sqrt ( theDphi*theDphi + theDeta*theDeta);
}

Double_t getDPHI( Double_t phi1, Double_t phi2) {
  Double_t dphi = phi1 - phi2;

  if ( dphi > 3.141592653589 )
    dphi = dphi - 2. * 3.141592653589;
  if ( dphi <= -3.141592653589 )
    dphi = dphi + 2. * 3.141592653589;

  if ( TMath::Abs(dphi) > 3.141592653589 ) {
    cout << " commonUtility::getDPHI error!!! dphi is bigger than 3.141592653589 " << endl;
  }

  return dphi;
}

Double_t getDETA(Double_t eta1, Double_t eta2){
    return eta1 - eta2;
}
