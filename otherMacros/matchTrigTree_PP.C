#include <TApplication.h>
#include "EventMatchingCMS.h"
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TString.h>
#include <TGraphAsymmErrors.h>
#include <TH1D.h>
#include <iostream>
#include <iomanip>
#include <TStyle.h>
#include "TH1F.h"
#include "TMath.h"
#include "TDirectory.h"

#include <string>
#include <vector>
#include <fstream>

#include "TLorentzVector.h"

#include "getLogBins.h"


const TString AnaPu3CaloTreename = "akPu3CaloJetAnalyzer/t";
const TString AnaPu4CaloTreename = "ak4CaloJetAnalyzer/t";
const TString AnaPu3PFTreename = "akPu3PFJetAnalyzer/t";
const TString AnaPu4PFTreename = "ak4PFJetAnalyzer/t";
const TString AnaVs3CaloTreename = "akVs3CaloJetAnalyzer/t";
const TString AnaVs4CaloTreename = "akVs4CaloJetAnalyzer/t";
const TString AnaPhotonTreename = "ggHiNtuplizer/EventTree";
const TString AnaHITreename = "hiEvtAnalyzer/HiTree";
const TString AnaSkimTreename = "skimanalysis/HltTree";
const TString AnaTrkTreename = "anaTrack/ppTrack";
const TString AnaGenTreename = "HiGenParticleAna/hi";

const TString HLTFilename = "openHLT_20150508_HIMinBias502_740F.root";

int matchTrigTree_HI(const std::string inHLTFile, const std::string inForestFile, const std::string inTrigFileName, const std::string inOptMaxName = "")
{

  std::string buffer;
  std::vector<std::string> listOfTrig;
  std::vector<Int_t> listOfThresh;
  Int_t nTrigTypeTemp = 0;
  int nLines = 0;
  std::ifstream* inTrigFile = new std::ifstream(inTrigFileName.data());
  std::cout << inTrigFileName << std::endl;

  if(!inTrigFile->is_open()){
    std::cout << "Error opening file. Exiting." <<std::endl;
    return 1;
  }
  else{
    std::cout << "Gettng Trig List" << std::endl;
    while(true){
      *inTrigFile >> buffer;
      if(inTrigFile->eof()) break;
      if(std::string::npos== buffer.find("HLT") && std::string::npos== buffer.find("L1")) nTrigTypeTemp++;
      listOfTrig.push_back(buffer);
      nLines++;
    }
  }
  delete inTrigFile;
  std::cout << "Trigger List Loaded" << std::endl;

  const Int_t nTrigType = nTrigTypeTemp;
  std::string trigType[nTrigType];
  Int_t trigTypeCount[nTrigType];

  std::cout << "nTrigType: " << nTrigType << std::endl;

  nLines = 0;

  for(Int_t iter = 0; iter < (Int_t)listOfTrig.size(); iter++){
    if(std::string::npos == listOfTrig[iter].find("HLT") && std::string::npos == listOfTrig[iter].find("L1")){

      std::size_t strIndex = 0;
      while(true){
        strIndex = listOfTrig[iter].find(",");
        if(strIndex == std::string::npos) break;
        listOfTrig[iter].replace(strIndex, std::string::npos, "");
      }

      trigType[nLines] = listOfTrig[iter];
      trigTypeCount[nLines] = 0;
      nLines++;
      listOfThresh.push_back(0);
    }
    else{
      trigTypeCount[nLines-1]++;

      std::size_t strIndex = 0;
      std::size_t strIndex2 = 0;
      while(true){
	strIndex = listOfTrig[iter].find(",");
	if(strIndex == std::string::npos) break;

        strIndex2 = listOfTrig[iter].substr(strIndex+1).find(",");
        listOfThresh.push_back(std::stoi(listOfTrig[iter].substr(strIndex+1, strIndex2)));
	listOfTrig[iter].replace(strIndex, std::string::npos, "");
      }      
    }
  }

  TFile *HLTFile_p =  new TFile(inHLTFile.c_str(), "READ");
  TTree *HLTTree = (TTree*)HLTFile_p->Get("HltTree");

  ULong64_t hlt_event;
  Int_t hlt_run, hlt_lumi;

  Int_t maxNTrig = -1;
  for(Int_t iter = 0; iter < nTrigType; iter++){
    if(trigTypeCount[iter] > maxNTrig) maxNTrig = trigTypeCount[iter];
  }
  const Int_t maxNTrig2 = maxNTrig;

  std::string trigName[nTrigType][maxNTrig2];
  Int_t trigThresh[nTrigType][maxNTrig2];
  Int_t trigVal[nTrigType][maxNTrig2];
  Int_t trigPrescale[nTrigType][maxNTrig2];
  Int_t nTrigFire[nTrigType][maxNTrig2];
  Float_t highPtMiss[nTrigType][maxNTrig2];
  Float_t highPtMiss_Eta[nTrigType][maxNTrig2];
  Float_t highPtMiss_Phi[nTrigType][maxNTrig2];
  Int_t highPtEventNum[nTrigType][maxNTrig2];
  Int_t highPtLumiNum[nTrigType][maxNTrig2];

  Float_t nZTrigFireNum[nTrigType][maxNTrig2];
  Float_t nZTrigFireDenom[nTrigType][maxNTrig2];

  nLines = 0;
  Int_t tempPosIter = 0;

  for(Int_t iter = 1; iter < (Int_t)(listOfTrig.size()); iter++){
    if(std::string::npos== listOfTrig[iter].find("HLT") && std::string::npos== listOfTrig[iter].find("L1")){
      nLines++;
      tempPosIter = 0;
    }
    else{
      trigName[nLines][tempPosIter] = listOfTrig[iter];
      trigThresh[nLines][tempPosIter] = listOfThresh[iter];
      trigVal[nLines][tempPosIter] = 0;
      trigPrescale[nLines][tempPosIter] = 0;
      nTrigFire[nLines][tempPosIter] = 0;
      highPtMiss[nLines][tempPosIter] = 0;
      highPtMiss_Eta[nLines][tempPosIter] = 0;
      highPtMiss_Phi[nLines][tempPosIter] = 0;
      highPtEventNum[nLines][tempPosIter] = 0;
      highPtLumiNum[nLines][tempPosIter] = 0;
      nZTrigFireNum[nLines][tempPosIter] = 0;
      nZTrigFireDenom[nLines][tempPosIter] = 0;
      tempPosIter++;
    }
  }

  HLTTree->SetBranchStatus("*", 0);

  HLTTree->SetBranchStatus("Event", 1);
  HLTTree->SetBranchStatus("Run", 1);
  HLTTree->SetBranchStatus("LumiBlock", 1);

  HLTTree->SetBranchAddress("Event", &hlt_event);
  HLTTree->SetBranchAddress("Run", &hlt_run);
  HLTTree->SetBranchAddress("LumiBlock", &hlt_lumi);

  //EDIT HERE FOR TRIGGER ANDING
  //prev -8
  for(Int_t iter = 0; iter < nTrigType; iter++){
    //    if(iter == 12 || iter == 13 || iter == 14) continue;
    for(Int_t iter2 = 0; iter2 < trigTypeCount[iter]; iter2++){
      HLTTree->SetBranchStatus(trigName[iter][iter2].c_str(), 1);
      HLTTree->SetBranchStatus(Form("%s_Prescl", trigName[iter][iter2].c_str()), 1);
      HLTTree->SetBranchAddress(trigName[iter][iter2].c_str(), &(trigVal[iter][iter2]));
      HLTTree->SetBranchAddress(Form("%s_Prescl", trigName[iter][iter2].c_str()), &(trigPrescale[iter][iter2]));
      std::cout << trigName[iter][iter2].c_str() << std::endl;
    }
  }

  TChain *AnaPu3CaloTree = new TChain(AnaPu3CaloTreename); 
  TChain *AnaPu4CaloTree = new TChain(AnaPu4CaloTreename); 
  TChain *AnaPu3PFTree = new TChain(AnaPu3PFTreename); 
  TChain *AnaPu4PFTree = new TChain(AnaPu4PFTreename); 
  TChain *AnaVs3CaloTree = new TChain(AnaVs3CaloTreename); 
  TChain *AnaVs4CaloTree = new TChain(AnaVs4CaloTreename); 
  TChain *AnaPhotonTree = new TChain(AnaPhotonTreename); 
  TChain *AnaHITree = new TChain(AnaHITreename); 
  TChain *AnaSkimTree = new TChain(AnaSkimTreename);
  TChain *AnaTrkTree = new TChain(AnaTrkTreename);
  TChain *AnaGenTree = new TChain(AnaGenTreename);
  
  for(Int_t iter = 0; iter < 5; iter++){
    //    if(iter == 1 || iter == 2 || iter == 3) continue;
   

    AnaPu3CaloTree->Add(Form("/mnt/hadoop/cms/store/user/dgulhan/mergedForest/PYTHIA_QCD_TuneCUET8M1_5020GeV_7_4_15/HiForest_PYTHIA_QCD_50_TuneCUETP8M1_5020GeV_7415_merged_forest_%d_0.root", iter));
    AnaPu4CaloTree->Add(Form("/mnt/hadoop/cms/store/user/dgulhan/mergedForest/PYTHIA_QCD_TuneCUET8M1_5020GeV_7_4_15/HiForest_PYTHIA_QCD_50_TuneCUETP8M1_5020GeV_7415_merged_forest_%d_0.root", iter));
    AnaPu3PFTree->Add(Form("/mnt/hadoop/cms/store/user/dgulhan/mergedForest/PYTHIA_QCD_TuneCUET8M1_5020GeV_7_4_15/HiForest_PYTHIA_QCD_50_TuneCUETP8M1_5020GeV_7415_merged_forest_%d_0.root", iter));
    AnaPu4PFTree->Add(Form("/mnt/hadoop/cms/store/user/dgulhan/mergedForest/PYTHIA_QCD_TuneCUET8M1_5020GeV_7_4_15/HiForest_PYTHIA_QCD_50_TuneCUETP8M1_5020GeV_7415_merged_forest_%d_0.root", iter));
    AnaVs3CaloTree->Add(Form("/mnt/hadoop/cms/store/user/dgulhan/mergedForest/PYTHIA_QCD_TuneCUET8M1_5020GeV_7_4_15/HiForest_PYTHIA_QCD_50_TuneCUETP8M1_5020GeV_7415_merged_forest_%d_0.root", iter));
    AnaVs4CaloTree->Add(Form("/mnt/hadoop/cms/store/user/dgulhan/mergedForest/PYTHIA_QCD_TuneCUET8M1_5020GeV_7_4_15/HiForest_PYTHIA_QCD_50_TuneCUETP8M1_5020GeV_7415_merged_forest_%d_0.root", iter));
    AnaPhotonTree->Add(Form("/mnt/hadoop/cms/store/user/dgulhan/mergedForest/PYTHIA_QCD_TuneCUET8M1_5020GeV_7_4_15/HiForest_PYTHIA_QCD_50_TuneCUETP8M1_5020GeV_7415_merged_forest_%d_0.root", iter));
    AnaHITree->Add(Form("/mnt/hadoop/cms/store/user/dgulhan/mergedForest/PYTHIA_QCD_TuneCUET8M1_5020GeV_7_4_15/HiForest_PYTHIA_QCD_50_TuneCUETP8M1_5020GeV_7415_merged_forest_%d_0.root", iter));
    AnaSkimTree->Add(Form("/mnt/hadoop/cms/store/user/dgulhan/mergedForest/PYTHIA_QCD_TuneCUET8M1_5020GeV_7_4_15/HiForest_PYTHIA_QCD_50_TuneCUETP8M1_5020GeV_7415_merged_forest_%d_0.root", iter));
    AnaTrkTree->Add(Form("/mnt/hadoop/cms/store/user/dgulhan/mergedForest/PYTHIA_QCD_TuneCUET8M1_5020GeV_7_4_15/HiForest_PYTHIA_QCD_50_TuneCUETP8M1_5020GeV_7415_merged_forest_%d_0.root", iter));
    AnaGenTree->Add(Form("/mnt/hadoop/cms/store/user/dgulhan/mergedForest/PYTHIA_QCD_TuneCUET8M1_5020GeV_7_4_15/HiForest_PYTHIA_QCD_50_TuneCUETP8M1_5020GeV_7415_merged_forest_%d_0.root", iter));
  }
  /* 
  const std::string totemString[19] = {"/mnt/hadoop/cms/store/user/cmcginn/hltJet/HI/20151105_TOTEMForest_753p1_3p0Thresh_GTv6_GRun_HIRCT_BeamSpot_2/HiForest_10_2_vHu.root",
"/mnt/hadoop/cms/store/user/cmcginn/hltJet/HI/20151105_TOTEMForest_753p1_3p0Thresh_GTv6_GRun_HIRCT_BeamSpot_2/HiForest_11_1_CCD.root",
"/mnt/hadoop/cms/store/user/cmcginn/hltJet/HI/20151105_TOTEMForest_753p1_3p0Thresh_GTv6_GRun_HIRCT_BeamSpot_2/HiForest_12_1_wKm.root",
"/mnt/hadoop/cms/store/user/cmcginn/hltJet/HI/20151105_TOTEMForest_753p1_3p0Thresh_GTv6_GRun_HIRCT_BeamSpot_2/HiForest_13_1_IE9.root",
"/mnt/hadoop/cms/store/user/cmcginn/hltJet/HI/20151105_TOTEMForest_753p1_3p0Thresh_GTv6_GRun_HIRCT_BeamSpot_2/HiForest_14_1_UKt.root",
"/mnt/hadoop/cms/store/user/cmcginn/hltJet/HI/20151105_TOTEMForest_753p1_3p0Thresh_GTv6_GRun_HIRCT_BeamSpot_2/HiForest_15_1_UR3.root",
"/mnt/hadoop/cms/store/user/cmcginn/hltJet/HI/20151105_TOTEMForest_753p1_3p0Thresh_GTv6_GRun_HIRCT_BeamSpot_2/HiForest_16_1_qAZ.root",
"/mnt/hadoop/cms/store/user/cmcginn/hltJet/HI/20151105_TOTEMForest_753p1_3p0Thresh_GTv6_GRun_HIRCT_BeamSpot_2/HiForest_17_1_6Dn.root",
"/mnt/hadoop/cms/store/user/cmcginn/hltJet/HI/20151105_TOTEMForest_753p1_3p0Thresh_GTv6_GRun_HIRCT_BeamSpot_2/HiForest_18_1_viw.root",
"/mnt/hadoop/cms/store/user/cmcginn/hltJet/HI/20151105_TOTEMForest_753p1_3p0Thresh_GTv6_GRun_HIRCT_BeamSpot_2/HiForest_19_1_7Qu.root",
"/mnt/hadoop/cms/store/user/cmcginn/hltJet/HI/20151105_TOTEMForest_753p1_3p0Thresh_GTv6_GRun_HIRCT_BeamSpot_2/HiForest_1_1_59n.root",
"/mnt/hadoop/cms/store/user/cmcginn/hltJet/HI/20151105_TOTEMForest_753p1_3p0Thresh_GTv6_GRun_HIRCT_BeamSpot_2/HiForest_2_1_6WX.root",
"/mnt/hadoop/cms/store/user/cmcginn/hltJet/HI/20151105_TOTEMForest_753p1_3p0Thresh_GTv6_GRun_HIRCT_BeamSpot_2/HiForest_3_1_NPE.root",
"/mnt/hadoop/cms/store/user/cmcginn/hltJet/HI/20151105_TOTEMForest_753p1_3p0Thresh_GTv6_GRun_HIRCT_BeamSpot_2/HiForest_4_2_iI9.root",
"/mnt/hadoop/cms/store/user/cmcginn/hltJet/HI/20151105_TOTEMForest_753p1_3p0Thresh_GTv6_GRun_HIRCT_BeamSpot_2/HiForest_5_1_FK0.root",
"/mnt/hadoop/cms/store/user/cmcginn/hltJet/HI/20151105_TOTEMForest_753p1_3p0Thresh_GTv6_GRun_HIRCT_BeamSpot_2/HiForest_6_1_jVN.root",
"/mnt/hadoop/cms/store/user/cmcginn/hltJet/HI/20151105_TOTEMForest_753p1_3p0Thresh_GTv6_GRun_HIRCT_BeamSpot_2/HiForest_7_1_6x9.root",
"/mnt/hadoop/cms/store/user/cmcginn/hltJet/HI/20151105_TOTEMForest_753p1_3p0Thresh_GTv6_GRun_HIRCT_BeamSpot_2/HiForest_8_1_rwa.root",
"/mnt/hadoop/cms/store/user/cmcginn/hltJet/HI/20151105_TOTEMForest_753p1_3p0Thresh_GTv6_GRun_HIRCT_BeamSpot_2/HiForest_9_1_mpS.root"};

   
  const std::string totemString[11] = {"HiForest_TOTEM_182BB6D4-FD74-E511-A418-02163E014655.root",
"HiForest_TOTEM_1CB375D0-0175-E511-8B65-02163E0119C3.root",
"HiForest_TOTEM_34C660CE-0175-E511-A6CF-02163E01343C.root",
"HiForest_TOTEM_46A2E8D4-FD74-E511-9A5F-02163E011B4D.root",
"HiForest_TOTEM_505062D5-FD74-E511-BADA-02163E011C72.root",
"HiForest_TOTEM_5A9E77D1-FD74-E511-A76C-02163E014223.root",
"HiForest_TOTEM_925C78E9-FD74-E511-B56A-02163E01437D.root",
"HiForest_TOTEM_94D211E2-FD74-E511-B8D2-02163E01434C.root",
"HiForest_TOTEM_AA2030D8-FD74-E511-BA1B-02163E014472.root",
"HiForest_TOTEM_BE2401D1-FD74-E511-907D-02163E0139C8.root",
"HiForest_TOTEM_EEF3E6FB-0175-E511-BED5-02163E0143E2.root"};
   
  for(Int_t iter =0 ; iter < 19; iter++){
    AnaPu3CaloTree->Add( totemString[iter].c_str());
    AnaPu4CaloTree->Add(totemString[iter].c_str());
    AnaPu3PFTree->Add(totemString[iter].c_str());
    AnaPu4PFTree->Add(totemString[iter].c_str());
    AnaVs3CaloTree->Add(totemString[iter].c_str());
    AnaVs4CaloTree->Add(totemString[iter].c_str());
    AnaPhotonTree->Add(totemString[iter].c_str());
    AnaHITree->Add(totemString[iter].c_str());
    AnaSkimTree->Add(totemString[iter].c_str());
    AnaTrkTree->Add(totemString[iter].c_str());
    AnaGenTree->Add(totemString[iter].c_str());
  }
  */
  Int_t ana_lumi, ana_run;
  Long64_t ana_event;
  Int_t hiBin;
  
  Int_t pcollisionEventSelection;

  const Int_t maxTrk = 50000;
  Int_t nTrk;
  Bool_t trkFake[maxTrk], highPurity[maxTrk];
  Float_t trkPt[maxTrk], trkPhi[maxTrk], trkEta[maxTrk];
  Float_t trkDz1[maxTrk], trkDzError1[maxTrk], trkDxy1[maxTrk], trkDxyError1[maxTrk], trkPtError[maxTrk];

  const Int_t maxGen = 100000;
  Int_t nGen;
  //  Float_t genPt[maxGen];
  //  Float_t genPhi[maxGen];
  //  Float_t genEta[maxGen];
  std::vector<float>* genPt=0;
  std::vector<float>* genEta=0;
  std::vector<float>* genPhi=0;
  std::vector<int>* genPDG=0;
  //  Int_t genPDG[maxGen];
  Int_t genStatus[maxGen];

  const Int_t maxJt = 500;
  Int_t nPu3Caloref;
  Float_t jtPu3Calopt[maxJt], jtPu3Caloeta[maxJt], jtPu3Calophi[maxJt];
  Float_t refPu3Calopt[maxJt], refPu3Caloeta[maxJt], refPu3Calophi[maxJt];
  Int_t nPu3Calogen;
  Float_t genPu3Calopt[maxJt], genPu3Caloeta[maxJt], genPu3Calophi[maxJt];

  Int_t nPu4Caloref;
  Float_t jtPu4Calopt[maxJt], jtPu4Caloeta[maxJt], jtPu4Calophi[maxJt], jtPu4Calorawpt[maxJt];
  Float_t refPu4Calopt[maxJt], refPu4Caloeta[maxJt], refPu4Calophi[maxJt];
  Int_t nPu4Calogen;
  Float_t genPu4Calopt[maxJt], genPu4Caloeta[maxJt], genPu4Calophi[maxJt];

  Int_t nPu3PFref;
  Float_t jtPu3PFpt[maxJt], jtPu3PFeta[maxJt], jtPu3PFphi[maxJt];
  Float_t refPu3PFpt[maxJt], refPu3PFeta[maxJt], refPu3PFphi[maxJt];
  Int_t nPu3PFgen;
  Float_t genPu3PFpt[maxJt], genPu3PFeta[maxJt], genPu3PFphi[maxJt];

  Int_t nPu4PFref;
  Float_t jtPu4PFpt[maxJt], jtPu4PFeta[maxJt], jtPu4PFphi[maxJt];
  Float_t refPu4PFpt[maxJt], refPu4PFeta[maxJt], refPu4PFphi[maxJt];
  Int_t nPu4PFgen;
  Float_t genPu4PFpt[maxJt], genPu4PFeta[maxJt], genPu4PFphi[maxJt];

  Int_t nVs3Caloref;
  Float_t jtVs3Calopt[maxJt], jtVs3Caloeta[maxJt], jtVs3Calophi[maxJt];
  Float_t refVs3Calopt[maxJt], refVs3Caloeta[maxJt], refVs3Calophi[maxJt];
  Int_t nVs3Calogen;
  Float_t genVs3Calopt[maxJt], genVs3Caloeta[maxJt], genVs3Calophi[maxJt];

  Int_t nVs4Caloref;
  Float_t jtVs4Calopt[maxJt], jtVs4Caloeta[maxJt], jtVs4Calophi[maxJt];
  Float_t refVs4Calopt[maxJt], refVs4Caloeta[maxJt], refVs4Calophi[maxJt];
  Int_t nVs4Calogen;
  Float_t genVs4Calopt[maxJt], genVs4Caloeta[maxJt], genVs4Calophi[maxJt];

  const Int_t maxGamma = 50;
  Int_t nPho;
  std::vector<float>* phoEt=0;
  std::vector<float>* phoEta=0;
  std::vector<float>* phoPhi=0;
  std::vector<float>* phoR9=0;
  std::vector<float>* phoSigmaIEtaIEta=0;
  std::vector<float>* phoHoverE=0;

  //  Float_t phoEt[maxGamma], phoEta[maxGamma], phoPhi[maxGamma];   //[nPhotons]
  Float_t seedTime[maxGamma], swissCrx[maxGamma], sigmaIphiIphi[maxGamma], sigmaIetaIeta[maxGamma];
  Float_t r9[maxGamma];
  Int_t isEle[maxGamma];

  
  AnaPu4CaloTree->SetBranchStatus("*", 0);
  AnaPu4CaloTree->SetBranchStatus("nref", 1);
  AnaPu4CaloTree->SetBranchStatus("jtpt", 1);
  AnaPu4CaloTree->SetBranchStatus("rawpt", 1);
  AnaPu4CaloTree->SetBranchStatus("jteta", 1);
  AnaPu4CaloTree->SetBranchStatus("jtphi", 1);
  AnaPu4CaloTree->SetBranchStatus("refpt", 1);
  AnaPu4CaloTree->SetBranchStatus("refeta", 1);
  AnaPu4CaloTree->SetBranchStatus("refphi", 1);
  AnaPu4CaloTree->SetBranchStatus("ngen", 1);
  AnaPu4CaloTree->SetBranchStatus("genpt", 1);
  AnaPu4CaloTree->SetBranchStatus("genphi", 1);
  AnaPu4CaloTree->SetBranchStatus("geneta", 1);
  AnaPu4CaloTree->SetBranchAddress("nref", &nPu4Caloref);
  AnaPu4CaloTree->SetBranchAddress("jtpt", jtPu4Calopt);
  AnaPu4CaloTree->SetBranchAddress("rawpt", jtPu4Calorawpt);
  AnaPu4CaloTree->SetBranchAddress("jteta", jtPu4Caloeta);
  AnaPu4CaloTree->SetBranchAddress("jtphi", jtPu4Calophi);
  AnaPu4CaloTree->SetBranchAddress("refpt", refPu4Calopt);
  AnaPu4CaloTree->SetBranchAddress("refeta", refPu4Caloeta);
  AnaPu4CaloTree->SetBranchAddress("refphi", refPu4Calophi);
  AnaPu4CaloTree->SetBranchAddress("ngen", &nPu4Calogen);
  AnaPu4CaloTree->SetBranchAddress("genpt", genPu4Calopt);
  AnaPu4CaloTree->SetBranchAddress("genphi", genPu4Calophi);
  AnaPu4CaloTree->SetBranchAddress("geneta", genPu4Caloeta);
  /*
  AnaPu4CaloTree->SetBranchStatus("evt", 1);
  AnaPu4CaloTree->SetBranchStatus("run", 1);
  AnaPu4CaloTree->SetBranchStatus("lumi", 1);
  AnaPu4CaloTree->SetBranchAddress("evt", &ana_event);
  AnaPu4CaloTree->SetBranchAddress("lumi", &ana_lumi);
  AnaPu4CaloTree->SetBranchAddress("run", &ana_run);
  */

  AnaPu4PFTree->SetBranchStatus("*", 0);
  AnaPu4PFTree->SetBranchStatus("nref", 1);
  AnaPu4PFTree->SetBranchStatus("jtpt", 1);
  AnaPu4PFTree->SetBranchStatus("jteta", 1);
  AnaPu4PFTree->SetBranchStatus("jtphi", 1);
  AnaPu4PFTree->SetBranchStatus("refpt", 1);
  AnaPu4PFTree->SetBranchStatus("refeta", 1);
  AnaPu4PFTree->SetBranchStatus("refphi", 1);
  AnaPu4PFTree->SetBranchStatus("ngen", 1);
  AnaPu4PFTree->SetBranchStatus("genpt", 1);
  AnaPu4PFTree->SetBranchStatus("genphi", 1);
  AnaPu4PFTree->SetBranchStatus("geneta", 1);
  AnaPu4PFTree->SetBranchAddress("nref", &nPu4PFref);
  AnaPu4PFTree->SetBranchAddress("jtpt", jtPu4PFpt);
  AnaPu4PFTree->SetBranchAddress("jteta", jtPu4PFeta);
  AnaPu4PFTree->SetBranchAddress("jtphi", jtPu4PFphi);
  AnaPu4PFTree->SetBranchAddress("refpt", refPu4PFpt);
  AnaPu4PFTree->SetBranchAddress("refeta", refPu4PFeta);
  AnaPu4PFTree->SetBranchAddress("refphi", refPu4PFphi);
  AnaPu4PFTree->SetBranchAddress("ngen", &nPu4PFgen);
  AnaPu4PFTree->SetBranchAddress("genpt", genPu4PFpt);
  AnaPu4PFTree->SetBranchAddress("genphi", genPu4PFphi);
  AnaPu4PFTree->SetBranchAddress("geneta", genPu4PFeta);

  std::cout << "A3" << std::endl;
  /*  
  AnaGenTree->SetBranchStatus("*", 0);
  AnaGenTree->SetBranchStatus("mult", 1);
  //  AnaGenTree->SetBranchStatus("sta", 1);
  AnaGenTree->SetBranchStatus("pdg", 1);
  AnaGenTree->SetBranchStatus("pt", 1);
  AnaGenTree->SetBranchStatus("phi", 1);
  AnaGenTree->SetBranchStatus("eta", 1);

  AnaGenTree->SetBranchAddress("mult", &nGen);
  //  AnaGenTree->SetBranchAddress("sta", genStatus);
  AnaGenTree->SetBranchAddress("pdg", &genPDG);
  AnaGenTree->SetBranchAddress("pt", &genPt);
  AnaGenTree->SetBranchAddress("phi", &genPhi);
  AnaGenTree->SetBranchAddress("eta", &genEta);
  */

  std::cout << "A4a" << std::endl;

  AnaPhotonTree->SetBranchStatus("*", 0);
  
  AnaPhotonTree->SetBranchStatus("event", 1);
  AnaPhotonTree->SetBranchStatus("run", 1);
  AnaPhotonTree->SetBranchStatus("lumis", 1);
  AnaPhotonTree->SetBranchAddress("event", &ana_event);
  AnaPhotonTree->SetBranchAddress("lumis", &ana_lumi);
  AnaPhotonTree->SetBranchAddress("run", &ana_run);
  
  
  AnaPhotonTree->SetBranchStatus("nPho", 1);
  AnaPhotonTree->SetBranchStatus("phoEt", 1);
  AnaPhotonTree->SetBranchStatus("phoEta", 1);
  AnaPhotonTree->SetBranchStatus("phoPhi", 1);
  AnaPhotonTree->SetBranchStatus("phoR9", 1);
  AnaPhotonTree->SetBranchStatus("phoHoverE", 1);
  AnaPhotonTree->SetBranchStatus("phoSigmaIEtaIEta", 1);
  /*
  AnaPhotonTree->SetBranchStatus("seedTime", 1);
  AnaPhotonTree->SetBranchStatus("swissCrx", 1);
  AnaPhotonTree->SetBranchStatus("sigmaIphiIphi", 1);
  AnaPhotonTree->SetBranchStatus("sigmaIetaIeta", 1);
  AnaPhotonTree->SetBranchStatus("isEle", 1);
  AnaPhotonTree->SetBranchStatus("r9", 1);


  std::cout << "A4b" << std::endl;
  */
  AnaPhotonTree->SetBranchAddress("nPho", &nPho);
  AnaPhotonTree->SetBranchAddress("phoEt", &phoEt);
  AnaPhotonTree->SetBranchAddress("phoEta", &phoEta);
  AnaPhotonTree->SetBranchAddress("phoPhi", &phoPhi);
  AnaPhotonTree->SetBranchAddress("phoR9", &phoR9);
  AnaPhotonTree->SetBranchAddress("phoHoverE", &phoHoverE);
  AnaPhotonTree->SetBranchAddress("phoSigmaIEtaIEta", &phoSigmaIEtaIEta);

  /*
  AnaPhotonTree->SetBranchAddress("seedTime", seedTime);
  AnaPhotonTree->SetBranchAddress("swissCrx", swissCrx);
  AnaPhotonTree->SetBranchAddress("sigmaIphiIphi", sigmaIphiIphi);
  AnaPhotonTree->SetBranchAddress("sigmaIetaIeta", sigmaIetaIeta);
  AnaPhotonTree->SetBranchAddress("isEle", isEle);
  AnaPhotonTree->SetBranchAddress("r9", r9);
  */
  std::cout << "A5" << std::endl;
  
  //FOR ADDITIONAL OFFLINE OBJECT MATCHING, EDIT HERE (1 of 2)


  //  const Int_t nPtBins[nTrigType] = {200, 200, 200, 200};
  //  const Int_t maxPt[nTrigType] = {200, 200, 200, 200};

  //  const Int_t nPtBins[nTrigType] = {160};
  //  const Int_t maxPt[nTrigType] = {160};

  //  const Int_t nPtBins[nTrigType] = {200, 120, 120, 120, 60, 60, 60, 30, 120, 60, 160, 160, 160, 160, 160, 160, 160, 160, 160, 200, 200, 200, 160, 30, 160, 160, 100};
  //  const Int_t maxPt[nTrigType] = {200, 120, 120, 120, 60, 60, 60, 30, 120, 60, 160, 160, 160, 160, 160, 160, 160, 160, 160, 100, 100, 100, 100, 100, 100, 100, 100, 100, 160, 160, 160, 200, 200, 200, 200, 200, 200, 200, 200, 160, 30, 160, 160, 100};

  const Int_t nPtBins[nTrigType] = {
    100,   //LEV1JET
    70,    //LEV1GAMMA
    160,   //R4JET
    160,   //R4PFJET
    100,   //GAMMAETA1P5
    100,   //GAMMAETA3P1
  };


  const Int_t maxPt[nTrigType] = {
    100,   //LEV1JET
    70,    //LEV1GAMMA
    160,   //R4JET
    160,   //R4PFJET
    100,   //GAMMAETA1P5
    100,   //GAMMAETA3P1
  };


  const Int_t nEtaBins = 50;
  TH1F *histsPt_p[nTrigType][maxNTrig2+1], *histsEta_p[nTrigType][maxNTrig2], *histsEtaTrig_p[nTrigType][maxNTrig2];
  TH1F* histsSpectPt_p[nTrigType];
  TH1F *histsSpectBoundPt_p[nTrigType][maxNTrig2+1];
  TH1F *histsTrkSpectBoundPt_p[nTrigType][maxNTrig2+1];
  TH1F *histsJetSpectBoundPt_p[nTrigType][maxNTrig2+1];

  const Int_t nTrkBins = 50;
  const Float_t trkLow = 0.5;
  const Float_t trkHi = 500;
  Float_t trkBins[nTrkBins+1];
  getLogBins(trkLow, trkHi, nTrkBins, trkBins);


  TH1D* ratesMatched_p[nTrigType];
  TH1D* ratesUnmatched_p[nTrigType];

  TH1F* maxMissedPt99_p[nTrigType];
  TH1F* maxMissedPt100_p[nTrigType];
  TH1F* prescale_p[nTrigType];

  TH1F* e1R9_p[nTrigType];
  TH1F* e2R9_p[nTrigType];

  TH1F* zPt_p = new TH1F("zPt_p", ";Z p_{T};Events", 101, -2.0, 200.0);
  TH1F* zMass_p = new TH1F("zMass_p", ";Z Mass;Events", 101, -2.0, 200.0);
  TH1F* e1Pt_p = new TH1F("e1Pt_p", ";e1 p_{T};Events", 101, -2.0, 200.0);
  TH1F* e2Pt_p = new TH1F("e2Pt_p", ";e2 p_{T};Events", 101, -2.0, 200.0);

  zMass_p->GetXaxis()->CenterTitle();
  zMass_p->GetYaxis()->CenterTitle();

  std::cout << "A7" << std::endl;

  for(Int_t iter = 0; iter < nTrigType; iter++){
    histsPt_p[iter][0] = new TH1F(Form("leading%s_pt", trigType[iter].c_str()), Form("leading%s_pt", trigType[iter].c_str()), nPtBins[iter], 0.0, maxPt[iter]);

    histsSpectPt_p[iter] = new TH1F(Form("spect%s_pt", trigType[iter].c_str()), Form("spect%s_pt", trigType[iter].c_str()), nPtBins[iter], 0.0, maxPt[iter]);

    Int_t max = trigThresh[iter][trigTypeCount[iter]-1];
    Int_t min = trigThresh[iter][0];
    Float_t minDiff = trigThresh[iter][1] - trigThresh[iter][0];

    for(Int_t iter2 = 1; iter2 < trigTypeCount[iter]-1; iter2++){
      if(trigThresh[iter][iter2+1] - trigThresh[iter][iter2] < minDiff) minDiff = trigThresh[iter][iter2+1] - trigThresh[iter][iter2];
    }

    minDiff /= 2;
    if(min < minDiff) minDiff = min;

    Int_t nBins = (max+minDiff - (min-minDiff))/(minDiff*2);
    if(trigTypeCount[iter] == 1){
      nBins = 1;
      minDiff = 20;
      if(min < minDiff) minDiff = min;
      min = trigThresh[iter][0] - minDiff;
      max = trigThresh[iter][0] + minDiff;
    }
    else{
      min -= minDiff;
      max += minDiff;
    }

    ratesMatched_p[iter] = new TH1D(Form("ratesMatched_%s_%d_%d", trigType[iter].c_str(), trigThresh[iter][0], trigThresh[iter][trigTypeCount[iter]-1]), Form("ratesMatched_%s_%d_%d", trigType[iter].c_str(), trigThresh[iter][0], trigThresh[iter][trigTypeCount[iter]-1]), nBins, min, max);
    ratesUnmatched_p[iter] = new TH1D(Form("ratesUnmatched_%s_%d_%d", trigType[iter].c_str(), trigThresh[iter][0], trigThresh[iter][trigTypeCount[iter]-1]), Form("ratesUnmatched_%s_%d_%d", trigType[iter].c_str(), trigThresh[iter][0], trigThresh[iter][trigTypeCount[iter]-1]), nBins, min, max);
    maxMissedPt99_p[iter] = new TH1F(Form("maxMissedPt99_%s_%d_%d", trigType[iter].c_str(), trigThresh[iter][0], trigThresh[iter][trigTypeCount[iter]-1]), Form("maxMissedPt99_%s_%d_%d", trigType[iter].c_str(), trigThresh[iter][0], trigThresh[iter][trigTypeCount[iter]-1]), nBins, min, max);
    maxMissedPt100_p[iter] = new TH1F(Form("maxMissedPt100_%s_%d_%d", trigType[iter].c_str(), trigThresh[iter][0], trigThresh[iter][trigTypeCount[iter]-1]), Form("maxMissedPt100_%s_%d_%d", trigType[iter].c_str(), trigThresh[iter][0], trigThresh[iter][trigTypeCount[iter]-1]), nBins, min, max);
    prescale_p[iter] = new TH1F(Form("prescale_%s_%d_%d", trigType[iter].c_str(), trigThresh[iter][0], trigThresh[iter][trigTypeCount[iter]-1]), Form("prescale_%s_%d_%d", trigType[iter].c_str(), trigThresh[iter][0], trigThresh[iter][trigTypeCount[iter]-1]), nBins, min, max);

    e1R9_p[iter] = new TH1F(Form("e1R9_%s_%d", trigType[iter].c_str(), trigThresh[iter][0]), Form("e1R9_%s_%d", trigType[iter].c_str(), trigThresh[iter][0]), 101, -0.05, 1.05);

    e2R9_p[iter] = new TH1F(Form("e2R9_%s_%d", trigType[iter].c_str(), trigThresh[iter][0]), Form("e2R9_%s_%d", trigType[iter].c_str(), trigThresh[iter][0]), 101, -0.05, 1.05);

    std::cout << "Min, max: " << min << ", " << max << std::endl;

    for(Int_t iter2 = 0; iter2 < trigTypeCount[iter]; iter2++){
      histsPt_p[iter][iter2+1] = (TH1F*)histsPt_p[iter][0]->Clone(Form("%s_%s_%d_pt", trigName[iter][iter2].c_str(), trigType[iter].c_str(), trigThresh[iter][iter2]));
      histsEta_p[iter][iter2] = new TH1F(Form("leading_%s_%d_eta", trigType[iter].c_str(), trigThresh[iter][iter2]), Form("leading_%s_%d_eta", trigType[iter].c_str(), trigThresh[iter][iter2]), nEtaBins, -5.0, 5.0);
      histsEtaTrig_p[iter][iter2] = new TH1F(Form("%s_%s_%d_eta", trigName[iter][iter2].c_str(), trigType[iter].c_str(), trigThresh[iter][iter2]), Form("%s_%s_%d_eta", trigName[iter][iter2].c_str(), trigType[iter].c_str(), trigThresh[iter][iter2]), nEtaBins, -5.0, 5.0);
      histsSpectBoundPt_p[iter][iter2] = new TH1F(Form("leadingPt_%s_%d", trigType[iter].c_str(), trigThresh[iter][iter2]), Form(";%s p_{T}^{reco};Events/GeV", trigType[iter].c_str()), nPtBins[iter], 0.0, maxPt[iter]);

      if(!strcmp("TRACK", trigType[iter].substr(0,5).c_str())){
	histsTrkSpectBoundPt_p[iter][iter2] = new TH1F(Form("trkPt_%s_%d", trigType[iter].c_str(), trigThresh[iter][iter2]), Form(";%s p_{T,trk}^{reco};#frac{1}{N_{evt}} N_{trk}", trigType[iter].c_str()), nTrkBins, trkBins);
	histsJetSpectBoundPt_p[iter][iter2] = new TH1F(Form("jetPt_%s_%d", trigType[iter].c_str(), trigThresh[iter][iter2]), Form(";%s p_{T,trk}^{reco};#frac{1}{N_{evt}} N_{trk}", trigType[iter].c_str()), nTrkBins, trkBins);
	/*
	histsTrkSpectBoundPt_p[iter][iter2] = new TH1F(Form("trkPt_%s_%d", trigType[iter].c_str(), trigThresh[iter][iter2]), Form(";%s p_{T,trk}^{reco};#frac{1}{N_{evt}} #frac{dN_{trk}}{dp_{T}}", trigType[iter].c_str()), nTrkBins, trkBins);
	histsJetSpectBoundPt_p[iter][iter2] = new TH1F(Form("jetPt_%s_%d", trigType[iter].c_str(), trigThresh[iter][iter2]), Form(";%s p_{T,trk}^{reco};#frac{1}{N_{evt}} #frac{dN_{trk}}{dp_{T}}", trigType[iter].c_str()), nTrkBins, trkBins);
	*/
      }
    }
  }
  

  //book histos

  std::cout << "Events in HLT file: " << HLTTree->GetEntries() << std::endl;
  std::cout << "Events in Ana file: " << AnaPu4CaloTree->GetEntries() << std::endl;

  //make map
  EventMatchingCMS *matcher = new EventMatchingCMS();

  for(Long64_t entry = 0; entry < HLTTree->GetEntries(); ++entry){

    HLTTree->GetEntry(entry);

    for(Int_t iter = 0; iter < nTrigType; iter++){
      for(Int_t iter2 = 0; iter2 < trigTypeCount[iter]; iter2++){
	if(trigVal[iter][iter2]) nTrigFire[iter][iter2]++;

	if(entry == 0) prescale_p[iter]->SetBinContent(prescale_p[iter]->FindBin(trigThresh[iter][iter2]), trigPrescale[iter][iter2]);
      }
    }
    
    matcher->addEvent(hlt_lumi, hlt_event, 0, entry);
  }

  // analysis loop
  std::cout << "COMMENCE LOOP" << std::endl;
  int matched = 0;
  for(Long64_t entry = 0; entry < AnaPu4CaloTree->GetEntries(); ++entry){
    if(entry%10000 == 0) std::cout << entry << std::endl;

    //    AnaTrkTree->GetEntry(entry);
    AnaPu4CaloTree->GetEntry(entry);
    AnaPu4PFTree->GetEntry(entry);
    AnaPhotonTree->GetEntry(entry);
    //    AnaGenTree->GetEntry(entry);
    
	long long hlt_entry = matcher->retrieveEvent(ana_lumi, ana_event, 0);
    if(hlt_entry == -1) continue;

    matched++;

    HLTTree->GetEntry(hlt_entry);

    Double_t maxTrkPt = -1;
    //    Double_t maxTrkEta = -100;
  
    for(int i = 0; i < nTrk; i++){
      if(fabs(trkEta[i]) > 2.4) continue;
      if(trkFake[i]) continue;

      
      if(!highPurity[i]) continue;
      if(TMath::Abs(trkDz1[i]/trkDzError1[i]) > 3) continue;
      if(TMath::Abs(trkDxy1[i]/trkDxyError1[i]) > 3) continue;
      if(trkPtError[i]/trkPt[i] > 0.1) continue;

      if(trkPt[i] > maxTrkPt){
	maxTrkPt = trkPt[i];
	//	maxTrkEta = trkEta[i];
      }
    }

    

    Double_t maxMuPt = -1;
    Double_t maxMuEta = -100;
    Double_t maxMuPhi = -100;
    /*
    for(int i = 0; i < nGen; i++){
      if(TMath::Abs(genPDG[i]) != 13) continue;
      if(fabs(genEta[i]) > 2.4) continue;
      if(genPt[i] > maxMuPt){
	maxMuPt = genPt[i];
	maxMuEta = genEta[i];
	maxMuPhi = genPhi[i];
      }
    }
    */
    /*
    Double_t maxDPt = -1;
    Double_t maxDEta = -100;

    for(int i = 0; i < nGen; i++){
      if(TMath::Abs(genPDG[i]) != 421) continue;
      if(fabs(genEta[i]) > 2.0) continue;
      if(genPt[i] > maxDPt){
	maxDPt = genPt[i];
	maxDEta = genEta[i];
      }
    }  
    */

    Double_t maxPu3CaloAnaPt = -1;
    Double_t maxPu3CaloAnaEta = -100;
    Double_t maxPu3CaloAnaPhi = -100;

    Double_t maxPu3CaloAnaPt_2p5 = -1;
    Double_t maxPu3CaloAnaEta_2p5 = -100;
    Double_t maxPu3CaloAnaPhi_2p5 = -100;
    
    for(int i = 0; i < nPu3Caloref; ++i){
      if(fabs(jtPu3Caloeta[i]) > 2.0) continue;
      if(jtPu3Calopt[i] > maxPu3CaloAnaPt){
	maxPu3CaloAnaPt = jtPu3Calopt[i];
	maxPu3CaloAnaEta = jtPu3Caloeta[i];
	maxPu3CaloAnaPhi = jtPu3Calophi[i];
      }
    }

    for(int i = 0; i < nPu3Caloref; ++i){
      if(fabs(jtPu3Caloeta[i]) > 2.4) continue;
      if(jtPu3Calopt[i] > maxPu3CaloAnaPt_2p5){
	maxPu3CaloAnaPt_2p5 = jtPu3Calopt[i];
	maxPu3CaloAnaEta_2p5 = jtPu3Caloeta[i];
	maxPu3CaloAnaPhi_2p5 = jtPu3Calophi[i];
      }
    }

    Double_t maxPu4CaloAnaPt_2p5 = -1;
    Double_t maxPu4CaloAnaEta_2p5 = -100;

    for(int i = 0; i < nPu4Caloref; ++i){
      if(fabs(jtPu4Caloeta[i]) > 2.4) continue;
      if(jtPu4Calopt[i] > maxPu4CaloAnaPt_2p5){
        maxPu4CaloAnaPt_2p5 = jtPu4Calopt[i];
        maxPu4CaloAnaEta_2p5 = jtPu4Caloeta[i];
      }
    }



    Double_t maxPu4CaloAnaPt = -1;
    Double_t maxPu4CaloAnaEta = -100;
    Double_t maxPu4CaloAnaPhi = -100;



    Double_t maxMidPu4CaloAnaPt = -1;
    //   Double_t maxMidPu4CaloAnaEta = -100;

    Double_t maxDiPu4CaloEta1p1AnaPt = -1;
    Double_t maxDiPu4CaloEta1p1AnaEta = -100;
    Double_t maxDiPu4CaloEta1p1AnaPhi = -100;

    Double_t twoDiPu4CaloEta1p1AnaPt = -1;
    //    Double_t twoDiPu4CaloAnaEta = -100;

    Double_t maxDiPu4CaloEta0p7AnaPt = -1;
    Double_t maxDiPu4CaloEta0p7AnaEta = -100;
    Double_t maxDiPu4CaloEta0p7AnaPhi = -100;

    Double_t twoDiPu4CaloEta0p7AnaPt = -1;
    //    Double_t twoDiPu4CaloAnaEta = -100;

    Double_t maxTriPu4CaloAnaPt = -1;
    Double_t maxTriPu4CaloAnaEta = -100;
    Double_t maxTriPu4CaloAnaPhi = -100;

    Double_t twoTriPu4CaloAnaPt = -1;
    //    Double_t twoTriPu4CaloAnaEta = -100;

    Double_t threeTriPu4CaloAnaPt = -1;
    //    Double_t threeTriPu4CaloAnaEta = -100;

    for(int i = 0; i < nPu4Caloref; ++i){
      if(fabs(jtPu4Caloeta[i]) > 2.0) continue;

      //      if(jtPu4Caloeta[i] < -1.52 && jtPu4Caloeta[i] > -1.53 && jtPu4Calophi[i] < 2.23 && jtPu4Calophi[i] > 2.22) continue;

      if(jtPu4Calopt[i] > maxPu4CaloAnaPt){
	maxPu4CaloAnaPt = jtPu4Calopt[i];
	maxPu4CaloAnaEta = jtPu4Caloeta[i];
	maxPu4CaloAnaPhi = jtPu4Calophi[i];
      }

      if(fabs(jtPu4Caloeta[i]) < 5.0){
	if(jtPu4Calopt[i] > maxMidPu4CaloAnaPt){
	  maxMidPu4CaloAnaPt = jtPu4Calopt[i];
	  //	  maxMidPu4CaloAnaEta = jtPu4Caloeta[i];
	}
      }

      if(fabs(jtPu4Caloeta[i]) < 0.6){
	if(jtPu4Calopt[i] > maxDiPu4CaloEta0p7AnaPt){
	  twoDiPu4CaloEta0p7AnaPt = maxDiPu4CaloEta0p7AnaPt;
	  //	  twoDiPu4CaloAnaEta = maxDiPu4CaloAnaEta

	  maxDiPu4CaloEta0p7AnaPt = jtPu4Calopt[i];
	  maxDiPu4CaloEta0p7AnaEta = jtPu4Caloeta[i];
	}
	else if(jtPu4Calopt[i] > twoDiPu4CaloEta0p7AnaPt){
          twoDiPu4CaloEta0p7AnaPt = jtPu4Calopt[i];
	  //	  twoDiPu4CaloAnaEta = jtPu4Caloeta[i];
	}
      }

      if(fabs(jtPu4Caloeta[i]) < 1.0){
	if(jtPu4Calopt[i] > maxDiPu4CaloEta1p1AnaPt){
	  twoDiPu4CaloEta1p1AnaPt = maxDiPu4CaloEta1p1AnaPt;
	  //	  twoDiPu4CaloAnaEta = maxDiPu4CaloAnaEta;

	  maxDiPu4CaloEta1p1AnaPt = jtPu4Calopt[i];
	  //	  maxDiPu4CaloEta1p1AnaEta = jtPu4Caloeta[i];
	}
	else if(jtPu4Calopt[i] > twoDiPu4CaloEta1p1AnaPt){
          twoDiPu4CaloEta1p1AnaPt = jtPu4Calopt[i];
	  //	  twoDiPu4CaloAnaEta = jtPu4Caloeta[i];
	}
      }

      if(jtPu4Calopt[i] > maxTriPu4CaloAnaPt){
	threeTriPu4CaloAnaPt = twoTriPu4CaloAnaPt;
	//	threeTriPu4CaloAnaEta = twoTriPu4CaloAnaEta;

	twoTriPu4CaloAnaPt = maxTriPu4CaloAnaPt;
	//  twoTriPu4CaloAnaEta = maxTriPu4CaloAnaEta;
	
	maxTriPu4CaloAnaPt = jtPu4Calopt[i];
	maxTriPu4CaloAnaEta = jtPu4Caloeta[i];
      }
      else if(jtPu4Calopt[i] > twoTriPu4CaloAnaPt){
        threeTriPu4CaloAnaPt = twoTriPu4CaloAnaPt;
	//	threeTriPu4CaloAnaEta = twoTriPu4CaloAnaEta;

	twoTriPu4CaloAnaPt = jtPu4Calopt[i];
	//	twoTriPu4CaloAnaEta = jtPu4Caloeta[i];
      }
      else if(jtPu4Calopt[i] > threeTriPu4CaloAnaPt){
        threeTriPu4CaloAnaPt = jtPu4Calopt[i];
	//	threeTriPu4CaloAnaEta = jtPu4Caloeta[i];                                                  
      }
    }


    Double_t maxPu3PFAnaPt = -1;
    //    Double_t maxPu3PFAnaEta = -100;
    
    for(int i = 0; i < nPu3PFref; ++i){
      if(fabs(jtPu3PFeta[i]) > 2.0) continue;
      if(jtPu3PFpt[i] > maxPu3PFAnaPt){
	maxPu3PFAnaPt = jtPu3PFpt[i];
	//	maxPu3PFAnaEta = jtPu3PFeta[i];
      }
    }

    Double_t maxPu4PFAnaPt = -1;
    Double_t maxPu4PFAnaEta = -100;
    Double_t maxPu4PFAnaPhi = -100;

    for(int i = 0; i < nPu4PFref; ++i){
      if(fabs(jtPu4PFeta[i]) > 2.0) continue;
      if(refPu4PFpt[i]/jtPu4PFpt[i] < .05) continue;
      if(jtPu4PFpt[i] > maxPu4PFAnaPt){
	maxPu4PFAnaPt = jtPu4PFpt[i];
	maxPu4PFAnaEta = jtPu4PFeta[i];
	maxPu4PFAnaPhi = jtPu4PFphi[i];
      }
    }


    /*
    Double_t maxVs3CaloAnaPt = -1;
    Double_t maxVs3CaloAnaEta = -100;
    
    for(int i = 0; i < nVs3Caloref; ++i){
      if(fabs(jtVs3Caloeta[i]) > 2.0) continue;
      if(jtVs3Calopt[i] > maxVs3CaloAnaPt){
	maxVs3CaloAnaPt = jtVs3Calopt[i];
	maxVs3CaloAnaEta = jtVs3Caloeta[i];
      }
    }

    Double_t maxVs4CaloAnaPt = -1;
    Double_t maxVs4CaloAnaEta = -100;

    for(int i = 0; i < nVs4Caloref; ++i){
      if(fabs(jtVs4Caloeta[i]) > 2.0) continue;
      if(jtVs4Calopt[i] > maxVs4CaloAnaPt){
	maxVs4CaloAnaPt = jtVs4Calopt[i];
	maxVs4CaloAnaEta = jtVs4Caloeta[i];
      }
    }
    */

    Double_t maxPhotonAnaPt = -1;
    Double_t maxPhotonAnaEta = -100;
    Double_t maxPhotonAnaPhi = -100;

    Double_t maxPhotonMidAnaPt = -1;
    Double_t maxPhotonMidAnaEta = -100;
    Double_t maxPhotonMidAnaPhi = -100;

    Double_t maxPhotonHiAnaPt = -1;
    Double_t maxPhotonHiAnaEta = -100;
    Double_t maxPhotonHiAnaPhi = -100;

    Double_t maxElectronAnaPt = -1;
    Double_t maxElectronAnaEta = -100;
    Double_t maxElectronAnaPhi = -100;
    Double_t maxElectronR9 = -1;

    Double_t twoElectronAnaPt = -1;
    Double_t twoElectronAnaEta = -100;
    Double_t twoElectronAnaPhi = -100;
    Double_t twoElectronR9 = -1;

    for(int i = 0; i < nPho; ++i){
      float tempEt = (*phoEt).at(i);
      float tempEta = (*phoEta).at(i);
      float tempPhi = (*phoPhi).at(i);
      float tempR9 = (*phoR9).at(i);
      float tempSigmaIEtaIEta = (*phoSigmaIEtaIEta).at(i);
      float tempHoverE = (*phoHoverE).at(i);


      //      if(tempR9 < 0.9) continue;
      //if(tempHoverE > 0.1) continue;
      //      if(tempSigmaIEtaIEta < 0.002) continue;
      if(fabs(tempEta) > 1.44 && fabs(tempEta) < 1.6) continue;

      //      if(fabs(tempEta) > 1.44) continue;
      //      if(TMath::Abs(tempPhoSeedTime) > 3) continue;
      //if(swissCrx[i] > 0.9) continue;
      //      if(isEle[i]) continue;

      if(tempEt > maxPhotonAnaPt && fabs(tempEta) < 1.44){	

	maxPhotonAnaPt = tempEt;
	maxPhotonAnaEta = tempEta;
	maxPhotonAnaPhi = tempPhi;

      }

      if(tempEt > maxPhotonHiAnaPt && fabs(tempEta) < 3.0){	

	maxPhotonHiAnaPt = tempEt;
	maxPhotonHiAnaEta = tempEta;
	maxPhotonHiAnaPhi = tempPhi;

      }


      /*
      if(fabs(phoEta->at(i)) > 2.0) continue;
      //      if(TMath::Abs(seedTime[i]) > 3) continue;
      //if(swissCrx[i] > 0.9) continue;
      //      if(isEle[i]) continue;

      if(phoEt->at(i) > maxPhotonAnaPt){	

	maxPhotonAnaPt = phoEt->at(i);
	maxPhotonAnaEta = phoEta->at(i);
	maxPhotonAnaPhi = phoPhi->at(i);

      }
      */
    }
    /*
    
    for(int i = 0; i < nPhotons; ++i){
      if(fabs(photonEta[i]) > 3.0) continue;
      if(TMath::Abs(seedTime[i]) > 3) continue;
      if(swissCrx[i] > 0.9) continue;
      //      if(isEle[i]) continue;

      if(photonEta[i] < 1.84 && photonEta[i] > 1.70 && photonPhi[i] < -2.00 && photonPhi[i] > -2.15) continue;

      if(photonEta[i] < -1.08 && photonEta[i] > -1.09 && photonPhi[i] < .70 && photonPhi[i] > .69) continue;

      if(photonPt[i] > maxPhotonHiAnaPt){	

	maxPhotonHiAnaPt = photonPt[i];
	maxPhotonHiAnaEta = photonEta[i];
	maxPhotonHiAnaPhi = photonPhi[i];
      }
      }*/
    /*
    for(int i = 0; i < nPhotons; ++i){
      if(fabs(photonEta[i]) > 1.44) continue;
      if(TMath::Abs(seedTime[i]) > 3) continue;
      if(swissCrx[i] > 0.9) continue;
      //      if(isEle[i]) continue;

      if(photonEta[i] < -1.08 && photonEta[i] > -1.09 && photonPhi[i] < .70 && photonPhi[i] > .69) continue;

      if(photonPt[i] > maxPhotonMidAnaPt){	

	maxPhotonMidAnaPt = photonPt[i];
	maxPhotonMidAnaEta = photonEta[i];
	maxPhotonMidAnaPhi = photonPhi[i];
      }
    }

    for(int i = 0; i < nPhotons; ++i){
      //      if(fabs(photonEta[i]) > 1.44) continue;
      //      if(!isEle[i]) continue;
      if(photonPt[i] > maxElectronAnaPt){	
	twoElectronAnaPt = maxElectronAnaPt;
	//	twoElectronAnaEta = maxElectronAnaEta;

	maxElectronAnaPt = photonPt[i];
	//	maxElectronAnaEta = photonEta[i];
	maxElectronR9 = r9[i];
      }
      else if(photonPt[i] > twoElectronAnaPt){
        twoElectronAnaPt = photonPt[i];
	//	twoElectronAnaEta = photonEta[i];

	twoElectronR9 = r9[i];
      }
    }
    */
    maxElectronAnaPt = -1;
    twoElectronAnaPt = -1;

    maxElectronAnaEta = -100;
    twoElectronAnaEta = -100;

    maxElectronAnaPhi = -100;
    twoElectronAnaPhi = -100;
    
    const Double_t electronMass = 0.000510998;
    TLorentzVector el1, el2, recombPhoton;
    Double_t trueZMass = -1;


    for(int i = 0; i < nGen; i++){
      float tempPt = (*genPt).at(i);
      float tempEta = (*genEta).at(i);
      float tempPhi = (*genPhi).at(i);
      float tempPDG = (*genPDG).at(i);


      //      if(genStatus[i] != 1) continue;
      if(TMath::Abs(tempPDG) != 22) continue;
      if(fabs(tempEta) > 3.0) continue;
      if(tempPt > maxElectronAnaPt){
	twoElectronAnaPt = maxElectronAnaPt;
	twoElectronAnaEta = maxElectronAnaEta;
	twoElectronAnaPhi = maxElectronAnaPhi;

	maxElectronAnaPt = tempPt;
	maxElectronAnaEta = tempEta;
	maxElectronAnaPhi = tempPhi;
      }
      else if(tempPt > twoElectronAnaPt){
        twoElectronAnaPt = tempPt;
	twoElectronAnaEta = tempEta;
	twoElectronAnaPhi = tempPhi;
      }
    }
    
    el1.SetPtEtaPhiM(maxElectronAnaPt, maxElectronAnaEta, maxElectronAnaPhi, electronMass);
    el2.SetPtEtaPhiM(twoElectronAnaPt, twoElectronAnaEta, twoElectronAnaPhi, electronMass);
    TLorentzVector z = el1+el2;

    trueZMass = z.M();

    zMass_p->Fill(trueZMass);
    e1Pt_p->Fill(maxElectronAnaPt);
    e2Pt_p->Fill(twoElectronAnaPt);
    


    Double_t trueZAnaPt = -1;
    Double_t trueZAnaEta = -100;
    Double_t trueZAnaPhi = -100;
    /*
    for(int i = 0; i < nGen; i++){
      if(TMath::Abs(genPDG[i]) == 23){
	trueZAnaPt = (Double_t)genPt[i];
	trueZAnaEta = (Double_t)genEta[i];
	trueZAnaPhi = (Double_t)genPhi[i];
      }
    }
    */
    zPt_p->Fill(trueZAnaPt);

    //FOR ADDITIONAL OFFLINE OBJECT MATCHING, EDIT HERE (2 of 2)
    //    const Int_t nOfflineObj = 4;
    //const Int_t nOfflineObj = 1;
    //current
        const Int_t nOfflineObj = 6;
    //    const Int_t nOfflineObj = 27;
    //    const Int_t nOfflineObj = 23;
    if(nOfflineObj != nTrigType){
      std::cout << "ERROR: OFFLINE OBJECT NUMBER MUST MATCH NUMBER OF TRIGGER 'TYPES' IN INPUT TEXT FILE; RETURN 1" << std::endl;
      std::cout << nTrigType << std::endl;
      return 1;
    }

    Double_t trigOfflinePt[nOfflineObj] = {
      maxPu4PFAnaPt,            //LEV1JET
      maxPhotonHiAnaPt,            //LEV1GAMMA
      maxPu4CaloAnaPt,            //R4JET
      maxPu4PFAnaPt,            //R4PFJET
      maxPhotonAnaPt,            //GAMMAETA1P5
      maxPhotonHiAnaPt,            //GAMMAETA3P1
    };

    Double_t trigOfflineEta[nOfflineObj] = {
      maxPu4PFAnaEta,            //LEV1JET
      maxPhotonHiAnaEta,            //LEV1GAMMA
      maxPu4CaloAnaEta,            //R4JET
      maxPu4PFAnaEta,            //R4PFJET
      maxPhotonAnaEta,            //GAMMAETA1P5
      maxPhotonHiAnaEta,            //GAMMAETA3P1
    };

    Double_t trigOfflinePhi[nOfflineObj] = {
      maxPu4PFAnaPhi,            //LEV1JET
      maxPhotonHiAnaPhi,            //LEV1GAMMA
      maxPu4CaloAnaPhi,            //R4JET
      maxPu4PFAnaPhi,            //R4PFJET
      maxPhotonAnaPhi,            //GAMMAETA1P5
      maxPhotonHiAnaPhi,            //GAMMAETA3P1
    };


    Bool_t trigCond[nOfflineObj] = {
      true,            //LEV1JET
      true,            //LEV1GAMMA
      true,            //R4JET
      true,            //R4PFJET
      true,            //GAMMAETA1P5
      true,            //GAMMAETA3P1
    };
    

    for(Int_t iter = 0; iter < nTrigType; iter++){
      if(trigOfflinePt[iter] > 0 && trigCond[iter]){
	histsPt_p[iter][0]->Fill(trigOfflinePt[iter]);

	if(trigVal[iter][0]){
	  e1R9_p[iter]->Fill(maxElectronR9);
	  e2R9_p[iter]->Fill(twoElectronR9);
	}

	for(Int_t iter2 = 0; iter2 < trigTypeCount[iter]; iter2++){
	  if(trigOfflinePt[iter] > trigThresh[iter][iter2]) histsEta_p[iter][iter2]->Fill(trigOfflineEta[iter]);
	  /*
	  if(iter2 != trigTypeCount[iter]-1){
	    if(trigOfflinePt[iter] > trigThresh[iter][iter2] && trigOfflinePt[iter] < trigThresh[iter][iter2+1]) histsSpectBoundPt_p[iter][iter2]->Fill(trigOfflinePt[iter]);
	  }
	  else{
	    if(trigOfflinePt[iter] > trigThresh[iter][iter2]) histsSpectBoundPt_p[iter][iter2]->Fill(trigOfflinePt[iter]);
	  }
	  */
	  if(trigOfflinePt[iter] > trigThresh[iter][iter2]*.8) histsSpectBoundPt_p[iter][iter2]->Fill(trigOfflinePt[iter]);

	  if(trigVal[iter][iter2]){
	    histsPt_p[iter][iter2+1]->Fill(trigOfflinePt[iter]);
	    if(trigOfflinePt[iter] > trigThresh[iter][iter2]) histsEtaTrig_p[iter][iter2]->Fill(trigOfflineEta[iter]);
	  }
	  else{
	    if(highPtMiss[iter][iter2] < trigOfflinePt[iter]){
	      highPtMiss[iter][iter2] = trigOfflinePt[iter];
	      highPtMiss_Eta[iter][iter2] = trigOfflineEta[iter];
	      highPtMiss_Phi[iter][iter2] = trigOfflinePhi[iter];
	      highPtEventNum[iter][iter2] = ana_event;
	      highPtLumiNum[iter][iter2] = ana_lumi;
	    }
	  }
	}
      }
      
      if(TMath::Abs(maxElectronAnaEta) < 2.4 && TMath::Abs(twoElectronAnaEta) < 2.4 && ((TMath::Abs(maxElectronAnaEta) > 2.0 && (TMath::Abs(twoElectronAnaEta) < 1.44 || TMath::Abs(twoElectronAnaEta) > 1.6)) || (TMath::Abs(twoElectronAnaEta) > 2.0 && (TMath::Abs(maxElectronAnaEta) < 1.44 || TMath::Abs(maxElectronAnaEta) > 1.6)))){
	for(Int_t iter2 = 0; iter2 < trigTypeCount[iter]; iter2++){
	  nZTrigFireDenom[iter][iter2]++;	  

	  if(trigVal[iter][iter2]) nZTrigFireNum[iter][iter2]++;
	}
      }
      
    }
  }

  std::cout << std::endl;
  std::cout << "Matched events: " << matched << std::endl;
  std::cout << "Trigger fires: " << std::endl;

  for(Int_t iter = 0; iter < nTrigType; iter++){
    std::cout << "  Trigger Type, raw #, rate matched (Hz), rate unmatched (Hz), highestMiss: " << trigType[iter] << std::endl;
    
    for(Int_t iter2 = 0; iter2 < trigTypeCount[iter]; iter2++){
      std::cout << "    " << trigName[iter][iter2] << ": " << nTrigFire[iter][iter2] << ", " << std::setprecision(5) << nTrigFire[iter][iter2]*1.7*1500*11246./matched << ", " << std::setprecision(5) << nTrigFire[iter][iter2]*1.7*1500*11246./HLTTree->GetEntries() << ", " << std::setprecision(5) << highPtMiss[iter][iter2] << ", " << highPtMiss_Eta[iter][iter2] << ", " << highPtMiss_Phi[iter][iter2] << ", " << highPtEventNum[iter][iter2] << ", " << highPtLumiNum[iter][iter2] << std::endl;

      Int_t bin = ratesMatched_p[iter]->FindBin(trigThresh[iter][iter2]);
      Float_t num = nTrigFire[iter][iter2]*1.7*1500*11246;
      Float_t numErr = (nTrigFire[iter][iter2] + TMath::Sqrt(nTrigFire[iter][iter2]))*1.7*1500*11246;

      //      Float_t num = nTrigFire[iter][iter2]*30000a.;
      //      Float_t numErr = (nTrigFire[iter][iter2] + TMath::Sqrt(nTrigFire[iter][iter2]))*30000.;

      ratesMatched_p[iter]->SetBinContent(bin, num/matched);
      ratesMatched_p[iter]->SetBinError(bin, numErr/matched - num/matched);
      ratesUnmatched_p[iter]->SetBinContent(bin, num/HLTTree->GetEntries());
      ratesUnmatched_p[iter]->SetBinError(bin, numErr/HLTTree->GetEntries() - num/HLTTree->GetEntries());

      maxMissedPt100_p[iter]->SetBinContent(bin, highPtMiss[iter][iter2]);
      maxMissedPt100_p[iter]->SetBinError(bin, 0);
    }
  }

  std::cout << std::endl;
  std::cout << "Z eff: " << std::endl;
  for(Int_t iter = 0; iter < nTrigType; iter++){
    std::cout << "  Trigger Type, raw #, denom #, eff" << trigType[iter] << std::endl;

    for(Int_t iter2 = 0; iter2 < trigTypeCount[iter]; iter2++){
      std::cout << "    " << trigName[iter][iter2] << ": " << nZTrigFireNum[iter][iter2] << ", " << nZTrigFireDenom[iter][iter2] << ", " << std::setprecision(5) << nZTrigFireNum[iter][iter2]/nZTrigFireDenom[iter][iter2] << std::endl;
    }
  }

  TGraphAsymmErrors *aPt_p[nTrigType][maxNTrig2], *aEta_p[nTrigType][maxNTrig2];

  for(Int_t iter = 0; iter < nTrigType; iter++){
    for(Int_t iter2 = 0; iter2 < trigTypeCount[iter]; iter2++){
      //      if(iter2 == 0) histsSpectPt_p[iter] = (TH1F*)(histsPt_p[iter][iter2+1]->Clone());
      histsSpectPt_p[iter]->Add(histsPt_p[iter][iter2+1]);

      aPt_p[iter][iter2] = new TGraphAsymmErrors();
      aPt_p[iter][iter2]->BayesDivide(histsPt_p[iter][iter2+1],histsPt_p[iter][0]);
      aPt_p[iter][iter2]->SetName(Form("%s_%s_pt_asymm", trigName[iter][iter2].c_str(), trigType[iter].c_str()));

      aEta_p[iter][iter2] = new TGraphAsymmErrors();
      aEta_p[iter][iter2]->BayesDivide(histsEtaTrig_p[iter][iter2],histsEta_p[iter][iter2]);
      aEta_p[iter][iter2]->SetName(Form("%s_%s_eta_asymm", trigName[iter][iter2].c_str(), trigType[iter].c_str()));
      Float_t pt99Miss = -1;
      for(Int_t binIter = 0; binIter < aPt_p[iter][iter2]->GetN(); binIter++){
	Double_t tempX = -1;
	Double_t tempY = -1;
	aPt_p[iter][iter2]->GetPoint(binIter, tempX, tempY);

	if(!strcmp("TRACK", trigType[iter].substr(0,5).c_str())){
	  if(!strcmp("R2", trigType[iter].substr(5,2).c_str()) && tempY < .95) pt99Miss = tempX;
	  else if(tempY < .99 && strcmp("R2", trigType[iter].substr(5,2).c_str()) != 0) pt99Miss = tempX;
	}
	else if(tempY < .99) pt99Miss = tempX;
      }
      Int_t bin = maxMissedPt99_p[iter]->FindBin(trigThresh[iter][iter2]);
      maxMissedPt99_p[iter]->SetBinContent(bin, pt99Miss);
      maxMissedPt99_p[iter]->SetBinError(bin, 0);
    }
  }

  TFile* optMaxFile_p;
  TH1F* optMaxMissed_p[nTrigType];
  Bool_t isMax = true;
  Int_t nEvt[nTrigType][maxNTrig2];

  for(Int_t iter = 0; iter < nTrigType; iter++){
    for(Int_t iter2 = 0; iter2 < maxNTrig2; iter2++){
      nEvt[iter][iter2] = 0;
    }
  }

  if(!strcmp(inOptMaxName.c_str(), "")) isMax = false;
 
  for(Int_t iter = 0; iter < nTrigType; iter++){
    if(strcmp("TRACK", trigType[iter].substr(0,5).c_str()) != 0) continue;

    if(isMax){
      optMaxFile_p = new TFile(inOptMaxName.c_str(), "READ");
      optMaxMissed_p[iter] = (TH1F*)optMaxFile_p->Get(Form("%sDir/maxMissedPt99_%s_%d_%d", trigType[iter].c_str(), trigType[iter].c_str(), trigThresh[iter][0], trigThresh[iter][trigTypeCount[iter] - 1]));
    }
    else optMaxMissed_p[iter] = (TH1F*)maxMissedPt99_p[iter]->Clone("temp");
  }


  if(isMax) delete optMaxFile_p;


  std::string outName = inHLTFile;
  const std::string inString = ".root";
  const std::string outString = "_HIST.root";
  std::size_t strIndex = 0;

  strIndex = outName.find(inString);
  if(!(strIndex == std::string::npos)){
    outName.replace(strIndex, inString.length(), outString); 
  }

  TFile* outFile_p = new TFile(outName.c_str(), "UPDATE");

  zPt_p->Write("", TObject::kOverwrite);
  zMass_p->Write("", TObject::kOverwrite);
  e1Pt_p->Write("", TObject::kOverwrite);
  e2Pt_p->Write("", TObject::kOverwrite);

  for(Int_t iter = 0; iter < nTrigType; iter++){
    outFile_p->cd();
    const std::string dirName = Form("%sDir", trigType[iter].c_str());
    TDirectory* dir_p = outFile_p->GetDirectory(Form("%s", dirName.c_str()));
    if(dir_p){
      dir_p->cd();
    }
    else{
      dir_p = outFile_p->mkdir(Form("%s", dirName.c_str()));
      dir_p->cd();
    }

    histsPt_p[iter][0]->Write("", TObject::kOverwrite);
    histsSpectPt_p[iter]->Write("", TObject::kOverwrite);

    for(Int_t iter2 = 0; iter2 < trigTypeCount[iter]; iter2++){
      histsPt_p[iter][iter2+1]->Write("", TObject::kOverwrite);
      histsEta_p[iter][iter2]->Write("", TObject::kOverwrite);
      histsEtaTrig_p[iter][iter2]->Write("", TObject::kOverwrite);

      histsSpectBoundPt_p[iter][iter2]->Write("", TObject::kOverwrite);

      if(!strcmp("TRACK", trigType[iter].substr(0, 5).c_str())){

	std::cout << nEvt[iter][iter2] << std::endl;
	histsTrkSpectBoundPt_p[iter][iter2]->Scale(1./nEvt[iter][iter2]);
	histsTrkSpectBoundPt_p[iter][iter2]->Write("", TObject::kOverwrite);

	histsJetSpectBoundPt_p[iter][iter2]->Scale(1./nEvt[iter][iter2]);
	histsJetSpectBoundPt_p[iter][iter2]->Write("", TObject::kOverwrite);
      }

      aPt_p[iter][iter2]->Write("", TObject::kOverwrite);
      aEta_p[iter][iter2]->Write("", TObject::kOverwrite);
    }
    ratesMatched_p[iter]->Write("", TObject::kOverwrite);
    ratesUnmatched_p[iter]->Write("", TObject::kOverwrite);

    maxMissedPt99_p[iter]->Write("", TObject::kOverwrite);
    maxMissedPt100_p[iter]->Write("", TObject::kOverwrite);
    prescale_p[iter]->Write("", TObject::kOverwrite);
    e1R9_p[iter]->Write("", TObject::kOverwrite);
    e2R9_p[iter]->Write("", TObject::kOverwrite);
  }

  outFile_p->Close();
  delete outFile_p;

  delete zPt_p;
  delete zMass_p;
  delete e1Pt_p;
  delete e2Pt_p;

  for(Int_t iter = 0; iter < nTrigType; iter++){
    delete histsPt_p[iter][0];
    delete histsSpectPt_p[iter];

    for(Int_t iter2 = 0; iter2 < trigTypeCount[iter]; iter2++){
      delete histsEta_p[iter][iter2];
      delete histsEtaTrig_p[iter][iter2];

      delete aPt_p[iter][iter2];
      delete aEta_p[iter][iter2];
    }
    delete ratesMatched_p[iter];
    delete ratesUnmatched_p[iter];

    delete maxMissedPt99_p[iter];
    delete maxMissedPt100_p[iter];
    delete prescale_p[iter];
    delete e1R9_p[iter];
    delete e2R9_p[iter];
  }


  delete AnaPu3CaloTree;
  delete AnaPu4CaloTree;
  delete AnaPu3PFTree;
  delete AnaPu4PFTree;
  delete AnaVs3CaloTree;
  delete AnaVs4CaloTree;
  delete AnaPhotonTree;
  delete AnaHITree;
  delete AnaSkimTree;
  delete AnaTrkTree;
  delete AnaGenTree;



  HLTFile_p->Close();
  delete HLTFile_p;

  return 0;
}

int main(int argc, char *argv[])
{
  TApplication a("a", 0, 0);

  if(argc != 4 && argc != 5){
    std::cout << "Usage: matchTrigTree_HI <inHLTFile> <inForestFile> <inTrigFileName> <inOptMaxName>" << std::endl;
    std::cout << "argNum: " << argc << std::endl;
    for(Int_t iter = 0; iter < argc; iter++){
      std::cout << "arg " << iter << ": " << argv[iter] << std::endl;
    }

    return -1;
  }

  int rStatus = -1;

  if(argc == 4) rStatus = matchTrigTree_HI(argv[1], argv[2], argv[3]);
  else if(argc == 5) rStatus = matchTrigTree_HI(argv[1], argv[2], argv[3], argv[4]);

  return rStatus;
}
