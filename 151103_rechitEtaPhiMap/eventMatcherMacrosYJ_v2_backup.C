/*
 * eventMatcherMacros.C
 *
 *
 * This is getting more and more confusing to me. I've asked Lindsey for more help debugging on the github page.

 Kaya, can you try to re-run the photon producer step during foresting (with your patch)? You should just be able to add "process.photons" to ana_step before the ggHiNtuplizer, and make sure to set useValMapIso in ggHiNtuplizer to False.
 Then, run on the parent RECO for the 75X HiForest you made:
 /mnt/hadoop/cms/store/user/tatar/HIHighPt/HiForest_HIHighPt_photon30_HIRun2011-v1/0.root
 (I *think* that the parent RECO is this one :
 /HIHighPt/dgulhan-HIHighPt_photon20and30_HIRun2011-v1_RECO_753_patch1-fd44351629dd155a25de2b4c109c824c/USER could you please confirm?)
 Compare the 53X sigmaIetaIeta from this forest : /mnt/hadoop/cms/store/user/luck/L1Emulator/HiForest_PbPb_photon2030.root to the new results you get from full5x5_sigmaIetaIeta, and see if they match better than the new sigmaIetaIeta. We're trying to see if the full5x5 as calculated by your patch is a better match to the 53X sieie for the same photons. I think you'll only need to run on a few hundred events (so no crab needed), and then do the usual Event Matching.

 Please keep me updated via email, and I'll jump in to help when I can find some time tomorrow.

 Alex

 ***********
 modified by Yeonju 
v2 : not using getMatchingEntry() function.
: using ematcher->retrieveEvent(event_file2, lumi_file2, run_file2);

*/

#include <TFile.h>
#include <TTree.h>
#include <TCut.h>
#include <TH1.h>
#include <TH1D.h>
#include <TList.h>

#include <vector>
#include <string>
#include <iostream>
#include <iomanip>

#include "../HIUtils/treeUtil.h"
#include "../HIUtils/histoUtil.h"
#include "../HIUtils/EventMatchingCMS.h"

const int MAXPHOTONS = 500;
const int MAXHITS = 50;
const float cutPt = 10;
const float cutEta = 1.4791;

//void eventMatcherMacros(const char* filePath1, const char* filePath2);

//void eventMatcherMacrosYJ(const char* filePath1, const char* filePath2)
void eventMatcherMacrosYJ_v2()
{
    const char* filePath1= "/mnt/hadoop/cms/store/user/luck/L1Emulator/HiForest_PbPb_photon2030.root";//entry = 349248
    const char* filePath2= "/mnt/hadoop/cms/store/user/luck/HIMinBiasUPC/2011_MB_750_hiForest/0.root";//entry = 1662337
    TFile* inputFile1= new TFile(filePath1, "READ");
    TFile* inputFile2= new TFile(filePath2, "READ");

    TTree *t_ee1 = (TTree*)inputFile1->Get("rechitanalyzer/ee");
    TTree *t_eb1 = (TTree*)inputFile1->Get("rechitanalyzer/eb");
    TTree *t_pho1 = (TTree*)inputFile1->Get("ggHiNtuplizer/EventTree");
    TTree *t_evt1 = (TTree*)inputFile1->Get("hiEvtAnalyzer/HiTree");
    TTree *t_ee2 = (TTree*)inputFile2->Get("rechitanalyzer/ee");
    TTree *t_eb2 = (TTree*)inputFile2->Get("rechitanalyzer/eb");
    TTree *t_pho2 = (TTree*)inputFile1->Get("multiPhotonAnalyzer/photon");
    TTree *t_evt2 = (TTree*)inputFile2->Get("hiEvtAnalyzer/HiTree");

    // Evt
    Int_t run1, evt1, lumi1;
    t_evt1->SetBranchAddress("run", &run1);
    t_evt1->SetBranchAddress("evt", &evt1);
    t_evt1->SetBranchAddress("lumi", &lumi1);
    // EB hits
    Int_t eb_n1;
    Float_t eb_e1[MAXHITS];
    Float_t eb_et1[MAXHITS];
    Float_t eb_eta1[MAXHITS];
    Float_t eb_phi1[MAXHITS];
    Bool_t  eb_isjet1[MAXHITS];
    t_eb1->SetBranchAddress("n",&eb_n1);
    t_eb1->SetBranchAddress("e",eb_e1);
    t_eb1->SetBranchAddress("et",eb_et1);
    t_eb1->SetBranchAddress("eta",eb_eta1);
    t_eb1->SetBranchAddress("phi",eb_phi1);
    t_eb1->SetBranchAddress("isjet",eb_isjet1);
    // EE hits
    Int_t ee_n1;
    Float_t ee_e1[MAXHITS];
    Float_t ee_et1[MAXHITS];
    Float_t ee_eta1[MAXHITS];
    Float_t ee_phi1[MAXHITS];
    Bool_t  ee_isjet1[MAXHITS];
    t_ee1->SetBranchAddress("n",&ee_n1);
    t_ee1->SetBranchAddress("e",ee_e1);
    t_ee1->SetBranchAddress("et",ee_et1);
    t_ee1->SetBranchAddress("eta",ee_eta1);
    t_ee1->SetBranchAddress("phi",ee_phi1);
    t_ee1->SetBranchAddress("isjet",ee_isjet1);

    //Evt
    Int_t run2, evt2, lumi2;
    t_evt2->SetBranchAddress("run", &run2);
    t_evt2->SetBranchAddress("evt", &evt2);
    t_evt2->SetBranchAddress("lumi", &lumi2);
    // EB hits
    Int_t eb_n2;
    Float_t eb_e2[MAXHITS];
    Float_t eb_et2[MAXHITS];
    Float_t eb_eta2[MAXHITS];
    Float_t eb_phi2[MAXHITS];
    Bool_t  eb_isjet2[MAXHITS];
    t_eb2->SetBranchAddress("n",&eb_n2);
    t_eb2->SetBranchAddress("e",eb_e2);
    t_eb2->SetBranchAddress("et",eb_et2);
    t_eb2->SetBranchAddress("eta",eb_eta2);
    t_eb2->SetBranchAddress("phi",eb_phi2);
    t_eb2->SetBranchAddress("isjet",eb_isjet2);
    // EE hits
    Int_t ee_n2;
    Float_t ee_e2[MAXHITS];
    Float_t ee_et2[MAXHITS];
    Float_t ee_eta2[MAXHITS];
    Float_t ee_phi2[MAXHITS];
    Bool_t  ee_isjet2[MAXHITS];
    t_ee2->SetBranchAddress("n",&ee_n2);
    t_ee2->SetBranchAddress("e",ee_e2);
    t_ee2->SetBranchAddress("et",ee_et2);
    t_ee2->SetBranchAddress("eta",ee_eta2);
    t_ee2->SetBranchAddress("phi",ee_phi2);
    t_ee2->SetBranchAddress("isjet",ee_isjet2);
/*
    Int_t nPho;
    std::vector<float>* phoEt=0;
    std::vector<float>* phoEta=0;
    std::vector<float>* phoSigmaIEtaIEta=0;
    std::vector<float>* phoSigmaIEtaIEta_2012=0;

    tree1->SetBranchAddress("nPho",&nPho);
    tree1->SetBranchAddress("phoEt",&phoEt);
    tree1->SetBranchAddress("phoEta",&phoEta);
    tree1->SetBranchAddress("phoSigmaIEtaIEta",&phoSigmaIEtaIEta);
    tree1->SetBranchAddress("phoSigmaIEtaIEta_2012",&phoSigmaIEtaIEta_2012);

    Int_t nPhotons;
    Float_t pt[MAXPHOTONS];
    Float_t eta[MAXPHOTONS];
    Float_t sigmaIetaIeta[MAXPHOTONS];

    tree2->SetBranchAddress("nPhotons",&nPhotons);
    tree2->SetBranchAddress("pt",pt);
    tree2->SetBranchAddress("eta",eta);
    tree2->SetBranchAddress("sigmaIetaIeta",sigmaIetaIeta);
*/


    const int MAXHIT = 100;
    Int_t run_, evt_, lumi_;
    Int_t ebn1_, een1_, ebn2, een2;
    Float_t ebe1_[MAXHIT];
    Float_t ebet1_[MAXHIT];
    Float_t ebeta1_[MAXHIT];
    Float_t ebphi1_[MAXHIT];
    Bool_t ebisjet1_[MAXHIT];
    Float_t eee1_[MAXHIT];
    Float_t eeet1_[MAXHIT];
    Float_t eeeta1_[MAXHIT];
    Float_t eephi1_[MAXHIT];
    Bool_t eeisjet1_[MAXHIT];
    Float_t ebe2_[MAXHIT];
    Float_t ebet2_[MAXHIT];
    Float_t ebeta2_[MAXHIT];
    Float_t ebphi2_[MAXHIT];
    Bool_t ebisjet2_[MAXHIT];
    Float_t eee2_[MAXHIT];
    Float_t eeet2_[MAXHIT];
    Float_t eeeta2_[MAXHIT];
    Float_t eephi2_[MAXHIT];
    Bool_t eeisjet2_[MAXHIT];
/*    std::vector<float> *ebe1_ = new std::vector<float>();
    std::vector<float> *ebet1_ = new std::vector<float>();
    std::vector<float> *ebeta1_ = new std::vector<float>();
    std::vector<float> *ebphi1_ = new std::vector<float>();
    std::vector<bool> *ebisjet1_ = new std::vector<bool>();
    std::vector<float> *eee1_ = new std::vector<float>();
    std::vector<float> *eeet1_ = new std::vector<float>();
    std::vector<float> *eeeta1_ = new std::vector<float>();
    std::vector<float> *eephi1_ = new std::vector<float>();
    std::vector<bool> *eeisjet1_ = new std::vector<bool>();
    std::vector<float> *ebe2_ = new std::vector<float>();
    std::vector<float> *ebet2_ = new std::vector<float>();
    std::vector<float> *ebeta2_ = new std::vector<float>();
    std::vector<float> *ebphi2_ = new std::vector<float>();
    std::vector<bool> *ebisjet2_ = new std::vector<bool>();
    std::vector<float> *eee2_ = new std::vector<float>();
    std::vector<float> *eeet2_ = new std::vector<float>();
    std::vector<float> *eeeta2_ = new std::vector<float>();
    std::vector<float> *eephi2_ = new std::vector<float>();
    std::vector<bool> *eeisjet2_ = new std::vector<bool>();
*/

    TTree *newt = new TTree("t","event matched rechitTree");
    newt->SetMaxTreeSize(MAXTREESIZE);
    newt->Branch("run",&run_,"run_/I");
    newt->Branch("evt",&evt_,"evt_/I");
    newt->Branch("lumi",&lumi_,"lumi_/I");
    newt->Branch("eb_n1",&ebn1_,"ebn1_/I");
    newt->Branch("eb_e1",ebe1_,"ebe1[eb_n1]/F");
/*
    newt->Branch("eb_et1",
    newt->Branch("eb_eta1",
    newt->Branch("eb_phi1",
    newt->Branch("eb_isjet1",
    newt->Branch("ee_n1",&een1_,"een1_/I");
    newt->Branch("ee_e1",
    newt->Branch("ee_et1",
    newt->Branch("ee_eta1",
    newt->Branch("ee_phi1",
    newt->Branch("ee_isjet1",
    newt->Branch("eb_n2",&ebn2_,"ebn2_/I");
    newt->Branch("eb_e2",
    newt->Branch("eb_et2",
    newt->Branch("eb_eta2",
    newt->Branch("eb_phi2",
    newt->Branch("eb_isjet2",
    newt->Branch("ee_n2",&een2_,"een2_/I");
    newt->Branch("ee_e2",
    newt->Branch("ee_et2",
    newt->Branch("ee_eta2",
    newt->Branch("ee_phi2",
    newt->Branch("ee_isjet2",
*/
    /*
    newt->Branch("eb_e1","vector<float>",&ebe1_);
    newt->Branch("eb_et1","vector<float>",&ebet1_);
    newt->Branch("eb_eta1","vector<float>",&ebeta1_);
    newt->Branch("eb_phi1","vector<float>",&ebphi1_);
    newt->Branch("eb_isjet1","vector<float>",&ebisjet1_);
    newt->Branch("ee_n1",&een1_,"een1_/I");
    newt->Branch("ee_e1","vector<float>",&eee1_);
    newt->Branch("ee_et1","vector<float>",&eeet1_);
    newt->Branch("ee_eta1","vector<float>",&eeeta1_);
    newt->Branch("ee_phi1","vector<float>",&eephi1_);
    newt->Branch("ee_isjet1","vector<float>",&eeisjet1_);
    newt->Branch("eb_n2",&ebn2_,"ebn2_/I");
    newt->Branch("eb_e2","vector<float>",&ebe2_);
    newt->Branch("eb_et2","vector<float>",&ebet2_);
    newt->Branch("eb_eta2","vector<float>",&ebeta2_);
    newt->Branch("eb_phi2","vector<float>",&ebphi2_);
    newt->Branch("eb_isjet2","vector<float>",&ebisjet2_);
    newt->Branch("ee_n2",&een2_,"een2_/I");
    newt->Branch("ee_e2","vector<float>",&eee2_);
    newt->Branch("ee_et2","vector<float>",&eeet2_);
    newt->Branch("ee_eta2","vector<float>",&eeeta2_);
    newt->Branch("ee_phi2","vector<float>",&eephi2_);
    newt->Branch("ee_isjet2","vector<float>",&eeisjet2_);
*/
    const char* outputFileName = "skimFiles/eventMatched_v2_53X_75X.root";
    TFile* outputFile = new TFile(outputFileName, "RECREATE");

    Long64_t entries1 = t_evt1->GetEntries();
    Long64_t entries_notMatched = 0;
    const char* selection;

    EventMatchingCMS* ematcher = new EventMatchingCMS();
    int matched_entry = -1;
    for (int i=0; i<entries1; i++) {
        t_ee1->GetEntry(i);
        t_eb1->GetEntry(i);
        t_evt1->GetEntry(i);
        ematcher->addEvent(evt1, lumi1, run1, i);
    }

    Long64_t entries2 = t_evt2->GetEntries();
    for(int j=0;j<entries2;++j){
        if (j % 2000 == 0)  {
            std::cout << "current entry = " <<j<<" out of "<<entries2<<" : "<<std::setprecision(2    )<<(double)j/entries2*100<<" %"<<std::endl;
        }
        t_ee2->GetEntry(j);
        t_eb2->GetEntry(j);
        t_evt2->GetEntry(j);
        matched_entry = ematcher->retrieveEvent(evt2, lumi2, run2);
        if (matched_entry == -1){
            entries_notMatched++;
            continue;
        }

        t_ee1->GetEntry(matched_entry);
        t_eb1->GetEntry(matched_entry);
        t_evt1->GetEntry(matched_entry);
/*
        ebe1_->clear();
        ebet1_->clear();
        ebeta1_->clear();
        ebphi1_->clear();
        ebisjet1_->clear();
        eee1_->clear();
        eeet1_->clear();
        eeeta1_->clear();
        eephi1_->clear();
        eeisjet1_->clear();
        ebe2_->clear();
        ebet2_->clear();
        ebeta2_->clear();
        ebphi2_->clear();
        ebisjet2_->clear();
        eee2_->clear();
        eeet2_->clear();
        eeeta2_->clear();
        eephi2_->clear();
        eeisjet2_->clear();
*/
        run_=run1;
        evt_=evt1;
        lumi_=lumi1;
        ebn1_=eb_n1;
        een1_=ee_n1;
        ebn2_=eb_n2;
        een2_=ee_n2;

        for(int i=0;i<eb_n1;i++){
            ebe1_->push_back(eb_e1[i]);
            ebet1_->push_back(eb_et1[i]);
            ebeta1_->push_back(eb_eta1[i]);
            ebphi1_->push_back(eb_phi1[i]);
            ebisjet1_->push_back(eb_isjet1[i]);
        }
/*
        eee1_->push_back();
        eeet1_->push_back();
        eeeta1_->push_back();
        eephi1_->push_back();
        eeisjet1_->push_back();
        ebe2_->push_back();
        ebet2_->push_back();
        ebeta2_->push_back();
        ebphi2_->push_back();
        ebisjet2_->push_back();
        eee2_->push_back();
        eeet2_->push_back();
        eeeta2_->push_back();
        eephi2_->push_back();
        eeisjet2_->push_back();
*/

    }
    cout << "toal matched entries " << entries2 - entries_notMatched << endl;

    outputFile->cd();
    newt->Write(); 
    outputFile->Write();
    outputFile->Close();

    /*
       TFile* outFile = new TFile("out_eventMatcherMacros.root","RECREATE");

       TH1D* h[3];
       h[0]=new TH1D("h_phoSigmaIEtaIEta","phoSigmaIEtaIEta pT>10, abs(eta)>1.4791",100,0,0.1);
       h[1]=new TH1D("h_phoSigmaIEtaIEta_2012","phoSigmaIEtaIEta_2012 pT>10, abs(eta)>1.4791",100,0,0.1);
       h[2]=new TH1D("h_sigmaIetaIeta","sigmaIetaIeta",100,0,0.1);

       Long64_t entries = tree1Event->GetEntries();
       Long64_t entries_notMatched = 0;
       const char* selection;
       std::cout<< "Loop : ggHiNtuplizer/EventTree" <<std::endl;
       for (int j=0; j<entries; ++j)
       {
       if (j % 2000 == 0)  {
       std::cout << "current entry = " <<j<<" out of "<<entries<<" : "<<std::setprecision(2)<<(double)j/entries*100<<" %"<<std::endl;
       }

       tree1Event->GetEntry(j);
       tree1->GetEntry(j);

    //           selection = Form("run == %d && event == %d && lumis == %d", run1, event1, lumis1);
    selection = Form("Run == %d && Event == %d && LumiBlock == %d", run1, event1, lumis1);
    Long64_t j2 = getMatchingEntry(tree2Event,selection);
    if (j2 == -1)
    {
    std::cout<< "mismatched event selection : " << selection << std::endl;
    entries_notMatched++;
    continue;
    }
    tree2Event->GetEntry(j2);
    tree2->GetEntry(j2);

    for(int i=0; i<nPho; ++i)
    {
    if(!(phoEt->at(i) > cutPt && TMath::Abs(phoEta->at(i))>cutEta)) continue;
    h[0]->Fill(phoSigmaIEtaIEta->at(i));
    h[1]->Fill(phoSigmaIEtaIEta_2012->at(i));
    }

    for(int i=0; i<nPhotons; ++i)
    {
    if(!(pt[i] > cutPt && TMath::Abs(eta[i])>cutEta)) continue;
    h[2]->Fill(sigmaIetaIeta[i]);
    }
    }
    std::cout<< "Loop ENDED : ggHiNtuplizer/EventTree" <<std::endl;

    std::cout << "entries            = " << entries << std::endl;
    std::cout << "entries_notMatched = " << entries_notMatched << std::endl;

    TCanvas* c1 = drawSame(h[0],h[2]);
    TCanvas* c2 = drawSame(h[1],h[2]);

    c1->Write();
    c2->Write();

    outFile->Write();
    outFile->Close();

    inputFile1->Close();
    inputFile2->Close();
}

int main(int argc, char** argv)
{
    const char* fileName1="/net/hisrv0001/home/tatar/output/HiForestAOD_full5x5_53X.root";
    const char* fileName2="/mnt/hadoop/cms/store/user/luck/L1Emulator/HiForest_PbPb_photon2030.root";

    std::cout << "comparing files" <<std::endl;
    std::cout << "fileName1 = " << fileName1 <<std::endl;
    std::cout << "fileName2 = " << fileName2 <<std::endl;

    eventMatcherMacros(fileName1, fileName2);

    return 0;
}
*/
}
