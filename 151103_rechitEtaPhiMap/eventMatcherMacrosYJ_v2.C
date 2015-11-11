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


    const char* outputFileName = "skimFiles/eventMatched_v2_53X_75X.root";
    TFile* outputFile = new TFile(outputFileName, "RECREATE");


    const int MAXHIT = 100;
    Int_t run_, evt_, lumi_;
    Int_t eb_n1_, ee_n1_, eb_n2_, ee_n2_;
    Float_t eb_e1_[MAXHIT];
    Float_t eb_et1_[MAXHIT];
    Float_t eb_eta1_[MAXHIT];
    Float_t eb_phi1_[MAXHIT];
    Bool_t eb_isjet1_[MAXHIT];
    Float_t ee_e1_[MAXHIT];
    Float_t ee_et1_[MAXHIT];
    Float_t ee_eta1_[MAXHIT];
    Float_t ee_phi1_[MAXHIT];
    Bool_t ee_isjet1_[MAXHIT];
    Float_t eb_e2_[MAXHIT];
    Float_t eb_et2_[MAXHIT];
    Float_t eb_eta2_[MAXHIT];
    Float_t eb_phi2_[MAXHIT];
    Bool_t eb_isjet2_[MAXHIT];
    Float_t ee_e2_[MAXHIT];
    Float_t ee_et2_[MAXHIT];
    Float_t ee_eta2_[MAXHIT];
    Float_t ee_phi2_[MAXHIT];
    Bool_t ee_isjet2_[MAXHIT];

    TTree *newt_eb1 = t_eb1->CloneTree(0); 
    TTree *newt_ee1 = t_ee1->CloneTree(0); 
    TTree *newt_evt1 = t_evt1->CloneTree(0); 
    TTree *newt_eb2 = t_eb2->CloneTree(0); 
    TTree *newt_ee2 = t_ee2->CloneTree(0); 
    newt_eb1->SetName("eb53x");
    newt_ee1->SetName("ee53x");
    newt_evt1->SetName("HiTree");
    newt_eb2->SetName("eb75x");
    newt_ee2->SetName("ee75x");

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
        if (j % 50000 == 0)  {
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

        newt_eb1->Fill();
        newt_ee1->Fill();
        newt_evt1->Fill();
        newt_eb2->Fill();
        newt_ee2->Fill();

    }
    cout << "toal matched entries " << entries2 - entries_notMatched << endl;

    outputFile->cd();
    newt_eb1->Write();
    newt_ee1->Write();
    newt_evt1->Write();
    newt_eb2->Write();
    newt_ee2->Write();

    outputFile->Write();
    outputFile->Close();
}
