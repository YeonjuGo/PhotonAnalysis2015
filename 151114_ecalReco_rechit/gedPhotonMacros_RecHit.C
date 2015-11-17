/* modified by Yeonju 
 * to check rechit distribution (energy density) around photon
 * from 2015 Nov. 14
 *
 *
 * gedPhotonMacros_RecHit.C
 *
 Hi  Kaya.

Let me repeat my comments about the GED photon problem in the last meeting.  I think it might have been unclear due to my crappy English speaking :-)

Anyway, I think the first problem we have to understand is "why are the number of photons are so different between GED and old style photon"?    We suspects that fake photons are playing a role and we don’t fully understand it.  So, let’s make some basic map of Ecal RecHit energies around the photons and get some intuitions how the fake photons look like using our naked eyes.   The basic idea is to turn on the EcalRecHiAnalyzer in the forest configuration and save the energy information of all the ecal crystals.   But, it will take million years to process all the events, so let’s pick handful number of events to save our time.

I would suggest this process :

From the 2011 data re-recoed forest file, find the events having high pT photon candidates in 3 pools.
First pool is.. it’s GED photon but not matched to old photon,
Second pool is.. it’s old photon coordinate but not matched to GED photon
Third pool is… it’s old photon coordinate and matched to GED photon ( or the order of algorithm can be changed)

Can you maybe get run number, event number, LumiSection number of those events?   I would suggest to start by 10 events from each of 3 pools.

Bringing the list of Run:LS:Event number back to the RECO dataset, turn on the recHitAnalyzer in the forest configuration file and use eventsToProcess filter to select only the events of your interest.   Please see this example.

process.source = cms.Source ("PoolSource",
                             fileNames = cms.untracked.vstring(“input file”),
                             eventsToProcess = cms.untracked.VEventRange(
        '183013:923:37329384’,   '183013:923:37329672’        # Run:LS:Event
        )
                             )
This filter allows you to select only the events of your interest.
Anyway, once you get the forest file, you can draw 2D histogram of pT, eta, phi of rechits around photon candidates.  And let’s see how does the shapes of the receipt distribution look like for 3 pools.

You can also check if the isolation energy is correctly counted using this way.
If anything is not clear, please let me know.

Best,
Yongsun
 */

#include <TFile.h>
#include <TTree.h>
#include <TCut.h>
#include <TH2D.h>
#include <TMath.h>

#include <iostream>
#include <iomanip>      // std::setprecision()

//#include "../HIUtils/smallPhotonUtil.h"
#include "../yjUtility.h"

const int MAXHITS = 20;
const int MAXPFCANDS = 15000;
const float pTcut = 30;
const float cutDeltaR  = 0.2;
const float cutDeltapT = 5;
const int cutpfId = 4;  // 4=photon, https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideParticleFlow#ParticleTypes


void gedPhotonMacros_RecHit(const char* inputfileName, const char* outputFileName);
void gedPhotonMacros_RecHit2(const char* inputfileName, const char* outputFileName);
void gedPhotonMacros_PfCand(const char* inputfileName, const char* outputFileName);
void make2DRecHit(Long64_t entry, TTree* t);

void gedPhotonMacros_RecHit(const char* inputfileName, const char* outputFileName)
{
    
    TFile *file = TFile::Open(inputfileName);
    TTree* ggHiNtuplizerTree = (TTree*)file->Get("ggHiNtuplizer/EventTree");
    TTree* ee = (TTree*)file->Get("rechitanalyzer/ee");
    TTree* eb = (TTree*)file->Get("rechitanalyzer/eb");
    TTree* masterTree = (TTree*)ggHiNtuplizerTree->Clone("t1");
    masterTree->AddFriend(ee,"ee");
    masterTree->AddFriend(eb,"eb");

    float ptCut = 30;
    TCut selection_event = "Entry$<1000";
    TCut selection_event2 = "Entry$<5000";
    TCut selection_pt_t1 = Form("Max$(t1.phoEt)>%f", ptCut);
    //TCut selection_pt_t1 = Form("Max$(t1.phoEt)>%f && t1.nPho<3", ptCut);
    TCut selection_recHit_n = "eb.n>0 || ee.n>0";

    std::cout<<"Scan started"<<std::endl;

    masterTree->SetScanField(0);
    std::cout<<selection_pt_t1.GetTitle()<<std::endl;
    masterTree->Scan("run:event:lumis:Entry$:t1.nPho:t1.phoEt:t1.phoEta:t1.phoPhi",
                                                    (selection_event && selection_pt_t1).GetTitle());
    masterTree->Scan("eb.n:eb.e:eb.et:eb.eta:eb.phi:eb.perp:eb.isjet", (selection_event && selection_pt_t1).GetTitle());
    Long64_t selected_pt_t1 = masterTree->GetEntries((selection_pt_t1).GetTitle());
    std::cout<<"selected_pt_t1 = "<<selected_pt_t1<<std::endl;
    std::cout<<"Scan finished"<<std::endl;

    TFile* outputFile = new TFile(outputFileName, "RECREATE");

    // selection_pt_t1
    make2DRecHit(837,masterTree);
/*
 * ***********************************************************************************************************************
 * *    Row   * Instance *       run *     event *     lumis *    Entry$ *   t1.nPho *  t1.phoEt * t1.phoEta * t1.phoPhi *
 * ***********************************************************************************************************************
*/
    // selection_pt_t1_t2
    outputFile->Write();
    outputFile->Close();
    file->Close();
}

void gedPhotonMacros_RecHit2(const char* inputfileName, const char* outputFileName)
{
    TFile *file = new TFile(inputfileName, "READ");
    TTree* ggHiNtuplizerTree = (TTree*)file->Get("ggHiNtuplizer/EventTree");
    TTree* ee = (TTree*)file->Get("rechitanalyzer/ee");
    TTree* eb = (TTree*)file->Get("rechitanalyzer/eb");

    // OLD RECO photons
    Int_t nPho;
    std::vector<float>* phoEt=0;
    std::vector<float>* phoEta=0;
    std::vector<float>* phoPhi=0;
    std::vector<float>* pho_ecalClusterIsoR4=0;
    std::vector<float>* pho_hcalRechitIsoR4=0;
    std::vector<float>* pho_trackIsoR4PtCut20=0;

    ggHiNtuplizerTree->SetBranchAddress("nPho",&nPho);
    ggHiNtuplizerTree->SetBranchAddress("phoEt",&phoEt);
    ggHiNtuplizerTree->SetBranchAddress("phoEta",&phoEta);
    ggHiNtuplizerTree->SetBranchAddress("phoPhi",&phoPhi);
    ggHiNtuplizerTree->SetBranchAddress("pho_ecalClusterIsoR4" ,&pho_ecalClusterIsoR4);
    ggHiNtuplizerTree->SetBranchAddress("pho_hcalRechitIsoR4"  ,&pho_hcalRechitIsoR4);
    ggHiNtuplizerTree->SetBranchAddress("pho_trackIsoR4PtCut20",&pho_trackIsoR4PtCut20);
    // EB hits
    Int_t eb_n;
    Float_t eb_e[MAXHITS];
    Float_t eb_et[MAXHITS];
    Float_t eb_eta[MAXHITS];
    Float_t eb_phi[MAXHITS];
    Bool_t  eb_isjet[MAXHITS];

    eb->SetBranchAddress("n",&eb_n);
    eb->SetBranchAddress("e",eb_e);
    eb->SetBranchAddress("et",eb_et);
    eb->SetBranchAddress("eta",eb_eta);
    eb->SetBranchAddress("phi",eb_phi);
    eb->SetBranchAddress("isjet",eb_isjet);

    // EE hits
    Int_t ee_n;
    Float_t ee_e[MAXHITS];
    Float_t ee_et[MAXHITS];
    Float_t ee_eta[MAXHITS];
    Float_t ee_phi[MAXHITS];
    Bool_t  ee_isjet[MAXHITS];

    ee->SetBranchAddress("n",&ee_n);
    ee->SetBranchAddress("e",ee_e);
    ee->SetBranchAddress("et",ee_et);
    ee->SetBranchAddress("eta",ee_eta);
    ee->SetBranchAddress("phi",ee_phi);
    ee->SetBranchAddress("isjet",ee_isjet);


    TFile* outputFile = new TFile(outputFileName, "RECREATE");
    std::cout << "pTcut      = " << pTcut << std::endl;
    std::cout << "cutDeltaR  = " << cutDeltaR << std::endl;
    std::cout << "cutDeltapT = " << cutDeltapT << std::endl;

    std::cout << "entering event loop" << std::endl;
    Long64_t entries = ggHiNtuplizerTree->GetEntries();
    std::cout << "number of entries = " << entries << std::endl;

    for(Long64_t jj = 0; jj < entries; ++jj)
    {
        if (jj % 100000 == 0)  {
            std::cout << "current entry = " <<jj<<" out of "<<entries<<" : "<<std::setprecision(2)<<(double)jj/entries*100<<" %"<<std::endl;
        }

        ggHiNtuplizerTree->GetEntry(jj);
        eb->GetEntry(jj);
        ee->GetEntry(jj);
    } // exited event loop

    outputFile->Write();
    outputFile->Close();
    file->Close();
}

void gedPhotonMacros_PfCand(const char* inputfileName, const char* outputFileName)
{
    TFile *file = TFile::Open(inputfileName);

    TTree* ggHiNtuplizerTree = (TTree*)file->Get("ggHiNtuplizer/EventTree");
    TTree* pfTree = (TTree*)file->Get("pfcandAnalyzer/pfTree");

    // OLD RECO photons
    Int_t nPho;
    std::vector<float>* phoEt=0;
    std::vector<float>* phoEta=0;
    std::vector<float>* phoPhi=0;
    std::vector<float>* pho_ecalClusterIsoR4=0;
    std::vector<float>* pho_hcalRechitIsoR4=0;
    std::vector<float>* pho_trackIsoR4PtCut20=0;

    ggHiNtuplizerTree->SetBranchAddress("nPho",&nPho);
    ggHiNtuplizerTree->SetBranchAddress("phoEt",&phoEt);
    ggHiNtuplizerTree->SetBranchAddress("phoEta",&phoEta);
    ggHiNtuplizerTree->SetBranchAddress("phoPhi",&phoPhi);
    ggHiNtuplizerTree->SetBranchAddress("pho_ecalClusterIsoR4" ,&pho_ecalClusterIsoR4);
    ggHiNtuplizerTree->SetBranchAddress("pho_hcalRechitIsoR4"  ,&pho_hcalRechitIsoR4);
    ggHiNtuplizerTree->SetBranchAddress("pho_trackIsoR4PtCut20",&pho_trackIsoR4PtCut20);
    // EB hits
    Int_t nPFpart;
//    Float_t eb_e[MAXHITS];
    Int_t pfId[MAXPFCANDS];
    Float_t pfPt[MAXPFCANDS];
    Float_t pfEta[MAXPFCANDS];
    Float_t pfPhi[MAXPFCANDS];
//    Bool_t  eb_jet[MAXHITS];

    pfTree->SetBranchAddress("nPFpart",&nPFpart);
//    eb->SetBranchAddress("e",eb_e);
    pfTree->SetBranchAddress("pfId",pfId);
    pfTree->SetBranchAddress("pfPt",pfPt);
    pfTree->SetBranchAddress("pfEta",pfEta);
    pfTree->SetBranchAddress("pfPhi",pfPhi);
//    eb->SetBranchAddress("isjet",eb_isjet);


    TFile* outputFile = new TFile(outputFileName, "RECREATE");
    outputFile->Write();
    outputFile->Close();
    file->Close();
}

void make2DRecHit(Long64_t entry, TTree* t)
{
    const char* selection = Form("Entry$ == %d", (int)entry);
    t->Draw("t1.run:t1.event:t1.lumis", selection, "goff");  // "goff" = graphics off
    int run = t->GetV1()[0];
    int event = t->GetV2()[0];
    int lumi = t->GetV3()[0];

    const char* title = Form("run = %d, event = %d, lumi = %d", run, event, lumi);
    const char* h2D_recHit_name = Form("h2D_recHit_entry%d", (int)entry);
    TH2D* h2D_recHit = new TH2D(h2D_recHit_name, title, 45, -4.5, 4.5, 32, -TMath::Pi(), TMath::Pi());

    const char* h2D_t1_name = Form("h2D_t1_entry%d", (int)entry);
    TH2D* h2D_t1 = new TH2D(h2D_t1_name, title, 45, -4.5, 4.5, 32, -TMath::Pi(), TMath::Pi());

    int selectedRows;
    t->Draw("eb.et:eb.eta:eb.phi", selection, "goff");  // "goff" = graphics off
    selectedRows = t->GetSelectedRows();
    for(int i=0; i<selectedRows; ++i)
    {
        h2D_recHit->Fill(t->GetV2()[i], t->GetV3()[i], t->GetV1()[i]);
    }

    t->Draw("ee.et:ee.eta:ee.phi", selection, "goff");  // "goff" = graphics off
    selectedRows = t->GetSelectedRows();
    for(int i=0; i<selectedRows; ++i)
    {
        h2D_recHit->Fill(t->GetV2()[i], t->GetV3()[i], t->GetV1()[i]);
    }

    t->Draw("t1.phoEt:t1.phoEta:t1.phoPhi", selection, "goff");  // "goff" = graphics off
    selectedRows = t->GetSelectedRows();
    for(int i=0; i<selectedRows; ++i)
    {
        h2D_t1->Fill(t->GetV2()[i], t->GetV3()[i], t->GetV1()[i]);
    }
}

int main(int argc, char** argv)
{
    const char* inputFileName;
    const char* outputFileName;

    if(argc == 1)
    {
        inputFileName = "/afs/cern.ch/work/y/ygo/public/ecalLocalReco/forest_AllQCDPhoton30_ecalGlobal_755p1.root";
        outputFileName = "RecHit1_global.root";
    }
    else if(argc == 2)
    {
        inputFileName = argv[1];
        outputFileName = "RecHit1_global.root";
    }
    else if(argc == 3)
    {
        inputFileName = argv[1];
        outputFileName = argv[2];
    }
    else
    {
        std::cout<<"wrong input"<<std::endl;
        return 0;
    }

    std::cout << "inputFileName  = " << inputFileName <<std::endl;
    std::cout << "outputFileName = " << outputFileName <<std::endl;
    std::cout<< "running gedPhotonMacros_RecHit()" << std::endl;
    //std::cout<< "running gedPhotonMacros_PfCand()" << std::endl;
    gedPhotonMacros_RecHit(inputFileName, outputFileName);
    //gedPhotonMacros_PfCand(inputFileName, outputFileName);
    return 0;
}


