/*
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

#include "/net/hisrv0001/home/tatar/code/HIUtils/smallPhotonUtil.h"

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
    TFile *file = new TFile(inputfileName, "READ");

    TTree* ggHiNtuplizerTree = (TTree*)file->Get("ggHiNtuplizer/EventTree");
    TTree* ggHiNtuplizerGEDTree = (TTree*)file->Get("ggHiNtuplizerGED/EventTree");

    TTree* ee = (TTree*)file->Get("rechitanalyzer/ee");
    TTree* eb = (TTree*)file->Get("rechitanalyzer/eb");

    TTree* masterTree = (TTree*)ggHiNtuplizerTree->Clone("t1");
    masterTree->AddFriend(ggHiNtuplizerGEDTree,"t2");
    masterTree->AddFriend(ee,"ee");
    masterTree->AddFriend(eb,"eb");

//    TCut phoCandCut   = "sigmaIetaIeta<0.010";
//    TCut phoDecayCut  = "(sigmaIetaIeta>0.011) && (sigmaIetaIeta<0.017)";
//    TCut jetCut       =  Form("abs(eta)<%f && pt>%f", (float)1.6, (float)30 );

    float ptCut = 30;
    TCut selection_event = "Entry$<1000";
    TCut selection_event_pt_t2 = "Entry$<10000";
    TCut selection_event2 = "Entry$<5000";
    // First pool is.. it’s GED photon but not matched to old photon
    TCut selection_pt_t2 = Form("Max$(t2.phoEt)>%f", ptCut);
    TCut selection_pt_t2_no_t1 = Form("Max$(t2.phoEt)>%f && Max$(t1.phoEt)<(%f)", ptCut, ptCut);
    // Second pool is.. it’s old photon coordinate but not matched to GED photon
    TCut selection_pt_t1 = Form("Max$(t1.phoEt)>%f", ptCut);
//    TCut selection_pt_t1_no_t2 = Form("Max$(t2.phoEt)>%f && Max$(t1.phoEt)<(%f-2)", ptCut, ptCut);
    TCut selection_pt_t1_no_t2 = Form("Max$(t1.phoEt)>%f && Max$(t2.phoEt)<(%f)", ptCut, ptCut);
    // Third pool is… it’s old photon coordinate and matched to GED photon ( or the order of algorithm can be changed)
    TCut selection_pt_t1_t2 = Form("Max$(t1.phoEt)>%f && Max$(t2.phoEt)>%f", ptCut, ptCut);

    TCut selection_recHit_n = "eb.n>0 || ee.n>0";

    std::cout<<"Scan started"<<std::endl;

    masterTree->SetScanField(0);
    std::cout<<"First pool is.. it's GED photon but not matched to old photon"<<std::endl;
    std::cout<<selection_pt_t2.GetTitle()<<std::endl;
    masterTree->Scan("run:event:lumis:Entry$:t1.nPho:t1.phoEt:t1.phoEta:t1.phoPhi:t2.nPho:t2.phoEt:t2.phoEta:t2.phoPhi",
                                                    (selection_event_pt_t2 && selection_pt_t2).GetTitle());
    masterTree->Scan("eb.n:eb.e:eb.et:eb.eta:eb.phi:eb.perp:eb.isjet", (selection_event_pt_t2 && selection_pt_t2).GetTitle());

    std::cout<<"Second pool is.. it's old photon coordinate but not matched to GED photon"<<std::endl;
    std::cout<<selection_pt_t1.GetTitle()<<std::endl;
    masterTree->Scan("run:event:lumis:Entry$:t1.nPho:t1.phoEt:t1.phoEta:t1.phoPhi:t2.nPho:t2.phoEt:t2.phoEta:t2.phoPhi",
                                                    (selection_event && selection_pt_t1).GetTitle());
    masterTree->Scan("eb.n:eb.e:eb.et:eb.eta:eb.phi:eb.perp:eb.isjet", (selection_event && selection_pt_t1).GetTitle());

    std::cout<<"Third pool is.. it's old photon coordinate and matched to GED photon ( or the order of algorithm can be changed)"<<std::endl;
    std::cout<<selection_pt_t1_t2.GetTitle()<<std::endl;
    masterTree->Scan("run:event:lumis:Entry$:t1.nPho:t1.phoEt:t1.phoEta:t1.phoPhi:t2.nPho:t2.phoEt:t2.phoEta:t2.phoPhi",
                                                    (selection_event2 && selection_pt_t1_t2).GetTitle());
    masterTree->Scan("eb.n:eb.e:eb.et:eb.eta:eb.phi:eb.perp:eb.isjet", (selection_event2 && selection_pt_t1_t2).GetTitle());


    std::cout<<(selection_pt_t2).GetTitle()<<std::endl;
    Long64_t selected_pt_t2 = masterTree->GetEntries((selection_pt_t2).GetTitle());
    std::cout<<"selected_pt_t2 = "<<selected_pt_t2<<std::endl;

    std::cout<<(selection_pt_t2 && selection_recHit_n).GetTitle()<<std::endl;
    Long64_t selected_t2 = masterTree->GetEntries((selection_pt_t2 && selection_recHit_n).GetTitle());
    std::cout<<"selected_t2 = "<<selected_t2<<std::endl;

    std::cout<<(selection_pt_t1).GetTitle()<<std::endl;
    Long64_t selected_pt_t1 = masterTree->GetEntries((selection_pt_t1).GetTitle());
    std::cout<<"selected_pt_t1 = "<<selected_pt_t1<<std::endl;

    std::cout<<(selection_pt_t1 && selection_recHit_n).GetTitle()<<std::endl;
    Long64_t selected_t1 = masterTree->GetEntries((selection_pt_t1 && selection_recHit_n).GetTitle());
    std::cout<<"selected_t1 = "<<selected_t1<<std::endl;

    std::cout<<(selection_pt_t2_no_t1 && selection_recHit_n).GetTitle()<<std::endl;
    Long64_t selected_t2_no_t1 = masterTree->GetEntries((selection_pt_t2_no_t1 && selection_recHit_n).GetTitle());
    std::cout<<"selected_t2_not_t1 = "<<selected_t2_no_t1<<std::endl;

    std::cout<<(selection_pt_t1_no_t2 && selection_recHit_n).GetTitle()<<std::endl;
    Long64_t selected_t1_no_t2 = masterTree->GetEntries((selection_pt_t1_no_t2 && selection_recHit_n).GetTitle());
    std::cout<<"selected_t1_no_t2 = "<<selected_t1_no_t2<<std::endl;

    std::cout<<(selection_pt_t1_t2 && selection_recHit_n).GetTitle()<<std::endl;
    Long64_t selected_t1_t2 = masterTree->GetEntries((selection_pt_t1_t2 && selection_recHit_n).GetTitle());
    std::cout<<"selected_t1_t2 = "<<selected_t1_t2<<std::endl;

    std::cout<<"Scan finished"<<std::endl;

    TFile* outputFile = new TFile(outputFileName, "RECREATE");

    // selection_pt_t1_no_t2
    make2DRecHit(6,masterTree);
    make2DRecHit(838,masterTree);
/*
    ***********************************************************************************************************
    *    Row   * Instance *      eb.n *      eb.e *     eb.et *    eb.eta *    eb.phi *   eb.perp *  eb.isjet *
    ***********************************************************************************************************
    *        6 *        0 *         1 * 108.00002 * 107.95961 * 0.0273604 * -1.423725 * 129.30691 *         0 *
    *      142 *        0 *         1 * 335.88659 * 150.39498 * -1.442281 * -2.976555 * 130.32800 *         0 *
    *      304 *        0 *         2 * 34.614543 * 25.794479 * 0.8050490 * -2.695347 * 129.97431 *         0 *
    *      327 *        0 *         1 * 44.568309 * 41.792812 * 0.3624596 * -3.061554 * 129.71817 *         1 *
    *      334 *        0 *         1 * 2007.2428 * 1615.7272 * 0.6828121 * -0.166210 * 129.84701 *         0 *
    *      483 *        0 *         2 * 59.067058 * 55.891361 * -0.335526 * -0.554044 * 129.34542 *         0 *
    *      483 *        1 *         2 * 1482.0382 * 812.57971 * 1.2087075 * -1.771346 * 130.34278 *         0 *
    *      838 *        0 *         1 * 69.748603 * 33.150783 * 1.3750121 * -2.884834 * 130.30023 *         1 *
    *      957 *        0 *         1 * 94.040130 * 83.812927 * 0.4891222 * 0.5488857 * 129.30520 *         0 *
    ***********************************************************************************************************
*/
    // selection_pt_t1_t2
    make2DRecHit(327,masterTree);
    make2DRecHit(1867,masterTree);


    outputFile->Write();
    outputFile->Close();
    file->Close();
}

void gedPhotonMacros_RecHit2(const char* inputfileName, const char* outputFileName)
{
    TFile *file = new TFile(inputfileName, "READ");

    TTree* ggHiNtuplizerTree = (TTree*)file->Get("ggHiNtuplizer/EventTree");
    TTree* ggHiNtuplizerGEDTree = (TTree*)file->Get("ggHiNtuplizerGED/EventTree");

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

    // GED RECO photons
    Int_t nPhoGED;
    std::vector<float>* phoEtGED=0;
    std::vector<float>* phoEtaGED=0;
    std::vector<float>* phoPhiGED=0;
    std::vector<float>* pho_ecalClusterIsoR4GED=0;
    std::vector<float>* pho_hcalRechitIsoR4GED=0;
    std::vector<float>* pho_trackIsoR4PtCut20GED=0;

    ggHiNtuplizerGEDTree->SetBranchAddress("nPho",&nPhoGED);
    ggHiNtuplizerGEDTree->SetBranchAddress("phoEt",&phoEtGED);
    ggHiNtuplizerGEDTree->SetBranchAddress("phoEta",&phoEtaGED);
    ggHiNtuplizerGEDTree->SetBranchAddress("phoPhi",&phoPhiGED);
    ggHiNtuplizerGEDTree->SetBranchAddress("pho_ecalClusterIsoR4" ,&pho_ecalClusterIsoR4GED);
    ggHiNtuplizerGEDTree->SetBranchAddress("pho_hcalRechitIsoR4"  ,&pho_hcalRechitIsoR4GED);
    ggHiNtuplizerGEDTree->SetBranchAddress("pho_trackIsoR4PtCut20",&pho_trackIsoR4PtCut20GED);

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
//    TH1D* h[10];
//    h[0] = new TH1D("h_calIso_OLD", "cal Iso - OLDmatched2RecHit;cal Iso",100,-50,150);
//    h[1] = new TH1D("h_trkIso_OLD", "trk Iso - OLDmatched2RecHit;trk Iso",100,-50,150);
//
//    h[2] = new TH1D("h_calIso_OLD", "cal Iso - OLDnoMatch2RecHit;cal Iso",100,-50,150);
//    h[3] = new TH1D("h_trkIso_OLD", "trk Iso - OLDnoMatch2RecHit;trk Iso",100,-50,150);
//
//    h[4] = new TH1D("h_calIso_GED", "cal Iso - GEDmatched2RecHit;cal Iso",100,-50,150);
//    h[5] = new TH1D("h_trkIso_GED", "trk Iso - GEDmatched2RecHit;trk Iso",100,-50,150);
//
//    h[6] = new TH1D("h_calIso_GED", "cal Iso - GEDnoMatch2RecHit;cal Iso",100,-50,150);
//    h[7] = new TH1D("h_trkIso_GED", "trk Iso - GEDnoMatch2RecHit;trk Iso",100,-50,150);

    TH1D* h_calIso_matched2RecHit_OLD = new TH1D("h_calIso_matched2RecHit_OLD", "cal Iso - OLDmatched2RecHit;cal Iso",100,-50,150);
    TH1D* h_trkIso_matched2RecHit_OLD = new TH1D("h_trkIso_matched2RecHit_OLD", "trk Iso - OLDmatched2RecHit;trk Iso",100,-50,150);

    TH1D* h_calIso_noMatch2RecHit_OLD = new TH1D("h_calIso_noMatch2RecHit_OLD", "cal Iso - OLDnoMatch2RecHit;cal Iso",100,-50,150);
    TH1D* h_trkIso_noMatch2RecHit_OLD = new TH1D("h_trkIso_noMatch2RecHit_OLD", "trk Iso - OLDnoMatch2RecHit;trk Iso",100,-50,150);

    TH1D* h_calIso_matched2RecHit_GED = new TH1D("h_calIso_matched2RecHit_GED", "cal Iso - GEDmatched2RecHit;cal Iso",100,-50,150);
    TH1D* h_trkIso_matched2RecHit_GED = new TH1D("h_trkIso_matched2RecHit_GED", "trk Iso - GEDmatched2RecHit;trk Iso",100,-50,150);

    TH1D* h_calIso_noMatch2RecHit_GED = new TH1D("h_calIso_noMatch2RecHit_GED", "cal Iso - GEDnoMatch2RecHit;cal Iso",100,-50,150);
    TH1D* h_trkIso_noMatch2RecHit_GED = new TH1D("h_trkIso_noMatch2RecHit_GED", "trk Iso - GEDnoMatch2RecHit;trk Iso",100,-50,150);

    // match GED to OLD
    Long64_t nGED = 0;
    Long64_t nGED_matched2OLD = 0;
    Long64_t nGED_noMatch2OLD = 0;

    Long64_t nGED_noMatch2OLD_matched2RecHit = 0;
    Long64_t nGED_noMatch2OLD_noMatch2RecHit = 0;

    Long64_t nGED_matched2OLD_matched2RecHit = 0;
    Long64_t nGED_matched2OLD_noMatch2RecHit = 0;

    // match OLD to GED
    Long64_t nOLD = 0;
    Long64_t nOLD_matched2GED = 0;
    Long64_t nOLD_noMatch2GED = 0;

    Long64_t nOLD_noMatch2GED_matched2RecHit = 0;
    Long64_t nOLD_noMatch2GED_noMatch2RecHit = 0;

    Long64_t nOLD_matched2GED_matched2RecHit = 0;
    Long64_t nOLD_matched2GED_noMatch2RecHit = 0;

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
        ggHiNtuplizerGEDTree->GetEntry(jj);
        eb->GetEntry(jj);
        ee->GetEntry(jj);

        // First pool is.. it's GED photon but not matched to old photon
        for (int i=0; i < nPhoGED; ++i)
        {
            bool passed = (phoEtGED->at(i) > pTcut);
            if (!passed) continue;

            nGED++;
            bool matchedGED2OLD = false;
            bool noMatchGED2OLD = true;
            for (int j=0; j < nPho; ++j)
            {
                bool passed2 = (phoEt->at(j) > pTcut);
                if (!passed2) continue;

                matchedGED2OLD = matchParticlePair(phoEtaGED->at(i), phoPhiGED->at(i), phoEtGED->at(i),
                                                   phoEta->at(j),    phoPhi->at(j),    phoEt->at(j),
                                                   cutDeltaR, cutDeltapT);
//                bool matchParticlePair(Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2, Double_t pt1, Double_t  pt2,
//                                       double thresholdDeltaR = 0.2, double thresholdDeltapT = 5);
                if (matchedGED2OLD)
                {
                    nGED_matched2OLD++;
                    noMatchGED2OLD = false;
/*
                    std::cout<< "jj = " << jj <<std::endl;
                    std::cout<< "nGED_matched2OLD = " << nGED_matched2OLD <<std::endl;
                    std::cout<< "GED index = " << i <<std::endl;
                    std::cout<< "OLD index = " << j <<std::endl;
*/
                    break;
                }
            }

            if (noMatchGED2OLD)
            {
                nGED_noMatch2OLD++;

                bool matchedGED2RecHit = false;
                bool noMatchGED2RecHit = true;
                for (int j_eb=0; j_eb < eb_n; ++j_eb)
                {
                    matchedGED2RecHit = matchParticlePair(phoEtaGED->at(i), phoPhiGED->at(i), phoEtGED->at(i),
                                                          eb_eta[j_eb],    eb_phi[j_eb],      eb_et[j_eb],
                                                          cutDeltaR, cutDeltapT);
                    if (matchedGED2RecHit)
                    {
                        nGED_noMatch2OLD_matched2RecHit++;
                        noMatchGED2RecHit = false;
                        h_calIso_matched2RecHit_GED->Fill((pho_ecalClusterIsoR4GED->at(i)+pho_hcalRechitIsoR4GED->at(i)));
                        h_trkIso_matched2RecHit_GED->Fill(pho_trackIsoR4PtCut20GED->at(i));
                        break;
                    }
                }
                for (int j_ee=0; j_ee < ee_n; ++j_ee)
                {
                    matchedGED2RecHit = matchParticlePair(phoEtaGED->at(i), phoPhiGED->at(i), phoEtGED->at(i),
                                                          ee_eta[j_ee],     ee_phi[j_ee],     ee_et[j_ee],
                                                          cutDeltaR, cutDeltapT);
                    if (matchedGED2RecHit)
                    {
                        nGED_noMatch2OLD_matched2RecHit++;
                        noMatchGED2RecHit = false;
                        h_calIso_matched2RecHit_GED->Fill((pho_ecalClusterIsoR4GED->at(i)+pho_hcalRechitIsoR4GED->at(i)));
                        h_trkIso_matched2RecHit_GED->Fill(pho_trackIsoR4PtCut20GED->at(i));
                        break;
                    }
                }

                if (noMatchGED2RecHit)
                {
                    nGED_noMatch2OLD_noMatch2RecHit++;
                    h_calIso_noMatch2RecHit_GED->Fill((pho_ecalClusterIsoR4GED->at(i)+pho_hcalRechitIsoR4GED->at(i)));
                    h_trkIso_noMatch2RecHit_GED->Fill(pho_trackIsoR4PtCut20GED->at(i));
                }
            }
            else { // ! noMatchGED2OLD

                bool matchedGED2RecHit = false;
                bool noMatchGED2RecHit = true;
                for (int j_eb=0; j_eb < eb_n; ++j_eb)
                {
                    matchedGED2RecHit = matchParticlePair(phoEtaGED->at(i), phoPhiGED->at(i), phoEtGED->at(i),
                                                          eb_eta[j_eb],     eb_phi[j_eb],     eb_et[j_eb],
                                                          cutDeltaR, cutDeltapT);
                    if (matchedGED2RecHit)
                    {
                        nGED_matched2OLD_matched2RecHit++;
                        noMatchGED2RecHit = false;
                        h_calIso_matched2RecHit_GED->Fill((pho_ecalClusterIsoR4GED->at(i)+pho_hcalRechitIsoR4GED->at(i)));
                        h_trkIso_matched2RecHit_GED->Fill(pho_trackIsoR4PtCut20GED->at(i));
                        break;
                    }
                }
                for (int j_ee=0; j_ee < ee_n; ++j_ee)
                {
                    matchedGED2RecHit = matchParticlePair(phoEtaGED->at(i), phoPhiGED->at(i), phoEtGED->at(i),
                                                          ee_eta[j_ee],    ee_phi[j_ee],      ee_et[j_ee],
                                                          cutDeltaR, cutDeltapT);
                    if (matchedGED2RecHit)
                    {
                        nGED_noMatch2OLD_matched2RecHit++;
                        noMatchGED2RecHit = false;
                        h_calIso_matched2RecHit_GED->Fill((pho_ecalClusterIsoR4GED->at(i)+pho_hcalRechitIsoR4GED->at(i)));
                        h_trkIso_matched2RecHit_GED->Fill(pho_trackIsoR4PtCut20GED->at(i));
                        break;
                    }
                }

                if (noMatchGED2RecHit)
                {
                    nGED_matched2OLD_noMatch2RecHit++;
                    h_calIso_noMatch2RecHit_GED->Fill((pho_ecalClusterIsoR4GED->at(i)+pho_hcalRechitIsoR4GED->at(i)));
                    h_trkIso_noMatch2RecHit_GED->Fill(pho_trackIsoR4PtCut20GED->at(i));
                }
            }

        }

        // Second pool is.. it's old photon coordinate but not matched to GED photon
        for (int i=0; i < nPho; ++i)
        {
            bool passed = (phoEt->at(i) > pTcut);
            if (!passed) continue;

            nOLD++;
            bool matchedOLD2GED = false;
            bool noMatchOLD2GED = true;
            for (int j=0; j < nPhoGED; ++j)
            {
                bool passed2 = (phoEtGED->at(j) > pTcut);
                if (!passed2) continue;

                matchedOLD2GED = matchParticlePair(phoEta->at(i),    phoPhi->at(i),    phoEt->at(i),
                                                   phoEtaGED->at(j), phoPhiGED->at(j), phoEtGED->at(j),
                                                   cutDeltaR, cutDeltapT);
                if (matchedOLD2GED)
                {
                    nOLD_matched2GED++;
                    noMatchOLD2GED = false;
/*
                    std::cout<< "jj = " << jj <<std::endl;
                    std::cout<< "nOLD_matched2GED = " << nOLD_matched2GED <<std::endl;
                    std::cout<< "GED index = " << j <<std::endl;
                    std::cout<< "OLD index = " << i <<std::endl;
*/
                    break;
                }
            }

            if (noMatchOLD2GED)
            {
                nOLD_noMatch2GED++;

                bool matchedOLD2RecHit = false;
                bool noMatchOLD2RecHit = true;
                for (int j_eb=0; j_eb < eb_n; ++j_eb)
                {
                    matchedOLD2RecHit = matchParticlePair(phoEta->at(i), phoPhi->at(i), phoEt->at(i),
                                                          eb_eta[j_eb],  eb_phi[j_eb],  eb_et[j_eb],
                                                          cutDeltaR, cutDeltapT);
                    if (matchedOLD2RecHit)
                    {
                        nOLD_noMatch2GED_matched2RecHit++;
                        noMatchOLD2RecHit = false;
                        h_calIso_matched2RecHit_OLD->Fill((pho_ecalClusterIsoR4->at(i)+pho_hcalRechitIsoR4->at(i)));
                        h_trkIso_matched2RecHit_OLD->Fill(pho_trackIsoR4PtCut20->at(i));
                        break;
                    }
                }
                for (int j_ee=0; j_ee < ee_n; ++j_ee)
                {
                    matchedOLD2RecHit = matchParticlePair(phoEta->at(i), phoPhi->at(i), phoEt->at(i),
                                                          ee_eta[j_ee],  ee_phi[j_ee],  ee_et[j_ee],
                                                          cutDeltaR, cutDeltapT);
                    if (matchedOLD2RecHit)
                    {
                        nOLD_noMatch2GED_matched2RecHit++;
                        noMatchOLD2RecHit = false;
                        h_calIso_matched2RecHit_OLD->Fill((pho_ecalClusterIsoR4->at(i)+pho_hcalRechitIsoR4->at(i)));
                        h_trkIso_matched2RecHit_OLD->Fill(pho_trackIsoR4PtCut20->at(i));
                        break;
                    }
                }

                if (noMatchOLD2RecHit)
                {
                    nOLD_noMatch2GED_noMatch2RecHit++;
                    h_calIso_noMatch2RecHit_OLD->Fill((pho_ecalClusterIsoR4->at(i)+pho_hcalRechitIsoR4->at(i)));
                    h_trkIso_noMatch2RecHit_OLD->Fill(pho_trackIsoR4PtCut20->at(i));
                }
            }
            else { // ! noMatchOLD2GED

                bool matchedOLD2RecHit = false;
                bool noMatchOLD2RecHit = true;
                for (int j_eb=0; j_eb < eb_n; ++j_eb)
                {
                    matchedOLD2RecHit = matchParticlePair(phoEta->at(i), phoPhi->at(i), phoEt->at(i),
                                                          eb_eta[j_eb],  eb_phi[j_eb],  eb_et[j_eb],
                                                          cutDeltaR, cutDeltapT);
                    if (matchedOLD2RecHit)
                    {
                        nOLD_matched2GED_matched2RecHit++;
                        noMatchOLD2RecHit = false;
                        h_calIso_matched2RecHit_OLD->Fill((pho_ecalClusterIsoR4->at(i)+pho_hcalRechitIsoR4->at(i)));
                        h_trkIso_matched2RecHit_OLD->Fill(pho_trackIsoR4PtCut20->at(i));
                        break;
                    }
                }
                for (int j_ee=0; j_ee < ee_n; ++j_ee)
                {
                    matchedOLD2RecHit = matchParticlePair(phoEta->at(i), phoPhi->at(i), phoEt->at(i),
                                                          ee_eta[j_ee],  ee_phi[j_ee],  ee_et[j_ee],
                                                          cutDeltaR, cutDeltapT);
                    if (matchedOLD2RecHit)
                    {
                        nOLD_noMatch2GED_matched2RecHit++;
                        noMatchOLD2RecHit = false;
                        h_calIso_matched2RecHit_OLD->Fill((pho_ecalClusterIsoR4->at(i)+pho_hcalRechitIsoR4->at(i)));
                        h_trkIso_matched2RecHit_OLD->Fill(pho_trackIsoR4PtCut20->at(i));
                        break;
                    }
                }

                if (noMatchOLD2RecHit)
                {
                    nOLD_matched2GED_noMatch2RecHit++;
                    h_calIso_noMatch2RecHit_OLD->Fill((pho_ecalClusterIsoR4->at(i)+pho_hcalRechitIsoR4->at(i)));
                    h_trkIso_noMatch2RecHit_OLD->Fill(pho_trackIsoR4PtCut20->at(i));
                }
            }
        }
    } // exited event loop

    // match results for GED to OLD
    std::cout << "nGED      = " << nGED << std::endl;
    std::cout << "nGED_matched2OLD = " << nGED_matched2OLD << std::endl;
    std::cout << "nGED_noMatch2OLD = " << nGED_noMatch2OLD << std::endl;

    std::cout << "nGED_noMatch2OLD_matched2RecHit = " << nGED_noMatch2OLD_matched2RecHit << std::endl;
    std::cout << "nGED_noMatch2OLD_noMatch2RecHit = " << nGED_noMatch2OLD_noMatch2RecHit << std::endl;

    std::cout << "nGED_matched2OLD_matched2RecHit = " << nGED_matched2OLD_matched2RecHit << std::endl;
    std::cout << "nGED_matched2OLD_noMatch2RecHit = " << nGED_matched2OLD_noMatch2RecHit << std::endl;

    // match results for OLD to GED
    std::cout << "nOLD      = " << nOLD << std::endl;
    std::cout << "nOLD_matched2GED = " << nOLD_matched2GED << std::endl;
    std::cout << "nOLD_noMatch2GED = " << nOLD_noMatch2GED << std::endl;

    std::cout << "nOLD_noMatch2GED_matched2RecHit = " << nOLD_noMatch2GED_matched2RecHit << std::endl;
    std::cout << "nOLD_noMatch2GED_noMatch2RecHit = " << nOLD_noMatch2GED_noMatch2RecHit << std::endl;

    std::cout << "nOLD_matched2GED_matched2RecHit = " << nOLD_matched2GED_matched2RecHit << std::endl;
    std::cout << "nOLD_matched2GED_noMatch2RecHit = " << nOLD_matched2GED_noMatch2RecHit << std::endl;

    std::cout << "nGED_matched2OLD and nOLD_matched2GED do not have to be same. Because there may be "
            "more than one GED photon that match to the same OLD photon or vice versa." << std::endl;

    outputFile->Write();
    outputFile->Close();
    file->Close();
}

void gedPhotonMacros_PfCand(const char* inputfileName, const char* outputFileName)
{
    TFile *file = new TFile(inputfileName, "READ");

    TTree* ggHiNtuplizerTree = (TTree*)file->Get("ggHiNtuplizer/EventTree");
    TTree* ggHiNtuplizerGEDTree = (TTree*)file->Get("ggHiNtuplizerGED/EventTree");

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

    // GED RECO photons
    Int_t nPhoGED;
    std::vector<float>* phoEtGED=0;
    std::vector<float>* phoEtaGED=0;
    std::vector<float>* phoPhiGED=0;
    std::vector<float>* pho_ecalClusterIsoR4GED=0;
    std::vector<float>* pho_hcalRechitIsoR4GED=0;
    std::vector<float>* pho_trackIsoR4PtCut20GED=0;

    ggHiNtuplizerGEDTree->SetBranchAddress("nPho",&nPhoGED);
    ggHiNtuplizerGEDTree->SetBranchAddress("phoEt",&phoEtGED);
    ggHiNtuplizerGEDTree->SetBranchAddress("phoEta",&phoEtaGED);
    ggHiNtuplizerGEDTree->SetBranchAddress("phoPhi",&phoPhiGED);
    ggHiNtuplizerGEDTree->SetBranchAddress("pho_ecalClusterIsoR4" ,&pho_ecalClusterIsoR4GED);
    ggHiNtuplizerGEDTree->SetBranchAddress("pho_hcalRechitIsoR4"  ,&pho_hcalRechitIsoR4GED);
    ggHiNtuplizerGEDTree->SetBranchAddress("pho_trackIsoR4PtCut20",&pho_trackIsoR4PtCut20GED);

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

    TH1D* h_calIso_matched2pfCand_OLD = new TH1D("h_calIso_matched2pfCand_OLD", "cal Iso - OLDmatched2pfCand;cal Iso",100,-50,150);
    TH1D* h_trkIso_matched2pfCand_OLD = new TH1D("h_trkIso_matched2pfCand_OLD", "trk Iso - OLDmatched2pfCand;trk Iso",100,-50,150);

    TH1D* h_calIso_noMatch2pfCand_OLD = new TH1D("h_calIso_noMatch2pfCand_OLD", "cal Iso - OLDnoMatch2pfCand;cal Iso",100,-50,150);
    TH1D* h_trkIso_noMatch2pfCand_OLD = new TH1D("h_trkIso_noMatch2pfCand_OLD", "trk Iso - OLDnoMatch2pfCand;trk Iso",100,-50,150);

    TH1D* h_calIso_matched2pfCand_GED = new TH1D("h_calIso_matched2pfCand_GED", "cal Iso - GEDmatched2pfCand;cal Iso",100,-50,150);
    TH1D* h_trkIso_matched2pfCand_GED = new TH1D("h_trkIso_matched2pfCand_GED", "trk Iso - GEDmatched2pfCand;trk Iso",100,-50,150);

    TH1D* h_calIso_noMatch2pfCand_GED = new TH1D("h_calIso_noMatch2pfCand_GED", "cal Iso - GEDnoMatch2pfCand;cal Iso",100,-50,150);
    TH1D* h_trkIso_noMatch2pfCand_GED = new TH1D("h_trkIso_noMatch2pfCand_GED", "trk Iso - GEDnoMatch2pfCand;trk Iso",100,-50,150);

    // match GED to OLD
    Long64_t nGED = 0;
    Long64_t nGED_matched2OLD = 0;
    Long64_t nGED_noMatch2OLD = 0;

    Long64_t nGED_noMatch2OLD_matched2pfCand = 0;
    Long64_t nGED_noMatch2OLD_noMatch2pfCand = 0;

    Long64_t nGED_matched2OLD_matched2pfCand = 0;
    Long64_t nGED_matched2OLD_noMatch2pfCand = 0;

    // match OLD to GED
    Long64_t nOLD = 0;
    Long64_t nOLD_matched2GED = 0;
    Long64_t nOLD_noMatch2GED = 0;

    Long64_t nOLD_noMatch2GED_matched2pfCand = 0;
    Long64_t nOLD_noMatch2GED_noMatch2pfCand = 0;

    Long64_t nOLD_matched2GED_matched2pfCand = 0;
    Long64_t nOLD_matched2GED_noMatch2pfCand = 0;

    std::cout << "pTcut      = " << pTcut << std::endl;
    std::cout << "cutDeltaR  = " << cutDeltaR << std::endl;
    std::cout << "cutDeltapT = " << cutDeltapT << std::endl;
    std::cout << "cutpfId    = " << cutpfId << std::endl;   // 4=photon

    std::cout << "entering event loop" << std::endl;
    Long64_t entries = ggHiNtuplizerTree->GetEntries();
    std::cout << "number of entries = " << entries << std::endl;

    for(Long64_t jj = 0; jj < entries; ++jj)
    {
        if (jj % 100000 == 0)  {
            std::cout << "current entry = " <<jj<<" out of "<<entries<<" : "<<std::setprecision(2)<<(double)jj/entries*100<<" %"<<std::endl;
        }

        ggHiNtuplizerTree->GetEntry(jj);
        ggHiNtuplizerGEDTree->GetEntry(jj);
        pfTree->GetEntry(jj);

        // First pool is.. it's GED photon but not matched to old photon
        for (int i=0; i < nPhoGED; ++i)
        {
            bool passed = (phoEtGED->at(i) > pTcut);
            if (!passed) continue;

            nGED++;
            bool matchedGED2OLD = false;
            bool noMatchGED2OLD = true;
            for (int j=0; j < nPho; ++j)
            {
                bool passed2 = (phoEt->at(j) > pTcut);
                if (!passed2) continue;

                matchedGED2OLD = matchParticlePair(phoEtaGED->at(i), phoPhiGED->at(i), phoEtGED->at(i),
                                                   phoEta->at(j),    phoPhi->at(j),    phoEt->at(j),
                                                   cutDeltaR, cutDeltapT);
                if (matchedGED2OLD)
                {
                    nGED_matched2OLD++;
                    noMatchGED2OLD = false;
/*
                    std::cout<< "jj = " << jj <<std::endl;
                    std::cout<< "nGED_matched2OLD = " << nGED_matched2OLD <<std::endl;
                    std::cout<< "GED index = " << i <<std::endl;
                    std::cout<< "OLD index = " << j <<std::endl;
*/
                    break;
                }
            }

            if (noMatchGED2OLD)
            {
                nGED_noMatch2OLD++;

                bool matchedGED2pfCand = false;
                bool noMatchGED2pfCand = true;
                for (int j_pfCand=0; j_pfCand < nPFpart; ++j_pfCand)
                {
                    bool passed3 = (pfId[j_pfCand] == cutpfId);
                    if(!passed3) continue;

                    matchedGED2pfCand = matchParticlePair(phoEtaGED->at(i), phoPhiGED->at(i),  phoEtGED->at(i),
                                                          pfEta[j_pfCand],    pfPhi[j_pfCand], pfPt[j_pfCand],
                                                          cutDeltaR, cutDeltapT);
                    if (matchedGED2pfCand)
                    {
                        nGED_noMatch2OLD_matched2pfCand++;
                        noMatchGED2pfCand = false;
                        h_calIso_matched2pfCand_GED->Fill((pho_ecalClusterIsoR4GED->at(i)+pho_hcalRechitIsoR4GED->at(i)));
                        h_trkIso_matched2pfCand_GED->Fill(pho_trackIsoR4PtCut20GED->at(i));
                        break;
                    }
                }

                if (noMatchGED2pfCand)
                {
                    nGED_noMatch2OLD_noMatch2pfCand++;
                    h_calIso_noMatch2pfCand_GED->Fill((pho_ecalClusterIsoR4GED->at(i)+pho_hcalRechitIsoR4GED->at(i)));
                    h_trkIso_noMatch2pfCand_GED->Fill(pho_trackIsoR4PtCut20GED->at(i));
                }
            }
            else { // ! noMatchGED2OLD

                bool matchedGED2pfCand = false;
                bool noMatchGED2pfCand = true;
                for (int j_pfCand=0; j_pfCand < nPFpart; ++j_pfCand)
                {
                    bool passed3 = (pfId[j_pfCand] == cutpfId);
                    if(!passed3) continue;

                    matchedGED2pfCand = matchParticlePair(phoEtaGED->at(i), phoPhiGED->at(i),   phoEtGED->at(i),
                                                          pfEta[j_pfCand],     pfPhi[j_pfCand], pfPt[j_pfCand],
                                                          cutDeltaR, cutDeltapT);
                    if (matchedGED2pfCand)
                    {
                        nGED_matched2OLD_matched2pfCand++;
                        noMatchGED2pfCand = false;
                        h_calIso_matched2pfCand_GED->Fill((pho_ecalClusterIsoR4GED->at(i)+pho_hcalRechitIsoR4GED->at(i)));
                        h_trkIso_matched2pfCand_GED->Fill(pho_trackIsoR4PtCut20GED->at(i));
                        break;
                    }
                }

                if (noMatchGED2pfCand)
                {
                    nGED_matched2OLD_noMatch2pfCand++;
                    h_calIso_noMatch2pfCand_GED->Fill((pho_ecalClusterIsoR4GED->at(i)+pho_hcalRechitIsoR4GED->at(i)));
                    h_trkIso_noMatch2pfCand_GED->Fill(pho_trackIsoR4PtCut20GED->at(i));
                }
            }

        }

        // Second pool is.. it's old photon coordinate but not matched to GED photon
        for (int i=0; i < nPho; ++i)
        {
            bool passed = (phoEt->at(i) > pTcut);
            if (!passed) continue;

            nOLD++;
            bool matchedOLD2GED = false;
            bool noMatchOLD2GED = true;
            for (int j=0; j < nPhoGED; ++j)
            {
                bool passed2 = (phoEtGED->at(j) > pTcut);
                if (!passed2) continue;

                matchedOLD2GED = matchParticlePair(phoEta->at(i),    phoPhi->at(i),    phoEt->at(i),
                                                   phoEtaGED->at(j), phoPhiGED->at(j), phoEtGED->at(j),
                                                   cutDeltaR, cutDeltapT);
                if (matchedOLD2GED)
                {
                    nOLD_matched2GED++;
                    noMatchOLD2GED = false;
/*
                    std::cout<< "jj = " << jj <<std::endl;
                    std::cout<< "nOLD_matched2GED = " << nOLD_matched2GED <<std::endl;
                    std::cout<< "GED index = " << j <<std::endl;
                    std::cout<< "OLD index = " << i <<std::endl;
*/
                    break;
                }
            }

            if (noMatchOLD2GED)
            {
                nOLD_noMatch2GED++;

                bool matchedOLD2pfCand = false;
                bool noMatchOLD2pfCand = true;
                for (int j_pfCand=0; j_pfCand < nPFpart; ++j_pfCand)
                {
                    bool passed3 = (pfId[j_pfCand] == cutpfId);
                    if(!passed3) continue;

                    matchedOLD2pfCand = matchParticlePair(phoEta->at(i), phoPhi->at(i),      phoEt->at(i),
                                                          pfEta[j_pfCand],  pfPhi[j_pfCand], pfPt[j_pfCand],
                                                          cutDeltaR, cutDeltapT);
                    if (matchedOLD2pfCand)
                    {
                        nOLD_noMatch2GED_matched2pfCand++;
                        noMatchOLD2pfCand = false;
                        h_calIso_matched2pfCand_OLD->Fill((pho_ecalClusterIsoR4->at(i)+pho_hcalRechitIsoR4->at(i)));
                        h_trkIso_matched2pfCand_OLD->Fill(pho_trackIsoR4PtCut20->at(i));
                        break;
                    }
                }

                if (noMatchOLD2pfCand)
                {
                    nOLD_noMatch2GED_noMatch2pfCand++;
                    h_calIso_noMatch2pfCand_OLD->Fill((pho_ecalClusterIsoR4->at(i)+pho_hcalRechitIsoR4->at(i)));
                    h_trkIso_noMatch2pfCand_OLD->Fill(pho_trackIsoR4PtCut20->at(i));
                }
            }
            else { // ! noMatchOLD2GED

                bool matchedOLD2pfCand = false;
                bool noMatchOLD2pfCand = true;
                for (int j_pfCand=0; j_pfCand < nPFpart; ++j_pfCand)
                {
                    bool passed3 = (pfId[j_pfCand] == cutpfId);
                    if(!passed3) continue;

                    matchedOLD2pfCand = matchParticlePair(phoEta->at(i), phoPhi->at(i),      phoEt->at(i),
                                                          pfEta[j_pfCand],  pfPhi[j_pfCand], pfPt[j_pfCand],
                                                          cutDeltaR, cutDeltapT);
                    if (matchedOLD2pfCand)
                    {
                        nOLD_matched2GED_matched2pfCand++;
                        noMatchOLD2pfCand = false;
                        h_calIso_matched2pfCand_OLD->Fill((pho_ecalClusterIsoR4->at(i)+pho_hcalRechitIsoR4->at(i)));
                        h_trkIso_matched2pfCand_OLD->Fill(pho_trackIsoR4PtCut20->at(i));
                        break;
                    }
                }

                if (noMatchOLD2pfCand)
                {
                    nOLD_matched2GED_noMatch2pfCand++;
                    h_calIso_noMatch2pfCand_OLD->Fill((pho_ecalClusterIsoR4->at(i)+pho_hcalRechitIsoR4->at(i)));
                    h_trkIso_noMatch2pfCand_OLD->Fill(pho_trackIsoR4PtCut20->at(i));
                }
            }
        }
    } // exited event loop

    // match results for GED to OLD
    std::cout << "nGED      = " << nGED << std::endl;
    std::cout << "nGED_matched2OLD = " << nGED_matched2OLD << std::endl;
    std::cout << "nGED_noMatch2OLD = " << nGED_noMatch2OLD << std::endl;

    std::cout << "nGED_noMatch2OLD_matched2pfCand = " << nGED_noMatch2OLD_matched2pfCand << std::endl;
    std::cout << "nGED_noMatch2OLD_noMatch2pfCand = " << nGED_noMatch2OLD_noMatch2pfCand << std::endl;

    std::cout << "nGED_matched2OLD_matched2pfCand = " << nGED_matched2OLD_matched2pfCand << std::endl;
    std::cout << "nGED_matched2OLD_noMatch2pfCand = " << nGED_matched2OLD_noMatch2pfCand << std::endl;

    // match results for OLD to GED
    std::cout << "nOLD      = " << nOLD << std::endl;
    std::cout << "nOLD_matched2GED = " << nOLD_matched2GED << std::endl;
    std::cout << "nOLD_noMatch2GED = " << nOLD_noMatch2GED << std::endl;

    std::cout << "nOLD_noMatch2GED_matched2pfCand = " << nOLD_noMatch2GED_matched2pfCand << std::endl;
    std::cout << "nOLD_noMatch2GED_noMatch2pfCand = " << nOLD_noMatch2GED_noMatch2pfCand << std::endl;

    std::cout << "nOLD_matched2GED_matched2pfCand = " << nOLD_matched2GED_matched2pfCand << std::endl;
    std::cout << "nOLD_matched2GED_noMatch2pfCand = " << nOLD_matched2GED_noMatch2pfCand << std::endl;

    std::cout << "nGED_matched2OLD and nOLD_matched2GED do not have to be same. Because there may be "
            "more than one GED photon that match to the same OLD photon or vice versa." << std::endl;

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

    const char* h2D_t2_name = Form("h2D_t2_entry%d", (int)entry);
    TH2D* h2D_t2 = new TH2D(h2D_t2_name, title, 45, -4.5, 4.5, 32, -TMath::Pi(), TMath::Pi());

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

    t->Draw("t2.phoEt:t2.phoEta:t2.phoPhi", selection, "goff");  // "goff" = graphics off
    selectedRows = t->GetSelectedRows();
    for(int i=0; i<selectedRows; ++i)
    {
        h2D_t2->Fill(t->GetV2()[i], t->GetV3()[i], t->GetV1()[i]);
    }
}

int main(int argc, char** argv)
{
    const char* inputFileName;
    const char* outputFileName;

    if(argc == 1)
    {
        inputFileName = "/mnt/hadoop/cms/store/user/luck/HIMinBiasUPC/2011_MB_750_hiForest/0.root";
        outputFileName = "gedPhotonMacros_RecHit2_2011_MB_750_hiForest.root";
    }
    else if(argc == 2)
    {
        inputFileName = argv[1];
        outputFileName = "gedPhotonMacros_RecHit2_2011_MB_750_hiForest.root";
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

    //    gedPhotonMacros_RecHit(inputfileName);
//        std::cout<< "running gedPhotonMacros_RecHit2()" << std::endl;
//        gedPhotonMacros_RecHit2(inputFileName, outputFileName);

        std::cout<< "running gedPhotonMacros_PfCand()" << std::endl;
        gedPhotonMacros_PfCand(inputFileName, outputFileName);
        return 0;
}


