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
const float cutDeltaR  = 0.05;
const float cutDeltapT = 5;
const int cutpfId = 4;  // 4=photon, https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideParticleFlow#ParticleTypes

class phoRecTree
{
   public:
   int nPhoton;
   vector<float> phoEt;
   vector<float> phoPhi;
   vector<float> phoEta;
   vector<float> phoSumIso;
   vector<float> phoEcalIso;
   vector<float> phoHcalIso;
   vector<float> phoTrkIso;
   vector<float> phoSigmaIEtaIEta;
   vector<float> phoHoverE;
   int nRechit;
   vector<float> ecalE;
   vector<float> ecalEt;
   vector<float> ecalEta;
   vector<float> ecalPhi;
   vector<float> deltaR;

   TTree *t;
   void init(){
      t = new TTree("t","");
      t->Branch("nPhoton",&nPhoton, "nPhoton/I");
      t->Branch("phoEt",&phoEt);
      t->Branch("phoEta",&phoEta);
      t->Branch("phoPhi",&phoPhi);
      t->Branch("phoSumIso",&phoSumIso);
      t->Branch("phoHcalIso",&phoHcalIso);
      t->Branch("phoEcalIso",&phoEcalIso);
      t->Branch("phoTrkIso",&phoTrkIso);
      t->Branch("phoSigmaIEtaIEta",&phoSigmaIEtaIEta);
      t->Branch("phoHoverE",&phoHoverE);
      t->Branch("nRechit",&nRechit, "nRechit/I");
      t->Branch("ecalE",&ecalE);
      t->Branch("ecalEt",&ecalEt);
      t->Branch("ecalEta",&ecalEta);
      t->Branch("ecalPhi",&ecalPhi);
      t->Branch("deltaR",&deltaR);
   }
   
   void clear(){
      nPhoton=0; 
      phoEt.clear();
      phoEta.clear();
      phoPhi.clear();
      phoSumIso.clear();
      phoEcalIso.clear();
      phoTrkIso.clear();
      phoHoverE.clear();
      phoSigmaIEtaIEta.clear();
      nRechit=0; 
      ecalE.clear();
      ecalEt.clear();
      ecalEta.clear();
      ecalPhi.clear();
      deltaR.clear();
   }
   
};
void rechitAroudPhoton(const char* inputfileName="", const char* outputFileName="")
{
    inputfileName="/afs/cern.ch/user/y/ygo/workspace/public/ecalLocalReco/forest_allqcdphoton30_multifit_ch2eError_755p1.root";
    //inputfileName="root://eoscm//eos/cms//store/group/phys_heavyions/velicanu/forest/Run2015E/HIExpressPhysics/Merged//HIForestExpress_run262548-v2.root";
    outputFileName="test.root";
    TFile *file = TFile::Open(inputfileName);
    TTree* ggHiNtuplizerTree = (TTree*)file->Get("ggHiNtuplizer/EventTree");
    TTree* ee = (TTree*)file->Get("rechitanalyzer/ee");
    TTree* eb = (TTree*)file->Get("rechitanalyzer/eb");
    TTree* masterTree = (TTree*)ggHiNtuplizerTree->Clone("t1");
    masterTree->AddFriend(ee,"ee");
    masterTree->AddFriend(eb,"eb");

    float ptCut = 30;
    float rad = 0.05; //0.1, 0.2

    // OLD RECO photons
    Int_t nPho;
    std::vector<float>* phoEt=0;
    std::vector<float>* phoEta=0;
    std::vector<float>* phoPhi=0;
    std::vector<float>* pho_ecalClusterIsoR4=0;
    std::vector<float>* pho_hcalRechitIsoR4=0;
    std::vector<float>* pho_trackIsoR4PtCut20=0;
    std::vector<float>* phoHoverE=0;
    std::vector<float>* phoSigmaIEtaIEta=0;
    std::vector<float>* pho_swissCrx=0;
    std::vector<float>* pho_seedTime=0;
    ggHiNtuplizerTree->SetBranchAddress("nPho",&nPho);
    ggHiNtuplizerTree->SetBranchAddress("phoEt",&phoEt);
    ggHiNtuplizerTree->SetBranchAddress("phoEta",&phoEta);
    ggHiNtuplizerTree->SetBranchAddress("phoPhi",&phoPhi);
    ggHiNtuplizerTree->SetBranchAddress("pho_ecalClusterIsoR4" ,&pho_ecalClusterIsoR4);
    ggHiNtuplizerTree->SetBranchAddress("pho_hcalRechitIsoR4"  ,&pho_hcalRechitIsoR4);
    ggHiNtuplizerTree->SetBranchAddress("pho_trackIsoR4PtCut20",&pho_trackIsoR4PtCut20);
    ggHiNtuplizerTree->SetBranchAddress("phoHoverE", &phoHoverE);
    ggHiNtuplizerTree->SetBranchAddress("phoSigmaIEtaIEta", &phoSigmaIEtaIEta);
    ggHiNtuplizerTree->SetBranchAddress("pho_swissCrx", &pho_swissCrx);
    ggHiNtuplizerTree->SetBranchAddress("pho_seedTime", &pho_seedTime);

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

    phoRecTree phorec;
    phorec.init();

    TFile* outputFile = new TFile(outputFileName, "RECREATE");
    std::cout << "pTcut      = " << ptCut << std::endl;
    std::cout << "cutDeltaR  = " << rad << std::endl;

    std::cout << "entering event loop" << std::endl;
    Long64_t entries = ggHiNtuplizerTree->GetEntries();
    std::cout << "number of entries = " << entries << std::endl;

    for(Long64_t jj = 0; jj < entries; ++jj)
    {
        if (jj % 10000 == 0)  {
            std::cout << "current entry = " <<jj<<" out of "<<entries<<" : "<<std::setprecision(2)<<(double)jj/entries*100<<" %"<<std::endl;
        }
        ggHiNtuplizerTree->GetEntry(jj);
        eb->GetEntry(jj);
        ee->GetEntry(jj);
        phorec.clear();
        int tmp_nPhoton=0;
        int tmp_nRechit=0;
        for(int ipho=0; ipho<nPho;ipho++){
            cout << "ss" << endl;
            if(phoEt->at(ipho)<ptCut) continue;
            if(phoHoverE->at(ipho)>0.1) continue;
            double sumIso = (pho_ecalClusterIsoR4->at(ipho)+pho_hcalRechitIsoR4->at(ipho)+pho_trackIsoR4PtCut20->at(ipho));
            //if( sumIso >1) continue;
            cout << "bs" << endl;
            //if(phoSigmaIEtaIEta->at(ipho)>0.01) continue;
            cout << "cs" << endl;
            if(abs(phoEta->at(ipho))<1.44 && (pho_swissCrx->at(ipho))<0.9 && abs(pho_seedTime->at(ipho)) < 3) continue;
            cout << "ds" << endl;

            cout << "es" << endl;
            bool passedDRee = 0; 
            cout << "ss" << endl;
            bool passedDReb = 0; 
            for(int iee=0;iee<ee_n;iee++){
                double deltaRtmp = getDR(phoEta->at(ipho), phoPhi->at(ipho), ee_eta[iee], ee_phi[iee] );
                passedDRee = (deltaRtmp < rad);
                if(!passedDRee) continue;
                phorec.ecalE.push_back(ee_e[iee]);
                phorec.ecalEt.push_back(ee_et[iee]);
                phorec.ecalEta.push_back(ee_eta[iee]);
                phorec.ecalPhi.push_back(ee_phi[iee]);
                phorec.deltaR.push_back(deltaRtmp);
                tmp_nRechit++;
            }//rechit loop
            for(int ieb=0;ieb<eb_n;ieb++){
                double deltaRtmp = getDR(phoEta->at(ipho), phoPhi->at(ipho), eb_eta[ieb], eb_phi[ieb] );
                passedDReb = (deltaRtmp < rad);
                if(!passedDReb) continue;
                phorec.ecalE.push_back(eb_e[ieb]);
                phorec.ecalEt.push_back(eb_et[ieb]);
                phorec.ecalEta.push_back(eb_eta[ieb]);
                phorec.ecalPhi.push_back(eb_phi[ieb]);
                phorec.deltaR.push_back(deltaRtmp);               
                tmp_nRechit++;
            }//rechit loop

            tmp_nPhoton++;
            phorec.phoEt.push_back(phoEt->at(ipho));
            phorec.phoEta.push_back(phoEta->at(ipho));
            phorec.phoPhi.push_back(phoPhi->at(ipho));
            phorec.phoSumIso.push_back(sumIso);
            phorec.phoEcalIso.push_back(pho_ecalClusterIsoR4->at(ipho));
            phorec.phoHcalIso.push_back(pho_hcalRechitIsoR4->at(ipho));
            phorec.phoTrkIso.push_back(pho_trackIsoR4PtCut20->at(ipho));
            phorec.phoSigmaIEtaIEta.push_back(phoSigmaIEtaIEta->at(ipho));
            phorec.phoHoverE.push_back(phoHoverE->at(ipho));
        }//photon loop 
        if(tmp_nPhoton>0 && tmp_){
            phorec.nPhoton=tmp_nPhoton;
            phorec.nRechit=tmp_nRechit;
            phorec.t->Fill();
        }
    } // exited event loop

    phorec.t->Write();
    outputFile->Close();
    file->Close();
}

