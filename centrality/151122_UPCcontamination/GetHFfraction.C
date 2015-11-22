// Author Yeonju Go
#include "TFile.h"
#include "TTree.h"
#include "TLine.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TMath.h"
#include "TDirectory.h"
#include "TProfile.h"
#include "TPad.h"
#include "TLorentzVector.h"
#include "stdio.h"
#include <iostream>
#include "../../yjUtility.h"
#include "mbemulationanalysis/L1UpgradeTree.h"
#include "mbemulationanalysis/L1UpgradeTree.C"


const double HFin = 2.87;
const double HFout = 5.02;

void GetHFfraction(const char* fname="root://eoscms//eos/cms/store/group/phys_heavyions/chflores/Foresting_RunPrep2015/STARLIGHTProd/HiForest_Starlight_Merge19112015.root",
        //const char* fname="/mnt/hadoop/cms/store/user/dgulhan/mergedForest/HiForest_Centrality_Unpacker_Hydjet_Quenched_MinBias_5020GeV_750_RECODEBUG_v0/HiForest_Centrality_Unpacker_Hydjet_Quenched_MinBias_5020GeV_750_RECODEBUG_v0_merged_forest_1.root",
        float eThr = 3,
        int nTower= 1,
        bool doL1 = 0,
        TString cap="")
{
    TH1::SetDefaultSumw2();
    gStyle->SetOptStat(0);
    SetHistTitleStyle();
    SetyjPadStyle();
    TFile *fin = TFile::Open(fname);
    cout << "open file : " << fname << endl;
    TTree *t = (TTree*) fin -> Get("skimanalysis/HltTree");
    TTree *t_gen = (TTree*) fin -> Get("HiGenParticleAna/hi");
    //TTree *t_l1 = (TTree*) fin -> Get("EmulatorResults/L1UpgradeTree");
    //TTree *t_evt = (TTree*) fin -> Get("hiEvtAnalyzer/HiTree");
    //TTree *t = (TTree*) fin -> Get("hltanalysis/HltTree");
    //t-> AddFriend(t_evt);
    //t-> AddFriend(t_skim);
    t-> AddFriend(t_gen);
    //t-> AddFriend(t_l1);
    cout << "Energh Threshold : " << eThr << " GeV, # of towers for firing HF : " << nTower << endl;
    //////////////////////////////////////////////////////////////// 
    //////////////////////////////////////////////////////////////// 
    // GEN FRACTION
    cout << " :::::::::::::::::: GEN FRACTION :::::::::::::::::: " << endl;
    Int_t mult;
    std::vector<float>* pt=0;
    std::vector<float>* eta=0;
    std::vector<float>* phi=0;
    std::vector<int>*   pdg=0;
    t->SetBranchAddress("mult",&mult);
    t->SetBranchAddress("pt",&pt);
    t->SetBranchAddress("eta",&eta);
    t->SetBranchAddress("phi",&phi);
    t->SetBranchAddress("pdg",&pdg);

    int totEvt = t->GetEntries();
    int hfANDevt = 0;
    int hfORevt = 0;
    int hfXORevt = 0;
    for(int jentry=0; jentry<totEvt;jentry++){
        t->GetEntry(jentry);
        int nPosTower = 0;
        int nNegTower = 0;
        for(int i=0; i<mult;i++){
            TLorentzVector * v = new TLorentzVector();  
            float m = -1.;
            if(fabs(pdg->at(i)) == 111) m = 0.13498; //pion 0
            else if(fabs(pdg->at(i)) == 211) m = 0.13957; //pion +-
            else if(fabs(pdg->at(i)) == 321) m = 0.49368; //kaon +-
            else if(fabs(pdg->at(i)) == 130) m = 0.49761; //kaon 0 long
            else if(fabs(pdg->at(i)) == 2212) m = 0.93827; //proton
            else if(fabs(pdg->at(i)) == 2112) m = 0.93957; //neutron
            else if(fabs(pdg->at(i)) == 22) m = 0.; //photon
            else if(fabs(pdg->at(i)) == 11) m = 0.000511; //electron
            else if(fabs(pdg->at(i)) == 13) m = 0.10566; //muon
            else continue;
            v->SetPtEtaPhiM(pt->at(i),eta->at(i),phi->at(i),m);
            double e = v->Energy();
            if(e<eThr) continue;
            if( (eta->at(i)>HFin ) && ( (eta->at(i))<HFout ) ) nPosTower++;
            else if( (eta->at(i)>-HFout ) && ( (eta->at(i))<-HFin ) ) nNegTower++;
            else continue;
        }
        if(nPosTower>=nTower && nNegTower>=nTower) { hfANDevt++; hfORevt++; }
        else if( ( nPosTower>=nTower && nNegTower<nTower ) || ( nPosTower<nTower && nNegTower>=nTower )) { hfORevt++; hfXORevt++; }
    }

    cout << " HF AND Efficiency : \t" << hfANDevt << "\t/\t" << totEvt << "\t=\t"<< (double)hfANDevt/totEvt*100 << "\t%"<< endl;
    cout << " HF OR Efficiency : \t" << hfORevt << "\t/\t" << totEvt << "\t=\t"<< (double)hfORevt/totEvt*100 << "\t%"<< endl;
    cout << " HF XOR Efficiency : \t" << hfXORevt << "\t/\t" << totEvt << "\t=\t"<< (double)hfXORevt/totEvt*100 << "\t%"<< endl;

    //////////////////////////////////////////////////////////////// 
    //////////////////////////////////////////////////////////////// 
    // RECO FRACTION
    cout << " :::::::::::::::::: RECO FRACTION ::::::::::::::::: " << endl;
    hfANDevt = 0;
    hfORevt = 0;
    hfXORevt = 0;
   
    Int_t phfPosTowers, phfNegTowers;
    Int_t phfPosFilter, phfPosFilter2, phfPosFilter3, phfPosFilter4, phfPosFilter5;
    Int_t phfNegFilter, phfNegFilter2, phfNegFilter3, phfNegFilter4, phfNegFilter5;
    Int_t phfCoincFilter, phfCoincFilter2, phfCoincFilter3, phfCoincFilter4, phfCoincFilter5;
    t->SetBranchAddress("phfPosFilter",&phfPosFilter);
    t->SetBranchAddress("phfPosFilter2",&phfPosFilter2);
    t->SetBranchAddress("phfPosFilter3",&phfPosFilter3);
    t->SetBranchAddress("phfPosFilter4",&phfPosFilter4);
    t->SetBranchAddress("phfPosFilter5",&phfPosFilter5);
    t->SetBranchAddress("phfNegFilter",&phfNegFilter);
    t->SetBranchAddress("phfNegFilter2",&phfNegFilter2);
    t->SetBranchAddress("phfNegFilter3",&phfNegFilter3);
    t->SetBranchAddress("phfNegFilter4",&phfNegFilter4);
    t->SetBranchAddress("phfNegFilter5",&phfNegFilter5);
    t->SetBranchAddress("phfCoincFilter",&phfCoincFilter);
    t->SetBranchAddress("phfCoincFilter2",&phfCoincFilter2);
    t->SetBranchAddress("phfCoincFilter3",&phfCoincFilter3);
    t->SetBranchAddress("phfCoincFilter4",&phfCoincFilter4);
    t->SetBranchAddress("phfCoincFilter5",&phfCoincFilter5);

    if(nTower==1){
        hfANDevt = t->GetEntries(Form("phfCoincFilter==1"));
        hfORevt = t->GetEntries(Form("phfPosFilter==1 || phfNegFilter==1"));
        hfXORevt = t->GetEntries(Form("(phfPosFilter==1 || phfNegFilter==1) && phfCoincFilter==0"));
    } else {
        hfANDevt = t->GetEntries(Form("phfCoincFilter%d==1",nTower));
        hfORevt = t->GetEntries(Form("phfPosFilter%d==1 || phfNegFilter%d==1",nTower,nTower));
        hfXORevt = t->GetEntries(Form("(phfPosFilter%d==1 || phfNegFilter%d==1) && phfCoincFilter%d==0",nTower,nTower,nTower));
    } 

    cout << " HF AND Efficiency : \t" << hfANDevt << "\t/\t" << totEvt << "\t=\t"<< (double)hfANDevt/totEvt*100 << "\t%"<< endl;
    cout << " HF OR Efficiency : \t" << hfORevt << "\t/\t" << totEvt << "\t=\t"<< (double)hfORevt/totEvt*100 << "\t%"<< endl;
    cout << " HF XOR Efficiency : \t" << hfXORevt << "\t/\t" << totEvt << "\t=\t"<< (double)hfXORevt/totEvt*100 << "\t%"<< endl;

    //////////////////////////////////////////////////////////////// 
    //////////////////////////////////////////////////////////////// 
    // L1 FRACTION
    hfANDevt = 0;
    hfORevt = 0;
    hfXORevt = 0;
    if(doL1){
        cout << " :::::::::::::::::: L1 FRACTION :::::::::::::::::: " << endl;
        L1UpgradeTree * tl1 = new L1UpgradeTree(fname);
        tl1->Loop(99999999);
        tl1->GetEntries();
    }
}
