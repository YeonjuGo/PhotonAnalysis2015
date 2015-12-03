#include "TFile.h"
#include "TTree.h"
#include "TCut.h"
#include "TH1D.h"
#include "TString.h"
#include "TStyle.h"
#include "TROOT.h"
#include "stdio.h"
#include "../yjUtility.h"

using namespace std;

void rechitFailFraction(const char* fname="/afs/cern.ch/work/y/ygo/public/ecalLocalReco/forest_jet80_2011data_755p1_global.root", bool isMC=0); 
void GetRechitFailFraction(){
    rechitFailFraction("root://eoscms//eos/cms/store/group/phys_heavyions/velicanu/forest/Run2015E/HIExpressPhysics/Merged/HIForestExpress_run262620-v6.root");
    rechitFailFraction("root://eoscms//eos/cms/store/group/phys_heavyions/velicanu/forest/Run2015E/HIExpressPhysics/Merged/HIForestExpress_run262837.root");
    //rechitFailFraction("/mnt/hadoop/cms/store/user/dgulhan/mergedForest/HiForest_Centrality_Unpacker_Hydjet_Quenched_MinBias_5020GeV_750_RECODEBUG_v0/HiForest_Centrality_Unpacker_Hydjet_Quenched_MinBias_5020GeV_750_RECODEBUG_v0_merged_forest_1.root",1);

}

void rechitFailFraction(const char* fname, bool isMC)
{
    gStyle->SetOptStat(0);  //0: don’t show statistic
    TFile* f = TFile::Open(fname);
    cout << "Open File : " << f->GetName() << endl;
    if(f->IsZombie()) return;
    TTree* tee;
    TTree* teb;
    TTree* tskim;
    TTree* thlt;
    tee = (TTree*)f->Get(Form("rechitanalyzer/ee"));
    teb = (TTree*)f->Get(Form("rechitanalyzer/eb"));
    tskim = (TTree*)f->Get(Form("skimanalysis/HltTree"));
    thlt= (TTree*)f->Get(Form("hltanalysis/HltTree"));
    teb->AddFriend(tskim);
    tee->AddFriend(tskim);
    teb->AddFriend(thlt);
    tee->AddFriend(thlt);

    Long64_t total_eb, failed_eb, frac_eb;  
    Long64_t total_ee, failed_ee, frac_ee;

    TCut evtCut = "HLT_HIL1MinimumBiasHF1ANDExpress_v1 && phfCoincFilter3 && pprimaryVertexFilter";
    if(isMC) evtCut = "phfCoincFilter3 && pprimaryVertexFilter"; 
    cout << "Event Selection Cut :: " << evtCut << endl;
    total_eb = teb->Draw("chi2",evtCut, "goff");
    failed_eb = teb->Draw("chi2",evtCut && "chi2>63", "goff");
    frac_eb = (double)failed_eb/total_eb;     
    cout << "<<<<Barrel>>>> total : "<< total_eb << ", failed : " << failed_eb << ", fraction : " << frac_eb << endl;

    total_ee = tee->Draw("chi2",evtCut, "goff");
    failed_ee = tee->Draw("chi2",evtCut && "chi2>63", "goff");
    frac_ee = (double)failed_ee/total_ee; 
    cout << "<<<<Endcap>>>> total : "<< total_ee << ", failed : " << failed_ee << ", fraction : " << frac_ee << endl;

}
