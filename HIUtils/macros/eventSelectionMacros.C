/*
 * eventSelectionMacros.C
 *
 * macro to study number of events passing event selection and trigger selection.
 */

#include <TFile.h>
#include <TTree.h>
#include <TString.h>

#include <iostream>
#include <string>
#include <iomanip>

#include "../EventMatchingCMS.h"

void eventSelectionMacros();
void eventSelectionMacros1(const char* fileName);
void eventSelectionMacros2(const char* fileName);

void eventSelectionMacros1(const char* fileName)
{
    TFile *file  = new TFile(fileName, "READ");

    TTree* t1 = (TTree*)file->Get("ggHiNtuplizer/EventTree");
    TTree* t2 = (TTree*)file->Get("hltanalysis/HltTree");
    TTree* t3 = (TTree*)file->Get("hiEvtAnalyzer/HiTree");
    TTree* t4 = (TTree*)file->Get("skimanalysis/HltTree");

    Long64_t entries_hltanalysis_HltTree      = t2->GetEntries();
    Long64_t HLT_HISinglePhoton30_v2_0        = t2->GetEntries("HLT_HISinglePhoton30_v2 == 0");
    Long64_t HLT_HISinglePhoton30_v2_1        = t2->GetEntries("HLT_HISinglePhoton30_v2 == 1");
    Long64_t HLT_HISinglePhoton30_v2_Prescl_0 = t2->GetEntries("HLT_HISinglePhoton30_v2_Prescl == 0");
    Long64_t HLT_HISinglePhoton30_v2_Prescl_1 = t2->GetEntries("HLT_HISinglePhoton30_v2_Prescl == 1");

    Long64_t entries_hiEvtAnalyzer_HiTree = t3->GetEntries();
    Long64_t vz_gt_15                     = t3->GetEntries("vz >= 15");
    Long64_t vz_lt_15                     = t3->GetEntries("vz <  15");

    Long64_t entries_skimanalysis_HltTree     = t4->GetEntries();
    Long64_t pcollisionEventSelection_0       = t4->GetEntries("pcollisionEventSelection == 0");
    Long64_t pcollisionEventSelection_1       = t4->GetEntries("pcollisionEventSelection == 1");
    Long64_t pHBHENoiseFilterResultProducer_0 = t4->GetEntries("pHBHENoiseFilterResultProducer == 0");
    Long64_t pHBHENoiseFilterResultProducer_1 = t4->GetEntries("pHBHENoiseFilterResultProducer == 1");
    Long64_t pprimaryVertexFilter_0           = t4->GetEntries("pprimaryVertexFilter == 0");
    Long64_t pprimaryVertexFilter_1           = t4->GetEntries("pprimaryVertexFilter == 1");

    std::cout << "entries_hltanalysis_HltTree      = " << entries_hltanalysis_HltTree << std::endl;
    std::cout << "HLT_HISinglePhoton30_v2_0        = " << HLT_HISinglePhoton30_v2_0 << std::endl;
    std::cout << "HLT_HISinglePhoton30_v2_1        = " << HLT_HISinglePhoton30_v2_1 << std::endl;
    std::cout << "HLT_HISinglePhoton30_v2_Prescl_0 = " << HLT_HISinglePhoton30_v2_Prescl_0 << std::endl;
    std::cout << "HLT_HISinglePhoton30_v2_Prescl_1 = " << HLT_HISinglePhoton30_v2_Prescl_1 << std::endl;

    std::cout << "entries_hiEvtAnalyzer_HiTree = " << entries_hiEvtAnalyzer_HiTree << std::endl;
    std::cout << "vz_gt_15                     = " << vz_gt_15 << std::endl;
    std::cout << "vz_lt_15                     = " << vz_lt_15 << std::endl;

    std::cout << "entries_skimanalysis_HltTree     = " << entries_skimanalysis_HltTree << std::endl;
    std::cout << "pcollisionEventSelection_0       = " << pcollisionEventSelection_0 << std::endl;
    std::cout << "pcollisionEventSelection_1       = " << pcollisionEventSelection_1 << std::endl;
    std::cout << "pHBHENoiseFilterResultProducer_0 = " << pHBHENoiseFilterResultProducer_0 << std::endl;
    std::cout << "pHBHENoiseFilterResultProducer_1 = " << pHBHENoiseFilterResultProducer_1 << std::endl;
    std::cout << "pprimaryVertexFilter_0           = " << pprimaryVertexFilter_0 << std::endl;
    std::cout << "pprimaryVertexFilter_1           = " << pprimaryVertexFilter_1 << std::endl;

    file->Close();
}

/*
 * look for duplicate events
 */
void eventSelectionMacros2(const char* fileName)
{
    TFile *file  = new TFile(fileName, "READ");

    TTree* t1 = (TTree*)file->Get("ggHiNtuplizer/EventTree");
    TTree* t2 = (TTree*)file->Get("hiEvtAnalyzer/HiTree");

    Int_t run, lumis;
    Long64_t event;
    t1->SetBranchAddress("run", &run);
    t1->SetBranchAddress("event", &event);
    t1->SetBranchAddress("lumis", &lumis);

    EventMatchingCMS* eventMatcher=new EventMatchingCMS();
    bool eventAdded;
    Long64_t duplicateEvents = 0;

    Long64_t entries = t1->GetEntries();
    const char* selection;
    std::cout<< "Loop : ggHiNtuplizer/EventTree" <<std::endl;
    for (int j=0; j<entries; ++j)
    {
        if (j % 100000 == 0)  {
            std::cout << "current entry = " <<j<<" out of "<<entries<<" : "<<std::setprecision(2)<<(double)j/entries*100<<" %"<<std::endl;
        }

        t1->GetEntry(j);

        eventAdded = eventMatcher->addEvent(event, lumis, run, j);
        if(!eventAdded) // this event is duplicate, skip this one.
        {
            duplicateEvents++;
            continue;
        }
    }
    std::cout<< "Loop ENDED : ggHiNtuplizer/EventTree" <<std::endl;

    std::cout << "entries         = " << entries << std::endl;
    std::cout << "duplicateEvents = " << duplicateEvents << std::endl;

    file->Close();
}


void eventSelectionMacros()
{
    const char* fileName = "/mnt/hadoop/cms/store/user/luck/HIMinBiasUPC/2011_MB_750_hiForest/0.root"; // NOTE : this is a big file. GetEntries() = 1662337
//    const char* fileName = "/mnt/hadoop/cms/store/user/tatar/HIHighPt/HiForest_HIHighPt_photon30_HIRun2011-v1.root"; // GetEntries() = 166993

    std::cout << "fileName  = " << fileName <<std::endl;
    eventSelectionMacros1(fileName);
    eventSelectionMacros2(fileName);
}

int main(int argc, char** argv)
{
    eventSelectionMacros();
}
