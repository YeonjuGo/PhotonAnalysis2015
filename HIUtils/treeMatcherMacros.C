/*
 * treeMathcerMacros.C
 *
 */

#include <TFile.h>
#include <TTree.h>
#include <TCut.h>

#include <vector>
#include <string>
#include <iostream>
#include <iomanip>

#include "systemUtil.h"
#include "treeUtil.h"

void sortTrees(const char* filePath1, const char* filePath2);
int compareTreesBitWise_v3(const char* filePath1, const char* filePath2, const char* treePath1, const char* treePath2);
int compareTreesBitWise_v2(const char* filePath1, const char* filePath2, const char* treePath);
int compareTreesBitWise(const char* filePath1, const char* filePath2, const char* treePath);
void checkEvents(const char* filePath1, const char* filePath2);

/*
 * sort trees of the second file in the order of the first tree
 * */
void sortTrees(const char* filePath1, const char* filePath2)
{
    TFile* inputFile1= new TFile(filePath1, "READ");
    TFile* inputFile2= new TFile(filePath2, "READ");
    TFile* outputFile=new TFile("file2_sorted.root","RECREATE");

    TTree* tree2_1 = (TTree*)inputFile2->Get("ggHiNtuplizer/EventTree");
    TTree* tree2_2 = (TTree*)inputFile2->Get("ggHiNtuplizerGED/EventTree");
    TTree* tree2_3 = (TTree*)inputFile2->Get("akPu3CaloJetAnalyzer/t");
    TTree* tree2_4 = (TTree*)inputFile2->Get("akVs3CaloJetAnalyzer/t");
    TTree* tree2_5 = (TTree*)inputFile2->Get("akVs3PFJetAnalyzer/t");
    TTree* tree2_6 = (TTree*)inputFile2->Get("akPu3PFJetAnalyzer/t");
    TTree* tree2_7 = (TTree*)inputFile2->Get("akPu4CaloJetAnalyzer/t");
    TTree* tree2_8 = (TTree*)inputFile2->Get("akVs4CaloJetAnalyzer/t");
    TTree* tree2_9 = (TTree*)inputFile2->Get("akVs4PFJetAnalyzer/t");
    TTree* tree2_10 = (TTree*)inputFile2->Get("akPu4PFJetAnalyzer/t");

    TTree* tree2_1_sorted = tree2_1->CloneTree(0);
    TTree* tree2_2_sorted = tree2_2->CloneTree(0);
    TTree* tree2_3_sorted = tree2_3->CloneTree(0);
    TTree* tree2_4_sorted = tree2_4->CloneTree(0);
    TTree* tree2_5_sorted = tree2_5->CloneTree(0);
    TTree* tree2_6_sorted = tree2_6->CloneTree(0);
    TTree* tree2_7_sorted = tree2_7->CloneTree(0);
    TTree* tree2_8_sorted = tree2_8->CloneTree(0);
    TTree* tree2_9_sorted = tree2_9->CloneTree(0);
    TTree* tree2_10_sorted = tree2_10->CloneTree(0);

    tree2_1_sorted->SetName("ggHiNtuplizer_EventTree");
    tree2_2_sorted->SetName("ggHiNtuplizerGED_EventTree");
    tree2_3_sorted->SetName("akPu3CaloJetAnalyzer_t");
    tree2_4_sorted->SetName("akVs3CaloJetAnalyzer_t");
    tree2_5_sorted->SetName("akVs3PFJetAnalyzer_t");
    tree2_6_sorted->SetName("akPu3PFJetAnalyzer_t");
    tree2_7_sorted->SetName("akPu4CaloJetAnalyzer_t");
    tree2_8_sorted->SetName("akVs4CaloJetAnalyzer_t");
    tree2_9_sorted->SetName("akVs4PFJetAnalyzer_t");
    tree2_10_sorted->SetName("akPu4PFJetAnalyzer_t");

    // event info is stored in that tree
    TTree *tree1Event = (TTree*)inputFile1->Get("ggHiNtuplizer/EventTree");
    TTree *tree2Event = (TTree*)inputFile2->Get("ggHiNtuplizer/EventTree");

    Int_t run1, lumis1;
    Long64_t event1;
    tree1Event->SetBranchAddress("run", &run1);
    tree1Event->SetBranchAddress("event", &event1);
    tree1Event->SetBranchAddress("lumis", &lumis1);

    Long64_t entries = tree1Event->GetEntries();
    const char* selection;
    std::cout<< "Loop : tree1Event" <<std::endl;
    for (int j=0; j<entries; ++j)
    {
        if (j % 2000 == 0)  {
          std::cout << "current entry = " <<j<<" out of "<<entries<<" : "<<std::setprecision(2)<<(double)j/entries*100<<" %"<<std::endl;
        }

        tree1Event->GetEntry(j);

        selection = Form("run == %d && event == %d && lumis == %d", run1, event1, lumis1);
        /*
        int duplicates = tree1Event->GetEntries(selection);
        int duplicatesOther = tree2Event->GetEntries(selection);
        if (duplicates > 1)
        {
            std::cout<< "selection = " << selection << std::endl;
            std::cout<< "duplicates = " << duplicates << std::endl;
        }
        else if  (duplicates != 1)
        {
            std::cout<< "something is wrong : duplicates = " << duplicates << std::endl;
        }
        if (duplicatesOther != 1)
        {
            std::cout<< "selection = " << selection << std::endl;
            std::cout<< "something is wrong : duplicatesOther = " << duplicatesOther << std::endl;
        }
        */

        tree2Event->Draw("Entry$", selection);

        Long64_t selected = tree2Event->GetSelectedRows();
        Double_t* values  = tree2Event->GetV1();
        int j2 = -1;
        if (selected != 1)
        {
            std::cout<< "something is wrong : selected = " << selected << std::endl;
            std::cout<< "selection = " << selection << std::endl;
        }
        else
        {
            j2=(int)values[0];
        }

        // get this event from the second file
        tree2_1->GetEntry(j2);
        tree2_2->GetEntry(j2);
        tree2_3->GetEntry(j2);
        tree2_4->GetEntry(j2);
        tree2_5->GetEntry(j2);
        tree2_6->GetEntry(j2);
        tree2_7->GetEntry(j2);
        tree2_8->GetEntry(j2);
        tree2_9->GetEntry(j2);
        tree2_10->GetEntry(j2);

        tree2_1_sorted->Fill();
        tree2_2_sorted->Fill();
        tree2_3_sorted->Fill();
        tree2_4_sorted->Fill();
        tree2_5_sorted->Fill();
        tree2_6_sorted->Fill();
        tree2_7_sorted->Fill();
        tree2_8_sorted->Fill();
        tree2_9_sorted->Fill();
        tree2_10_sorted->Fill();
    }
    std::cout<< "Loop ENDED : tree1Event" <<std::endl;

    outputFile->Write();

    outputFile->Close();
    inputFile1->Close();
    inputFile2->Close();
}

int compareTreesBitWise_v3(const char* filePath1, const char* filePath2, const char* treePath1, const char* treePath2)
{
//    TFile* inputFile1= new TFile(filePath1, "READ");
//    TFile* inputFile2= new TFile(filePath2, "READ");
//
//    TTree *tree1 = (TTree*)inputFile1->Get(treePath1);
//    TTree *tree2 = (TTree*)inputFile2->Get(treePath2);

    std::string filepath1_str = filePath1;
    std::string filepath2_str = filePath2;

    int lenBranchNames = 7;
    const char* branchNames[] = {"nMC",
                                 "mcPt",
                                 "mcEta",
                                 "mcPID",
                                 "mcMomPID",
                                 "mcStatus",
                                 "mcCalIsoDR04"};

//    bool treesAreSame = compareTrees(tree1, tree2, 0);

    bool treesAreSame = compareTrees(new TFile(filepath1_str.c_str()), treePath1,
                                     new TFile(filepath2_str.c_str()), treePath2,
                                     0);

    return treesAreSame;
}


/*
 * compare by
 * sorting the second tree in the order of the first tree
 * */
int compareTreesBitWise_v2(const char* filePath1, const char* filePath2, const char* treePath)
{

    TFile* inputFile1= new TFile(filePath1, "READ");
    TFile* inputFile2= new TFile(filePath2, "READ");

    TTree *tree1 = (TTree*)inputFile1->Get(treePath);
    TTree *tree2 = (TTree*)inputFile2->Get(treePath);
    TFile* a=new TFile("out.root","recreate");
    TTree *tree2_sorted = tree2->CloneTree(0);

    // event info is stored in that tree
    TTree *tree1Event = (TTree*)inputFile1->Get("ggHiNtuplizer/EventTree");
    TTree *tree2Event = (TTree*)inputFile2->Get("ggHiNtuplizer/EventTree");

    // add as friend to retrieve event info
    tree1->AddFriend(tree1Event);
    tree2->AddFriend(tree2Event);

    Int_t run1, lumis1;
    Long64_t event1;
    tree1Event->SetBranchAddress("run", &run1);
    tree1Event->SetBranchAddress("event", &event1);
    tree1Event->SetBranchAddress("lumis", &lumis1);

    Long64_t entries = tree1Event->GetEntries();
    const char* selection;
    std::cout<< " tree1" <<std::endl;
    for (int j=0; j<entries; ++j)
    {
        if (j % 2000 == 0)  {
          std::cout << "current entry = " <<j<<" out of "<<entries<<" : "<<std::setprecision(2)<<(double)j/entries*100<<" %"<<std::endl;
        }

        tree1Event->GetEntry(j);

        selection = Form("run == %d && event == %d && lumis == %d", run1, event1, lumis1);
        int duplicates = tree1Event->GetEntries(selection);
        int duplicatesOther = tree2Event->GetEntries(selection);
        if (duplicates > 1)
        {
            std::cout<< "selection = " << selection << std::endl;
            std::cout<< "duplicates = " << duplicates << std::endl;
        }
        else if  (duplicates != 1)
        {
            std::cout<< "something is wrong : duplicates = " << duplicates << std::endl;
        }
        if (duplicatesOther != 1)
        {
            std::cout<< "selection = " << selection << std::endl;
            std::cout<< "something is wrong : duplicatesOther = " << duplicatesOther << std::endl;
        }

        tree2Event->Draw("Entry$", selection);

        Long64_t selected = tree2Event->GetSelectedRows();
        Double_t* values  = tree2Event->GetV1();
        int j2 = -1;
        if (selected != 1)
        {
            std::cout<< "selection = " << selection << std::endl;
            std::cout<< "something is wrong : selected = " << selected << std::endl;
        }
        else
        {
            j2=(int)values[0];
        }

        tree2->GetEntry(j2);
        tree2_sorted->Fill();
    }

    tree2_sorted->Write();
    a->Close();

    inputFile1->Close();
    inputFile2->Close();

    return 0;
}

int compareTreesBitWise(const char* filePath1, const char* filePath2, const char* treePath)
{
    int notSame = 0;

    TFile* inputFile1= new TFile(filePath1, "READ");
    TFile* inputFile2= new TFile(filePath2, "READ");

    TTree *tree1 = (TTree*)inputFile1->Get(treePath);
    TTree *tree2 = (TTree*)inputFile2->Get(treePath);

    // event info is stored in that tree
    TTree *tree1Event = (TTree*)inputFile1->Get("ggHiNtuplizer/EventTree");
    TTree *tree2Event = (TTree*)inputFile2->Get("ggHiNtuplizer/EventTree");

    // add as friend to retrieve event info
    tree1->AddFriend(tree1Event);
    tree2->AddFriend(tree2Event);

    Int_t run1, lumis1;
    Long64_t event1;
    tree1Event->SetBranchAddress("run", &run1);
    tree1Event->SetBranchAddress("event", &event1);
    tree1Event->SetBranchAddress("lumis", &lumis1);

    Long64_t entries = tree1Event->GetEntries();
    const char* selection_entry1;
    const char* selection_entry2;
    for (int j=0; j<entries; ++j)
    {
        if (j % 2000 == 0)  {
          std::cout << "current entry = " <<j<<" out of "<<entries<<" : "<<std::setprecision(2)<<(double)j/entries*100<<" %"<<std::endl;
        }

        tree1Event->GetEntry(j);

        selection_entry1 = Form("run == %d && event == %d && lumis == %d", run1, event1, lumis1);
        selection_entry2 = Form("run == %d && event == %d && lumis == %d", run1, event1, lumis1);
        // the same exists in the other tree. had checked that.

        bool eventHasSameContent = compareTrees(tree1, tree2, selection_entry1, selection_entry2);
        if (!eventHasSameContent)
        {
            std::cout<< "this event does not have same content" << std::endl;
            std::cout<< "selection = " << selection_entry1 << std::endl;

            std::cout << "current entry = " << j << std::endl;
            tree1->Scan("Entry$", selection_entry1);
            tree2->Scan("Entry$", selection_entry2);

            notSame++;
        }
    }

    inputFile1->Close();
    inputFile2->Close();

    return notSame;
}

/*
 * check duplicate events inside a file : duplicates
 * check whether event in one file is present also in the other file : duplicatesOther
 * */
void checkEvents(const char* filePath1, const char* filePath2)
{
       TFile* inputFile1= new TFile(filePath1, "READ");
       TFile* inputFile2= new TFile(filePath2, "READ");

       TTree *tree1 = (TTree*)inputFile1->Get("ggHiNtuplizer/EventTree");
       TTree *tree2 = (TTree*)inputFile2->Get("ggHiNtuplizer/EventTree");

       Int_t run1, lumis1;
       Long64_t event1;
       tree1->SetBranchAddress("run", &run1);
       tree1->SetBranchAddress("event", &event1);
       tree1->SetBranchAddress("lumis", &lumis1);

       Int_t run2, lumis2;
       Long64_t event2;
       tree2->SetBranchAddress("run", &run2);
       tree2->SetBranchAddress("event", &event2);
       tree2->SetBranchAddress("lumis", &lumis2);

       Long64_t entries = tree1->GetEntries();
       const char* selection;
       std::cout<< " tree1" <<std::endl;
       for (int j=0; j<entries; ++j)
       {
           if (j % 2000 == 0)  {
             std::cout << "current entry = " <<j<<" out of "<<entries<<" : "<<std::setprecision(2)<<(double)j/entries*100<<" %"<<std::endl;
           }

           tree1->GetEntry(j);

           selection = Form("run == %d && event == %d && lumis == %d", run1, event1, lumis1);
           int duplicates = tree1->GetEntries(selection);
           int duplicatesOther = tree2->GetEntries(selection);
           if (duplicates > 1)
           {
               std::cout<< "selection = " << selection << std::endl;
               std::cout<< "duplicates = " << duplicates << std::endl;
           }
           else if  (duplicates != 1)
           {
               std::cout<< "something is wrong : duplicates = " << duplicates << std::endl;
           }
           if (duplicatesOther != 1)
           {
               std::cout<< "selection = " << selection << std::endl;
               std::cout<< "something is wrong : duplicatesOther = " << duplicatesOther << std::endl;
           }
       }

       entries = tree2->GetEntries();
       std::cout<< " tree2" <<std::endl;
       for (int j=0; j<entries; ++j)
       {
           if (j % 2000 == 0)  {
             std::cout << "current entry = " <<j<<" out of "<<entries<<" : "<<std::setprecision(2)<<(double)j/entries*100<<" %"<<std::endl;
           }

           tree2->GetEntry(j);

           selection = Form("run == %d && event == %d && lumis == %d", run2, event2, lumis2);
           int duplicates = tree2->GetEntries(selection);
           int duplicatesOther = tree1->GetEntries(selection);
           if (duplicates > 1)
           {
               std::cout<< "selection = " << selection << std::endl;
               std::cout<< "duplicates = " << duplicates << std::endl;
           }
           else if  (duplicates != 1)
           {
               std::cout<< "something is wrong : duplicates = " << duplicates << std::endl;
           }
           if (duplicatesOther != 1)
           {
               std::cout<< "selection = " << selection << std::endl;
               std::cout<< "something is wrong : duplicatesOther = " << duplicatesOther << std::endl;
           }
       }

       inputFile1->Close();
       inputFile2->Close();
}

int main(int argc, char** argv)
{
//    const char* dir1="/mnt/hadoop/cms/store/user/richard/reLinkingTest/RelValZEEMM_13_HI/RelValZEEMM_13_HI_CMSSW_7_6_0_pre4-76X_mcRun2_HeavyIon_v1-v1_HiForest/150924_200129/0000/";
//    const char* dir2="/mnt/hadoop/cms/store/user/richard/reLinkingTest/RelValZEEMM_13_HI/RelValZEEMM_13_HI_mnguyen-PFup_RelValZEEMM_13_HI_CMSSW_7_6_0_pre4-76X_mcRun2_HeavyIon_v1_HiForest/150924_204237/0000/";

    const char* dir1="/mnt/hadoop/cms/store/user/tatar/reLinkingTest/RelValZEEMM_13_HI/";
    const char* filePath1 = Form("%s/HiForest_merged_CMSSW_7_6_0_pre4.root",dir1);
    const char* filePath2 = Form("%s/HiForest_merged_mnguyen.root",dir1);


////    checkEvents(filePath1, filePath2);
//    const char* treePath = "ggHiNtuplizer/EventTree";
//    int numDifferent = compareTreesBitWise(filePath1, filePath2, treePath);
//    if(numDifferent == 0)
//    {
//        std::cout << "files have the same tree content" << std::endl;
//        std::cout << "tree = " << treePath << std::endl;
//    }
//    else
//    {
//        std::cout << "files have different tree content" << std::endl;
//        std::cout << "tree = "         << treePath << std::endl;
//        std::cout << "numDifferent = " << numDifferent << std::endl;
//    }

//    const char* treePath = "ggHiNtuplizer/EventTree";
//    compareTreesBitWise_v2(filePath1, filePath2, treePath);

//    sortTrees(filePath1, filePath2);

    int treesAreSame = compareTreesBitWise_v3(filePath1, "file2_sorted.root", "ggHiNtuplizer/EventTree", "ggHiNtuplizer_EventTree");
    std::cout << "treesAreSame = " << treesAreSame << std::endl;

    const char* filePath11 = Form("%s/HiForest_merged_CMSSW_7_6_0_pre4.root",dir1);
    int treesAreSame2 = compareTreesBitWise_v3(filePath11, "file2_sorted.root", "akPu3CaloJetAnalyzer/t", "akPu3CaloJetAnalyzer_t");
    std::cout << "treesAreSame2 = " << treesAreSame2 << std::endl;

    return 0;
}
