#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TMath.h"

#include "TH1F.h"
#include "TCanvas.h"

void derivePtHatWeights(std::string fList = "")
{
  std::string buffer;
  std::vector<std::string> listOfFiles;
  int nLines = 0;
  ifstream inFile(fList.data());

  std::cout << fList << std::endl;
  std::cout << inFile.is_open() << std::endl;

  if(!inFile.is_open()){
    std::cout << "Error opening file. Exiting." <<std::endl;
    return;
  }
  else{
    while(true){
      inFile >> buffer;
      if(inFile.eof()) break;
      listOfFiles.push_back(buffer);
      nLines++;
    }
  }

  std::cout << "FileList Loaded" << std::endl;

  TChain* evtChain_p = new TChain("HltTree");

  for(Int_t iter = 0; iter < (Int_t)(listOfFiles.size()); iter++){
    evtChain_p->Add(listOfFiles[iter].c_str());

    std::cout << listOfFiles[iter] << std::endl;
  }

  Int_t evt_ = 0;
  Int_t nref = 0;

  evtChain_p->SetBranchStatus("*", 0);
  evtChain_p->SetBranchStatus("Event", 1);
  evtChain_p->SetBranchAddress("Event", &evt_);

  evtChain_p->SetBranchStatus("HLT_AK5CaloJet15_v1", 1);
  evtChain_p->SetBranchAddress("HLT_AK5CaloJet15_v1", &nref);

  Int_t nEntries = evtChain_p->GetEntries();
  std::cout << nEntries << std::endl;

  Int_t nDup_ = 0;

  std::vector<Int_t>* evtVect_p = new std::vector<Int_t>;

  for(Int_t evtIter = 0; evtIter < nEntries; evtIter++){
    evtChain_p->GetEntry(evtIter);

    if(evtIter%10000 == 0) std::cout << evtIter << nDup_ << std::endl;

    for(Int_t iter = 0; iter < (Int_t)evtVect_p->size(); iter++){
      if(evt_ == evtVect_p->at(iter)){
	if(nref != 0) std::cout << "A MATCH: " << evt_ << ", " << nref << std::endl;
	nDup_++;
	break;
      }
    }

    evtVect_p->push_back(evt_);
  }


  return;
}
