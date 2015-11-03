/*
 * macro to study different photon Reconstruction algorithms
 * modifiedy by Yeonju
 * */
#include <vector>
#ifdef __MAKECINT__
#pragma link C++ class vector<float>+;
#endif

#include "../gedPhotonUtility.h" 
static const long MAXTREESIZE = 10000000000;

void avgGenPtAroundGamma(const char* hiForestfileName="/afs/cern.ch/work/y/ygo/public/PFphoton/merged_EmEnrichedDijet30_PbPb_5TeV_forest_pfisoAdded_genThr500MeV.root", 
        TString outName="EmEnrichedDijet30", TString treePath="ggHiNtuplizer", float ptThr=20)
{

    TFile* inputFile = new TFile(hiForestfileName, "READ");
    std::cout << "input HiForest : " << inputFile->GetName() << std::endl;

    TString outputFileName = Form("skimFiles/jskim_%s_%s.root",outName.Data(),treePath.Data());
    TFile* fout = new TFile(outputFileName, "RECREATE");
    fout->cd();


    TTree* hiEvtAnalyzerTree = (TTree*)inputFile->Get("hiEvtAnalyzer/HiTree");
    TTree* HiGenParticleAnaTree = (TTree*)inputFile->Get("HiGenParticleAna/hi");
    TTree* ggHiNtuplizerTree = (TTree*)inputFile->Get(Form("%s/EventTree",treePath.Data()));

    Int_t hiBin;
    hiEvtAnalyzerTree->SetBranchAddress("hiBin", &hiBin);

    //hiGenTree    
    Int_t mult;
    std::vector<float>* pt=0;
    std::vector<float>* eta=0;
    std::vector<float>* phi=0;
    std::vector<int>*   pdg=0;
    HiGenParticleAnaTree->SetBranchAddress("mult", &mult);
    HiGenParticleAnaTree->SetBranchAddress("pt", &pt);
    HiGenParticleAnaTree->SetBranchAddress("eta", &eta);
    HiGenParticleAnaTree->SetBranchAddress("phi", &phi);
    HiGenParticleAnaTree->SetBranchAddress("pdg", &pdg);

    // mc Info.
    Int_t nMC;
    std::vector<float>* mcPt=0;
    std::vector<float>* mcEta=0;
    std::vector<float>* mcPhi=0;
    std::vector<int>*   mcPID=0;
    std::vector<int>*   mcMomPID=0;
    std::vector<int>*   mcStatus=0;
    ggHiNtuplizerTree->SetBranchAddress("nMC",&nMC);
    ggHiNtuplizerTree->SetBranchAddress("mcPt",&mcPt);
    ggHiNtuplizerTree->SetBranchAddress("mcEta",&mcEta);
    ggHiNtuplizerTree->SetBranchAddress("mcPhi",&mcPhi);
    ggHiNtuplizerTree->SetBranchAddress("mcPID",&mcPID);
    ggHiNtuplizerTree->SetBranchAddress("mcMomPID",&mcMomPID);
    ggHiNtuplizerTree->SetBranchAddress("mcStatus",&mcStatus);

    // RECO photons
    Int_t nPho;
    std::vector<float>* phoEt=0;
    std::vector<float>* phoEta=0;
    std::vector<float>* phoPhi=0;
    std::vector<float>* pho_ecalClusterIsoR2=0;
    std::vector<float>* pho_ecalClusterIsoR3=0;
    std::vector<float>* pho_ecalClusterIsoR4=0;
    std::vector<float>* pho_ecalClusterIsoR5=0;
    std::vector<float>* pho_hcalRechitIsoR2=0;
    std::vector<float>* pho_hcalRechitIsoR3=0;
    std::vector<float>* pho_hcalRechitIsoR4=0;
    std::vector<float>* pho_hcalRechitIsoR5=0;
    std::vector<float>* pho_trackIsoR2PtCut20=0;
    std::vector<float>* pho_trackIsoR3PtCut20=0;
    std::vector<float>* pho_trackIsoR4PtCut20=0;
    std::vector<float>* pho_trackIsoR5PtCut20=0;
    std::vector<float>* phoR9=0;
    std::vector<float>* phoHoverE=0;
    std::vector<float>* phoSigmaIEtaIEta=0;

    ggHiNtuplizerTree->SetBranchAddress("nPho",&nPho);
    ggHiNtuplizerTree->SetBranchAddress("phoEt",&phoEt);
    ggHiNtuplizerTree->SetBranchAddress("phoEta",&phoEta);
    ggHiNtuplizerTree->SetBranchAddress("phoPhi",&phoPhi);
    ggHiNtuplizerTree->SetBranchAddress("pho_ecalClusterIsoR2",&pho_ecalClusterIsoR2);
    ggHiNtuplizerTree->SetBranchAddress("pho_ecalClusterIsoR3",&pho_ecalClusterIsoR3);
    ggHiNtuplizerTree->SetBranchAddress("pho_ecalClusterIsoR4",&pho_ecalClusterIsoR4);
    ggHiNtuplizerTree->SetBranchAddress("pho_ecalClusterIsoR5",&pho_ecalClusterIsoR5);
    ggHiNtuplizerTree->SetBranchAddress("pho_hcalRechitIsoR2",&pho_hcalRechitIsoR2);
    ggHiNtuplizerTree->SetBranchAddress("pho_hcalRechitIsoR3",&pho_hcalRechitIsoR3);
    ggHiNtuplizerTree->SetBranchAddress("pho_hcalRechitIsoR4",&pho_hcalRechitIsoR4);
    ggHiNtuplizerTree->SetBranchAddress("pho_hcalRechitIsoR5",&pho_hcalRechitIsoR5);
    ggHiNtuplizerTree->SetBranchAddress("pho_trackIsoR2PtCut20",&pho_trackIsoR2PtCut20);
    ggHiNtuplizerTree->SetBranchAddress("pho_trackIsoR3PtCut20",&pho_trackIsoR3PtCut20);
    ggHiNtuplizerTree->SetBranchAddress("pho_trackIsoR4PtCut20",&pho_trackIsoR4PtCut20);
    ggHiNtuplizerTree->SetBranchAddress("pho_trackIsoR5PtCut20",&pho_trackIsoR5PtCut20);
    ggHiNtuplizerTree->SetBranchAddress("phoR9",&phoR9);
    ggHiNtuplizerTree->SetBranchAddress("phoHoverE",&phoHoverE);
    ggHiNtuplizerTree->SetBranchAddress("phoSigmaIEtaIEta",&phoSigmaIEtaIEta);

    // RECO photons (pfIso)
    std::vector<float>* pfcIso1=0; 
    std::vector<float>* pfcIso2=0; 
    std::vector<float>* pfcIso3=0; 
    std::vector<float>* pfcIso4=0; 
    std::vector<float>* pfcIso5=0; 
    std::vector<float>* pfpIso1=0; 
    std::vector<float>* pfpIso2=0; 
    std::vector<float>* pfpIso3=0; 
    std::vector<float>* pfpIso4=0; 
    std::vector<float>* pfpIso5=0; 
    std::vector<float>* pfnIso1=0; 
    std::vector<float>* pfnIso2=0; 
    std::vector<float>* pfnIso3=0; 
    std::vector<float>* pfnIso4=0; 
    std::vector<float>* pfnIso5=0; 
    std::vector<float>* pfcVsIso1=0; 
    std::vector<float>* pfcVsIso2=0; 
    std::vector<float>* pfcVsIso3=0; 
    std::vector<float>* pfcVsIso4=0; 
    std::vector<float>* pfcVsIso5=0; 
    std::vector<float>* pfpVsIso1=0; 
    std::vector<float>* pfpVsIso2=0; 
    std::vector<float>* pfpVsIso3=0; 
    std::vector<float>* pfpVsIso4=0; 
    std::vector<float>* pfpVsIso5=0; 
    std::vector<float>* pfnVsIso1=0; 
    std::vector<float>* pfnVsIso2=0; 
    std::vector<float>* pfnVsIso3=0; 
    std::vector<float>* pfnVsIso4=0; 
    std::vector<float>* pfnVsIso5=0;
    ggHiNtuplizerTree->SetBranchAddress("pfcIso1", &pfcIso1); 
    ggHiNtuplizerTree->SetBranchAddress("pfcIso2", &pfcIso2); 
    ggHiNtuplizerTree->SetBranchAddress("pfcIso3", &pfcIso3); 
    ggHiNtuplizerTree->SetBranchAddress("pfcIso4", &pfcIso4); 
    ggHiNtuplizerTree->SetBranchAddress("pfcIso5", &pfcIso5); 
    ggHiNtuplizerTree->SetBranchAddress("pfpIso1", &pfpIso1); 
    ggHiNtuplizerTree->SetBranchAddress("pfpIso2", &pfpIso2); 
    ggHiNtuplizerTree->SetBranchAddress("pfpIso3", &pfpIso3); 
    ggHiNtuplizerTree->SetBranchAddress("pfpIso4", &pfpIso4); 
    ggHiNtuplizerTree->SetBranchAddress("pfpIso5", &pfpIso5); 	
    ggHiNtuplizerTree->SetBranchAddress("pfnIso1", &pfnIso1); 
    ggHiNtuplizerTree->SetBranchAddress("pfnIso2", &pfnIso2); 
    ggHiNtuplizerTree->SetBranchAddress("pfnIso3", &pfnIso3); 
    ggHiNtuplizerTree->SetBranchAddress("pfnIso4", &pfnIso4); 
    ggHiNtuplizerTree->SetBranchAddress("pfnIso5", &pfnIso5); 

    ggHiNtuplizerTree->SetBranchAddress("pfcVsIso1", &pfcVsIso1); 
    ggHiNtuplizerTree->SetBranchAddress("pfcVsIso2", &pfcVsIso2); 
    ggHiNtuplizerTree->SetBranchAddress("pfcVsIso3", &pfcVsIso3); 
    ggHiNtuplizerTree->SetBranchAddress("pfcVsIso4", &pfcVsIso4); 
    ggHiNtuplizerTree->SetBranchAddress("pfcVsIso5", &pfcVsIso5); 
    ggHiNtuplizerTree->SetBranchAddress("pfpVsIso1", &pfpVsIso1); 
    ggHiNtuplizerTree->SetBranchAddress("pfpVsIso2", &pfpVsIso2); 
    ggHiNtuplizerTree->SetBranchAddress("pfpVsIso3", &pfpVsIso3); 
    ggHiNtuplizerTree->SetBranchAddress("pfpVsIso4", &pfpVsIso4); 
    ggHiNtuplizerTree->SetBranchAddress("pfpVsIso5", &pfpVsIso5); 	
    ggHiNtuplizerTree->SetBranchAddress("pfnVsIso1", &pfnVsIso1); 
    ggHiNtuplizerTree->SetBranchAddress("pfnVsIso2", &pfnVsIso2); 
    ggHiNtuplizerTree->SetBranchAddress("pfnVsIso3", &pfnVsIso3); 
    ggHiNtuplizerTree->SetBranchAddress("pfnVsIso4", &pfnVsIso4); 
    ggHiNtuplizerTree->SetBranchAddress("pfnVsIso5", &pfnVsIso5); 

    std::vector<float> *genSumIso2 = new std::vector<float>();
    std::vector<float> *genSumIso3 = new std::vector<float>();
    std::vector<float> *genSumIso4 = new std::vector<float>();
    std::vector<float> *genSumIso5 = new std::vector<float>();

    TTree* t_phogen;
    t_phogen= ggHiNtuplizerTree->CloneTree(0);
    t_phogen->SetName("t_phogen");
    t_phogen->SetMaxTreeSize(MAXTREESIZE);
    t_phogen->Branch("genSumIso2","vector<float>",&genSumIso2);
    t_phogen->Branch("genSumIso3","vector<float>",&genSumIso3);
    t_phogen->Branch("genSumIso4","vector<float>",&genSumIso4);
    t_phogen->Branch("genSumIso5","vector<float>",&genSumIso5);
    t_phogen->Branch("pt","vector<float>",&pt);
    t_phogen->Branch("hiBin",&hiBin,"hiBin/I");

    TH1::SetDefaultSumw2();
    std::cout << "entering event loop" << std::endl;
    Long64_t entries = ggHiNtuplizerTree->GetEntries();
    std::cout << "number of entries = " << entries << std::endl;

    std::clock_t    start_loop, end_loop;
    start_loop = std::clock();
    for(Long64_t jj = 0; jj < entries; ++jj)
    {
        if (jj % 100000 == 0)  {
            std::cout << "current entry = " <<jj<<" out of "<<entries<<" : "<<std::setprecision(2)<<(double)jj/entries*100<<" %"<<std::endl;
        }

        hiEvtAnalyzerTree->GetEntry(jj);
        HiGenParticleAnaTree->GetEntry(jj);
        ggHiNtuplizerTree->GetEntry(jj);

        for (int i=0; i < nPho; ++i)
        {
            bool passedPtThr = (phoEt->at(i) > ptThr);
            if(!passedPtThr ) continue;

            genSumIso2->clear();
            genSumIso3->clear();
            genSumIso4->clear();
            genSumIso5->clear();

            double tmp2=0;
            double tmp3=0;
            double tmp4=0;
            double tmp5=0;
          
           // find GEN photons that match to RECO photons
            bool passedGENselection;      // selections for GEN photon
            bool passedDR;
            double deltaRtmp;
            double deltaRMin = 999;
            int matchedIndex=-1;           // index of the matched GEN photon in  to this RECO photon

            for (int j=0; j<mult; ++j){
                if(TMath::Abs(pdg->at(j))!= PDG_PHOTON) continue;
//                bool passedGenPtThr = ( pt->at(j) > (ptThr-10.0) );
//                if(!passedGenPtThr) continue;
                deltaRtmp = getDR(phoEta->at(i), phoPhi->at(i), eta->at(j), phi->at(j));
                passedGENselection = (pdg->at(j) == PDG_PHOTON) && (deltaRtmp < cutdeltaR);
                if (!passedGENselection) continue;
                if (deltaRtmp < deltaRMin){
                    deltaRMin  = deltaRtmp;
                    matchedIndex = j;
                    //cout << "deltaRMin : " << deltaRMin << ", matchedIndex : " << matchedIndex << endl;
                }
            }
            if(matchedIndex==-1) continue;
            
            for (int j=0; j<mult; ++j){
                deltaRtmp = getDR(eta->at(matchedIndex), phi->at(matchedIndex), eta->at(j), phi->at(j));
#if 1 
                if(deltaRtmp >= 0.0 && deltaRtmp < 0.2) {
                    tmp2  += pt->at(j);
                    tmp3  += pt->at(j);
                    tmp4  += pt->at(j);
                    tmp5  += pt->at(j);
                } else if(deltaRtmp >= 0.2 && deltaRtmp < 0.3) {
                    tmp3  += pt->at(j);
                    tmp4  += pt->at(j);
                    tmp5  += pt->at(j);
                } else if(deltaRtmp >= 0.3 && deltaRtmp < 0.4) {
                    tmp4  += pt->at(j);
                    tmp5  += pt->at(j);
                } else if(deltaRtmp >= 0.4 && deltaRtmp < 0.5) {
                    tmp5  += pt->at(j);
                }
#endif
#if 0
                if(deltaRtmp >= 0.0 && deltaRtmp < 0.2) {
                    genSumIso2 += pt->at(j);
                    genSumIso3 += pt->at(j);
                    genSumIso4 += pt->at(j);
                    genSumIso5 += pt->at(j);
                } else if(deltaRtmp >= 0.2 && deltaRtmp < 0.3) {
                    genSumIso3 += pt->at(j);
                    genSumIso4 += pt->at(j);
                    genSumIso5 += pt->at(j);
                } else if(deltaRtmp >= 0.3 && deltaRtmp < 0.4) {
                    genSumIso4 += pt->at(j);
                    genSumIso5 += pt->at(j);
                } else if(deltaRtmp >= 0.4 && deltaRtmp < 0.5) {
                    genSumIso5 += pt->at(j);
                }
#endif

            }//mc loop
            genSumIso2->push_back(tmp2);
            genSumIso3->push_back(tmp3);
            genSumIso4->push_back(tmp4);
            genSumIso5->push_back(tmp5);
            
        }// photon loop
        t_phogen->Fill();
    } // exited event loop

    end_loop = std::clock();
    std::cout.precision(6);      // get back to default precision
    std::cout << "LOOP finished in             : " << (end_loop - start_loop) / (double)CLOCKS_PER_SEC << " seconds" << std::endl;
    std::cout << "exited event loop" << std::endl;

    t_phogen->Write();
    //=============================================================================================
    //=============================================================================================
    // Save the histograms 
    fout->Close();
    inputFile->Close();
}
