/*
 * macro to study different photon Reconstruction algorithms
 * modifiedy by Yeonju
 * */

#include "gedPhotonUtility.h" 

void gedPhoton_treeProducer(const char* hiForestfileName="HiForest.root", TString outName="AllQCDPhoton30", TString treePath="ggHiNtuplizer", float ptThr=30)
{
	//TString outSuffix = getCondSuffix(cond_);

	TFile* inputFile = new TFile(hiForestfileName, "READ");
	std::cout << "input HiForest : " << inputFile->GetName() << std::endl;

	TTree* hiEvtAnalyzerTree = (TTree*)inputFile->Get("hiEvtAnalyzer/HiTree");
	TTree* ggHiNtuplizerTree = (TTree*)inputFile->Get(Form("%s/EventTree",treePath.Data()));

	Int_t hiBin;
	hiEvtAnalyzerTree->SetBranchAddress("hiBin", &hiBin);

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


	TString outputFileName = Form("skimFiles/skim_gedPhoton_%s_%s_ptThr%d.root",outName.Data(),treePath.Data(),(int)ptThr);
	TFile* fout = new TFile(outputFileName, "RECREATE");
	fout->cd();
	TTree* genMatchedTree = ggHiNtuplizerTree->CloneTree(0);
	TTree* fakeTree = ggHiNtuplizerTree->CloneTree(0);

	Int_t hiBin_genMatched, hiBin_fake;	
	genMatchedTree->Branch("hiBin", hiBin_genMatched, "hiBin/I");	
	fakeTree->Branch("hiBin", hiBin_fake, "hiBin/I");

	// histos[4][5][5][3];
	// [4][][][] : eta cuts : no cut, Barrel, Endcap, Endcap2
	// [][5][][] : hiBin cuts : all centralities, 0-10%, 10-30%, 30-50%, 50-100%
	// [][][5][] : mcMomPID cuts : all gen matched photons, prompt photon, fragmentation photon, decay photon (bkg.), decay photon from pi0

	const int nEtaCut = 4;
	const int nCentCut = 5;
	const int nMomIdCut = 5;
	const int nRadius = 6;

	float eta_gt[nEtaCut] = {  0,            0, cutetaBarrel, cutetaEndCap};
	float eta_lt[nEtaCut] = {999, cutetaBarrel, cutetaEndCap,       999};
	int hiBin_gt[nCentCut] = {-999,  0, 20,  60, 100};
	int hiBin_lt[nCentCut] = { 999, 20, 60, 100, 999};
	int mcMomPID_gt[nMomIdCut] = {-999, 21, -999,  22, 110};
	int mcMomPID_lt[nMomIdCut] = { 999, 23,   22, 999, 112};

	TH1::SetDefaultSumw2();

	std::cout << "entering event loop" << std::endl;
	Long64_t entries = ggHiNtuplizerTree->GetEntries();
	std::cout << "number of entries = " << entries << std::endl;

	std::clock_t    start_loop, end_loop;
	start_loop = std::clock();
	for(Long64_t jj = 0; jj < 100; ++jj)
		//for(Long64_t jj = 0; jj < entries; ++jj)
	{
		if (jj % 100000 == 0)  {
			std::cout << "current entry = " <<jj<<" out of "<<entries<<" : "<<std::setprecision(2)<<(double)jj/entries*100<<" %"<<std::endl;
		}

		hiEvtAnalyzerTree->GetEntry(jj);
		ggHiNtuplizerTree->GetEntry(jj);

		int centFlag= -99; // 0: all centralities, 1: 0-10%, 2: 10-30%, 3: 30-50%, 4: 50-100%
		for(int icent=1; icent<nCentCut ; ++icent){
			if( (hiBin>= hiBin_gt[icent]) && (hiBin < hiBin_lt[icent]) ) {
				centFlag = icent;
				//cout << "hiBin : " << hiBin << ", centFlag : " << centFlag << endl;
				break;
			}	
		}

		for (int i=0; i < nPho; ++i)
		{
			//cout << "nPho : " << nPho << ", pt : " << phoEt->at(i) << endl;
			bool passedPtThr = (phoEt->at(i) > ptThr);
			if(!passedPtThr ) continue;

			bool passedHoe = ( phoHoverE->at(i) < 0.1 );
			bool passedSigma= ( phoSigmaIEtaIEta->at(i) < 0.01 );

			if( (cond_ == hoeC) && !passedHoe ) {
				//cout << "I'm hoeC and don't pass hoe cut with h/e : " << phoHoverE->at(i) << endl; 
				continue;
			}
			if( (cond_ == sigmaC) && !passedSigma ) continue;
			if( (cond_ == hoeAndSigmaC) && !(passedHoe&&passedSigma) ) continue;

			int etaFlag = -99;// 0: no cut, 1: Barrel, 2: Endcap(1.5-3), 3: Endcap2(>3)
			int momFlag = -99;// 0: all photons, 1:prompt photon, 2:fragmentation photon, 3:decay photon (bkg.), 4:decay photon from pi0
			int fakeFlag = -99;// 0: all, 1: matched, 2: fake(not matched to MC)

			for(int ieta=1; ieta<nEtaCut; ++ieta){
				if( (TMath::Abs(phoEta->at(i))>= eta_gt[ieta]) && (TMath::Abs(phoEta->at(i)) < eta_lt[ieta]) ) {
					etaFlag = ieta;
					//cout << "eta : " << phoEta->at(i) << ", etaFlag : " << etaFlag << endl;
				}
			}

			// find GEN photons that match to RECO photons
			bool passedGENselection;      // selections for GEN photon
			bool passedDR;
			double deltaRtmp;
			double deltaRMin = 999;
			int matchedIndex=-1;           // index of the matched GEN photon in  to this RECO photon

			for (int j=0; j<nMC; ++j)
			{
				if(TMath::Abs(mcPID->at(j))!= PDG_PHOTON) continue;
				deltaRtmp = getDR(phoEta->at(i), phoPhi->at(i), mcEta->at(j), mcPhi->at(j));
				//cout << "nMC : " << nMC << ", deltaRtmp : " << deltaRtmp << endl;
				passedGENselection = (mcPID->at(j) == PDG_PHOTON) && (mcStatus->at(j) == cutmcStatus) && (deltaRtmp < cutdeltaR);
				if (!passedGENselection) continue;
				//cout << "gen selection pass" << endl;

				// matched GEN photon is the one closest to the RECO photon
				if (deltaRtmp < deltaRMin)
				{
					deltaRMin  = deltaRtmp;
					matchedIndex = j;
					//cout << "deltaRMin : " << deltaRMin << ", matchedIndex : " << matchedIndex << endl;
				}
			}

			if (matchedIndex==-1) { // if no matched GEN photon, then this RECO photon is fake.
				hiBin_fake = hiBin;
				fakeTree->Fill();
			}
			else if(matchedIndex!=-1){ //matched
				hiBin_genMatched = hiBin;
				genMatchedTree->Fill();
			}

			// if you want to get No gen mathcing all photons, then use this. 
			// (No Gen Mathcing) all photons(matched+fake)
		}// photon loop
	} // exited event loop

	end_loop = std::clock();
	std::cout.precision(6);      // get back to default precision
	std::cout << "LOOP finished in             : " << (end_loop - start_loop) / (double)CLOCKS_PER_SEC << " seconds" << std::endl;
	std::cout << "exited event loop" << std::endl;

	//=============================================================================================
	//=============================================================================================
	// Save the histograms 
	genMatchedTree->Write();
	fakeTree->Write();
	fout->Write();
	fout->Close();
	inputFile->Close();
}
