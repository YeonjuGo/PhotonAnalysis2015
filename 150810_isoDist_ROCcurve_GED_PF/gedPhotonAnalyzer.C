/*
 * macro to study different photon Reconstruction algorithms
 * modifiedy by Yeonju
 * */

#include "gedPhotonUtility.h" 

void gedPhotonAnalyzer(const char* hiForestfileName="HiForest.root", TString outName="AllQCDPhoton30", TString treePath="ggHiNtuplizer", float ptThr=30, condition cond_ = noC, int isoBin = 200)
{
	TString outSuffix = getCondSuffix(cond_);

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

	// histograms for photons
	TH1D* h_phoEt[nEtaCut][nCentCut][nMomIdCut];
	TH1D* h_phoPhi[nEtaCut][nCentCut][nMomIdCut];
	TH1D* h_phoEta[nEtaCut][nCentCut][nMomIdCut];
	TH1D* h_phoR9[nEtaCut][nCentCut][nMomIdCut];
	TH1D* h_phoHoverE[nEtaCut][nCentCut][nMomIdCut];
	TH1D* h_phoSigmaIEtaIEta[nEtaCut][nCentCut][nMomIdCut];
	TH1D* h_pho_ecalClusterIso[nEtaCut][nCentCut][nMomIdCut][nRadius];
	TH1D* h_pho_hcalRechitIso[nEtaCut][nCentCut][nMomIdCut][nRadius];
	TH1D* h_pho_trackIsoPtCut20[nEtaCut][nCentCut][nMomIdCut][nRadius];
	TH1D* h_pho_sumIso[nEtaCut][nCentCut][nMomIdCut][nRadius];

	TH1D* h_phoEt_fake[nEtaCut][nCentCut];
	TH1D* h_phoPhi_fake[nEtaCut][nCentCut];
	TH1D* h_phoEta_fake[nEtaCut][nCentCut];
	TH1D* h_phoR9_fake[nEtaCut][nCentCut];
	TH1D* h_phoHoverE_fake[nEtaCut][nCentCut];
	TH1D* h_phoSigmaIEtaIEta_fake[nEtaCut][nCentCut];
	TH1D* h_pho_ecalClusterIso_fake[nEtaCut][nCentCut][nRadius];
	TH1D* h_pho_hcalRechitIso_fake[nEtaCut][nCentCut][nRadius];
	TH1D* h_pho_trackIsoPtCut20_fake[nEtaCut][nCentCut][nRadius];
	TH1D* h_pho_sumIso_fake[nEtaCut][nCentCut][nRadius];

	TH1D* h_pfcIso[nEtaCut][nCentCut][nMomIdCut][nRadius];
	TH1D* h_pfnIso[nEtaCut][nCentCut][nMomIdCut][nRadius];
	TH1D* h_pfpIso[nEtaCut][nCentCut][nMomIdCut][nRadius];
	TH1D* h_sum_pfIso[nEtaCut][nCentCut][nMomIdCut][nRadius];
	TH1D* h_pfcVsIso[nEtaCut][nCentCut][nMomIdCut][nRadius];
	TH1D* h_pfnVsIso[nEtaCut][nCentCut][nMomIdCut][nRadius];
	TH1D* h_pfpVsIso[nEtaCut][nCentCut][nMomIdCut][nRadius];
	TH1D* h_sum_pfVsIso[nEtaCut][nCentCut][nMomIdCut][nRadius];
	TH1D* h_pfcBKGIso[nEtaCut][nCentCut][nMomIdCut][nRadius];
	TH1D* h_pfnBKGIso[nEtaCut][nCentCut][nMomIdCut][nRadius];
	TH1D* h_pfpBKGIso[nEtaCut][nCentCut][nMomIdCut][nRadius];
	TH1D* h_sum_pfBKGIso[nEtaCut][nCentCut][nMomIdCut][nRadius];

	int nBin_iso=isoBin;
	for(int j=0; j<nEtaCut; ++j){
		if(j==1) nBin_iso=isoBin;
		else if(j==2) nBin_iso=(int)isoBin*0.7;
		else if(j==3) nBin_iso=(int)isoBin*0.5;
		for(int k=0; k<nCentCut; ++k){
			for(int m=0; m<nMomIdCut; ++m){
				TString name = Form("eta%d_cent%d_momId%d",j,k,m);
				h_phoEt[j][k][m] = new TH1D(Form("phoEt_%s",name.Data()), Form("%s;p_{T} (GeV)",name.Data()),75,0,150);
				h_phoPhi[j][k][m] = new TH1D(Form("phoPhi_%s",name.Data()), Form("%s;#phi",name.Data()),30,-TMath::Pi(),TMath::Pi());
				h_phoEta[j][k][m] = new TH1D(Form("phoEta_%s",name.Data()), Form("%s;#eta",name.Data()),30,-5,5);
				h_phoR9[j][k][m] = new TH1D(Form("phoR9_%s",name.Data()), Form("%s;R9",name.Data()),100,-5,5);
				h_phoHoverE[j][k][m] = new TH1D(Form("phoHoverE_%s",name.Data()), Form("%s;h/e",name.Data()),100,-5,5);
				h_phoSigmaIEtaIEta[j][k][m] = new TH1D(Form("phoSigmaIEtaIEta_%s",name.Data()), Form("%s;#sigma_{I#etaI#eta}",name.Data()),50,0,0.05);
				for(int iR=2; iR<nRadius; iR++){
					TString nameR = Form("R%d_eta%d_cent%d_momId%d",iR,j,k,m);
					h_pho_ecalClusterIso[j][k][m][iR] = new TH1D(Form("pho_ecalClusterIso%s",nameR.Data()), Form("%s;ecalClusterIsoR%d",nameR.Data(),iR),nBin_iso,-40,90);
					h_pho_hcalRechitIso[j][k][m][iR] = new TH1D(Form("pho_hcalRechitIso%s",nameR.Data()), Form("%s;hcalRechitIsoR%d",nameR.Data(),iR),nBin_iso,-40,90);
					h_pho_trackIsoPtCut20[j][k][m][iR] = new TH1D(Form("pho_trackIsoPtCut20%s",nameR.Data()), Form("%s;trackIsoR%dPtCut20",nameR.Data(),iR),nBin_iso,-20,100);
					h_pho_sumIso[j][k][m][iR] = new TH1D(Form("pho_sumIso%s",nameR.Data()), Form("%s;sumIsoR%d(ecal+hcal+tracker Iso)",nameR.Data(),iR),nBin_iso,-100,250);
					h_pfcIso[j][k][m][iR] = new TH1D(Form("pfcIso%s",nameR.Data()), Form("%s;pfcIso%d",nameR.Data(),iR),nBin_iso,0,120);
					h_pfnIso[j][k][m][iR] = new TH1D(Form("pfnIso%s",nameR.Data()), Form("%s;pfnIso%d",nameR.Data(),iR),nBin_iso,0,300);
					h_pfpIso[j][k][m][iR] = new TH1D(Form("pfpIso%s",nameR.Data()), Form("%s;pfpIso%d",nameR.Data(),iR),nBin_iso,0,350);
					h_sum_pfIso[j][k][m][iR] = new TH1D(Form("sum_pfIso%s",nameR.Data()), Form("%s;sum_pfIso%d(pfcIso+pfpIso+pfnIso)",nameR.Data(),iR),nBin_iso,-100,250);
					h_pfcVsIso[j][k][m][iR] = new TH1D(Form("pfcVsIso%s",nameR.Data()), Form("%s;pfcVsIso%d",nameR.Data(),iR),nBin_iso,-100,150);
					h_pfnVsIso[j][k][m][iR] = new TH1D(Form("pfnVsIso%s",nameR.Data()), Form("%s;pfnVsIso%d",nameR.Data(),iR),nBin_iso,-30,260);
					h_pfpVsIso[j][k][m][iR] = new TH1D(Form("pfpVsIso%s",nameR.Data()), Form("%s;pfpVsIso%d",nameR.Data(),iR),nBin_iso,-50,250);
					h_sum_pfVsIso[j][k][m][iR] = new TH1D(Form("sum_pfVsIso%s",nameR.Data()), Form("%s;sum_pfVsIso%d(pfcVsIso+pfpVsIso+pfnVsIso)",nameR.Data(),iR),nBin_iso,-100,250);
					h_pfcBKGIso[j][k][m][iR] = new TH1D(Form("pfcBKGIso%s",nameR.Data()), Form("%s;pfcIso%d-pfcVsIso%d",nameR.Data(),iR,iR),nBin_iso,-50,120);
					h_pfnBKGIso[j][k][m][iR] = new TH1D(Form("pfnBKGIso%s",nameR.Data()), Form("%s;pfnIso%d-pfnVsIso%d",nameR.Data(),iR,iR),nBin_iso,-50,300);
					h_pfpBKGIso[j][k][m][iR] = new TH1D(Form("pfpBKGIso%s",nameR.Data()), Form("%s;pfpIso%d-pfpVsIso%d",nameR.Data(),iR,iR),nBin_iso,-50,350);
					h_sum_pfBKGIso[j][k][m][iR] = new TH1D(Form("sum_pfBKGIso%s",nameR.Data()), Form("%s;(sum_pfIso%d) - (sum_pfVsIso%d)",nameR.Data(),iR,iR),nBin_iso,-100,250);	
				}
			}
			TString name = Form("eta%d_cent%d_fake",j,k);
			h_phoEt_fake[j][k] = new TH1D(Form("phoEt_%s",name.Data()), Form("%s;p_{T} (GeV)",name.Data()),75,0,150);
			h_phoPhi_fake[j][k] = new TH1D(Form("phoPhi_%s",name.Data()), Form("%s;#phi",name.Data()),30,-TMath::Pi(),TMath::Pi());
			h_phoEta_fake[j][k] = new TH1D(Form("phoEta_%s",name.Data()), Form("%s;#eta",name.Data()),30,-5,5);
			h_phoR9_fake[j][k] = new TH1D(Form("phoR9_%s",name.Data()), Form("%s;R9",name.Data()),100,-5,5);
			h_phoHoverE_fake[j][k] = new TH1D(Form("phoHoverE_%s",name.Data()), Form("%s;h/e",name.Data()),100,-5,5);
			h_phoSigmaIEtaIEta_fake[j][k] = new TH1D(Form("phoSigmaIEtaIEta_%s",name.Data()), Form("%s;#sigma_{I#etaI#eta}",name.Data()),50,0,0.05);
			for(int iR=2; iR<nRadius; iR++){
				TString nameR = Form("R%d_eta%d_cent%d_fake",iR,j,k);
				h_pho_ecalClusterIso_fake[j][k][iR] = new TH1D(Form("pho_ecalClusterIso%s",nameR.Data()), Form("%s;ecalClusterIsoR%d",nameR.Data(),iR),50,-40,90);
				h_pho_hcalRechitIso_fake[j][k][iR] = new TH1D(Form("pho_hcalRechitIso%s",nameR.Data()), Form("%s;hcalRechitIsoR%d",nameR.Data(),iR),50,-40,90);
				h_pho_trackIsoPtCut20_fake[j][k][iR] = new TH1D(Form("pho_trackIsoPtCut20%s",nameR.Data()), Form("%s;trackIsoR%dPtCut20",nameR.Data(),iR),50,-20,100);
				h_pho_sumIso_fake[j][k][iR] = new TH1D(Form("pho_sumIso%s",nameR.Data()), Form("%s;sumIsoR%d(ecal+hcal+tracker Iso)",nameR.Data(),iR),50,-20,300);
			}
		}
	}

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
#if 0
				h_phoEt_fake[etaFlag][centFlag]->Fill(phoEt->at(i));
				h_phoPhi_fake[etaFlag][centFlag]->Fill(phoPhi->at(i)); 
				h_phoEta_fake[etaFlag][centFlag]->Fill(phoEta->at(i));
				h_phoR9_fake[etaFlag][centFlag]->Fill(phoR9->at(i));
				h_phoHoverE_fake[etaFlag][centFlag]->Fill(phoHoverE->at(i));
				h_phoSigmaIEtaIEta_fake[etaFlag][centFlag]->Fill(phoSigmaIEtaIEta->at(i));
				h_pho_ecalClusterIso_fake[etaFlag][centFlag][2]->Fill(pho_ecalClusterIsoR2->at(i));
				h_pho_ecalClusterIso_fake[etaFlag][centFlag][3]->Fill(pho_ecalClusterIsoR3->at(i));
				h_pho_ecalClusterIso_fake[etaFlag][centFlag][4]->Fill(pho_ecalClusterIsoR4->at(i));
				h_pho_ecalClusterIso_fake[etaFlag][centFlag][5]->Fill(pho_ecalClusterIsoR5->at(i));
				h_pho_hcalRechitIso_fake[etaFlag][centFlag][2]->Fill(pho_hcalRechitIsoR2->at(i));
				h_pho_hcalRechitIso_fake[etaFlag][centFlag][3]->Fill(pho_hcalRechitIsoR3->at(i));
				h_pho_hcalRechitIso_fake[etaFlag][centFlag][4]->Fill(pho_hcalRechitIsoR4->at(i));
				h_pho_hcalRechitIso_fake[etaFlag][centFlag][5]->Fill(pho_hcalRechitIsoR5->at(i));
				h_pho_trackIsoPtCut20_fake[etaFlag][centFlag][2]->Fill(pho_trackIsoR2PtCut20->at(i));
				h_pho_trackIsoPtCut20_fake[etaFlag][centFlag][3]->Fill(pho_trackIsoR3PtCut20->at(i));
				h_pho_trackIsoPtCut20_fake[etaFlag][centFlag][4]->Fill(pho_trackIsoR4PtCut20->at(i));
				h_pho_trackIsoPtCut20_fake[etaFlag][centFlag][5]->Fill(pho_trackIsoR5PtCut20->at(i));
				h_pho_sumIso_fake[etaFlag][centFlag][2]->Fill(pho_ecalClusterIsoR2->at(i)+pho_hcalRechitIsoR2->at(i)+pho_trackIsoR2PtCut20->at(i));
				h_pho_sumIso_fake[etaFlag][centFlag][3]->Fill(pho_ecalClusterIsoR3->at(i)+pho_hcalRechitIsoR3->at(i)+pho_trackIsoR3PtCut20->at(i));
				h_pho_sumIso_fake[etaFlag][centFlag][4]->Fill(pho_ecalClusterIsoR4->at(i)+pho_hcalRechitIsoR4->at(i)+pho_trackIsoR4PtCut20->at(i));
				h_pho_sumIso_fake[etaFlag][centFlag][5]->Fill(pho_ecalClusterIsoR5->at(i)+pho_hcalRechitIsoR5->at(i)+pho_trackIsoR5PtCut20->at(i));
#endif
			}
			else if(matchedIndex!=-1){ //matched
				if( TMath::Abs(mcMomPID->at(matchedIndex)) == PDG_PHOTON ) momFlag = 1; //prompt photon
				else if ( TMath::Abs(mcMomPID->at(matchedIndex)) < 22) momFlag = 2; //fragmentation photon
				else if ( TMath::Abs(mcMomPID->at(matchedIndex)) > 22) momFlag = 3; //decay photon 
				bool isPi0 = 0;
				if (TMath::Abs(mcMomPID->at(matchedIndex)) == cutmcMomPID_pi0) isPi0 =1; 
				//cout << "mcMomPID : " << mcMomPID->at(matchedIndex) << ", momFlag : "<< momFlag << ", isPi0 : " << isPi0 << endl;

				for(int kk=0; kk<=(int)isPi0; kk++){
					int tmpMomFlag = momFlag+kk;
					h_phoEt[etaFlag][centFlag][momFlag+kk]->Fill(phoEt->at(i));
					h_phoPhi[etaFlag][centFlag][momFlag+kk]->Fill(phoPhi->at(i));
					h_phoEta[etaFlag][centFlag][momFlag+kk]->Fill(phoEta->at(i));
					h_phoR9[etaFlag][centFlag][momFlag+kk]->Fill(phoR9->at(i));
					h_phoHoverE[etaFlag][centFlag][momFlag+kk]->Fill(phoHoverE->at(i));
					h_phoSigmaIEtaIEta[etaFlag][centFlag][momFlag+kk]->Fill(phoSigmaIEtaIEta->at(i));
					h_pho_ecalClusterIso[etaFlag][centFlag][momFlag+kk][2]->Fill(pho_ecalClusterIsoR2->at(i));
					h_pho_ecalClusterIso[etaFlag][centFlag][momFlag+kk][3]->Fill(pho_ecalClusterIsoR3->at(i));
					h_pho_ecalClusterIso[etaFlag][centFlag][momFlag+kk][4]->Fill(pho_ecalClusterIsoR4->at(i));
					h_pho_ecalClusterIso[etaFlag][centFlag][momFlag+kk][5]->Fill(pho_ecalClusterIsoR5->at(i));
					h_pho_hcalRechitIso[etaFlag][centFlag][momFlag+kk][2]->Fill(pho_hcalRechitIsoR2->at(i));
					h_pho_hcalRechitIso[etaFlag][centFlag][momFlag+kk][3]->Fill(pho_hcalRechitIsoR3->at(i));
					h_pho_hcalRechitIso[etaFlag][centFlag][momFlag+kk][4]->Fill(pho_hcalRechitIsoR4->at(i));
					h_pho_hcalRechitIso[etaFlag][centFlag][momFlag+kk][5]->Fill(pho_hcalRechitIsoR5->at(i));
					h_pho_trackIsoPtCut20[etaFlag][centFlag][momFlag+kk][2]->Fill(pho_trackIsoR2PtCut20->at(i));
					h_pho_trackIsoPtCut20[etaFlag][centFlag][momFlag+kk][3]->Fill(pho_trackIsoR3PtCut20->at(i));
					h_pho_trackIsoPtCut20[etaFlag][centFlag][momFlag+kk][4]->Fill(pho_trackIsoR4PtCut20->at(i));
					h_pho_trackIsoPtCut20[etaFlag][centFlag][momFlag+kk][5]->Fill(pho_trackIsoR5PtCut20->at(i));
					h_pho_sumIso[etaFlag][centFlag][momFlag+kk][2]->Fill(pho_ecalClusterIsoR2->at(i)+pho_hcalRechitIsoR2->at(i)+pho_trackIsoR2PtCut20->at(i));
					h_pho_sumIso[etaFlag][centFlag][momFlag+kk][3]->Fill(pho_ecalClusterIsoR3->at(i)+pho_hcalRechitIsoR3->at(i)+pho_trackIsoR3PtCut20->at(i));
					h_pho_sumIso[etaFlag][centFlag][momFlag+kk][4]->Fill(pho_ecalClusterIsoR4->at(i)+pho_hcalRechitIsoR4->at(i)+pho_trackIsoR4PtCut20->at(i));
					h_pho_sumIso[etaFlag][centFlag][momFlag+kk][5]->Fill(pho_ecalClusterIsoR5->at(i)+pho_hcalRechitIsoR5->at(i)+pho_trackIsoR5PtCut20->at(i));
				}//for loop pi0
#if 1
	
				// if you want to get gen mathced all photons, then use this. 
				// (Gen Mathced) all photons(matched), not separated by momId.
				h_phoEt[etaFlag][centFlag][0]->Fill(phoEt->at(i));
				h_phoPhi[etaFlag][centFlag][0]->Fill(phoPhi->at(i));
				h_phoEta[etaFlag][centFlag][0]->Fill(phoEta->at(i));
				h_phoR9[etaFlag][centFlag][0]->Fill(phoR9->at(i));
				h_phoHoverE[etaFlag][centFlag][0]->Fill(phoHoverE->at(i));
				h_phoSigmaIEtaIEta[etaFlag][centFlag][0]->Fill(phoSigmaIEtaIEta->at(i));
				h_pho_ecalClusterIso[etaFlag][centFlag][0][2]->Fill(pho_ecalClusterIsoR2->at(i));
				h_pho_ecalClusterIso[etaFlag][centFlag][0][3]->Fill(pho_ecalClusterIsoR3->at(i));
				h_pho_ecalClusterIso[etaFlag][centFlag][0][4]->Fill(pho_ecalClusterIsoR4->at(i));
				h_pho_ecalClusterIso[etaFlag][centFlag][0][5]->Fill(pho_ecalClusterIsoR5->at(i));
				h_pho_hcalRechitIso[etaFlag][centFlag][0][2]->Fill(pho_hcalRechitIsoR2->at(i));
				h_pho_hcalRechitIso[etaFlag][centFlag][0][3]->Fill(pho_hcalRechitIsoR3->at(i));
				h_pho_hcalRechitIso[etaFlag][centFlag][0][4]->Fill(pho_hcalRechitIsoR4->at(i));
				h_pho_hcalRechitIso[etaFlag][centFlag][0][5]->Fill(pho_hcalRechitIsoR5->at(i));
				h_pho_trackIsoPtCut20[etaFlag][centFlag][0][2]->Fill(pho_trackIsoR2PtCut20->at(i));
				h_pho_trackIsoPtCut20[etaFlag][centFlag][0][3]->Fill(pho_trackIsoR3PtCut20->at(i));
				h_pho_trackIsoPtCut20[etaFlag][centFlag][0][4]->Fill(pho_trackIsoR4PtCut20->at(i));
				h_pho_trackIsoPtCut20[etaFlag][centFlag][0][5]->Fill(pho_trackIsoR5PtCut20->at(i));
				h_pho_sumIso[etaFlag][centFlag][0][2]->Fill(pho_ecalClusterIsoR2->at(i)+pho_hcalRechitIsoR2->at(i)+pho_trackIsoR2PtCut20->at(i));
				h_pho_sumIso[etaFlag][centFlag][0][3]->Fill(pho_ecalClusterIsoR3->at(i)+pho_hcalRechitIsoR3->at(i)+pho_trackIsoR3PtCut20->at(i));
				h_pho_sumIso[etaFlag][centFlag][0][4]->Fill(pho_ecalClusterIsoR4->at(i)+pho_hcalRechitIsoR4->at(i)+pho_trackIsoR4PtCut20->at(i));
				h_pho_sumIso[etaFlag][centFlag][0][5]->Fill(pho_ecalClusterIsoR5->at(i)+pho_hcalRechitIsoR5->at(i)+pho_trackIsoR5PtCut20->at(i));
				
				h_pfcIso[etaFlag][centFlag][0][2]->Fill(pfcIso2->at(i));
				h_pfcIso[etaFlag][centFlag][0][3]->Fill(pfcIso3->at(i));
				h_pfcIso[etaFlag][centFlag][0][4]->Fill(pfcIso4->at(i));
				h_pfcIso[etaFlag][centFlag][0][5]->Fill(pfcIso5->at(i));
				h_pfnIso[etaFlag][centFlag][0][2]->Fill(pfnIso2->at(i));
				h_pfnIso[etaFlag][centFlag][0][3]->Fill(pfnIso3->at(i));
				h_pfnIso[etaFlag][centFlag][0][4]->Fill(pfnIso4->at(i));
				h_pfnIso[etaFlag][centFlag][0][5]->Fill(pfnIso5->at(i));
				h_pfpIso[etaFlag][centFlag][0][2]->Fill(pfpIso2->at(i));
				h_pfpIso[etaFlag][centFlag][0][3]->Fill(pfpIso3->at(i));
				h_pfpIso[etaFlag][centFlag][0][4]->Fill(pfpIso4->at(i));
				h_pfpIso[etaFlag][centFlag][0][5]->Fill(pfpIso5->at(i));
				h_sum_pfIso[etaFlag][centFlag][0][2]->Fill(pfcIso2->at(i)+pfnIso2->at(i)+pfpIso2->at(i));
				h_sum_pfIso[etaFlag][centFlag][0][3]->Fill(pfcIso3->at(i)+pfnIso3->at(i)+pfpIso3->at(i));
				h_sum_pfIso[etaFlag][centFlag][0][4]->Fill(pfcIso4->at(i)+pfnIso4->at(i)+pfpIso4->at(i));
				h_sum_pfIso[etaFlag][centFlag][0][5]->Fill(pfcIso5->at(i)+pfnIso5->at(i)+pfpIso5->at(i));

				h_pfcVsIso[etaFlag][centFlag][0][2]->Fill(pfcVsIso2->at(i));
				h_pfcVsIso[etaFlag][centFlag][0][3]->Fill(pfcVsIso3->at(i));
				h_pfcVsIso[etaFlag][centFlag][0][4]->Fill(pfcVsIso4->at(i));
				h_pfcVsIso[etaFlag][centFlag][0][5]->Fill(pfcVsIso5->at(i));
				h_pfnVsIso[etaFlag][centFlag][0][2]->Fill(pfnVsIso2->at(i));
				h_pfnVsIso[etaFlag][centFlag][0][3]->Fill(pfnVsIso3->at(i));
				h_pfnVsIso[etaFlag][centFlag][0][4]->Fill(pfnVsIso4->at(i));
				h_pfnVsIso[etaFlag][centFlag][0][5]->Fill(pfnVsIso5->at(i));
				h_pfpVsIso[etaFlag][centFlag][0][2]->Fill(pfpVsIso2->at(i));
				h_pfpVsIso[etaFlag][centFlag][0][3]->Fill(pfpVsIso3->at(i));
				h_pfpVsIso[etaFlag][centFlag][0][4]->Fill(pfpVsIso4->at(i));
				h_pfpVsIso[etaFlag][centFlag][0][5]->Fill(pfpVsIso5->at(i));
				h_sum_pfVsIso[etaFlag][centFlag][0][2]->Fill(pfcVsIso2->at(i)+pfnVsIso2->at(i)+pfpVsIso2->at(i));
				h_sum_pfVsIso[etaFlag][centFlag][0][3]->Fill(pfcVsIso3->at(i)+pfnVsIso3->at(i)+pfpVsIso3->at(i));
				h_sum_pfVsIso[etaFlag][centFlag][0][4]->Fill(pfcVsIso4->at(i)+pfnVsIso4->at(i)+pfpVsIso4->at(i));
				h_sum_pfVsIso[etaFlag][centFlag][0][5]->Fill(pfcVsIso5->at(i)+pfnVsIso5->at(i)+pfpVsIso5->at(i));

				h_pfcBKGIso[etaFlag][centFlag][0][2]->Fill(pfcIso2->at(i) - pfcVsIso2->at(i));
				h_pfcBKGIso[etaFlag][centFlag][0][3]->Fill(pfcIso3->at(i) - pfcVsIso3->at(i));
				h_pfcBKGIso[etaFlag][centFlag][0][4]->Fill(pfcIso4->at(i) - pfcVsIso4->at(i));
				h_pfcBKGIso[etaFlag][centFlag][0][5]->Fill(pfcIso5->at(i) - pfcVsIso5->at(i));
				h_pfnBKGIso[etaFlag][centFlag][0][2]->Fill(pfnIso2->at(i) - pfnVsIso2->at(i));
				h_pfnBKGIso[etaFlag][centFlag][0][3]->Fill(pfnIso3->at(i) - pfnVsIso3->at(i));
				h_pfnBKGIso[etaFlag][centFlag][0][4]->Fill(pfnIso4->at(i) - pfnVsIso4->at(i));
				h_pfnBKGIso[etaFlag][centFlag][0][5]->Fill(pfnIso5->at(i) - pfnVsIso5->at(i));
				h_pfpBKGIso[etaFlag][centFlag][0][2]->Fill(pfpIso2->at(i) - pfpVsIso2->at(i));
				h_pfpBKGIso[etaFlag][centFlag][0][3]->Fill(pfpIso3->at(i) - pfpVsIso3->at(i));
				h_pfpBKGIso[etaFlag][centFlag][0][4]->Fill(pfpIso4->at(i) - pfpVsIso4->at(i));
				h_pfpBKGIso[etaFlag][centFlag][0][5]->Fill(pfpIso5->at(i) - pfpVsIso5->at(i));
				h_sum_pfBKGIso[etaFlag][centFlag][0][2]->Fill( (pfcIso2->at(i)+pfnIso2->at(i)+pfpIso2->at(i)) - (pfcVsIso2->at(i)+pfnVsIso2->at(i)+pfpVsIso2->at(i)) );
				h_sum_pfBKGIso[etaFlag][centFlag][0][3]->Fill( (pfcIso3->at(i)+pfnIso3->at(i)+pfpIso3->at(i)) - (pfcVsIso3->at(i)+pfnVsIso3->at(i)+pfpVsIso3->at(i)) );
				h_sum_pfBKGIso[etaFlag][centFlag][0][4]->Fill( (pfcIso4->at(i)+pfnIso4->at(i)+pfpIso4->at(i)) - (pfcVsIso4->at(i)+pfnVsIso4->at(i)+pfpVsIso4->at(i)) );
				h_sum_pfBKGIso[etaFlag][centFlag][0][5]->Fill( (pfcIso5->at(i)+pfnIso5->at(i)+pfpIso5->at(i)) - (pfcVsIso5->at(i)+pfnVsIso5->at(i)+pfpVsIso5->at(i)) );	
#endif

			}//if matched
#if 0 
			// if you want to get No gen mathcing all photons, then use this. 
			// (No Gen Mathcing) all photons(matched+fake)
			h_phoEt[etaFlag][centFlag][0]->Fill(phoEt->at(i));
			h_phoPhi[etaFlag][centFlag][0]->Fill(phoPhi->at(i));
			h_phoEta[etaFlag][centFlag][0]->Fill(phoEta->at(i));
			h_phoR9[etaFlag][centFlag][0]->Fill(phoR9->at(i));
			h_phoHoverE[etaFlag][centFlag][0]->Fill(phoHoverE->at(i));
			h_phoSigmaIEtaIEta[etaFlag][centFlag][0]->Fill(phoSigmaIEtaIEta->at(i));
			h_pho_ecalClusterIso[etaFlag][centFlag][0][2]->Fill(pho_ecalClusterIsoR2->at(i));
			h_pho_ecalClusterIso[etaFlag][centFlag][0][3]->Fill(pho_ecalClusterIsoR3->at(i));
			h_pho_ecalClusterIso[etaFlag][centFlag][0][4]->Fill(pho_ecalClusterIsoR4->at(i));
			h_pho_ecalClusterIso[etaFlag][centFlag][0][5]->Fill(pho_ecalClusterIsoR5->at(i));
			h_pho_hcalRechitIso[etaFlag][centFlag][0][2]->Fill(pho_hcalRechitIsoR2->at(i));
			h_pho_hcalRechitIso[etaFlag][centFlag][0][3]->Fill(pho_hcalRechitIsoR3->at(i));
			h_pho_hcalRechitIso[etaFlag][centFlag][0][4]->Fill(pho_hcalRechitIsoR4->at(i));
			h_pho_hcalRechitIso[etaFlag][centFlag][0][5]->Fill(pho_hcalRechitIsoR5->at(i));
			h_pho_trackIsoPtCut20[etaFlag][centFlag][0][2]->Fill(pho_trackIsoR2PtCut20->at(i));
			h_pho_trackIsoPtCut20[etaFlag][centFlag][0][3]->Fill(pho_trackIsoR3PtCut20->at(i));
			h_pho_trackIsoPtCut20[etaFlag][centFlag][0][4]->Fill(pho_trackIsoR4PtCut20->at(i));
			h_pho_trackIsoPtCut20[etaFlag][centFlag][0][5]->Fill(pho_trackIsoR5PtCut20->at(i));
			h_pho_sumIso[etaFlag][centFlag][0][2]->Fill(pho_ecalClusterIsoR2->at(i)+pho_hcalRechitIsoR2->at(i)+pho_trackIsoR2PtCut20->at(i));
			h_pho_sumIso[etaFlag][centFlag][0][3]->Fill(pho_ecalClusterIsoR3->at(i)+pho_hcalRechitIsoR3->at(i)+pho_trackIsoR3PtCut20->at(i));
			h_pho_sumIso[etaFlag][centFlag][0][4]->Fill(pho_ecalClusterIsoR4->at(i)+pho_hcalRechitIsoR4->at(i)+pho_trackIsoR4PtCut20->at(i));
			h_pho_sumIso[etaFlag][centFlag][0][5]->Fill(pho_ecalClusterIsoR5->at(i)+pho_hcalRechitIsoR5->at(i)+pho_trackIsoR5PtCut20->at(i));
#endif
		}// photon loop
	} // exited event loop

	end_loop = std::clock();
	std::cout.precision(6);      // get back to default precision
	std::cout << "LOOP finished in             : " << (end_loop - start_loop) / (double)CLOCKS_PER_SEC << " seconds" << std::endl;
	std::cout << "exited event loop" << std::endl;

	//=============================================================================================
	//=============================================================================================
	// Save the histograms 
	TString outputFileName = Form("histFiles/gedPhotonHist_%s_%s_ptThr%d%s.root",outName.Data(),treePath.Data(),(int)ptThr,outSuffix.Data());
	TFile* fout = new TFile(outputFileName, "RECREATE");
	fout->cd();

	for(int j=1; j<nEtaCut; ++j){
		for(int k=0; k<nCentCut; ++k){
			for(int m=0; m<nMomIdCut; ++m){
				h_phoEt[j][k][m]->Scale(1./entries); 
				h_phoPhi[j][k][m]->Scale(1./entries); 
				h_phoEta[j][k][m]->Scale(1./entries);
				h_phoR9[j][k][m]->Scale(1./entries);
				h_phoHoverE[j][k][m]->Scale(1./entries);
				h_phoSigmaIEtaIEta[j][k][m] ->Scale(1./entries);
				h_phoEt[j][k][m]->Write(); 
				h_phoPhi[j][k][m]->Write(); 
				h_phoEta[j][k][m]->Write();
				h_phoR9[j][k][m]->Write();
				h_phoHoverE[j][k][m]->Write();
				h_phoSigmaIEtaIEta[j][k][m]->Write();			
				for(int iR=2; iR<nRadius; iR++){
					h_pho_ecalClusterIso[j][k][m][iR]->Scale(1./entries);
					h_pho_hcalRechitIso[j][k][m][iR]->Scale(1./entries);
					h_pho_trackIsoPtCut20[j][k][m][iR] ->Scale(1./entries);
					h_pho_sumIso[j][k][m][iR] ->Scale(1./entries);
					h_pho_ecalClusterIso[j][k][m][iR]->Write();
					h_pho_hcalRechitIso[j][k][m][iR]->Write();
					h_pho_trackIsoPtCut20[j][k][m][iR] ->Write();
					h_pho_sumIso[j][k][m][iR] ->Write();
				
					h_pfcIso[j][k][m][iR]->Scale(1./entries);
					h_pfnIso[j][k][m][iR]->Scale(1./entries);
					h_pfpIso[j][k][m][iR]->Scale(1./entries);					
					h_sum_pfIso[j][k][m][iR]->Scale(1./entries);					
					h_pfcVsIso[j][k][m][iR]->Scale(1./entries);
					h_pfnVsIso[j][k][m][iR]->Scale(1./entries);
					h_pfpVsIso[j][k][m][iR]->Scale(1./entries);
					h_sum_pfVsIso[j][k][m][iR]->Scale(1./entries);	
					h_pfcBKGIso[j][k][m][iR]->Scale(1./entries);
					h_pfnBKGIso[j][k][m][iR]->Scale(1./entries);
					h_pfpBKGIso[j][k][m][iR]->Scale(1./entries);					
					h_sum_pfBKGIso[j][k][m][iR]->Scale(1./entries);					
					h_pfcIso[j][k][m][iR]->Write();
					h_pfnIso[j][k][m][iR]->Write();
					h_pfpIso[j][k][m][iR]->Write();
					h_sum_pfIso[j][k][m][iR]->Write();
					h_pfcVsIso[j][k][m][iR]->Write();
					h_pfnVsIso[j][k][m][iR]->Write();
					h_pfpVsIso[j][k][m][iR]->Write();
					h_sum_pfVsIso[j][k][m][iR]->Write();
					h_pfcBKGIso[j][k][m][iR]->Write();
					h_pfnBKGIso[j][k][m][iR]->Write();
					h_pfpBKGIso[j][k][m][iR]->Write();
					h_sum_pfBKGIso[j][k][m][iR]->Write();	
				}	
			}
			h_phoEt_fake[j][k] ->Scale(1./entries);
			h_phoPhi_fake[j][k] ->Scale(1./entries);
			h_phoEta_fake[j][k]->Scale(1./entries);
			h_phoR9_fake[j][k] ->Scale(1./entries);
			h_phoHoverE_fake[j][k] ->Scale(1./entries);
			h_phoSigmaIEtaIEta_fake[j][k] ->Scale(1./entries);
			h_phoEt_fake[j][k] ->Write();
			h_phoPhi_fake[j][k] ->Write();
			h_phoEta_fake[j][k]->Write();
			h_phoR9_fake[j][k] ->Write();
			h_phoHoverE_fake[j][k] ->Write();
			h_phoSigmaIEtaIEta_fake[j][k] ->Write();
			for(int iR=2; iR<nRadius; iR++){
				h_pho_ecalClusterIso_fake[j][k][iR] ->Scale(1./entries);
				h_pho_hcalRechitIso_fake[j][k][iR] ->Scale(1./entries);
				h_pho_trackIsoPtCut20_fake[j][k][iR] ->Scale(1./entries);
				h_pho_sumIso_fake[j][k][iR] ->Scale(1./entries);
				h_pho_ecalClusterIso_fake[j][k][iR] ->Write();
				h_pho_hcalRechitIso_fake[j][k][iR] ->Write();
				h_pho_trackIsoPtCut20_fake[j][k][iR] ->Write();
				h_pho_sumIso_fake[j][k][iR] ->Write();
			}
		}
	}
	fout->Close();
	inputFile->Close();
}
