#include "../HiForestAnalysis/hiForest.h"

class photonTree
{
   public:
   int nPho;
   vector<bool>* isGenMatched;
   vector<int> matchedMomPID;
   vector<float> matchedCalIsoDR04;
   vector<float> matchedTrkIsoDR04;
   vector<float> matchedEt;
   vector<float> phoEt;
   vector<float> phoSCRawE;
   vector<float> phoE5x5;
   vector<float> phoEta;
   vector<float> phoPhi;
   vector<float> phoSCEtaWidth;
   vector<float> phoSCPhiWidth;
   vector<float> phoSumIso;
   vector<float> phoEcalIso;
   vector<float> phoTrkIso;
   vector<float> phoSigmaIEtaIEta;
   vector<float> phoHoverE;
   float hiBin;
   TTree *t;
   void init(){
      t = new TTree("t","");
      t->Branch("isGenMatched",&isGenMatched);
      t->Branch("phoEt",&phoEt);
      t->Branch("phoSCRawE",&phoSCRawE);
      t->Branch("phoE5x5",&phoE5x5);
      t->Branch("phoEta",&phoEta);
      t->Branch("phoPhi",&phoPhi);
      t->Branch("phoSCEtaWidth",&phoSCEtaWidth);
      t->Branch("phoSCPhiWidth",&phoSCPhiWidth);
      t->Branch("phoSigmaIEtaIEta",&phoSigmaIEtaIEta);
      t->Branch("phoSumIso",&phoSumIso);
      t->Branch("phoEcalIso",&phoEcalIso);
      t->Branch("phoTrkIso",&phoTrkIso);
      t->Branch("matchedEt",&matchedEt);
      t->Branch("matchedMomPID",&matchedMomPID);
      t->Branch("matchedCalIsoDR04",&matchedCalIsoDR04);
      t->Branch("matchedTrkIsoDR04",&matchedTrkIsoDR04);
      t->Branch("phoHoverE",&phoHoverE);
      t->Branch("hiBin",&hiBin);
   }
   
   void clear(){
      isGenMatched.clear();
      matchedMomPID.clear();
      matchedCalIsoDR04.clear();
      matchedTrkIsoDR04.clear();
      matchedEt.clear();
      phoEt.clear();
      phoE5x5.clear();
      phoSCRawE.clear();
      phoEta.clear();
      phoPhi.clear();
      phoSCEtaWidth.clear();
      phoSCPhiWidth.clear();
      phoSumIso.clear();
      phoEcalIso.clear();
      phoTrkIso.clear();
      phoHoverE.clear();
      phoSigmaIEtaIEta.clear();
      
   }
   
};

void photonAna(const char *infName)
{
   HiForest h(infName,"",cPbPb);
   h.verbose=0;
   h.InitTree();
   cout << "s" << endl; 
   TFile *outf = new TFile("photonTree.root","recreate");
   cout << "s" << endl; 
   photonTree photon;   
   photon.init();

   for (int i=0;i<h.GetEntries();i++)
   {
 //     cout <<i<<endl;
      h.GetEntry(i);
      if (i%1000==0) cout <<i<<" / "<<h.GetEntries()<<endl;
      h.MatchGenPhoton(); 
      photon.clear();
      photon.hiBin=h.evt.hiBin;
      for (int j=0;j<h.photon.nPho;j++){
//         cout <<j<<endl;	
         photon.isGenMatched->push_back(h.photon.isGenMatched->at(j));
	// photon.phoSumIso.push_back((h.photon.pho_ecalClusterIsoR4->at(j)+h.photon.pho_hcalRechitIsoR4->at(j)+h.photon.pho_trackIsoR4PtCut20->at(j))/0.9);
	// photon.phoEcalIso.push_back((h.photon.pho_ecalClusterIsoR4->at(j))/0.9);
	// photon.phoTrkIso.push_back((h.photon.pho_trackIsoR4PtCut20->at(j))/0.9);
         photon.phoEt.push_back(h.photon.phoEt->at(j));
	 photon.phoE5x5.push_back(h.photon.phoE5x5->at(j));
	 photon.phoSCRawE.push_back(h.photon.phoSCRawE->at(j));
	 photon.phoEta.push_back(h.photon.phoEta->at(j));
	 photon.phoPhi.push_back(h.photon.phoPhi->at(j));
	 photon.phoSCEtaWidth.push_back(h.photon.phoSCEtaWidth->at(j));
	 photon.phoSCPhiWidth.push_back(h.photon.phoSCPhiWidth->at(j));
	 photon.phoSigmaIEtaIEta.push_back(h.photon.phoSigmaIEtaIEta->at(j));
	 photon.phoHoverE.push_back(h.photon.phoHoverE->at(j));
         if (h.photon.isGenMatched->at(j)!=0) {
 	    photon.matchedMomPID.push_back(h.photon.mcMomPID->at(h.photon.genMatchedIdx->at(j)));
            photon.matchedEt.push_back(h.photon.mcPt->at(h.photon.genMatchedIdx->at(j)));
            photon.matchedCalIsoDR04.push_back(h.photon.mcCalIsoDR04->at(h.photon.genMatchedIdx->at(j)));
            photon.matchedTrkIsoDR04.push_back(h.photon.mcTrkIsoDR04->at(h.photon.genMatchedIdx->at(j)));
         } else {
	    photon.matchedMomPID.push_back(0);
            photon.matchedEt.push_back(0);
            photon.matchedCalIsoDR04.push_back(0);
            photon.matchedTrkIsoDR04.push_back(0);
         }
      }
      photon.t->Fill();
   }
   photon.t->Write();
   outf->Close();

}

int main(int argc, char** argv)
{
    photonAna(argv[0]);
    return 0;
}
