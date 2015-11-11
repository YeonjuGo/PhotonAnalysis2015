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

const int nCentBin = 2;
const int centBin[nCentBin+1] = {0,60,200};

double col[] = {12,8,9,41,46};
void makeBkgTree(const char* hiForestfileName="/afs/cern.ch/work/y/ygo/public/PFphoton/Pyquen_Unquenched_AllQCDPhoton30_PhotonFilter20GeV_eta24_TuneZ2_PbPb_5020GeV.root",
        TString outName="AllQCDPhoton30", Float_t etThr=0.0, TString det="tower")
{
    TH1::SetDefaultSumw2();
    gStyle->SetOptStat(0);

    const int nEP = 4;
    const double epmax = TMath::Pi();
    double epBin[nEP+1];
    for(int i=0; i<nEP+1;i++){
        double binWidth = (2*epmax)/(double)nEP;
        epBin[i] = -epmax + binWidth*i;
    }
    TFile* inputFile = new TFile(hiForestfileName, "READ");
    std::cout << "input HiForest : " << inputFile->GetName() << std::endl;
    TString outputFileName = Form("skimFiles/jskim_%s_towerHist.root",outName.Data());
    TFile* fout = new TFile(outputFileName, "RECREATE");
    fout->cd();

    TTree* hiEvtAnalyzerTree = (TTree*)inputFile->Get("hiEvtAnalyzer/HiTree");
    TTree* towerTree= (TTree*)inputFile->Get("rechitanalyzer/tower");
    TTree* hbheTree= (TTree*)inputFile->Get("rechitanalyzer/hbhe");
    TTree* hfTree= (TTree*)inputFile->Get("rechitanalyzer/hf");
    TTree* eeTree= (TTree*)inputFile->Get("rechitanalyzer/ee");
    TTree* ebTree= (TTree*)inputFile->Get("rechitanalyzer/eb");

    towerTree->AddFriend(hiEvtAnalyzerTree);
    hbheTree->AddFriend(hiEvtAnalyzerTree);
    hfTree->AddFriend(hiEvtAnalyzerTree);
    eeTree->AddFriend(hiEvtAnalyzerTree);
    ebTree->AddFriend(hiEvtAnalyzerTree);

    Int_t hiBin, hiNevtPlane;
    Float_t hiEvtPlanes[50];
    hiEvtAnalyzerTree->SetBranchAddress("hiBin", &hiBin);
    hiEvtAnalyzerTree->SetBranchAddress("hiNevtPlane", &hiNevtPlane);
    hiEvtAnalyzerTree->SetBranchAddress("hiEvtPlanes", hiEvtPlanes);

    // tower Info.
    Int_t n_tower;
    Float_t e_tower[6000];
    Float_t et_tower[6000];
    Float_t eta_tower[6000];
    Float_t phi_tower[6000];
    Float_t emEt_tower[6000];
    Float_t hadEt_tower[6000];

    // hbhe Info.
    Int_t n_hbhe;
    Float_t e_hbhe[10];
    Float_t et_hbhe[10];
    Float_t eta_hbhe[10];
    Float_t phi_hbhe[10];
    Float_t perp_hbhe[10];

    // hf Info.
    Int_t n_hf;
    Float_t e_hf[10];
    Float_t et_hf[10];
    Float_t eta_hf[10];
    Float_t phi_hf[10];
    Float_t perp_hf[10];

    // ee Info.
    Int_t n_ee;
    Float_t e_ee[10];
    Float_t et_ee[10];
    Float_t eta_ee[10];
    Float_t phi_ee[10];
    Float_t perp_ee[10];

    // eb Info.
    Int_t n_eb;
    Float_t e_eb[10];
    Float_t et_eb[10];
    Float_t eta_eb[10];
    Float_t phi_eb[10];
    Float_t perp_eb[10];

    towerTree->SetBranchAddress("n",&n_tower);
    towerTree->SetBranchAddress("e",e_tower);
    towerTree->SetBranchAddress("et",et_tower);
    towerTree->SetBranchAddress("eta",eta_tower);
    towerTree->SetBranchAddress("phi",phi_tower);
    towerTree->SetBranchAddress("emEt",emEt_tower);
    towerTree->SetBranchAddress("hadEt",hadEt_tower);
    ebTree->SetBranchAddress("n",&n_eb);
    ebTree->SetBranchAddress("e",e_eb);
    ebTree->SetBranchAddress("et",et_eb);
    ebTree->SetBranchAddress("eta",eta_eb);
    ebTree->SetBranchAddress("phi",phi_eb);
    eeTree->SetBranchAddress("n",&n_ee);
    eeTree->SetBranchAddress("e",e_ee);
    eeTree->SetBranchAddress("et",et_ee);
    eeTree->SetBranchAddress("eta",eta_ee);
    eeTree->SetBranchAddress("phi",phi_ee);
    hfTree->SetBranchAddress("n",&n_hf);
    hfTree->SetBranchAddress("e",e_hf);
    hfTree->SetBranchAddress("et",et_hf);
    hfTree->SetBranchAddress("eta",eta_hf);
    hfTree->SetBranchAddress("phi",phi_hf);
    hbheTree->SetBranchAddress("n",&n_hbhe);
    hbheTree->SetBranchAddress("e",e_hbhe);
    hbheTree->SetBranchAddress("et",et_hbhe);
    hbheTree->SetBranchAddress("eta",eta_hbhe);
    hbheTree->SetBranchAddress("phi",phi_hbhe);

    TH1::SetDefaultSumw2();
    std::cout << "entering event loop" << std::endl;
    Long64_t entries = towerTree->GetEntries();
    std::cout << "number of entries = " << entries << std::endl;

    std::clock_t    start_loop, end_loop;
    start_loop = std::clock();

    TH1D* h1D_dphi[5][5][5];//[eta][cent][evtpl_order]
    //TH1D* h1D_dphi[5][5][5][10];//[eta][cent][evtpl_order][evp_dphi]
    for(int ieta=1; ieta<nEtaCut; ieta++){
        for(int icent=0;icent<nCentBin;icent++){
            for(int iep=0;iep<5;iep++){
                //for(int iep_dphi=0;iep_dphi<nEP;iep_dphi++){
                    h1D_dphi[ieta][icent][iep] = new TH1D(Form("h1D_dphi_eta%d_cent%d_epOrder%d",ieta,icent,iep),";;",20,-TMath::Pi(),TMath::Pi());
                    //h1D_dphi[ieta][icent][iep][iep_dphi] = new TH1D(Form("h1D_dphi_eta%d_cent%d_epOrder%d_epDphi%d",ieta,icent,iep,iep_dphi),";;",100,-TMath::Pi(),TMath::Pi());
               // }
            }
        }
    }

    //for(Long64_t jj = 0; jj < 10; ++jj){}
    for(Long64_t jj = 0; jj < entries; ++jj)
    {
        if (jj % 100000 == 0)  {
            std::cout << "current entry = " <<jj<<" out of "<<entries<<" : "<<std::setprecision(2)<<(double)jj/entries*100<<" %"<<std::endl;
        }

        hiEvtAnalyzerTree->GetEntry(jj);
        towerTree->GetEntry(jj);
        hbheTree->GetEntry(jj);
        hfTree->GetEntry(jj);
        eeTree->GetEntry(jj);
        ebTree->GetEntry(jj);
        //cout << "event : " << jj << "/// nTower : " << n_tower << ", n_hbhe : " << n_hbhe << ", hiBin : "<< hiBin<< endl;

        int centNum = -9;
        int etaNum = -9;
        int epDphiNum = -9;
        if(hiBin>=0 && hiBin < 60) centNum = 0;
        else if(hiBin>=60 && hiBin < 200) centNum = 1;

        double dphi_evp[5];

        for (int i=0; i < n_tower; ++i){
            if(et_tower[i] < etThr) continue;
            dphi_evp[1] = getDPHI(phi_tower[i],hiEvtPlanes[2]);
            dphi_evp[2] = getDPHI(phi_tower[i],hiEvtPlanes[8]);
            dphi_evp[3] = getDPHI(phi_tower[i],hiEvtPlanes[15]);
            dphi_evp[4] = getDPHI(phi_tower[i],hiEvtPlanes[21]);

            if(eta_tower[i]>=0 && eta_tower[i]<1.4442) etaNum = 1;
            else if(eta_tower[i]>=1.566 && eta_tower[i] < 2.00) etaNum = 2;
            else if(eta_tower[i]>=2.00) etaNum = 3;
            else etaNum=-9;
            if(etaNum==-9) continue;

            for(int iep=0;iep<nEP;iep++){
                h1D_dphi[etaNum][centNum][iep]->Fill(dphi_evp[iep+1]);
            }
        }//tower loop
    } // exited event loop

    end_loop = std::clock();
    std::cout.precision(6);      // get back to default precision
    std::cout << "LOOP finished in             : " << (end_loop - start_loop) / (double)CLOCKS_PER_SEC << " seconds" << std::endl;
    std::cout << "exited event loop" << std::endl;

    TCanvas* can2 = new TCanvas("can2","", 1000,500);
    makeMultiPanelCanvas(can2,3,2,0.0,0.0,0.2,0.15,0.02);

    TCanvas* c1 = new TCanvas("c1", "", 100,100);
    TLegend* leg2 = new TLegend(0.6,0.03,0.95,0.4);
    legStyle(leg2);
    double xpos = 0.5;
    double ypos = 0.8;
    double dy = 0.06;
    for(int ieta=1; ieta<nEtaCut; ieta++){
       for(int icent=0;icent<nCentBin;icent++){
            TString centSt;
            if(icent==0) centSt = "0-30%";
            else centSt = "30-100%";

            for(int iep=0;iep<nEP;iep++){
                can2->cd(ieta+3*icent);
                if(ieta==1 && icent==0) leg2->AddEntry(h1D_dphi[ieta][icent][iep],Form("order %d",iep+1));
                h1D_dphi[ieta][icent][iep]->SetMarkerStyle(20);
                h1D_dphi[ieta][icent][iep]->SetMarkerColor(col[iep]);
                h1D_dphi[ieta][icent][iep]->SetAxisRange(0.04,0.06,"Y");
                h1D_dphi[ieta][icent][iep]->GetYaxis()->SetRangeUser(0.04,0.06);
                h1D_dphi[ieta][icent][iep]->Scale(1./(h1D_dphi[ieta][icent][iep]->GetEntries()));
                if(iep==1) h1D_dphi[ieta][icent][1]->DrawCopy("p");
            }
            
            for(int iep=0;iep<nEP;iep++){
                can2->cd(ieta+3*icent);
                if(iep!=1) h1D_dphi[ieta][icent][iep]->DrawCopy("p same");
            }

            if(ieta==1 && icent==0) {leg2->Draw("same"); drawText(centSt,xpos,ypos);}
            if(ieta==1 && icent==1) drawText(centSt,xpos,ypos);
        }
        can2->cd(ieta);
        drawText(Form("%.2f<|eta|<%.2f",eta_gt[ieta],eta_lt[ieta]),xpos,ypos+2*dy); 
    }
    can2->SaveAs(Form("png/towerDphiDist_etThr%d.png",(int)etThr));
    //=============================================================================================
    //=============================================================================================
    // Save the histograms 
/*
    for(int ieta=1; ieta<nEtaCut; ieta++){
        for(int icent=0;icent<nCentBin;icent++){
            for(int iep=0;iep<5;iep++){
                h1D_dphi[ieta][icent][iep]->Write();
            }
        }
    }
    fout->Close();
*/
    inputFile->Close();
}
