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

    TTree* hiEvtAnalyzerTree = (TTree*)inputFile->Get("hiEvtAnalyzer/HiTree");
    TTree* tree= (TTree*)inputFile->Get(Form("rechitanalyzer/%s",det.Data()));
    tree->AddFriend(hiEvtAnalyzerTree);

    Int_t hiBin, hiNevtPlane;
    Float_t hiEvtPlanes[50];
    hiEvtAnalyzerTree->SetBranchAddress("hiBin", &hiBin);
    hiEvtAnalyzerTree->SetBranchAddress("hiNevtPlane", &hiNevtPlane);
    hiEvtAnalyzerTree->SetBranchAddress("hiEvtPlanes", hiEvtPlanes);

    Int_t n;
    Float_t e[6000];
    Float_t et[6000];
    Float_t eta[6000];
    Float_t phi[6000];
    tree->SetBranchAddress("n",&n);
    tree->SetBranchAddress("e",e);
    tree->SetBranchAddress("et",et);
    tree->SetBranchAddress("eta",eta);
    tree->SetBranchAddress("phi",phi);
/*
    if(det=="tower"){
        Float_t emEt[6000];
        Float_t hadEt[6000];
        tree->SetBranchAddress("emEt",emEt);
        tree->SetBranchAddress("hadEt",hadEt);
    }
*/
    TH1::SetDefaultSumw2();
    std::cout << "entering event loop" << std::endl;
    Long64_t entries = tree->GetEntries();
    std::cout << "number of entries = " << entries << std::endl;

    std::clock_t    start_loop, end_loop;
    start_loop = std::clock();

    TH1D* h1D_dphi[5][5][5];//[eta][cent][evtpl_order]
    for(int ieta=1; ieta<nEtaCut; ieta++){
        for(int icent=0;icent<nCentBin;icent++){
            for(int iep=0;iep<5;iep++){
                    h1D_dphi[ieta][icent][iep] = new TH1D(Form("h1D_dphi_eta%d_cent%d_epOrder%d",ieta,icent,iep),";;",20,-TMath::Pi(),TMath::Pi());
                    //h1D_dphi[ieta][icent][iep][iep_dphi] = new TH1D(Form("h1D_dphi_eta%d_cent%d_epOrder%d_epDphi%d",ieta,icent,iep,iep_dphi),";;",100,-TMath::Pi(),TMath::Pi());
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
        tree->GetEntry(jj);
       //cout << "event : " << jj << "/// nTower : " << n_tower << ", n_hbhe : " << n_hbhe << ", hiBin : "<< hiBin<< endl;

        int centNum = -9;
        int etaNum = -9;
        int epDphiNum = -9;
        if(hiBin>=0 && hiBin < 60) centNum = 0;
        else if(hiBin>=60 && hiBin < 200) centNum = 1;

        double dphi_evp[5];

        for (int i=0; i < n; ++i){
            if(et[i] < etThr) continue;
            dphi_evp[1] = getDPHI(phi[i],hiEvtPlanes[2]);
            dphi_evp[2] = getDPHI(phi[i],hiEvtPlanes[8]);
            dphi_evp[3] = getDPHI(phi[i],hiEvtPlanes[15]);
            dphi_evp[4] = getDPHI(phi[i],hiEvtPlanes[21]);

            if(eta[i]>=0 && eta[i]<1.4442) etaNum = 1;
            else if(eta[i]>=1.566 && eta[i] < 2.00) etaNum = 2;
            else if(eta[i]>=2.00) etaNum = 3;
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

    cout << "ssssssssssssssss" << endl;
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
    cout << "ssssssssssssssss" << endl;
    can2->SaveAs(Form("png/%sDphiDist_etThr%d.png",det.Data(),(int)etThr));
    //inputFile->Close();
}
