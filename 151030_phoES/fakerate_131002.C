// yeonju
// pt fine bin 
#include "../HiForestAnalysis/hiForest.h"
#include "../gammaJetAnalysis/CutAndBinCollection2012.h"
#include <TStyle.h>
#include <TH1D.h>
#include <TTree.h>
#include <TBranch.h>
#include <TFile.h>
#include <iostream>
#include <TCanvas.h>
#include <TF1.h>
#include <TGraphAsymmErrors.h>
#include <TLegend.h>
#include <TCut.h>

double myFunc(double *x, double *par){
   return par[0]+ par[1]/sqrt(x[0]) + par[2]/x[0];
}

double myFunc2(double *x, double *par){
    return sqrt(0.0561326*0.0561326+ 0.819547*0.819547/x[0] + par[0]*par[0]/(x[0]*x[0]));
}  

double myFunc3(double *x, double *par){
    return sqrt(par[0]*par[0]+ par[1]*par[1]/x[0] + par[2]*par[2]/(x[0]*x[0]));
}


void fakerate_131002(int flvOpt = 0, int genOpt = 0){
    /*  const int kHIcentral = 0; // 0-30%
        const int kHIperipheral = 1;//30-100%
        const int kPP = 2;
        const int kPA = 3;
        const int kHI010 = 4; //0-10%
        const int kHI1030 = 5; //10-30%
        const int kHI3050 = 6;//30-50%
        const int kHI50100 = 7;//50-100%

*/
    const int collision = 2;
    TLegend *l1 = new TLegend(0.4365615,0.6445304,0.9577623,0.846736,NULL,"brNDC");

    TCut centCut = "";
    if ( (collision ==0) )   {   
        centCut = "cBin > 0 && cBin< 12";
        easyLeg(l1,"Pb+Pb 0-30%");
    }
    else if (  (collision ==1) ){  // if it's pp 
        centCut = "cBin >=12";
        easyLeg(l1,"Pb+Pb 30-100%");
    }
    else if (collision == 2 || collision == 3){ // pPb
        centCut = "";
        if (collision == 2) easyLeg(l1,"p+p");
        else easyLeg(l1, "p+Pb");
    }
    else if (  (collision ==4) ){  //HI
        centCut = "cBin > 0 && cBin < 4";
        easyLeg(l1,"Pb+Pb 0-10%");
    } else if (  (collision ==5) ){  //HI
        centCut = "cBin >= 4 && cBin < 12";
        easyLeg(l1,"Pb+Pb 10-30%");
    } else if (  (collision ==6) ){  //HI
        centCut = "cBin >= 12 && cBin < 20";
        easyLeg(l1,"Pb+Pb 30-50%");
    } else if (  (collision ==7) ){  //HI
        centCut = "cBin >= 20";
        easyLeg(l1,"Pb+Pb 50-100%");
    }    

    //        
    TH1::SetDefaultSumw2();

    // gStyle->SetOptFit(0);
    gStyle -> SetOptStat(0);
    //  gStyle -> SetTitleYOffset(2.35);
    gStyle -> SetTitleYSize(0.04);





    TCut partonCut = "";
    if (flvOpt ==0 ) partonCut = "";
    else if (flvOpt == 1 ) partonCut = "refPartonFlv == 21";
    else if (flvOpt == 2 ) partonCut = "abs(refPartonFlv) < 21";
    else partonCut = "refPartonFlv < -200";


    const int nPtBin = 12;
    //const int nPtBin = 17;
    //double ptBin[nPtBin+1] = {30.0, 40.0, 50.0, 80.0, 9999.0}; 
    //double AvePtBin[nPtBin+1] = { 35.0, 45.0, 65.0, 100.0};
    double ptBin[nPtBin+1];
    double AvePtBin[nPtBin];
   
        ptBin[0] = 30.0;
        ptBin[1] = 40.0;
        ptBin[2] = 50.0;
        ptBin[3] = 60.0;
        ptBin[4] = 70.0;
        ptBin[5] = 80.0;
        ptBin[6] = 90.0;
        ptBin[7] = 100.0;
        ptBin[8] = 120.0;
        ptBin[9] = 140.0;
        ptBin[10] = 160.0;
        ptBin[11] = 180.0;
        ptBin[12] = 200.0;

/*
    for(int i=0;i<=nPtBin;i++){
        ptBin[i] = 30.0+i*10.0;
    }
  */  for(int i=0;i<nPtBin;i++){
        AvePtBin[i] = (ptBin[i+1]+ptBin[i])/2.0;
    }

    //    const int nCentBinPa = 3;
    //   double centBinPa[nCentBinPa+1] = {0,20,30,100};
    //
    // fake rate
    //
    int nJetmax = 100;
    float refPt[nJetmax], pt[nJetmax], eta[nJetmax], dphi[nJetmax];
    int nJet, cBin, refPartonFlv[nJetmax];
    EvtSel evtImb;
    TBranch *b_evt;

    multiTreeUtil* yJet = new multiTreeUtil();
    //:PbPb
    if(collision ==1 || collision == 0  || collision == 4 || collision == 5 || collision == 6|| collision == 7){
        yJet -> addFile("/home/jazzitup/forestFiles/yskimmedFiles/yskim_qcdAllPhoton30to50_genPhotonPtCut40_allCent.root", "yJet", "",32796./32796);
        yJet -> addFile("/home/jazzitup/forestFiles/yskimmedFiles/yskim_qcdAllPhoton50to80_genPhotonPtCut40_allCent.root", "yJet", "",21470./53876);
        yJet -> addFile("/home/jazzitup/forestFiles/yskimmedFiles/yskim_qcdAllPhoton80to9999_genPhotonPtCut40_allCent.root", "yJet", "",6462./58781);  
    } else if (collision ==3){
        // pA
        yJet->addFile("/home/jazzitup/forestFiles/yskimmedFiles/yskim_PA2013_pyquen_allQCDPhoton30to50_forestv85.root", "yJet","", 56669./50385);
        yJet->addFile("/home/jazzitup/forestFiles/yskimmedFiles/yskim_PA2013_pyquen_allQCDPhoton50to80_forestv85.root", "yJet","", 41906./114136);
        yJet->addFile("/home/jazzitup/forestFiles/yskimmedFiles/yskim_PA2013_pyquen_allQCDPhoton80to120_forestv85.root", "yJet","",12044./103562);
        yJet->addFile("/home/jazzitup/forestFiles/yskimmedFiles/yskim_PA2013_pyquen_allQCDPhoton120to9999_forestv85.root", "yJet","", 4481./151511);
    } else if (collision ==2){
        //pp 
        yJet->addFile("/home/jazzitup/forestFiles/yskimmedFiles/yskim_merged_allQCDPhoton30to50_genPhotonPtCut30_CMSSW538HIp2.root", "yJet", "",29329. / 29329.);
        yJet->addFile("/home/jazzitup/forestFiles/yskimmedFiles/yskim_merged_allQCDPhoton50to80_genPhotonPtCut30_CMSSW538HIp2.root", "yJet", "",8098. / 87988.);
        yJet->addFile("/home/jazzitup/forestFiles/yskimmedFiles/yskim_merged_allQCDPhoton80to120_genPhotonPtCut30_CMSSW538HIp2.root", "yJet", "",1680. / 96756.);
        yJet->addFile("/home/jazzitup/forestFiles/yskimmedFiles/yskim_merged_allQCDPhoton120to9999_genPhotonPtCut30_CMSSW538HIp2.root", "yJet", "", 438.   / 90972.);
    }

    yJet -> AddFriend("tgj");
    yJet->AddFriend("yPhotonTree");


    TCanvas* c4 = new TCanvas("c4","",400,800);
    c4->Divide(1,2);
    c4->cd(1);
    TH1D* hpt1 = new TH1D("hpt1",";p_{T}^{RECO};Entries",20,0,200);
    TH1D* hpt0 = (TH1D*)hpt1->Clone("hpt2");
    TH1D* hpt2 = (TH1D*)hpt1->Clone("hpt2");
    TH1D* hpt3 = (TH1D*)hpt1->Clone("hpt3");

    yJet -> Draw2(hpt0, "refPt", Form(" photonEt>40 && genPhotonEt> 30 && abs(genMomId)<=22 && (abs(eta) < 1.6) && (dphi > 7*3.141592/8.0) "),"");
    yJet -> Draw2(hpt1, "refPt", Form("  photonEt>40 &&  genPhotonEt> 30 && abs(genMomId)<=22 && (abs(eta) < 1.6) && (dphi > 7*3.141592/8.0) && refPartonFlv == 21"),"");
    yJet -> Draw2(hpt2, "refPt", Form(" photonEt>40 &&  genPhotonEt> 30 && abs(genMomId)<=22 && (abs(eta) < 1.6) && (dphi > 7*3.141592/8.0) && abs(refPartonFlv)<21 "),"");
    yJet -> Draw2(hpt3, "refPt", Form(" photonEt>40 &&  genPhotonEt> 30 && abs(genMomId)<=22 && (abs(eta) < 1.6) && (dphi > 7*3.141592/8.0) && refPartonFlv < -200"),"");

    handsomeTH1(hpt0,1);
    handsomeTH1(hpt1,1);
    handsomeTH1(hpt2,2);
    handsomeTH1(hpt3,4);

    hpt0->GetYaxis()->SetTitleOffset(1.8);

    hpt0->DrawCopy("hist");
    hpt1->DrawCopy("same");
    hpt2->DrawCopy("same");
    hpt3->DrawCopy("same");
    jumSun(30,0,30,7400,2);

    c4->cd(2);
    hpt1->Divide(hpt0);
    hpt2->Divide(hpt0);
    hpt3->Divide(hpt0);
    hpt1->SetAxisRange(0,1,"Y");
    hpt1->SetYTitle("Ratio");
    hpt1->DrawCopy();
    hpt2->DrawCopy("same");
    hpt3->DrawCopy("same");
    jumSun(30,0,30,1,2);

    //    return;


    // pt spectrum
    //    TCanvas* c1=  new TCanvas("c1", "", 500,500);
    //   TH1D* hptHat = new TH1D("hptHat",";pt hat;Entries",200,0,200);
    //   yJet->Draw3( hptHat, "yPhotonTree.ptHat"," photonEt>40 && genPhotonEt> 30 && abs(genMomId)<=22","");
    //yJet->Draw3( hptHat, "yPhotonTree.pt","","");
    //    return;


    // Energy Scale
    TCanvas* c2 = new TCanvas("c2", "pt/refPt distribution", 1200, 900); 
    TCanvas* ccc = new TCanvas("ccc", "pt/refpt 30-40GeV", 400, 400); 
    makeMultiPanelCanvas(c2,5,4,0.0,0.0,0.2,0.15,0.02);


    TH1D* Escale[nCentBinPa+5][nPtBin];
    double mean[nCentBinPa+5][nPtBin], var[nCentBinPa+5][nPtBin], resol[nCentBinPa+5][nPtBin], resolVar[nCentBinPa+5][nPtBin];
    int icent =1;
    // for(int icent=1; icent <= nCentBinPa ; icent++){
        for(int i=0; i < nPtBin ; i++){
          //  c2 -> cd((icent-1)*4+i+1);
            c2 -> cd(i+1);
            Escale[icent][i] = new TH1D(Form("Escale%d_%d",icent, i) , " ;  p_{T}^{RECO}/p_{T}^{GEN}", 50, 0, 2);

            if ( genOpt == 0 )  {
                yJet -> Draw2(Escale[icent][i], "pt/refPt", centCut && partonCut && Form(" (abs(eta) < 1.6) && (dphi > 7*3.141592/8.0) && (refPt >= %d && refPt < %d)", (int)ptBin[i], (int)ptBin[i+1]),""); 
            }
            else if (genOpt == 1)  {
                yJet -> Draw2(Escale[icent][i], "pt/refPt", centCut && partonCut && Form(" (abs(eta) < 1.6) && (dphi > 7*3.141592/8.0) && (pt >= %d && pt < %d) ", (int)ptBin[i], (int)ptBin[i+1]),"");
            }
            //         Escale[icent][i] -> Draw(); 
            TF1* ff = cleverGaus(Escale[icent][i]);
            gPad->SetLogy();
            mean[icent][i] = ff->GetParameter(1);
            var[icent][i] = ff->GetParError(1);
            resol[icent][i] = ff->GetParameter(2);
            resolVar[icent][i] = ff->GetParError(2);

            float dx1;    
            // ((icent==1)||(icent==4))? dx1=0.15 : dx1=0 ;
            dx1=0;
            //            if ( icent == nCentBinPa )
            //  drawText(Form("E_{T}^{HF|#eta|>4} > %dGeV, ", (int)centBinPa[icent-1]), 0.12+dx1,0.929118,1,15);//yeonju 130805
            //else
            //    drawText(Form("%dGeV < E_{T}^{HF|#eta|>4} < %dGeV, ", (int)centBinPa[icent-1], (int)centBinPa[icent]), 0.12+dx1,0.929118,1,15);

            if ( i+1 == nPtBin ) 
                drawText(Form("p_{T}^{GEN Jet} > %dGeV, ", (int)ptBin[i]), 0.17+dx1,0.84,1,15);//yeonju 130823
            else
                drawText(Form("%dGeV < p_{T}^{GEN Jet} < %dGeV, ", (int)ptBin[i], (int)ptBin[i+1]), 0.17+dx1,0.84,1,12);//yeonju 130823

            //            TLegend *l1 = new TLegend(0.6365615,0.6445304,0.9577623,0.846736,NULL,"brNDC");
            //            easyLeg(l1,"p+Pb 5.02TeV");
            //    l1->AddEntry(hxjgNorm[kPADATA][icent + kPADATA*50][j],"pPb ","p");
            //                  if ( icent==1 && j==1)   l1->Draw();



        }
   // }
    c2 -> Update();
    

    ccc -> Divide(2,1);
    ccc -> cd(1);
Escale[1][0]->Draw();
    ccc -> cd(2);
Escale[1][1]->Draw();





    //
    // Energy Scale
    //

    TCanvas *c3 = new TCanvas("c3", "Energy Scale vs. jet pt",147,37,930,465);
    c3->Divide(2,1);
    c3->cd(1);

    TH1D* h = new TH1D("h",";JetP_{T}^{GEN} (GeV); Energy Scale ",1000,20,200);
    handsomeTH1(h,1);
    h -> SetAxisRange(0.9, 1.15, "Y");
    h -> Draw();
    jumSun(20,1,200,1);
    jumSun(30,0.9,30,1.15,2);

    l1 -> Draw();
    TGraphAsymmErrors* g[nCentBinPa+5];
    g[1] = new TGraphAsymmErrors();
    g[2] = new TGraphAsymmErrors();
    g[3] = new TGraphAsymmErrors();
    TLegend *l3=new TLegend(0,0,0.4490239,0.08695652,NULL,"brNDC");
    l3->SetTextFont(42);
    l3->SetTextSize(0.04);
    l3->SetFillColor(0);
    l3->SetLineColor(0);

    TF1 *fit0 = new TF1("fit0","pol1(0)",30,100);
    TF1 *fit1 = new TF1("fit1","pol1(0)",100,200);
    for(int icent=1; icent <= nCentBinPa ; icent++){
        for(int i=0; i < nPtBin ; i++){ 
            g[icent] -> SetPoint(i+1, AvePtBin[i], mean[icent][i]);
            g[icent] -> SetPointError(i+1, 0.001,0.001, var[icent][i], var[icent][i]);
            g[icent] -> SetMarkerStyle(20);
            g[icent] -> SetMarkerSize(1);
            g[icent] -> SetMarkerColor(icent);
            //          g[icent] -> Fit("fit0","R");
            //          g[icent] -> Fit("fit1","R");
            g[icent] -> Draw("same p");
        }
        //   l3->AddEntry(g[icent], Form("%d GeV < E_{T}^{HF|#eta|>4} < %d GeV", (int)centBinPa[icent-1], (int)centBinPa[icent]), "p"); 
    }
    //     l3 -> Draw();
   TF1 *fs = new TF1("fs", myFunc, 30, 200, 3);
   /*  f -> SetParameters(0.7);
    f -> SetParNames("N");
    */
    fs -> SetParameters(0.04,0.8,0.001);
    fs -> SetParNames("C", "S", "N");

    double params[3];
   
    g[1] -> Fit("fs", "RLL");
//    double a = f->GetParameter(0);
//    double b = f->GetParameter(1);
//    double c = f->GetParameter(2);
    fs->GetParameters(params);
    cout << " scale C : " << params[0] << ", S : " << params[1] << ", N : : " << params[2] << endl;
    //cout << " C : " << a << ", S : " << b << ", N : : " << c << endl;
   // cout << " N : " << param[0] << endl;

    fs->Draw("same");
 

    //
    // Energy Resolution
    //
    c3->cd(2);
    TH1D* h_resol = new TH1D("h_resol",";JetP_{T} (GeV); Energy Resolution ",1000,20,200);
    h_resol -> SetAxisRange(0.0, 0.3, "Y");
    h_resol->GetYaxis()->CenterTitle();

    h_resol -> Draw();

    jumSun(30,0.0,30,0.3,2);

    TGraphAsymmErrors* gResol[nCentBinPa+5];
    gResol[1] = new TGraphAsymmErrors();
    gResol[2] = new TGraphAsymmErrors();
    gResol[3] = new TGraphAsymmErrors();
    
    TLegend *l_resol=new TLegend(0.40,0.20,0.85,0.42);
    l_resol->SetTextFont(42);
    l_resol->SetTextSize(0.04);
    l_resol->SetFillColor(0);
    l_resol->SetLineColor(0);

    for(int icent=1; icent <= nCentBinPa ; icent++){
        for(int i=0; i < nPtBin ; i++){ 
            //            cout << "resol : " << resol[icent][i] << ", resolVar : " << resolVar[icent][i]<< endl;
            gResol[icent] -> SetPoint(i+1, AvePtBin[i], resol[icent][i]);
            gResol[icent] -> SetPointError(i+1, 0.001,0.001, resolVar[icent][i], resolVar[icent][i]);
            gResol[icent] -> SetMarkerStyle(20);
            gResol[icent] -> SetMarkerSize(1);
            gResol[icent] -> SetMarkerColor(icent);
            gResol[icent] -> Draw("same p");
        }
        l_resol->AddEntry(gResol[icent], Form("%d GeV < E_{T}^{HF|#eta|>4} < %d GeV", (int)centBinPa[icent-1], (int)centBinPa[icent]), "p"); 
    }
    TF1 *f = new TF1("f", myFunc3, 30, 200, 3);
   /*  f -> SetParameters(0.7);
    f -> SetParNames("N");
    */
    f -> SetParameters(0.04,0.8,0.001);
    f -> SetParNames("C", "S", "N");

    double param[3];
   
    gResol[1] -> Fit("f", "RLL");
    f->GetParameters(param);
    cout << " resol C : " << param[0] << ", S : " << param[1] << ", N : : " << param[2] << endl;
    //cout << " C : " << a << ", S : " << b << ", N : : " << c << endl;
    //cout << " N : " << param[0] << endl;

    f->Draw("same");
    l1 -> Draw();
    //  l_resol->Draw();
    //flvOpt = 0, int genOpt
    // c1 -> SaveAs(".gif");
    c2 -> SaveAs(Form("c2_flvOpt%d_genOpt%d.pdf",flvOpt,genOpt));
    c3 -> SaveAs(Form("c3_flvOpt%d_genOpt%d.pdf",flvOpt,genOpt));
    //    c_fake -> SaveAs(".gif");


}
