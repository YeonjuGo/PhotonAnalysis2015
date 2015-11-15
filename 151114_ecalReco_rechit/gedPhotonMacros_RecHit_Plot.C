/*
 * macro to study different photon Reconstruction algorithms
 * version 6 - work on data sample, there will be no GEN particles
 *           - disable any plot related to GEN particles.
 *           - this code is a subset of previous versions. commented out a lot of areas.
 *               so the code is more complicated and less efficient than it should be.
 *               there are dummy operations
 * */

#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TList.h>
#include <TString.h>
#include <TLegend.h>

#include <iostream>     // std::cout
#include <ctime>        // std::clock()
#include <algorithm>    // std::find()
#include <iomanip>      // std::setprecision()
#include <vector>

//#include "/net/hisrv0001/home/tatar/code/HIUtils/histoUtil.h"
#include "../HIUtils/histoUtil.h"

void   gedPhotonMacros_RecHit_Plot();
void   gedPhotonMacros_RecHit2_Plot();
void   gedPhotonMacros_pfCand_Plot(const char* inputFileName, const char* outputDirName);
TList* drawSame(TList* histos1, TList* histos2, TList* histos3);
TList* draw2D_vRecHit(TList* histos);
TList* getListOfALLHistograms2D(TDirectoryFile* dir);

void gedPhotonMacros_RecHit_Plot()
{
    const char* fileName = "RecHit1_Dijet_NcollFilt_pthat80_740pre8_hiForest.root";
    TFile* file = new TFile(fileName, "READ");
    std::cout << "input file for Histograms : " << file->GetName() << std::endl;

    // TH2D
    TList* histos2D = getListOfALLHistograms2D(file);

    TList* histos_RecHit = new TList();
    TList* histos_t1 = new TList();
    //TList* histos_t2 = new TList();
    TH2D* h;
    TString h_name;
    TIter* iter = new TIter(histos2D);
    while ((h=(TH2D*)iter->Next())) {

        h_name=h->GetName();
        if(h_name.Contains("recHit"))
        {
            histos_RecHit->Add(h);
        }
        else if (h_name.Contains("t1")){
            histos_t1->Add(h);
        }
        else if (h_name.Contains("t2")){
            histos_t2->Add(h);
        }
    }

    TList* canvases2DSame = drawSame(histos_RecHit, histos_t1);
    //TList* canvases2DSame = drawSame(histos_RecHit, histos_t1, histos_t2);
    TList* canvases2D = draw2D(histos2D);

    const char* dirName = "gif";
    saveAllCanvasesToPicture(canvases2DSame, "gif", dirName);
    saveAllCanvasesToPicture(canvases2D, "gif", dirName);
}

void  gedPhotonMacros_RecHit2_Plot()
{
    const char* fileName = "/home/kaya/Desktop/gedPhotonResults_v6/gedPhotonMacros_RecHit2_2011_MB_750_hiForest.root";
    TFile* file = new TFile(fileName, "READ");
    std::cout << "input file for Histograms : " << file->GetName() << std::endl;

    const char* dirName = "/home/kaya/Desktop/gedPhotonResults_v6/gedPhotonMacros_RecHit2_2011_MB_750_hiForest";

    saveAllHistogramsToPicture(file, "gif", dirName, 1);

//    h->SetMarkerStyle(21);
//    h->SetMarkerColor(kBlack);
}

void  gedPhotonMacros_pfCand_Plot(const char* inputFileName, const char* outputDirName)
{
    TFile* file = new TFile(inputFileName, "READ");
    std::cout << "input file for Histograms : " << file->GetName() << std::endl;

    saveAllHistogramsToPicture(file, "gif", outputDirName, 1);

//    h->SetMarkerStyle(21);
//    h->SetMarkerColor(kBlack);
}

TList* drawSame(TList* histos1, TList* histos2, TList* histos3)
{
    TList* canvasList = new TList();

    TH2D* h1;
    TH2D* h2;
    TH2D* h3;
    TCanvas* c;
    TIter* iter1 = new TIter(histos1);
    TIter* iter2 = new TIter(histos2);
    TIter* iter3 = new TIter(histos3);
    while((h1=(TH2D*)iter1->Next())) {
        h2=(TH2D*)iter2->Next();
        h3=(TH2D*)iter3->Next();

        h1->SetMarkerStyle(34);
        h2->SetMarkerStyle(24);
        h3->SetMarkerStyle(24);

        h1->SetMarkerColor(kBlue);
        h2->SetMarkerColor(kBlack);
        h3->SetMarkerColor(kRed);

        h1->SetMarkerSize(5);
        h2->SetMarkerSize(2);
        h3->SetMarkerSize(2);

        c=new TCanvas(Form("%s_t1_t2",h1->GetName()));
        c->Clear();

        h1->Draw();
        h2->Draw("SAME");
        h3->Draw("SAME");

        // special cases
        {
            h1->SetStats(0);
            h2->SetStats(0);
            h3->SetStats(0);
        }

        canvasList->Add(c);
    }

    return canvasList;
}

TList* draw2D_vRecHit(TList* histos)
{
    TList* canvasList = new TList();

    TH2D* h;
    TCanvas* c;
    TIter* iter = new TIter(histos);
    while((h=(TH2D*)iter->Next())) {

        c=new TCanvas(h->GetName());

        h->Draw("colz");
//        c->SetLogz();

        canvasList->Add(c);
    }

    return canvasList;
}

/*
 * get list of all histograms under a directory "dir" for objects of a given "type"
 */
TList* getListOfALLHistograms2D(TDirectoryFile* dir)
{
    TList* histos=new TList();
    TList* keysHisto = getListOfALLKeys(dir, "TH2D");

    TIter* iter = new TIter(keysHisto);
    TKey*  key;
    while ((key=(TKey*)iter->Next()))
    {
        histos->Add((TH2D*)key->ReadObj());
    }

    return histos;
}

int main(int argc, char** argv)
{
//    gedPhotonMacros_RecHit_Plot();
//    gedPhotonMacros_RecHit2_Plot();

    const char* inputFileName = "/home/kaya/Desktop/gedPhotonResults_v6/gedPhotonMacros_PfCand_2011_MB_750_hiForest.root";
    const char* outputDirName = "/home/kaya/Desktop/gedPhotonResults_v6/gedPhotonMacros_PfCand_2011_MB_750_hiForest";
//    const char* inputFileName = "/home/kaya/Desktop/gedPhotonResults_v6/gedPhotonMacros_PfCand_HiForest_HIHighPt_photon30_HIRun2011-v1.root";
//    const char* outputDirName = "/home/kaya/Desktop/gedPhotonResults_v6/gedPhotonMacros_PfCand_HiForest_HIHighPt_photon30_HIRun2011-v1";
    gedPhotonMacros_pfCand_Plot(inputFileName, outputDirName);

    return 0;
}
