
#include <stdio.h>
//#include <iostream>
#include <math.h>
#include <TTree.h>
#include <TH1.h>
#include <TH1D.h>
#include <TROOT.h>
#include <TFile.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TRandom.h>

double gaussianRandom(double average, double stdev) {
        double v1, v2, s, temp;
        do {
                v1 =  2 * ((double) rand() / RAND_MAX) - 1;      // -1.0 ~ 1.0 ?~L?~@?~]~X ?~R
                v2 =  2 * ((double) rand() / RAND_MAX) - 1;      // -1.0 ~ 1.0 ?~L?~@?~]~X ?~R
                s = v1 * v1 + v2 * v2;
        } while (s >= 1 || s == 0);

        s = sqrt( (-2 * log(s)) / s );
        temp = v1*s;
        temp = (stdev*temp)+average;
        return temp;
}

void GausRandTest(){
	TH1D* h1 = new TH1D("h1", "h1", 100, -10, 10);
	srand((unsigned int)time(NULL));
	for(int i=0;i<1000;i++){
		
		h1 -> Fill(gaussianRandom(0,2));
	//	h1 -> Fill(gRandom->Gaus(0,2));
	
	}
	cout << "dlsfjsdk" << endl;
	TCanvas* c1 = new TCanvas();
	h1 -> Draw();
}

