#include <math.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <functional>
#include <unistd.h>
#include <dirent.h>
#include <stdint.h>
#include <sys/stat.h>
#include <cstdlib>
#include <fstream>

#include <TROOT.h>
#include <TApplication.h>
#include <TFile.h>
#include <TChain.h>
#include <TVector.h>
#include <TLine.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>

#include "filereader.h"
#include "utilities.h"
#include "stylesheet.h"

using namespace std;

bool save_output = 1;

string dir, filename;

//////////////////////////////////////////////////////////////////////////////
// main loop
//////////////////////////////////////////////////////////////////////////////

int main() {

  TROOT _MakePad("MakePad", "Display data"); // initialise ROOT for access inside compiled c++ program
  TApplication theApp("App", NULL, NULL);

  int N = 200; // n bins

  dir = "/home/christoph/Documents/MasterThesis/Data/CoMPASS/";
  // dir = "/home/christoph/Documents/MasterThesis/Data/CoMPASS/Stolberg/";
  filename = "onboard_PSD_test_1";
  // filename = "ABside549_CDside539";
  // filename = "Messung3";
  // filename = "onboard_PSD_AmBe";
  // cout << "Reading: " << dir+filename << endl;
  // TFile *f = new TFile((dir+filename+(string)"_PSD.root").c_str());


  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  // load data from algorithm
  TFile *f = new TFile((dir+filename+"/Plots/"+filename+"_25_270_PSD_adc.root").c_str());
  if ( f->IsOpen() ) printf("File opened successfully\n");

  vector<double> *vPeak;
  vector<double> *vPSD;
  vector<int> *vQS;
  vector<int> *vQL;
  vector<int> *vValid;
  f->GetObject("vPeak", vPeak);
  f->GetObject("vPSD", vPSD);
  f->GetObject("vQS", vQS);
  f->GetObject("vQL", vQL);
  f->GetObject("vValid", vValid);

  TH2F *h2 = new TH2F("h2d","h2d", N, 0, 2e5, N, 0, 1);
  // TH2F *h2 = new TH2F("h2d","h2d", N, 0, 2e3, N, 0, 1);
  TH1F *hE = new TH1F("h1E","h1E", 300, 0, 1e4);
  TH1F *hPSD = new TH1F("h1PSD","h1PSD", 300, 0, 1);
  double mVToADC = 16.384/2.;

  for (int i=0; i<vPeak->size(); i++) {
    // if (vValid->at(i)==1) {
      // h2->Fill(vPeak->at(i), vPSD->at(i));
      h2->Fill(vQL->at(i), vPSD->at(i));
      hE->Fill(vQL->at(i));
      hPSD->Fill(vPSD->at(i));
    // }
  }

  // // run pulse analyser
  // output data = filereader(filename, DeltaS, DeltaL, false);
  // vector<double> vPeak, vPSD, vQS, vQL;
  // vector<int> vValid;
  // vPeak = data.vPeak;
  // vPSD = data.vPSD;
  // vQS = data.vQS;
  // vQL = data.vQL;
  // vValid = data.vValid;
  //
  // // fill histogram
  // for (int i=0; i<vPSD.size(); i++) {
  //   h2d->Fill(vPeak[i], vPSD[i]);
  //   h1->Fill(vPeak[i]);
  // }

  TCanvas *cAlg = new TCanvas("cAlg", "cAlg", 900, 600);
  h2->SetTitle("");
  h2->Draw("colz");
  h2->GetXaxis()->SetTitle("long Int / ADC");
  h2->GetYaxis()->SetTitle("PSD");
  cAlg->Print(("Plots/"+filename+"psd2d_alg.png").c_str());


  // plot 1d spectrum
  TCanvas *cE = new TCanvas("cEnergy", "cEnergy", 900, 600);
  cE->SetLogy();
  gStyle->SetOptStat(10);
  hE->SetTitle("Energy Spectrum");
  // h1->GetXaxis()->SetTitle("PH / mV");
  hE->GetXaxis()->SetTitle("long Int / ADC");
  hE->GetYaxis()->SetTitle("Counts");
  hE->Draw();
  cE->Print(("Plots/"+filename+"E_alg.png").c_str());


  // plot 1d spectrum
  TCanvas *cPSD = new TCanvas("cPSD", "cPSD", 900, 600);
  cPSD->SetLogy();
  gStyle->SetOptStat(10);
  hPSD->SetTitle("PSD Spectrum");
  // h1->GetXaxis()->SetTitle("PH / mV");
  hPSD->GetXaxis()->SetTitle("PSD");
  hPSD->GetYaxis()->SetTitle("Counts");
  hPSD->Draw();
  cPSD->Print(("Plots/"+filename+"psd_alg.png").c_str());



  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  // load data from CAEN board
  TChain *tree = new TChain("Events");
	tree->Add((dir+filename+"/"+filename+"_0_0.root").c_str());
	if(tree == nullptr) {
		cout << "Tree could not be loaded. Closing program." << endl;
		exit(EXIT_FAILURE);
	}
  uint16_t energy;
  uint16_t energyShort;
  tree->SetBranchAddress("energy", &energy);
  tree->SetBranchAddress("energyShort", &energyShort);
  int treeSize;
	treeSize = tree->GetEntries();

  TH2F *h2_caen = new TH2F("h2d_caen","h2d_caen", 200, 0, 2e5, 200, 0, 1);
  TH1F *hE_caen = new TH1F("h1E_caen","h1E_caen", 2000, 0, 3e3);
  TH1F *hPSD_caen = new TH1F("h1PSD_caen","h1PSD_caen", 200, 0, 1);

	for(int iEvent=0; iEvent < treeSize; ++iEvent) {
		tree->GetEntry(iEvent);
    h2_caen->Fill((double)energy, ((double)energy-(double)energyShort)/(double)energy);
    hE_caen->Fill(energy);
    // hPSD_caen->Fill(energyShort);
    hPSD_caen->Fill(((double)energy-(double)energyShort)/(double)energy);
  }

  TCanvas *cCAEN = new TCanvas("cCAEN", "cCAEN", 900, 600);
  h2_caen->SetTitle("");
  h2_caen->Draw("colz");
  h2_caen->GetXaxis()->SetTitle("long Int / ADC");
  h2_caen->GetYaxis()->SetTitle("PSD");
  cCAEN->Print(("Plots/"+filename+"psd2d_caen.png").c_str());


  // plot 1d spectrum
  TCanvas *cE_caen = new TCanvas("cEnergy_caen", "cEnergy_caen", 900, 600);
  cE_caen->SetLogy();
  gStyle->SetOptStat(10);
  hE_caen->SetTitle("Energy Spectrum");
  // h1->GetXaxis()->SetTitle("PH / mV");
  hE_caen->GetXaxis()->SetTitle("long Int / ADC");
  hE_caen->GetYaxis()->SetTitle("Counts");
  hE_caen->Draw();
  cE_caen->Print(("Plots/"+filename+"E_caen.png").c_str());


  // plot 1d spectrum
  TCanvas *cPSD_caen = new TCanvas("cPSD_caen", "cPSD_caen", 900, 600);
  cPSD_caen->SetLogy();
  gStyle->SetOptStat(10);
  hPSD_caen->SetTitle("PSD Spectrum");
  hPSD_caen->GetXaxis()->SetTitle("PSD");
  hPSD_caen->GetYaxis()->SetTitle("Counts");
  hPSD_caen->Draw();
  cPSD_caen->Print(("Plots/"+filename+"psd_caen.png").c_str());




  theApp.Run();

  return 0;
}
