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
#include <stdint.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TPaveStats.h>


#include "/home/christoph/Documents/MasterThesis/Analyse/Utlities/utilities.h"
#include "/home/christoph/Documents/MasterThesis/Analyse/Utlities/dataReader.cpp"

using namespace std;

double mean(vector<double> x) {
  double result = 0;
  for (int i=0; i<x.size(); i++) {
    result+=x[i];
  }
  return result / x.size();
}

double std_vec(vector<double> x) {
  double result=0;
  double mu = mean(x);
  for (int i=0; i<x.size(); i++) {
    result += pow(x[i]-mu,2);
  }
  return sqrt(result/(x.size()-1));
}

double exponential(double *x, double *par) {
  return par[0] * exp(- x[0] / (par[1]) * log(2));
}

double exponential_off(double *x, double *par) {
  return par[0] * exp(- x[0] / (par[1]) * log(2)) + par[2];
}

double fConstant(double *x, double *par) {
  return par[0];
}

bool save_output = 1;
bool fitWithOffset = 1;
int startAtZero = 1;
int iStart = 0;
int pxSelect = 1;
bool singlePixel = false;
// bool mean_pixels = false;
int doAbsolute = 0;
vector<double> neutronRate0, gammaRate0;
double systErr = 0.03;

string dir, measurement, calibration_file;
vector<string> measurements;
vector<string> channelMapping = {"A1", "A2", "A3", "A4", "B1", "B2", "B3", "B4", "C1", "C2", "C3", "C4", "D1", "D2", "D3", "D4"};
vector<double> vPx = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
vector<double> vPx_err = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

vector<vector<double>> neutronRates, neutronRates_err, gammaRates, gammaRates_err;
vector<double> mean_neutronRates, mean_neutronRates_err, mean_gammaRates, mean_gammaRates_err;
vector<double> *nRate, *nRate_err, *gRate, *gRate_err;
double thickness_err;
vector<double> vThickness, vThickness_err = {};
vector<double> nResiduals, gResiduals = {};

vector<double> x12n, x12g, x0n, x0g, En, Eg;
vector<double> x12n_err, x12g_err, x0n_err, x0g_err, En_err, Eg_err;
double mean_x12n, mean_x12g, mean_x0n, mean_x0g, mean_En, mean_Eg;
double mean_x12n_err, mean_x12g_err, mean_x0n_err, mean_x0g_err, mean_En_err, mean_Eg_err;

vector<double> An, Ag, A0n, A0g;
vector<double> An_err, Ag_err, A0n_err, A0g_err;
double mean_An, mean_Ag, mean_A0n, mean_A0g;
double mean_An_err, mean_Ag_err, mean_A0n_err, mean_A0g_err;




// ./analyseAbsorption measurement fitWithOffset startAtZero
int main (int argc, char **argv) {

  TROOT _MakePad("MakePad", "Display data"); // initialise ROOT for access inside compiled c++ program
  TApplication theApp("App", NULL, NULL);

  // fit function
  TF1 *fitFunction = new TF1("fit1", exponential, -0.2, 11, 2);
  TF1 *fitFunction_off = new TF1("fit1_off", exponential_off, -0.2, 11, 3);
  fitFunction->SetLineWidth(1);
  fitFunction_off->SetLineWidth(1);
  char *fitSet = "+Q";

  // data directory
  dir = "/home/christoph/Documents/MasterThesis/Data/Absorption/";

  if ((string)argv[argc-1] == (string)"--help") {
    cout << "usage: ./analyseAbsorption measurement fitWithOffset startAtZero doAbsolute" << endl;
    cout << "or: ./analyseAbsorption measurement fitWithOffset startAtZero doAbsolute iPixel" << endl;
    return 0;
  }
  if (argc == 5) {
    measurement = (string) argv[argc-4];
    fitWithOffset = atoi(argv[argc-3]);
    startAtZero = atoi(argv[argc-2]);
    doAbsolute = atoi(argv[argc-1]);
  }
  else if (argc == 6) {
    // cout << "argc = " << argc << endl;
    cout << "fitting single pixel" << endl;
    measurement = (string) argv[argc-5];
    fitWithOffset = atoi(argv[argc-4]);
    startAtZero = atoi(argv[argc-3]);
    doAbsolute = atoi(argv[argc-2]);
    pxSelect = atoi(argv[argc-1]);
    singlePixel = true;
  }
  cout << "doAbsolute " << doAbsolute << endl;
  dir = dir + measurement + "/";
  cout << "fitWithOffset: " << fitWithOffset << endl;
  cout << "startAtZero: " << startAtZero << endl;
  cout << "measurement: " << measurement << endl;

  if (measurement=="H2O") {
    measurements = {"H2O1cm_20190617", "H2O2cm_20190617", "H2O3cm_20190617", "H2O4cm_20190617", "H2O5cm_20190617", "H2O6cm_20190618", "H2O7cm_20190618", "H2O8cm_20190618"};
    calibration_file = "Calibration_20190617";
    thickness_err = 0.1;
  }
  if (measurement=="Al") {
    measurements = {"Al1cm_20190611", "Al2cm_20190611", "Al3cm_20190611", "Al4cm_20190611", "Al5cm_20190611", "Al6cm_20190611",
      "Al7cm_20190611", "Al8cm_20190611", "Al9cm_20190611", "Al10cm_20190611"};
    calibration_file = "Calibration_20190611";
    // calibration_file = "Calibration_20190604_1";
    thickness_err = 0.01;
  }
  if (measurement=="Fe") {
    measurements = {"Fe_1cm_20190624", "Fe_2cm_20190624", "Fe_3cm_20190624", "Fe_4cm_20190624", "Fe_5cm_20190624", "Fe_6cm_20190624", "Fe_7cm_20190624", "Fe_8cm_20190624", "Fe_9cm_20190624", "Fe_10cm_20190624"};
    calibration_file = "Calibration_20190624";
    thickness_err = 0.01;
  }
  if (measurement=="Cu") {
    measurements = {"Cu_1cm_20190701", "Cu_2cm_20190701", "Cu_3cm_20190701", "Cu_4cm_20190701", "Cu_5cm_20190701", "Cu_6cm_20190701", "Cu_7cm_20190701", "Cu_8cm_20190701", "Cu_9cm_20190701", "Cu_10cm_20190701"};
    calibration_file = "Calibration_20190701";
    thickness_err = 0.01;
  }
  if (measurement=="Sand") {
    measurements = {"Sand_1cm_20190704", "Sand_2cm_20190704", "Sand_3cm_20190704", "Sand_4cm_20190704", "Sand_5cm_20190704", "Sand_6cm_20190704", "Sand_7cm_20190704", "Sand_8cm_20190704", "Sand_9cm_20190704"};
    calibration_file = "Calibration_20190703";
    thickness_err = 0.1;
  }
  if (measurement=="Paraffin") {
    measurements = {"Paraffin_1cm_20190705", "Paraffin_2cm_20190705", "Paraffin_3cm_20190705", "Paraffin_4cm_20190705", "Paraffin_5cm_20190705", "Paraffin_6cm_20190705", "Paraffin_7cm_20190705", "Paraffin_8cm_20190705", "Paraffin_9cm_20190705", "Paraffin_10cm_20190705", "Paraffin_11cm_20190705"};
    calibration_file = "Calibration_20190705";
    thickness_err = 0.03;
  }
  // calibration_file = "Calibration_20190617";

  // read calibration file data
  cout << "reading calibration file: " << dir+calibration_file+"/Results/"+calibration_file+"_Results.root" << endl;
  TFile *fCalib = new TFile((dir+calibration_file+"/Results/"+calibration_file+"_Results.root").c_str());
  neutronRates.push_back({});
  neutronRates_err.push_back({});
  gammaRates.push_back({});
  gammaRates_err.push_back({});
  fCalib->GetObject("neutronRates", nRate);
  fCalib->GetObject("neutronRates_err", nRate_err);
  fCalib->GetObject("gammaRates", gRate);
  fCalib->GetObject("gammaRates_err", gRate_err);
  for (int i=0; i<16; i++) {
    if (!doAbsolute) {
      neutronRates[0].push_back(1.);
      neutronRates_err[0].push_back(nRate_err->at(i) / nRate->at(i));
      gammaRates[0].push_back(1.);
      gammaRates_err[0].push_back(gRate_err->at(i) / gRate->at(i));
    }
    else {
      neutronRates[0].push_back(nRate->at(i));
      gammaRates[0].push_back(gRate->at(i));
      neutronRates_err[0].push_back(nRate_err->at(i));
      gammaRates_err[0].push_back(gRate_err->at(i));
    }
    neutronRate0.push_back(nRate->at(i));
    gammaRate0.push_back(gRate->at(i));
    // cout << "n err : " << neutronRates_err[0][0] << endl;
    // cout << "g err : " << gammaRates_err[0][0] << endl;
  }
  fCalib->Close();

  // read absorption file data
  for (int i=0; i<measurements.size(); i++) {
    cout << "reading absorption file: " << dir+measurements[i]+"/Images/"+measurements[i]+"_Results.root" << endl;
    neutronRates.push_back({});
    neutronRates_err.push_back({});
    gammaRates.push_back({});
    gammaRates_err.push_back({});

    TFile *fAbsorb = new TFile((dir+measurements[i]+"/Images/"+measurements[i]+"_Results.root").c_str());
    if (!doAbsolute) {
      fAbsorb->GetObject("neutronRates_subs", nRate);
      fAbsorb->GetObject("neutronRates_subs_err", nRate_err);
      fAbsorb->GetObject("gammaRates_subs", gRate);
      fAbsorb->GetObject("gammaRates_subs_err", gRate_err);
    }
    else {
      fAbsorb->GetObject("neutronRates_abs", nRate);
      fAbsorb->GetObject("neutronRates_abs_err", nRate_err);
      fAbsorb->GetObject("gammaRates_abs", gRate);
      fAbsorb->GetObject("gammaRates_abs_err", gRate_err);
    }
    // cout << i+1 << endl;
    for (int j=0; j<16; j++) {
      if (!doAbsolute) {
        neutronRates[i+1].push_back((nRate->at(j) + 100.) / 100.);
        neutronRates_err[i+1].push_back(sqrt( pow(nRate_err->at(j) / 100.,2) + pow(systErr,2) ));
        gammaRates[i+1].push_back((gRate->at(j) + 100.) / 100.);
        gammaRates_err[i+1].push_back(sqrt( pow(gRate_err->at(j) / 100.,2) + pow(systErr,2) ));
      }
      else {
        neutronRates[i+1].push_back(nRate->at(j));
        neutronRates_err[i+1].push_back(sqrt( pow(nRate_err->at(j),2) + pow(nRate->at(i)*systErr,2) ));
        gammaRates[i+1].push_back(gRate->at(j));
        gammaRates_err[i+1].push_back(sqrt( pow(gRate_err->at(j),2) + pow(gRate->at(i)*systErr,2) ));
      }
    }
    fAbsorb->Close();
  }

  // output plot and fitting
  // gStyle->SetOptFit(111);
  TCanvas *cFit = new TCanvas("cFit", "cFit", 900, 600);
  // TPad *pad1 = new TPad("pad1","pad1", 0, 0.3, 1, 1);
  // TPad *pad2 = new TPad("pad2","pad2", 0, 0.1, 1, 0.3);
  // pad1->SetBottomMargin(0.00001);
  // pad1->SetBorderMode(0);
  // pad2->SetTopMargin(0.00001);
  // pad2->SetBottomMargin(0.1);
  // pad2->SetBorderMode(0);
  // pad1->Draw();
  // pad2->Draw();
  // pad1->cd();

  TGraphErrors *grExpoN;
  TGraphErrors *grExpoG;

for (int iPixel=0; iPixel<16; iPixel++) {

  pxSelect = iPixel;


  vThickness.clear();
  vThickness_err.clear();
  mean_neutronRates.clear();
  mean_neutronRates_err.clear();
  mean_gammaRates.clear();
  mean_gammaRates_err.clear();


  // average rates and create thickness vectors
  if (!startAtZero)
  iStart = 1;
  for (int i=iStart; i<neutronRates.size(); i++) {
    vThickness.push_back(i);
    // if (iStart==0)
    // vThickness_err.push_back(0.);
    // else
    vThickness_err.push_back(thickness_err);

    // average over all pixels
    if (!singlePixel) {
      mean_neutronRates.push_back(mean(neutronRates[i]));
      mean_neutronRates_err.push_back(std_vec(neutronRates[i]));
      mean_gammaRates.push_back(mean(gammaRates[i]));
      mean_gammaRates_err.push_back(std_vec(gammaRates[i]));
    }
    // work with single pixel data
    else {
      mean_neutronRates.push_back(neutronRates[i][pxSelect]);
      mean_neutronRates_err.push_back(neutronRates_err[i][pxSelect] / neutronRate0[pxSelect]);
      mean_gammaRates.push_back(gammaRates[i][0]);
      mean_gammaRates_err.push_back(gammaRates_err[i][pxSelect] / gammaRate0[pxSelect]);
    }
  }


  // do the fitting
  grExpoN = new TGraphErrors(mean_neutronRates.size(), &vThickness[0], &mean_neutronRates[0], &vThickness_err[0], &mean_neutronRates_err[0]);
  grExpoG = new TGraphErrors(mean_neutronRates.size(), &vThickness[0], &mean_gammaRates[0], &vThickness_err[0], &mean_gammaRates_err[0]);
  // grExpoN->SetTitle(measurement.c_str());
  // if (singlePixel)
  grExpoN->SetTitle((measurement + " px" + to_string((int)iPixel)).c_str());
  grExpoG->SetMarkerColor(kRed);
  grExpoN->SetMarkerColor(kGreen);

  // cout << "\nfitting neutrons..." << endl;
  fitFunction->SetLineColor(kGreen);
  fitFunction_off->SetLineColor(kGreen);
  if (fitWithOffset) {
    if (!doAbsolute)
    fitFunction_off->SetParameters(1, 3.5, 0.1);
    else
    fitFunction_off->SetParameters(5, 3.5, 0.1);
    if (measurement == (string)"H2O")
    fitFunction->SetParameter(1,35);

    fitFunction_off->SetParNames("A^{n}", "x_{1/2}^{n} / cm", "A_{0}^{n}");
    grExpoN->Fit("fit1_off", fitSet, "", -0.2, vThickness.back()+0.2);

    x12n.push_back(fitFunction_off->GetParameter(1));
    x0n.push_back(x12n.back() / log(2));
    En.push_back(1 / x0n.back());
    x12n_err.push_back(fitFunction_off->GetParError(1));
    x0n_err.push_back(x12n_err.back() / log(2));
    En_err.push_back(1 / pow(x0n.back(),2) * x0n_err.back());

    An.push_back(fitFunction_off->GetParameter(0));
    A0n.push_back(fitFunction_off->GetParameter(2));
    An_err.push_back(fitFunction_off->GetParError(0));
    A0n_err.push_back(fitFunction_off->GetParError(2));

    for (int i=0; i<mean_neutronRates.size(); i++) {
      // nResiduals.push_back( (mean_neutronRates[i] - fitFunction_off->Eval(vThickness[i]) ) / mean_neutronRates_err[i]);
      nResiduals.push_back( (mean_neutronRates[i] - fitFunction_off->Eval(vThickness[i]) ));
    }
  }
  else {
    if (!doAbsolute)
    fitFunction->SetParameters(1, 3.5);
    else
    fitFunction->SetParameters(5, 3.5);

    fitFunction->SetParNames("A^{n}", "x_{1/2}^{n} / cm");
    grExpoN->Fit("fit1", fitSet, "", -0.2, vThickness.back()+0.2);

    x12n.push_back(fitFunction->GetParameter(1));
    x0n.push_back(x12n.back() / log(2));
    En.push_back(1 / x0n.back());
    x12n_err.push_back(fitFunction->GetParError(1));
    x0n_err.push_back(x12n_err.back() / log(2));
    En_err.push_back(1 / pow(x0n.back(),2) * x0n_err.back());

    An.push_back(fitFunction->GetParameter(0));
    An_err.push_back(fitFunction->GetParError(0));


    for (int i=0; i<mean_neutronRates.size(); i++)
    // nResiduals.push_back( (mean_neutronRates[i] - fitFunction->Eval(vThickness[i]) ) / mean_neutronRates_err[i]);
    nResiduals.push_back( (mean_neutronRates[i] - fitFunction->Eval(vThickness[i]) ));
  }

  // cout << "x_1/2=\t" << x12n << " +- " << x12n_err << endl;
  // cout << "x_0=\t" << x0n << " +- " << x0n_err << endl;
  // cout << "E=\t" << En << " +- " << En_err << endl;



  // cout << "\nfitting gammas..." << endl;
  fitFunction->SetLineColor(kRed);
  fitFunction_off->SetLineColor(kRed);
  if (measurement == (string)"H2O")
    fitFunction->SetParameter(1,35);
  if (fitWithOffset) {
    if (!doAbsolute)
    fitFunction_off->SetParameters(1, 3.5, 0.1);
    else
    fitFunction_off->SetParameters(5, 3.5, 0.1);
    if (measurement == (string)"H2O")
    fitFunction->SetParameter(1,35);

    fitFunction_off->SetParNames("A^{g}", "x_{1/2}^{g} / cm", "A_{0}^{g}");
    grExpoG->Fit("fit1_off", fitSet, "", -0.2, vThickness.back()+0.2);

    x12g.push_back(fitFunction_off->GetParameter(1));
    x0g.push_back(x12g.back() / log(2));
    Eg.push_back(1. / x0g.back());
    x12g_err.push_back(fitFunction_off->GetParError(1));
    x0g_err.push_back(x12g_err.back() / log(2));
    Eg_err.push_back(1 / pow(x0g.back(),2) * x0g_err.back());

    Ag.push_back(fitFunction_off->GetParameter(0));
    A0g.push_back(fitFunction_off->GetParameter(2));
    Ag_err.push_back(fitFunction_off->GetParError(0));
    A0g_err.push_back(fitFunction_off->GetParError(2));

    // cout << "size: " << mean_gammaRates.size() << endl;
    for (int i=0; i<mean_gammaRates.size(); i++)
    // gResiduals.push_back( (mean_gammaRates[i] - fitFunction_off->Eval(vThickness[i])) / mean_gammaRates_err[i] );
    gResiduals.push_back( (mean_gammaRates[i] - fitFunction_off->Eval(vThickness[i])));
  }
  else {
    if (!doAbsolute)
    fitFunction->SetParameters(1, 3.5);
    else
    fitFunction->SetParameters(5, 3.5);

    fitFunction->SetParNames("A^{g}", "x_{1/2}^{g} / cm");
    grExpoG->Fit("fit1", fitSet, "", -0.2, vThickness.back()+0.2);

    x12g.push_back(fitFunction->GetParameter(1));
    x0g.push_back(x12g.back() / log(2));
    Eg.push_back(1. / x0g.back());
    x12g_err.push_back(fitFunction->GetParError(1));
    x0g_err.push_back(x12g_err.back() / log(2));
    Eg_err.push_back(1 / pow(x0g.back(),2) * x0g_err.back());

    Ag.push_back(fitFunction_off->GetParameter(0));
    Ag_err.push_back(fitFunction_off->GetParError(0));




    for (int i=0; i<mean_gammaRates.size(); i++)
    // gResiduals.push_back( (mean_gammaRates[i] - fitFunction->Eval(vThickness[i])) / mean_gammaRates_err[i] );
    gResiduals.push_back( (mean_gammaRates[i] - fitFunction->Eval(vThickness[i])));
  }

  // cout << "x_1/2=\t" << x12g << " +- " << x12g_err << endl;
  // cout << "x_0=\t" << x0g << " +- " << x0g_err << endl;
  // cout << "E=\t" << Eg << " +- " << Eg_err << endl;


}



  // cosmetics
  grExpoN->GetXaxis()->SetNdivisions(20, 5, 0, kTRUE);
  if (measurement==(string)"H2O")
    grExpoN->GetXaxis()->SetNdivisions(10, 5, 0, kTRUE);
  grExpoG->SetMarkerStyle(20);
  grExpoG->SetMarkerSize(0.7);
  grExpoN->SetMarkerStyle(20);
  grExpoN->SetMarkerSize(0.7);
  grExpoN->GetXaxis()->SetTitle("thickness / cm");
  grExpoN->GetYaxis()->SetTitle("relative rate");
  if (doAbsolute)
    grExpoN->GetYaxis()->SetTitle("absolute rate / Hz");
  //
  grExpoN->Draw("AP");
  grExpoG->Draw("P same");
  grExpoN->GetXaxis()->SetRangeUser(-0.5, 10.5);
  if (doAbsolute)
    grExpoN->GetYaxis()->SetRangeUser(.2, 7);
  else
    grExpoN->GetYaxis()->SetRangeUser(0., 1.1);

  // TText *tn = new TText(.5, 1.,"Neutrons");
  // tn->SetTextColor(kGreen);
  // tn->SetTextFont(43);
  // tn->SetTextSize(15);
  // tn->Draw();
  //
  // TText *tg = new TText(.5, 1.5,"Gammas");
  // tg->SetTextColor(kRed);
  // tg->SetTextFont(43);
  // tg->SetTextSize(15);
  // tg->Draw();

  cFit->Modified();
  cFit->Update();
  //
  //
  //
  // // residual plot
  // pad2->cd();
  // gPad->SetFillStyle(0);
  // TGraphErrors *grResidualN = new TGraphErrors(nResiduals.size(), &vThickness[0], &nResiduals[0], &vThickness_err[0], &mean_neutronRates_err[0]);
  // TGraphErrors *grResidualG = new TGraphErrors(gResiduals.size(), &vThickness[0], &gResiduals[0], &vThickness_err[0], &mean_gammaRates_err[0]);
  // grResidualN->SetTitle("");
  // grResidualG->SetTitle("");
  // // grResidualN->GetYaxis()->SetTitle("residuals");
  // grResidualN->GetXaxis()->SetLabelSize(0.1);
  // grResidualN->GetYaxis()->SetLabelSize(0.1);
  // grResidualN->SetMarkerStyle(20);
  // grResidualG->SetMarkerStyle(20);
  // grResidualN->SetMarkerColor(kRed);
  // grResidualG->SetMarkerColor(kGreen);
  // grResidualN->SetMarkerSize(0.7);
  // grResidualG->SetMarkerSize(0.7);
  //
  // grResidualN->GetYaxis()->SetNdivisions(5, 3, 0, kTRUE);
  // grResidualN->GetXaxis()->SetNdivisions(20, 5, 0, kTRUE);
  // if (measurement==(string)"H2O")
  //   grResidualN->GetXaxis()->SetNdivisions(10, 5, 0, kTRUE);
  //
  //
  // grResidualN->Draw("AP");
  // grResidualG->Draw("P same");
  // // grResidualG->Draw("AP");
  //
  // grResidualN->GetXaxis()->SetRangeUser(-0.5, 10.5);
  // // grResidualN->GetYaxis()->SetRangeUser(-3, 3);
  // grResidualN->GetYaxis()->SetRangeUser(-0.3, 0.3);
  // cFit->Modified();
  // cFit->Update();
  //
  // // TLine *l1 = new TLine(0., 0, grResidualN->GetXaxis()->GetXmax(), 0);
  // TLine *l1 = new TLine(0., 0, 10, 0);
  // if (!startAtZero) {
  //   delete l1;
  //   TLine *l1 = new TLine(0., 0, grResidualN->GetXaxis()->GetXmax()-0.5, 0);
  // }
  // l1->SetLineStyle(7);
  // l1->SetLineColor(kRed);
  // l1->Draw("");
  //
  // // residuals axis label
  // cFit->cd();
  // TText *t = new TText(.8,.05,"thickness / cm");
  // t->SetNDC(true);
  // t->SetTextFont(43);
  // t->SetTextSize(15);
  // t->Draw();
  //
  // TText *t1 = new TText(.06, .14,"residuals / Hz");
  // t1->SetNDC(true);
  // t1->SetTextFont(43);
  // t1->SetTextSize(15);
  // t1->SetTextAngle(90);
  // t1->Draw();
  // cFit->Update();

  mean_x12n = mean(x12n);
  mean_x12g = mean(x12g);
  mean_x0n = mean(x0n);
  mean_x0g = mean(x0g);
  mean_En = mean(En);
  mean_Eg = mean(Eg);

  mean_An = mean(An);
  mean_Ag = mean(Ag);
  mean_A0n = mean(A0n);
  mean_A0g = mean(A0g);

  mean_x12n_err = std_vec(x12n_err);
  mean_x12g_err = std_vec(x12g_err);
  mean_x0n_err = std_vec(x0n_err);
  mean_x0g_err = std_vec(x0g_err);
  mean_En_err = std_vec(En_err);
  mean_Eg_err = std_vec(Eg_err);

  mean_An_err = std_vec(An);
  mean_Ag_err = std_vec(Ag);
  mean_A0n_err = std_vec(A0n);
  mean_A0g_err = std_vec(A0g);

  cout << "mean_x12n\t=\t" << mean_x12n << " +- " << mean_x12n_err << endl;
  cout << "mean_x12g\t=\t" << mean_x12g << " +- " << mean_x12g_err << endl;
  cout << "mean_x0n\t=\t" << mean_x0n << " +- " << mean_x0n_err << endl;
  cout << "mean_x0g\t=\t" << mean_x0g << " +- " << mean_x0g_err << endl;
  cout << "mean_En\t=\t" << mean_En << " +- " << mean_En_err << endl;
  cout << "mean_Eg\t=\t" << mean_Eg << " +- " << mean_Eg_err << endl;


  TF1 *fConst = new TF1("fitConst", fConstant, -0.2, 11, 1);
  fConst->SetLineWidth(1);


  TCanvas *cAverage = new TCanvas("cAverage", "cAverage", 900, 600);
  gStyle->SetOptFit(111);

  // TH1F *hResultsN = new TH1F("hResN", "hResN", 50, 2, 4);
  // TH1F *hResultsG = new TH1F("hResG", "hResG", 50, 2, 4);
  //
  // hResultsN->SetLineColor(kGreen);
  // hResultsG->SetLineColor(kRed);
  //
  // for (int i=0; i<16; i++) {
  //   hResultsN->Fill(x12n[i]);
  //   hResultsG->Fill(x12g[i]);
  // }
  // hResultsN->Draw();
  // hResultsG->Draw("same");

  TGraphErrors *grx12n = new TGraphErrors(x12n.size(), &vPx[0], &x12n[0], &vPx_err[0], &x12n_err[0]);
  TGraphErrors *grx12g = new TGraphErrors(x12g.size(), &vPx[0], &x12g[0], &vPx_err[0], &x12g_err[0]);

  // TGraphErrors *grx12n = new TGraphErrors(En.size(), &vPx[0], &En[0], &vPx_err[0], &En_err[0]);
  // TGraphErrors *grx12g = new TGraphErrors(Eg.size(), &vPx[0], &Eg[0], &vPx_err[0], &Eg_err[0]);

  // TGraphErrors *grx12n = new TGraphErrors(En.size(), &vPx[0], &An[0], &vPx_err[0], &An_err[0]);
  // TGraphErrors *grx12g = new TGraphErrors(Eg.size(), &vPx[0], &Ag[0], &vPx_err[0], &Ag_err[0]);

  // TGraphErrors *grx12n = new TGraphErrors(En.size(), &vPx[0], &A0n[0], &vPx_err[0], &A0n_err[0]);
  // TGraphErrors *grx12g = new TGraphErrors(Eg.size(), &vPx[0], &A0g[0], &vPx_err[0], &A0g_err[0]);



  fConst->SetLineColor(kGreen);
  grx12n->Fit("fitConst", "", "+Q", -0.5, 15.5);
  fConst->SetLineColor(kRed);
  grx12g->Fit("fitConst", "", "+Q", -0.5, 15.5);

  grx12n->SetTitle(("Fit Results " + measurement).c_str());
  grx12n->GetXaxis()->SetTitle("pixel");
  grx12n->GetYaxis()->SetTitle("x_{1/2}");
  grx12n->SetMarkerStyle(20);
  grx12g->SetMarkerStyle(20);
  grx12n->SetMarkerColor(kGreen);
  grx12g->SetMarkerColor(kRed);
  grx12n->Draw("AP");
  grx12g->Draw("P same");

  TText *tn = new TText(.5, 4.,"Neutrons");
  tn->SetTextColor(kGreen);
  tn->SetTextFont(43);
  tn->SetTextSize(15);
  tn->Draw();

  TText *tg = new TText(.5, 3,"Gammas");
  tg->SetTextColor(kRed);
  tg->SetTextFont(43);
  tg->SetTextSize(15);
  tg->Draw();

  grx12n->GetXaxis()->SetRangeUser(-0.5,15.5);
  grx12n->GetYaxis()->SetRangeUser(1, 8);
  if (measurement == (string)"H2O" || measurement == (string)"Sand")
    grx12n->GetYaxis()->SetRangeUser(1, 50);
  // grx12n->GetYaxis()->SetRangeUser(En[0]-0.05, En[0]+0.06);

  cAverage->Modified();
  cAverage->Update();




  theApp.Run();

  return 0;
}
