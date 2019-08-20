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

bool save_output = 1;
bool fitWithOffset = 1;
int startAtZero = 1;
int iStart = 0;
int pxSelect = 1;
bool singlePixel = false;
// bool mean_pixels = false;
int doAbsolute = 0;
vector<double> neutronRate0, gammaRate0;
double systErr = 0.02;

string dir, measurement, calibration_file;
vector<string> measurements;
vector<string> channelMapping = {"A1", "A2", "A3", "A4", "B1", "B2", "B3", "B4", "C1", "C2", "C3", "C4", "D1", "D2", "D3", "D4"};

vector<vector<double>> neutronRates, neutronRates_err, gammaRates, gammaRates_err;
vector<double> mean_neutronRates, mean_neutronRates_err, mean_gammaRates, mean_gammaRates_err;
vector<double> *nRate, *nRate_err, *gRate, *gRate_err;
double thickness_err;
vector<double> vThickness, vThickness_err = {};
vector<double> nResiduals, gResiduals = {};

double x12n, x12g, x0n, x0g, En, Eg;
double x12n_err, x12g_err, x0n_err, x0g_err, En_err, Eg_err;



// ./analyseAbsorption measurement fitWithOffset startAtZero
int main (int argc, char **argv) {

  TROOT _MakePad("MakePad", "Display data"); // initialise ROOT for access inside compiled c++ program
  TApplication theApp("App", NULL, NULL);

  // fit function
  TF1 *fitFunction = new TF1("fit1", exponential, -0.2, 11, 2);
  TF1 *fitFunction_off = new TF1("fit1_off", exponential_off, -0.2, 11, 3);
  fitFunction->SetLineWidth(1);
  fitFunction_off->SetLineWidth(1);
  char *fitSet = "+";

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
  // thickness_err = 0.5;

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
  gStyle->SetOptFit(111);
  // TStyle* mcStyle = new TStyle("mcStyle","Style 1");
  // mcStyle->SetTitleSize(0.035,"xyz");
  // mcStyle->SetOptFit(111);
  // gROOT->SetStyle("mcStyle");

  // gStyle->SetFitFormat("%.2f");
  TCanvas *cFit = new TCanvas("cFit", "cFit", 900, 600);
  TPad *pad1 = new TPad("pad1","pad1", 0, 0.3, 1, 1);
  TPad *pad2 = new TPad("pad2","pad2", 0, 0.1, 1, 0.3);
  pad1->SetBottomMargin(0.00001);
  pad1->SetBorderMode(0);
  pad2->SetTopMargin(0.00001);
  pad2->SetBottomMargin(0.1);
  pad2->SetBorderMode(0);
  pad1->Draw();
  pad2->Draw();
  pad1->cd();



  // average rates and create thickness vectors
  if (!startAtZero)
    iStart = 1;
  for (int i=iStart; i<neutronRates.size(); i++) {
    vThickness.push_back(i);
    // if (iStart==0)
    //   vThickness_err.push_back(0.);
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
  TGraphErrors *grExpoN = new TGraphErrors(mean_neutronRates.size(), &vThickness[0], &mean_neutronRates[0], &vThickness_err[0], &mean_neutronRates_err[0]);
  TGraphErrors *grExpoG = new TGraphErrors(mean_neutronRates.size(), &vThickness[0], &mean_gammaRates[0], &vThickness_err[0], &mean_gammaRates_err[0]);
  grExpoN->SetTitle(measurement.c_str());
  if (singlePixel)
    grExpoN->SetTitle((measurement + " px" + to_string((int)pxSelect)).c_str());
  grExpoG->SetMarkerColor(kRed);
  grExpoN->SetMarkerColor(kGreen);

  cout << "\nfitting neutrons..." << endl;
  fitFunction->SetLineColor(kGreen);
  fitFunction_off->SetLineColor(kGreen);
  if (fitWithOffset) {
    if (!doAbsolute)
    fitFunction_off->SetParameters(1, 3.5, 0.1);
    else
      fitFunction_off->SetParameters(5, 3.5, 0.1);

    fitFunction_off->SetParNames("A^{n}", "x_{1/2}^{n} / cm", "A_{0}^{n}");
    grExpoN->Fit("fit1_off", fitSet, "", -0.2, vThickness.back()+0.2);

    x12n = fitFunction_off->GetParameter(1);
    x0n = x12n / log(2);
    En = 1 / x0n;
    x12n_err = fitFunction_off->GetParError(1);
    x0n_err = x12n_err / log(2);
    En_err = 1 / pow(x0n,2) * x0n_err;

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

    x12n = fitFunction->GetParameter(1);
    x0n = x12n / log(2);
    En = 1 / x0n;
    x12n_err = fitFunction_off->GetParError(1);
    x0n_err = x12n_err / log(2);
    En_err = 1 / pow(x0n,2) * x0n_err;

    for (int i=0; i<mean_neutronRates.size(); i++)
    // nResiduals.push_back( (mean_neutronRates[i] - fitFunction->Eval(vThickness[i]) ) / mean_neutronRates_err[i]);
      nResiduals.push_back( (mean_neutronRates[i] - fitFunction->Eval(vThickness[i]) ));
  }

  cout << "x_1/2=\t" << x12n << " +- " << x12n_err << endl;
  cout << "x_0=\t" << x0n << " +- " << x0n_err << endl;
  cout << "E=\t" << En << " +- " << En_err << endl;



  cout << "\nfitting gammas..." << endl;
  fitFunction->SetLineColor(kRed);
  fitFunction_off->SetLineColor(kRed);
  if (fitWithOffset) {
    fitFunction_off->SetParNames("A^{g}", "x_{1/2}^{g} / cm", "A_{0}^{g}");
    grExpoG->Fit("fit1_off", fitSet, "", -0.2, vThickness.back()+0.2);

    x12g = fitFunction_off->GetParameter(1);
    x0g = x12g / log(2);
    Eg = 1. / x0g;
    x12g_err = fitFunction_off->GetParError(1);
    x0g_err = x12g_err / log(2);
    Eg_err = 1 / pow(x0g,2) * x0g_err;

    // cout << "size: " << mean_gammaRates.size() << endl;
    for (int i=0; i<mean_gammaRates.size(); i++)
    // gResiduals.push_back( (mean_gammaRates[i] - fitFunction_off->Eval(vThickness[i])) / mean_gammaRates_err[i] );
      gResiduals.push_back( (mean_gammaRates[i] - fitFunction_off->Eval(vThickness[i])));
  }
  else {
    fitFunction->SetParNames("A^{g}", "x_{1/2}^{g} / cm");
    grExpoG->Fit("fit1", fitSet, "", -0.2, vThickness.back()+0.2);

    x12g = fitFunction->GetParameter(1);
    x0g = x12g / log(2);
    Eg = 1. / x0g;
    x12g_err = fitFunction_off->GetParError(1);
    x0g_err = x12g_err / log(2);
    Eg_err = 1 / pow(x0g,2) * x0g_err;

    for (int i=0; i<mean_gammaRates.size(); i++)
    // gResiduals.push_back( (mean_gammaRates[i] - fitFunction->Eval(vThickness[i])) / mean_gammaRates_err[i] );
      gResiduals.push_back( (mean_gammaRates[i] - fitFunction->Eval(vThickness[i])));
  }

  cout << "x_1/2=\t" << x12g << " +- " << x12g_err << endl;
  cout << "x_0=\t" << x0g << " +- " << x0g_err << endl;
  cout << "E=\t" << Eg << " +- " << Eg_err << endl;


  // // ge = (TGraphErrors *)c1->GetListOfPrimitives()->FindObject(measurement.c_str());
  // // ge->Fit("gaus");
  // cFit->Modified();
  // cFit->Update();
  // // TPaveText *pt = (TPaveText *)cFit->GetListOfPrimitives()->FindObject(measurement.c_str());
  // // pt->SetX1NDC(0.4);
  // // pt->SetX2NDC(0.7);
  // // pt->SetY1NDC(0.15);
  // // pt->SetY2NDC(0.25);
  // TPaveStats *ps = (TPaveStats *)grExpoN->GetListOfFunctions()->FindObject("stats");
  // ps->SetX1NDC(0.15);
  // ps->SetX2NDC(0.55);
  // cFit->Modified();
  // cFit->Update();

  // gStyle->SetStatX(0.3);
  // gPad->Update();
  // auto stat = dynamic_cast<TPaveStats*>(grExpoN->GetHistogram()->FindObject("stats"));
  // if (stat) {
  //   std::cout << " X1NDC: " << stat->GetX1NDC() << " X2NDC: " << stat->GetX2NDC() << std::endl;
  //   stat->SetX1NDC(0.1); stat->SetX2NDC(0.3);
  //   stat->Draw();
  // } else {
  //   cerr << "No stats box found!\n";
  // }

  // cosmetics
  grExpoN->GetXaxis()->SetNdivisions(20, 5, 0, kTRUE);
  if (measurement==(string)"H2O" || measurement == (string)"Sand"  || measurement == (string)"Paraffin")
    grExpoN->GetXaxis()->SetNdivisions(10, 5, 0, kTRUE);
  grExpoG->SetMarkerStyle(20);
  grExpoG->SetMarkerSize(0.7);
  grExpoN->SetMarkerStyle(20);
  grExpoN->SetMarkerSize(0.7);
  grExpoN->GetXaxis()->SetTitle("thickness / cm");
  grExpoN->GetYaxis()->SetTitle("relative rate");
  if (doAbsolute)
    grExpoN->GetYaxis()->SetTitle("absolute rate / Hz");

  grExpoN->Draw("AP");
  grExpoG->Draw("P same");
  grExpoN->GetXaxis()->SetRangeUser(-0.5, 10.5);
  if (doAbsolute)
    grExpoN->GetYaxis()->SetRangeUser(.2, 7);
  else
    grExpoN->GetYaxis()->SetRangeUser(0., 1.1);

  TText *tn = new TText(.5, 1.,"Neutrons");
  tn->SetTextColor(kGreen);
  tn->SetTextFont(43);
  tn->SetTextSize(15);
  tn->Draw();

  TText *tg = new TText(.5, 1.5,"Gammas");
  tg->SetTextColor(kRed);
  tg->SetTextFont(43);
  tg->SetTextSize(15);
  tg->Draw();

  cFit->Modified();
  cFit->Update();



  // residual plot
  pad2->cd();
  gPad->SetFillStyle(0);
  TGraphErrors *grResidualN = new TGraphErrors(nResiduals.size(), &vThickness[0], &nResiduals[0], &vThickness_err[0], &mean_neutronRates_err[0]);
  TGraphErrors *grResidualG = new TGraphErrors(gResiduals.size(), &vThickness[0], &gResiduals[0], &vThickness_err[0], &mean_gammaRates_err[0]);
  grResidualN->SetTitle("");
  grResidualG->SetTitle("");
  // grResidualN->GetYaxis()->SetTitle("residuals");
  grResidualN->GetXaxis()->SetLabelSize(0.1);
  grResidualN->GetYaxis()->SetLabelSize(0.1);
  grResidualN->SetMarkerStyle(20);
  grResidualG->SetMarkerStyle(20);
  grResidualN->SetMarkerColor(kRed);
  grResidualG->SetMarkerColor(kGreen);
  grResidualN->SetMarkerSize(0.7);
  grResidualG->SetMarkerSize(0.7);

  grResidualN->GetYaxis()->SetNdivisions(5, 3, 0, kTRUE);
  grResidualN->GetXaxis()->SetNdivisions(20, 5, 0, kTRUE);
  if (measurement==(string)"H2O" || measurement == (string)"Sand" || measurement == (string)"Paraffin")
    grResidualN->GetXaxis()->SetNdivisions(10, 5, 0, kTRUE);


  grResidualN->Draw("AP");
  grResidualG->Draw("P same");
  // grResidualG->Draw("AP");

  grResidualN->GetXaxis()->SetRangeUser(-0.5, 10.5);
  // grResidualN->GetYaxis()->SetRangeUser(-3, 3);
  grResidualN->GetYaxis()->SetRangeUser(-0.3, 0.3);
  cFit->Modified();
  cFit->Update();

  // TLine *l1 = new TLine(0., 0, grResidualN->GetXaxis()->GetXmax(), 0);
  TLine *l1 = new TLine(0., 0, 10, 0);
  if (!startAtZero) {
    delete l1;
    TLine *l1 = new TLine(0., 0, grResidualN->GetXaxis()->GetXmax()-0.5, 0);
  }
  l1->SetLineStyle(7);
  l1->SetLineColor(kRed);
  l1->Draw("");

  // residuals axis label
  cFit->cd();
  TText *t = new TText(.8,.05,"thickness / cm");
  t->SetNDC(true);
  t->SetTextFont(43);
  t->SetTextSize(15);
  t->Draw();

  TText *t1 = new TText(.06, .14,"residuals / Hz");
  t1->SetNDC(true);
  t1->SetTextFont(43);
  t1->SetTextSize(15);
  t1->SetTextAngle(90);
  t1->Draw();
  cFit->Update();






  theApp.Run();

  return 0;
}
