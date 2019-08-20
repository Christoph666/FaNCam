// program that sums all neutron / gamma counts from images to create an expoenential
// fit. alternatively single pixels can be selected.
// usage: e.g.: ./absorptionPixelSum Fe 1 1 1 0 10



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

double max_vec(vector<double> x) {
  double result = x[0];
  for (int i=1; i<x.size(); i++) {
    if (x[i] > result) {
      result = x[i];
    }
  }
  return result;
}

double min_vec(vector<double> x) {
  double result = x[0];
  for (int i=1; i<x.size(); i++) {
    if (x[i] < result) {
      result = x[i];
    }
  }
  return result;
}

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

double sum_vec(vector<double> *x) {
  double result = 0;
  for (int i=0; i<x->size(); i++) {
    result +=x->at(i);
  }
  return result;
}

// double exponential(double *x, double *par) {
//   return par[0] * exp(- x[0] / (par[1]) * log(2));
// }
//
// double exponential_off(double *x, double *par) {
//   return par[0] * exp(- x[0] / (par[1]) * log(2)) + par[2];
// }

double exponential(double *x, double *par) {
  return par[0] * exp(- x[0] / (par[1]));
}

double exponential_off(double *x, double *par) {
  return par[0] * exp(- x[0] / (par[1])) + par[2];
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
int iPixel=0;
bool saveToTxt = false;
// bool mean_pixels = false;
int doAbsolute = 0;
vector<double> neutronCounts0, gammaRate0;
double systErr = 0.01;

string dir, measurement, calibration_file;
vector<string> measurements;
vector<string> channelMapping = {"A1", "A2", "A3", "A4", "B1", "B2", "B3", "B4", "C1", "C2", "C3", "C4", "D1", "D2", "D3", "D4"};
vector<double> vPx = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
vector<double> vPx_err = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

double measurementTime_calib;
double measurementTime_absorb;

vector<double> neutronCounts, neutronCounts_err, gammaCounts, gammaCounts_err;
vector<double> mean_neutronCounts, mean_neutronCounts_err, mean_gammaCounts, mean_gammaCounts_err;
vector<double> *nCounts, *nCounts_err, *gCounts, *gCounts_err;
double thickness_err;
vector<double> vThickness, vThickness_err = {};
vector<double> nResiduals, gResiduals = {};

double chin, chig, NDFn, NDFg, An, Ag, xn, xg, A0n, A0g, probn, probg, En, Eg;
double An_err, Ag_err, xn_err, xg_err, A0n_err, A0g_err, En_err, Eg_err;


// ./analyseAbsorption measurement fitWithOffset startAtZero doAbsolute pxSelect
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
    cout << "usage: ./analyseAbsorption measurement fitWithOffset startAtZero doAbsolute singlePixel iPixel saveToTxt" << endl;
    // cout << "or: ./analyseAbsorption measurement fitWithOffset startAtZero doAbsolute iPixel" << endl;
    return 0;
  }
  if (argc == 8) {
    measurement = (string) argv[argc-7];
    fitWithOffset = atoi(argv[argc-6]);
    startAtZero = atoi(argv[argc-5]);
    doAbsolute = atoi(argv[argc-4]);
    singlePixel = atoi(argv[argc-3]);
    iPixel = atoi(argv[argc-2]);
    saveToTxt = atoi(argv[argc-1]);
  }
  else {
    cout << "wrong arguments" << endl;
  }
  cout << "iPixel " << iPixel << endl;
  cout << "singlePixel " << singlePixel << endl;
  cout << "doAbsolute " << doAbsolute << endl;
  dir = dir + measurement + "/";
  cout << "fitWithOffset: " << fitWithOffset << endl;
  cout << "startAtZero: " << startAtZero << endl;
  cout << "measurement: " << measurement << endl;
  cout << "saveToTxt: " << saveToTxt << endl;

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
    thickness_err = 0.01;
  }
  // calibration_file = "Calibration_20190617";

  measurementTime_calib = readMeasurementTime(dir+calibration_file+"/run.info");
  measurementTime_absorb = readMeasurementTime(dir+measurements[0]+"/run.info");

  // read calibration file data
  cout << "reading calibration file: " << dir+calibration_file+"/Results/"+calibration_file+"_Results.root" << endl;
  TFile *fCalib = new TFile((dir+calibration_file+"/Results/"+calibration_file+"_Results.root").c_str());
  fCalib->GetObject("N_neutron", nCounts);
  fCalib->GetObject("N_gamma", gCounts);

  if (!singlePixel)
  neutronCounts.push_back(sum_vec(nCounts));
  else
  neutronCounts.push_back(nCounts->at(iPixel));

  neutronCounts_err.push_back( sqrt( neutronCounts[0] + pow(neutronCounts[0]*systErr, 2) ) );

  if (!singlePixel)
  gammaCounts.push_back(sum_vec(gCounts));
  else
  gammaCounts.push_back(gCounts->at(iPixel));

  gammaCounts_err.push_back( sqrt( gammaCounts[0] + pow(gammaCounts[0]*systErr, 2) ) );

  neutronCounts[0] *= (measurementTime_absorb/measurementTime_calib);
  neutronCounts_err[0] *= (measurementTime_absorb/measurementTime_calib);
  gammaCounts[0] *= (measurementTime_absorb/measurementTime_calib);
  gammaCounts_err[0] *= (measurementTime_absorb/measurementTime_calib);
  fCalib->Close();

  // read absorption file data
  for (int i=0; i<measurements.size(); i++) {
    cout << "reading absorption file: " << dir+measurements[i]+"/Images/"+measurements[i]+"_Results.root" << endl;
    TFile *fAbsorb = new TFile((dir+measurements[i]+"/Images/"+measurements[i]+"_Results.root").c_str());
    fAbsorb->GetObject("neutronCounts_abs", nCounts);
    fAbsorb->GetObject("gammaCounts_abs", gCounts);

    if (!singlePixel)
    neutronCounts.push_back(sum_vec(nCounts));
    else
    neutronCounts.push_back(nCounts->at(iPixel));

    neutronCounts_err.push_back(sqrt( neutronCounts.back() + pow(neutronCounts.back()*systErr, 2) ) );
    // cout << neutronCounts.back() << endl;
    // cout << neutronCounts_err.back() << endl;

    if (!singlePixel)
    gammaCounts.push_back(sum_vec(gCounts));
    else
    gammaCounts.push_back(gCounts->at(iPixel));

    gammaCounts_err.push_back(sqrt( gammaCounts.back() + pow(neutronCounts.back()*systErr, 2) ) );
    fAbsorb->Close();
  }


  for (int i=0; i<neutronCounts.size(); i++) {
    neutronCounts[i] /= measurementTime_absorb;
    neutronCounts_err[i] /= measurementTime_absorb;
    gammaCounts[i] /= measurementTime_absorb;
    gammaCounts_err[i] /= measurementTime_absorb;
  }
  double neutronNorm = neutronCounts[0];
  double gammaNorm = gammaCounts[0];

  if (!doAbsolute) {
    for (int i=0; i<neutronCounts.size(); i++) {
      neutronCounts[i] = neutronCounts[i] / neutronNorm;
      neutronCounts_err[i] = neutronCounts_err[i] / neutronNorm;

      gammaCounts[i] = gammaCounts[i] / gammaNorm;
      gammaCounts_err[i] = gammaCounts_err[i] / gammaNorm;
    }
  }
  // output plot and fitting
  // gStyle->SetOptFit(1111);
  gStyle->SetOptFit(0);
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

  TGraphErrors *grExpoN;
  TGraphErrors *grExpoG;

  // average rates and create thickness vectors
  if (!startAtZero)
    iStart = 1;
  for (int i=iStart; i<neutronCounts.size(); i++) {
    vThickness.push_back(i);
    vThickness_err.push_back(thickness_err);
  }

  if (!startAtZero) {
    neutronCounts.erase(neutronCounts.begin());
    neutronCounts_err.erase(neutronCounts_err.begin());
    gammaCounts.erase(gammaCounts.begin());
    gammaCounts_err.erase(gammaCounts_err.begin());
  }


  // do the fitting
  grExpoN = new TGraphErrors(neutronCounts.size(), &vThickness[0], &neutronCounts[0], &vThickness_err[0], &neutronCounts_err[0]);
  grExpoG = new TGraphErrors(gammaCounts.size(), &vThickness[0], &gammaCounts[0], &vThickness_err[0], &gammaCounts_err[0]);
  // grExpoN->SetTitle(measurement.c_str());
  grExpoN->SetTitle((measurement + " Sum").c_str());
  if (measurement == (string)"Cu")
  grExpoN->SetTitle(((string)"Co" + " Sum").c_str());
  if (singlePixel){
    grExpoN->SetTitle((measurement + " Pixel " + channelMapping[iPixel]).c_str());
    if (measurement == (string)"Cu")
    grExpoN->SetTitle(((string)"Co" + " Pixel " + channelMapping[iPixel]).c_str());

  }
  grExpoG->SetMarkerColor(kRed);
  grExpoN->SetMarkerColor(kGreen);



  // cout << "\nfitting neutrons..." << endl;
  fitFunction->SetLineColor(kGreen);
  fitFunction_off->SetLineColor(kGreen);

  // with offset
  if (fitWithOffset) {
    if (!doAbsolute)
      fitFunction_off->SetParameters(1, 3.5, 0.1);
    else
      fitFunction_off->SetParameters(neutronCounts[0], 3.5, 0.1*neutronCounts[0]);
    if (measurement == (string)"H2O")
      fitFunction->SetParameter(1,35);

    fitFunction_off->SetParNames("A^{n}", "x^{n} / cm", "A_{0}^{n}");
    grExpoN->Fit("fit1_off", fitSet, "", -0.2, vThickness.back()+0.2);

    for (int i=0; i<neutronCounts.size(); i++)
      nResiduals.push_back( (neutronCounts[i] - fitFunction_off->Eval(vThickness[i]) ));

    An=fitFunction_off->GetParameter(0);
    An_err=fitFunction_off->GetParError(0);
    xn=fitFunction_off->GetParameter(1);
    xn_err=fitFunction_off->GetParError(1);
    A0n=fitFunction_off->GetParameter(2);
    A0n_err=fitFunction_off->GetParError(2);
    chin=fitFunction_off->GetChisquare();
    NDFn=fitFunction_off->GetNDF();
    probn=fitFunction_off->GetProb();
  }
  // without offset
  else {
    if (!doAbsolute)
      fitFunction->SetParameters(1, 3.5);
    else
      fitFunction->SetParameters(gammaCounts[0], 3.5);

    fitFunction->SetParNames("A^{n}", "x^{n} / cm");
    grExpoN->Fit("fit1", fitSet, "", -0.2, vThickness.back()+0.2);

    for (int i=0; i<neutronCounts.size(); i++)
      nResiduals.push_back( (neutronCounts[i] - fitFunction->Eval(vThickness[i]) ));

    An=fitFunction->GetParameter(0);
    An_err=fitFunction->GetParError(0);
    xn=fitFunction->GetParameter(1);
    xn_err=fitFunction->GetParError(1);
    chin=fitFunction->GetChisquare();
    NDFn=fitFunction->GetNDF();
    probn=fitFunction->GetProb();
  }
  En=log(2)/xn;
  En_err=pow(log(2)/xn, 2) * xn_err;




  // cout << "\nfitting gammas..." << endl;
  fitFunction->SetLineColor(kRed);
  fitFunction_off->SetLineColor(kRed);
  if (measurement == (string)"H2O")
    fitFunction->SetParameter(1, 35);
  // with offset
  if (fitWithOffset) {
    if (!doAbsolute)
      fitFunction_off->SetParameters(1, 3.5, 0.1);
    else
      fitFunction_off->SetParameters(gammaCounts[0], 3.5, 0.1*gammaCounts[0]);
    if (measurement == (string)"H2O")
      fitFunction->SetParameter(1,35);

    fitFunction_off->SetParNames("A^{g}", "x^{g} / cm", "A_{0}^{g}");
    grExpoG->Fit("fit1_off", fitSet, "", -0.2, vThickness.back()+0.2);

    // cout << "size: " << mean_gammaCounts.size() << endl;
    for (int i=0; i<gammaCounts.size(); i++)
      gResiduals.push_back( (gammaCounts[i] - fitFunction_off->Eval(vThickness[i])));

    Ag=fitFunction_off->GetParameter(0);
    Ag_err=fitFunction_off->GetParError(0);
    xg=fitFunction_off->GetParameter(1);
    xg_err=fitFunction_off->GetParError(1);
    A0g=fitFunction_off->GetParameter(2);
    A0g_err=fitFunction_off->GetParError(2);
    chig=fitFunction_off->GetChisquare();
    NDFg=fitFunction_off->GetNDF();
    probg=fitFunction_off->GetProb();
  }
  // without offset
  else {
    if (!doAbsolute)
      fitFunction->SetParameters(1, 3.5);
    else
      fitFunction->SetParameters(gammaCounts[0], 3.5);

    fitFunction->SetParNames("A^{g}", "x^{g} / cm");
    grExpoG->Fit("fit1", fitSet, "", -0.2, vThickness.back()+0.2);

    for (int i=0; i<gammaCounts.size(); i++)
      gResiduals.push_back( (gammaCounts[i] - fitFunction->Eval(vThickness[i])));

    Ag=fitFunction->GetParameter(0);
    Ag_err=fitFunction->GetParError(0);
    xg=fitFunction->GetParameter(1);
    xg_err=fitFunction->GetParError(1);
    chig=fitFunction->GetChisquare();
    NDFg=fitFunction->GetNDF();
    probg=fitFunction->GetProb();
  }
  Eg=log(2)/xg;
  Eg_err=pow(log(2)/xg, 2) * xg_err;



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
  // if (doAbsolute)
  // grExpoN->GetYaxis()->SetRangeUser(min_vec(neutronCounts)*0.8, max_vec(gammaCounts)*1.1);
  if (doAbsolute) {
    if (singlePixel)
    grExpoN->GetYaxis()->SetRangeUser(0, 9);
    else
    grExpoN->GetYaxis()->SetRangeUser(0, 110);
  }
  else
  grExpoN->GetYaxis()->SetRangeUser(0., 1.1);
  // else
  //   grExpoN->GetYaxis()->SetRangeUser(0., 1.1);

  // TText *tn = new TText(.2, 0.05,"Neutrons");
  // tn->SetNDC(true);
  // tn->SetTextColor(kGreen);
  // tn->SetTextFont(43);
  // tn->SetTextSize(15);
  // tn->Draw();
  //
  // TText *tg = new TText(.2, 0.09,"Gammas");
  // tg->SetNDC(true);
  // tg->SetTextColor(kRed);
  // tg->SetTextFont(43);
  // tg->SetTextSize(15);
  // tg->Draw();

  TPaveText *rn=new TPaveText(0.65,0.65,0.9,0.9,"brNDC");
  TPaveText *rg=new TPaveText(0.65,0.4,0.9,0.65,"brNDC");
  rn->SetBorderSize(1);
  rg->SetBorderSize(1);

  rn->AddText("Neutrons");
  ((TText*)rn->GetListOfLines()->Last())->SetTextColor(kGreen);
  rn->AddText(Form("#chi^{2} / NDF = %3.2lf / %1.0lf ",chin,NDFn));
  // rn->AddText(Form("Prob = %3.4lf  ",probn));

  if (doAbsolute)
    rn->AddText(Form("A = %3.2lf #pm %3.2lf Hz",An,An_err));
  else
    rn->AddText(Form("A = %3.2lf #pm %3.2lf ",An,An_err));

  rn->AddText(Form("x = %3.2lf #pm %3.2lf cm",xn,xn_err));

  if (fitWithOffset) {
    if (doAbsolute)
      rn->AddText(Form("A_{0} = %3.2lf #pm %3.2lf Hz",A0n,A0n_err));
    else
      rn->AddText(Form("A_{0} = %3.2lf #pm %3.2lf ",A0n,A0n_err));
  }
  rn->AddText(Form("#Sigma = %3.4lf #pm %3.4lf cm^{-1}",En,En_err));
  rn->Draw();

  rg->AddText("Gammas");
  ((TText*)rg->GetListOfLines()->Last())->SetTextColor(kRed);
  rg->AddText(Form("#chi^{2} / NDF = %3.2lf / %1.0lf ",chig,NDFg));
  // rg->AddText(Form("Prob = %3.4lf  ",probg));

  if (doAbsolute)
  rg->AddText(Form("A = %3.2lf #pm %3.2lf Hz",Ag,Ag_err));
  else
  rg->AddText(Form("A = %3.2lf #pm %3.2lf ",Ag,Ag_err));

  rg->AddText(Form("x = %3.2lf #pm %3.2lf cm",xg,xg_err));
  if (fitWithOffset) {
    if (doAbsolute)
    rg->AddText(Form("A_{0} = %3.2lf #pm %3.2lf Hz",A0g,A0g_err));
    else
    rg->AddText(Form("A_{0} = %3.2lf #pm %3.2lf ",A0g,A0g_err));
  }

  rg->AddText(Form("#Sigma = %3.4lf #pm %3.4lf cm^{-1}",Eg,Eg_err));
  rg->Draw();

  cFit->Modified();
  cFit->Update();




  // residual plot
  pad2->cd();
  gPad->SetFillStyle(0);
  TGraphErrors *grResidualN = new TGraphErrors(nResiduals.size(), &vThickness[0], &nResiduals[0], &vThickness_err[0], &neutronCounts_err[0]);
  TGraphErrors *grResidualG = new TGraphErrors(gResiduals.size(), &vThickness[0], &gResiduals[0], &vThickness_err[0], &gammaCounts_err[0]);

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
  if (measurement==(string)"H2O")
    grResidualN->GetXaxis()->SetNdivisions(10, 5, 0, kTRUE);


  grResidualN->Draw("AP");
  grResidualG->Draw("P same");
  // grResidualG->Draw("AP");

  grResidualN->GetXaxis()->SetRangeUser(-0.5, 10.5);
  grResidualN->GetYaxis()->SetRangeUser(-0.06*neutronCounts[0], 0.06*neutronCounts[0]);
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

  mkdir((dir + "Fit_Results").c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  if (saveToTxt) {
    if (singlePixel) {
      cFit->Print((dir + "/Fit_Results/Fit_Result_px" + to_string((int)iPixel) + ".png").c_str());
      ofstream ofile;
      if (fitWithOffset) {
        if (iPixel==0)
        remove((dir + "Fit_Results/Results_off.txt").c_str());
        ofile.open((dir + "Fit_Results/Results_off.txt").c_str(), ios::app);
        cout << "saving results to file: " << dir + "Fit_Results/Results_off.txt" << endl;
        if (iPixel==0) {
          ofile << "iPixel" << "\t" << "An" << "\t" << "An_err" << "\t" << "xn" << "\t" << "xn_err" << "\t" << "A0n" << "\t" << "A0n_err" << "\t" << "En" << "\t" << "En_err";
          ofile << "\t" << "Ag" << "\t" << "Ag_err" << "\t" << "xg" << "\t" << "xg_err" << "\t" << "A0g" << "\t" << "A0g_err" << "\t" << "Eg" << "\t" << "Eg_err" << endl;
        }
      }
      else {
        if (iPixel==0)
        remove((dir + "Fit_Results/Results_no_off.txt").c_str());
        ofile.open((dir + "Fit_Results/Results_no_off.txt").c_str(), ios::app);
        cout << "saving results to file: " << dir + "Fit_Results/Results_no_off.txt" << endl;
        if (iPixel==0) {
          ofile << "iPixel" << "\t" << "An" << "\t" << "An_err" << "\t" << "xn" << "\t" << "xn_err" << "\t" << "En" << "\t" << "En_err";
          ofile << "\t" << "Ag" << "\t" << "Ag_err" << "\t" << "xg" << "\t" << "xg_err" << "\t" << "\t" << "Eg" << "\t" << "Eg_err" << endl;
        }
      }
      if (fitWithOffset) {
        ofile << iPixel << "\t" << An << "\t" << An_err << "\t" << xn << "\t" << xn_err << "\t" << A0n << "\t" << A0n_err << "\t" << En << "\t" << En_err;
        ofile << "\t" << Ag << "\t" << Ag_err << "\t" << xg << "\t" << xg_err << "\t" << A0g << "\t" << A0g_err << "\t" << Eg << "\t" << Eg_err << endl;
      }
      else {
        ofile << iPixel << "\t" << An << "\t" << An_err << "\t" << xn << "\t" << xn_err << "\t" << En << "\t" << En_err;
        ofile << "\t" << Ag << "\t" << Ag_err << "\t" << xg << "\t" << xg_err << "\t" << "\t" << Eg << "\t" << Eg_err << endl;
      }
      ofile.close();
    }
    else {
      ofstream ofile;
      if (fitWithOffset) {
        remove((dir + "Fit_Results/Results_sum_off.txt").c_str());
        ofile.open((dir + "Fit_Results/Results_sum_off.txt").c_str(), ios::app);
        cout << "saving results to file: " << dir + "Fit_Results/Results_sum_off.txt" << endl;
        ofile << "An" << "\t" << "An_err" << "\t" << "xn" << "\t" << "xn_err" << "\t" << "A0n" << "\t" << "A0n_err" << "\t" << "En" << "\t" << "En_err";
        ofile << "\t" << "Ag" << "\t" << "Ag_err" << "\t" << "xg" << "\t" << "xg_err" << "\t" << "A0g" << "\t" << "A0g_err" << "\t" << "Eg" << "\t" << "Eg_err" << endl;
      }
      else {
        remove((dir + "Fit_Results/Results_sum_no_off.txt").c_str());
        ofile.open((dir + "Fit_Results/Results_sum_no_off.txt").c_str(), ios::app);
        cout << "saving results to file: " << dir + "Fit_Results/Results_sum_no_off.txt" << endl;
        ofile << "An" << "\t" << "An_err" << "\t" << "xn" << "\t" << "xn_err" << "\t" << "En" << "\t" << "En_err";
        ofile << "\t" << "Ag" << "\t" << "Ag_err" << "\t" << "xg" << "\t" << "xg_err" << "\t" << "\t" << "Eg" << "\t" << "Eg_err" << endl;
      }
      if (fitWithOffset) {
        ofile << An << "\t" << An_err << "\t" << xn << "\t" << xn_err << "\t" << A0n << "\t" << A0n_err << "\t" << En << "\t" << En_err;
        ofile << "\t" << Ag << "\t" << Ag_err << "\t" << xg << "\t" << xg_err << "\t" << A0g << "\t" << A0g_err << "\t" << Eg << "\t" << Eg_err << endl;
      }
      else {
        ofile << An << "\t" << An_err << "\t" << xn << "\t" << xn_err << "\t" << En << "\t" << En_err;
        ofile << "\t" << Ag << "\t" << Ag_err << "\t" << xg << "\t" << xg_err << "\t" << "\t" << Eg << "\t" << Eg_err << endl;
      }
      ofile.close();
    }

    if (doAbsolute)
    cFit->Print((dir + "Fit_Results/Fit_Result_sum_abs.png").c_str());
    else
    cFit->Print((dir + "Fit_Results/Fit_Result_sum_rel.png").c_str());
  }



  if (!singlePixel)
  theApp.Run();

  return 0;
}
