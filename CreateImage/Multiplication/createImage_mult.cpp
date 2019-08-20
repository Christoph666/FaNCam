// program to compute image out of calibration data with
// the calibration rates and the bands. needs psd ph data
// of absorption image. Image is computed by multtraction.
// todo: instead of fraction implement different measurements times
// i.e. 10min, 20min,.. -> more anschaulich

// todo: check why plot doesnt work for multiple fractions!

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

// #include "filereader.h"
#include "utilities.h"
// #include "stylesheet.h"

using namespace std;

bool save_output = 1;

string filename_calib, filename_absorb, dir;

vector<string> channelMapping = {"A1", "A2", "A3", "A4", "B1", "B2", "B3", "B4", "C1", "C2", "C3", "C4", "D1", "D2", "D3", "D4"};
int labelColor = kWhite;

int iPixel = 0; // pixel counter
int NEvents; // events in one file
int iEvent = 0; // event counter
double measurementTime_absorb, measurementTime_calib, measurementTime_frac;
double fraction;
// vector<double> vFraction = {1., 0.5, 0.2, 0.1, 0.01}; // fraction of measurementTime
// vector<double> vFraction = {0.01, 0.02, 0.1, 0.5, 0.7, 1.}; // fraction of measurementTime
vector<double> vFraction = {1}; // fraction of measurementTime
int NEvents_cut;
int NEvents_total;

vector<double> *vPeak;
vector<double> *vPSD;
vector<int> *vValid;
vector<double> *vSlices_FOM_e;
vector<double> *vSlices_FOM;
vector<double> *vSlices_n;
vector<double> *vSlices_g;
vector<double> *mu_n;
vector<double> *sig_n;
vector<double> *mu_g;
vector<double> *sig_g;

vector<int> *neutronCounts_calib;
vector<int> *gammaCounts_calib;
vector<double> *neutronRates_calib;
vector<double> *gammaRates_calib;
vector<double> neutronRates_calib_err;
vector<double> gammaRates_calib_err;

vector<int> neutronCounts_abs;
vector<int> gammaCounts_abs;
vector<double> neutronRates_abs;
vector<double> gammaRates_abs;
vector<double> neutronRates_abs_err;
vector<double> gammaRates_abs_err;

vector<double> neutronRates_mult;
vector<double> gammaRates_mult;
vector<double> neutronRates_mult_err;
vector<double> gammaRates_mult_err;

vector<double> neutronRates_mult_slice;
vector<double> neutronRates_mult_err_slice;

vector<double> neutronCalibFactors;
vector<double> neutronCalibFactors_err;
vector<double> gammaCalibFactors;
vector<double> gammaCalibFactors_err;


int main () {

  dir = "/home/christoph/Documents/MasterThesis/Data/CoMPASS/";

  // calibration and absorption filenames
  filename_calib = "onlyAmBe_20190305";
  measurementTime_calib = 2.*3600. + 31.*60. + 38.;

  // filename_absorb = "SteelStick_side1_20190306";
  // measurementTime_absorb = 3.*3600. + 50.*60. + 10.;

  filename_absorb = "only_AmBe_20190314";
  measurementTime_absorb = 3.*3600. + 19.*60. + 42.;

  mkdir((dir+filename_absorb+"/Images").c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);



  TROOT _MakePad("MakePad", "Display data"); // initialise ROOT for access inside compiled c++ program
  TApplication theApp("App", NULL, NULL);

  // read calibration image rates
  TFile *file_calib_image = new TFile((dir+filename_calib+"/Plots/FOM/"+filename_calib+"_stats.root").c_str());
  // if ( file_calib_image->IsOpen() ) printf("Calibration image file opened successfully\n");
  file_calib_image->GetObject("vRate_neutron", neutronRates_calib);
  file_calib_image->GetObject("vRate_gamma", gammaRates_calib);
  file_calib_image->GetObject("vAccN", neutronCounts_calib);
  file_calib_image->GetObject("vAccG", gammaCounts_calib);

  for (int i=0; i<neutronRates_calib->size(); i++) {
    neutronRates_calib_err.push_back(sqrt(neutronCounts_calib->at(i)) / measurementTime_calib);
    gammaRates_calib_err.push_back(sqrt(gammaCounts_calib->at(i)) / measurementTime_calib);
    // cout << "calib error : " << neutronRates_calib_err[i] << " rel: " << neutronRates_calib_err[i] / neutronRates_calib->at(i) << endl;
    // cout << "counts : " << neutronCounts_calib->at(i) << " rel error: " << 1. / sqrt(neutronCounts_calib->at(i)) << endl;
  }

  // neutrons
  TH2F *h2dNCalib = new TH2F("Calibration Image","Calibration Image", 4, 0, 4, 4, 0, 4);
  h2dNCalib->GetZaxis()->SetTitle("Acc. neutron rate / Hz");
  h2dNCalib->GetZaxis()->SetTitleOffset(1.4);
  TH2F *h2dNAbsorb = new TH2F("Absorption Image","Absorption Image", 4, 0, 4, 4, 0, 4);
  h2dNAbsorb->GetZaxis()->SetTitle("Acc. neutron rate / Hz");
  h2dNAbsorb->GetZaxis()->SetTitleOffset(1.4);
  TH2F *h2dNmult = new TH2F("mult Image","mult Image", 4, 0, 4, 4, 0, 4);
  h2dNmult->GetZaxis()->SetTitle("#Delta neutron rate / Hz");
  h2dNmult->GetZaxis()->SetTitleOffset(1.4);

  // gammas
  TH2F *h2dGCalib = new TH2F("Gamma Calibration Image","Gamma Calibration Image", 4, 0, 4, 4, 0, 4);
  h2dGCalib->GetZaxis()->SetTitle("Acc. gamma rate / Hz");
  h2dGCalib->GetZaxis()->SetTitleOffset(1.4);
  TH2F *h2dGAbsorb = new TH2F("Gamma Absorption Image","Gamma Absorption Image", 4, 0, 4, 4, 0, 4);
  h2dGAbsorb->GetZaxis()->SetTitle("Acc. gamma rate / Hz");
  h2dGAbsorb->GetZaxis()->SetTitleOffset(1.4);
  TH2F *h2dGMult = new TH2F("Gamma Mult Image","Gamma Mult Image", 4, 0, 4, 4, 0, 4);
  h2dGMult->GetZaxis()->SetTitle("#Delta gamma rate / Hz");
  h2dGMult->GetZaxis()->SetTitleOffset(1.4);

  // loop over different fractions of measurement time
  for (int iFraction=0; iFraction<vFraction.size(); iFraction++) {

    // reset variables and compute fraction of measurement time
    neutronCounts_abs = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    gammaCounts_abs = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    NEvents_total=0;
    fraction = vFraction[iFraction];
    measurementTime_frac = measurementTime_absorb * fraction;
    cout << "Fraction : " << fraction << endl;
    cout << "measurementTime : " << measurementTime_frac << " s" << endl;

    // loop over pixels
    for (int iDEV=0; iDEV<2; iDEV++) {
      for (int iCHN=0; iCHN<8; iCHN++) {

        cout << '\r' << "Dev: " << iDEV << " Chn: " << iCHN << flush;

        // read bands from calibration image
        TFile *file_calib = new TFile((dir+filename_calib+"/Plots/FOM/"+filename_calib+"_"+to_string((int)iDEV)+"_"+to_string((int)iCHN)+"/"+filename_calib+"_"+to_string((int)iDEV)+"_"+to_string((int)iCHN)+"_FOM_results.root").c_str());
        // if ( file_calib->IsOpen() ) printf("Calibration file opened successfully\n");
        file_calib->GetObject("vSlices_n", vSlices_n);
        file_calib->GetObject("mu_n", mu_n);
        file_calib->GetObject("sig_n", sig_n);
        file_calib->GetObject("vSlices_g", vSlices_g);
        file_calib->GetObject("mu_g", mu_g);
        file_calib->GetObject("sig_g", sig_g);

        // read psd, ph data from absorption measurement
        TFile *file_absorb = new TFile((dir+filename_absorb+"/Plots/"+filename_absorb+"_"+to_string((int)iDEV)+"_"+to_string((int)iCHN)+"_25_270_PSD.root").c_str());
        // if ( file_absorb->IsOpen() ) printf("Absorb file opened successfully\n");
        file_absorb->GetObject("vPeak", vPeak);
        file_absorb->GetObject("vPSD", vPSD);
        file_absorb->GetObject("vValid", vValid);

        NEvents = vPSD->size();
        NEvents_cut = (int) NEvents * fraction;
        iEvent = 0;
        // cout <<  "nevents cut : " << NEvents_cut << endl;

        // counts gammas / neutrons
        for (int i=0; i<NEvents; i++) {

          // check for neutron
          if (vPeak->at(i)>vSlices_n->front() && vPeak->at(i)<vSlices_n->back()) {
            if (vPSD->at(i) < interpolate1d(vPeak->at(i), vSlices_n, mu_n) + 3*interpolate1d(vPeak->at(i), vSlices_n, sig_n)
            && vPSD->at(i) > interpolate1d(vPeak->at(i), vSlices_n, mu_n) - 3*interpolate1d(vPeak->at(i), vSlices_n, sig_n)
            && vPSD->at(i) > interpolate1d(vPeak->at(i), vSlices_g, mu_g) + 3*interpolate1d(vPeak->at(i), vSlices_g, sig_g))
            {
              neutronCounts_abs[iPixel]++;
            }
          }
          // check for gamma
          if (vPeak->at(i)>vSlices_g->front() && vPeak->at(i)<vSlices_g->back()) {
            if (vPSD->at(i) < interpolate1d(vPeak->at(i), vSlices_g, mu_g) + 3*interpolate1d(vPeak->at(i), vSlices_g, sig_g)
            && vPSD->at(i) > interpolate1d(vPeak->at(i), vSlices_g, mu_g) - 3*interpolate1d(vPeak->at(i), vSlices_g, sig_g)
            && vPSD->at(i) < interpolate1d(vPeak->at(i), vSlices_n, mu_n) - 3*interpolate1d(vPeak->at(i), vSlices_n, sig_n))
            {
              gammaCounts_abs[iPixel]++;
            }
          }
          iEvent++;
          NEvents_total++;
          if (NEvents_cut > 0 && iEvent >= NEvents_cut) {
            // cout << "breaking " << iEvent << endl;
            break;
          }
        }

        // cout << iDEV << " " <<  iCHN << "Neutrons : " << neutronCounts_abs[iPixel] << endl;

        iPixel++;

        // clear variables
        file_calib->Close();
        file_absorb->Close();
        delete file_calib, file_absorb;
        vSlices_n->clear();
        mu_n->clear();
        sig_n->clear();
        vSlices_g->clear();
        mu_g->clear();
        sig_g->clear();
        vPeak->clear();
        vPSD->clear();
        vValid->clear();
      }
    }

    cout << endl;
    cout << "Toatal number of events : " << NEvents_total << endl;

    // compute normalization factors
    float maxCalib = neutronRates_calib->at(0);
    for (int i=0; i<neutronRates_calib->size(); i++) {
      if (neutronRates_calib->at(i) > maxCalib) {
        maxCalib = neutronRates_calib->at(i);
      }
    }
    for (int i=0; i<neutronRates_calib->size(); i++) {
      neutronCalibFactors.push_back(neutronRates_calib->at(i) / maxCalib);
      neutronCalibFactors_err.push_back(neutronRates_calib_err[i] / maxCalib);
    }
    float maxCalibGamma = gammaRates_calib->at(0);
    for (int i=0; i<gammaRates_calib->size(); i++) {
      if (gammaRates_calib->at(i) > maxCalib) {
        maxCalib = gammaRates_calib->at(i);
      }
    }
    for (int i=0; i<gammaRates_calib->size(); i++) {
      gammaCalibFactors.push_back(gammaRates_calib->at(i) / maxCalib);
      gammaCalibFactors_err.push_back(gammaRates_calib_err[i] / maxCalib);
    }

    // compute rates and multtract calibration image
    for (int i=0; i<neutronCounts_abs.size(); i++) {
      neutronRates_abs.push_back((float)neutronCounts_abs[i] / measurementTime_frac);
      gammaRates_abs.push_back((float)gammaCounts_abs[i] / measurementTime_frac);
      neutronRates_abs_err.push_back(sqrt((float)neutronCounts_abs[i]) / measurementTime_frac);
      gammaRates_abs_err.push_back(sqrt((float)gammaCounts_abs[i]) / measurementTime_frac);
      // cout << "absorption error : " << neutronRates_abs_err[i] << " rel: " << neutronRates_abs_err[i] / neutronRates_abs[i] << endl;
      // cout << "counts : " << neutronCounts_abs[i] << " rel error: " << 1. / sqrt(neutronCounts_abs[i]) << endl << endl;

      // cout << (float)neutronCounts_abs[i] << endl;
      neutronRates_mult.push_back((float)neutronCounts_abs[i] / measurementTime_frac / neutronCalibFactors[i]);
      neutronRates_mult_err.push_back( neutronRates_mult[i] * sqrt(pow(neutronRates_abs_err[i]/neutronRates_abs[i],2) + pow(neutronCalibFactors_err[i]/neutronCalibFactors[i],2) ));
      gammaRates_mult.push_back((float)gammaCounts_abs[i] / measurementTime_frac / gammaCalibFactors[i]);
      gammaRates_mult_err.push_back( gammaRates_mult[i] * sqrt(pow(gammaRates_abs_err[i]/gammaRates_abs[i],2) + pow(gammaCalibFactors_err[i]/gammaCalibFactors[i],2) ));
      // cout << "multtraction error : " << neutronRates_mult_err[i] << " rel: " << neutronRates_mult_err[i] / neutronRates_mult[i] << endl << endl;
      // cout << "counts : " << neutronCounts_mult[i] << " rel error: " << 1. / sqrt(neutronCounts_mult[i]) << endl << endl;
    }

    // TCanvas *cImages = new TCanvas("cImages", "cImages", 1500, 600);

    // cImages->Divide(3,1);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kViridis);

    ////////////////////////////////////////////////////////////////////////////
    // neutron image
    ////////////////////////////////////////////////////////////////////////////
    // calibration image
    // cImages->cd(1);
    TCanvas *cCalibration = new TCanvas("cCalibration", "cCalibration", 600, 600);
    cCalibration->SetRightMargin(0.15);
    iPixel=0;
    for (int iBinx=1; iBinx<5; iBinx++) {
      for (int iBiny=4; iBiny>0; iBiny--) {
        h2dNCalib->SetBinContent(iBinx, iBiny, neutronRates_calib->at(iPixel));
        iPixel++;
      }
    }
    h2dNCalib->Draw("colz");
    iPixel=0;
    for (int iBinx=1; iBinx<5; iBinx++) {
      for (int iBiny=4; iBiny>0; iBiny--) {
        TText *t1 = new TText(iBinx-0.5, iBiny-0.5, (channelMapping[iPixel]).c_str());
        t1->SetTextColor(labelColor);
        t1->SetTextFont(43);
        t1->SetTextSize(20);
        t1->Draw("same");
        iPixel++;
      }
    }
    cCalibration->Modified();
    cCalibration->Update();

    // absorption image
    // cImages->cd(2);
    TCanvas *cAbsorbtion = new TCanvas("cAbsorbtion", "cAbsorbtion", 600, 600);
    cAbsorbtion->SetRightMargin(0.15);
    iPixel=0;
    for (int iBinx=1; iBinx<5; iBinx++) {
      for (int iBiny=4; iBiny>0; iBiny--) {
        h2dNAbsorb->SetBinContent(iBinx, iBiny, neutronRates_abs[iPixel]);
        iPixel++;
      }
    }
    h2dNAbsorb->Draw("colz");
    iPixel=0;
    for (int iBinx=1; iBinx<5; iBinx++) {
      for (int iBiny=4; iBiny>0; iBiny--) {
        TText *t1 = new TText(iBinx-0.5, iBiny-0.5, (channelMapping[iPixel]).c_str());
        t1->SetTextColor(labelColor);
        t1->SetTextFont(43);
        t1->SetTextSize(20);
        t1->Draw("same");
        iPixel++;
      }
    }
    cAbsorbtion->Modified();
    cAbsorbtion->Update();

    // multtracted image
    // cImages->cd(3);
    TCanvas *cmult = new TCanvas("cmult", "cmult", 600, 600);
    cmult->SetRightMargin(0.15);
    iPixel=0;
    for (int iBinx=1; iBinx<5; iBinx++) {
      for (int iBiny=4; iBiny>0; iBiny--) {
        h2dNmult->SetBinContent(iBinx, iBiny, neutronRates_mult[iPixel]);
        iPixel++;
      }
    }
    h2dNmult->Draw("colz");
    iPixel=0;
    for (int iBinx=1; iBinx<5; iBinx++) {
      for (int iBiny=4; iBiny>0; iBiny--) {
        TText *t1 = new TText(iBinx-0.5, iBiny-0.5, (channelMapping[iPixel]).c_str());
        t1->SetTextColor(labelColor);
        t1->SetTextFont(43);
        t1->SetTextSize(20);
        t1->Draw("same");
        iPixel++;
      }
    }
    cmult->Modified();
    cmult->Update();

    // plot the projection of the rates with errors
    TCanvas *c1d = new TCanvas("c1d", "c1d", 600, 600);
    c1d->SetLeftMargin(0.15);
    TMultiGraph *mg = new TMultiGraph();
    vector<int> colors = {kRed, kBlue, kGreen, kTeal};
    vector<string> names = {"A", "B", "C", "D"};
    for (int iRow=0; iRow<4; iRow++) {
      TGraphErrors *gr1d = new TGraphErrors();
      for (int i=0; i<4; i++) {
        gr1d->SetPoint(i, i+1, neutronRates_mult[i+iRow*4]);
        gr1d->SetPointError(i, 0, neutronRates_mult_err[i+iRow*4]);
        // cout << "error " << i << " " << neutronRates_mult_err[i] << endl;
      }
      gr1d->GetXaxis()->SetRangeUser(0.5, 4.5);
      gr1d->SetMarkerStyle(20);
      gr1d->SetMarkerColor(colors[iRow]);
      mg->Add(gr1d);
      // gr1d->Clear();

    }
    mg->Draw("AP");
    mg->GetXaxis()->SetTitle("i pixel");
    mg->GetYaxis()->SetTitle("#Delta neutron rate / Hz");

    float xpos=3.;
    float ypos=7.8;
    float xsep=0.25;
    TText *tA = new TText(xpos+xsep*0, ypos, (names[0]).c_str());
    tA->SetTextColor(colors[0]);
    tA->SetTextFont(43);
    tA->SetTextSize(20);
    tA->Draw("same");
    TText *tB = new TText(xpos+xsep*1, ypos, (names[1]).c_str());
    tB->SetTextColor(colors[1]);
    tB->SetTextFont(43);
    tB->SetTextSize(20);
    tB->Draw("same");
    TText *tC = new TText(xpos+xsep*2, ypos, (names[2]).c_str());
    tC->SetTextColor(colors[2]);
    tC->SetTextFont(43);
    tC->SetTextSize(20);
    tC->Draw("same");
    TText *tD = new TText(xpos+xsep*3, ypos, (names[3]).c_str());
    tD->SetTextColor(colors[3]);
    tD->SetTextFont(43);
    tD->SetTextSize(20);
    tD->Draw("same");

    c1d->Modified();
    c1d->Update();
    mg->GetXaxis()->SetNdivisions(4);

    ////////////////////////////////////////////////////////////////////////////
    // gamma image
    ////////////////////////////////////////////////////////////////////////////
    // calibration image
    // cImages->cd(1);
    TCanvas *cCalibrationGamma = new TCanvas("cCalibrationGamma", "cCalibrationGamma", 600, 600);
    cCalibrationGamma->SetRightMargin(0.15);
    iPixel=0;
    for (int iBinx=1; iBinx<5; iBinx++) {
      for (int iBiny=4; iBiny>0; iBiny--) {
        h2dGCalib->SetBinContent(iBinx, iBiny, gammaRates_calib->at(iPixel));
        iPixel++;
      }
    }
    h2dGCalib->Draw("colz");
    iPixel=0;
    for (int iBinx=1; iBinx<5; iBinx++) {
      for (int iBiny=4; iBiny>0; iBiny--) {
        TText *t1 = new TText(iBinx-0.5, iBiny-0.5, (channelMapping[iPixel]).c_str());
        t1->SetTextColor(labelColor);
        t1->SetTextFont(43);
        t1->SetTextSize(20);
        t1->Draw("same");
        iPixel++;
      }
    }
    cCalibration->Modified();
    cCalibration->Update();

    // absorption image
    // cImages->cd(2);
    TCanvas *cAbsorbtionGamma = new TCanvas("cAbsorbtionGamma", "cAbsorbtionGamma", 600, 600);
    cAbsorbtionGamma->SetRightMargin(0.15);
    iPixel=0;
    for (int iBinx=1; iBinx<5; iBinx++) {
      for (int iBiny=4; iBiny>0; iBiny--) {
        h2dGAbsorb->SetBinContent(iBinx, iBiny, gammaRates_abs[iPixel]);
        iPixel++;
      }
    }
    h2dGAbsorb->Draw("colz");
    iPixel=0;
    for (int iBinx=1; iBinx<5; iBinx++) {
      for (int iBiny=4; iBiny>0; iBiny--) {
        TText *t1 = new TText(iBinx-0.5, iBiny-0.5, (channelMapping[iPixel]).c_str());
        t1->SetTextColor(labelColor);
        t1->SetTextFont(43);
        t1->SetTextSize(20);
        t1->Draw("same");
        iPixel++;
      }
    }
    cAbsorbtionGamma->Modified();
    cAbsorbtionGamma->Update();

    // Mult image
    // cImages->cd(3);
    TCanvas *cMultGamma = new TCanvas("cMultGamma", "cMultGamma", 600, 600);
    cMultGamma->SetRightMargin(0.15);
    iPixel=0;
    for (int iBinx=1; iBinx<5; iBinx++) {
      for (int iBiny=4; iBiny>0; iBiny--) {
        h2dGMult->SetBinContent(iBinx, iBiny, gammaRates_mult[iPixel]);
        iPixel++;
      }
    }
    h2dGMult->Draw("colz");
    iPixel=0;
    for (int iBinx=1; iBinx<5; iBinx++) {
      for (int iBiny=4; iBiny>0; iBiny--) {
        TText *t1 = new TText(iBinx-0.5, iBiny-0.5, (channelMapping[iPixel]).c_str());
        t1->SetTextColor(labelColor);
        t1->SetTextFont(43);
        t1->SetTextSize(20);
        t1->Draw("same");
        iPixel++;
      }
    }
    cMultGamma->Modified();
    cMultGamma->Update();

    // plot the projection of the rates with errors
    TCanvas *c1dGamma = new TCanvas("c1dGamma", "c1dGamma", 600, 600);
    c1dGamma->SetLeftMargin(0.15);
    TMultiGraph *mgGamma = new TMultiGraph();
    for (int iRow=0; iRow<4; iRow++) {
      TGraphErrors *gr1dGamma = new TGraphErrors();
      for (int i=0; i<4; i++) {
        gr1dGamma->SetPoint(i, i+1, gammaRates_mult[i+iRow*4]);
        gr1dGamma->SetPointError(i, 0, gammaRates_mult_err[i+iRow*4]);
        // cout << "error " << i << " " << neutronRates_mult_err[i] << endl;
      }
      gr1dGamma->GetXaxis()->SetRangeUser(0.5, 4.5);
      gr1dGamma->SetMarkerStyle(20);
      gr1dGamma->SetMarkerColor(colors[iRow]);
      mgGamma->Add(gr1dGamma);
      // gr1d->Clear();

    }
    mgGamma->Draw("AP");
    mgGamma->GetXaxis()->SetTitle("i pixel");
    mgGamma->GetYaxis()->SetTitle("#Delta gamma rate / Hz");

    xpos=3;
    ypos=11.;
    xsep=0.25;
    TText *tAg = new TText(xpos+xsep*0, ypos, (names[0]).c_str());
    tAg->SetTextColor(colors[0]);
    tAg->SetTextFont(43);
    tAg->SetTextSize(20);
    tAg->Draw("same");
    TText *tBg = new TText(xpos+xsep*1, ypos, (names[1]).c_str());
    tBg->SetTextColor(colors[1]);
    tBg->SetTextFont(43);
    tBg->SetTextSize(20);
    tBg->Draw("same");
    TText *tCg = new TText(xpos+xsep*2, ypos, (names[2]).c_str());
    tCg->SetTextColor(colors[2]);
    tCg->SetTextFont(43);
    tCg->SetTextSize(20);
    tCg->Draw("same");
    TText *tDg = new TText(xpos+xsep*3, ypos, (names[3]).c_str());
    tDg->SetTextColor(colors[3]);
    tDg->SetTextFont(43);
    tDg->SetTextSize(20);
    tDg->Draw("same");

    c1dGamma->Modified();
    c1dGamma->Update();
    mg->GetXaxis()->SetNdivisions(4);


    if (save_output) {
      // cImages->Print(("Images/"+filename_absorb+"_"+to_string((int)(measurementTime_frac/60.))+"min.pdf").c_str());
      // cImages->Print(("Images/"+filename_absorb+"_"+to_string((int)(measurementTime_frac/60.))+"min.png").c_str());

      cCalibration->Print((dir+filename_absorb+"/Images/"+filename_absorb+"_mult_Calibration_"+to_string((int)(measurementTime_frac/60.))+"min.pdf").c_str());
      cCalibration->Print((dir+filename_absorb+"/Images/"+filename_absorb+"_mult_Calibration_"+to_string((int)(measurementTime_frac/60.))+"min.png").c_str());

      cAbsorbtion->Print((dir+filename_absorb+"/Images/"+filename_absorb+"_mult_Absorption_"+to_string((int)(measurementTime_frac/60.))+"min.pdf").c_str());
      cAbsorbtion->Print((dir+filename_absorb+"/Images/"+filename_absorb+"_mult_Absorption_"+to_string((int)(measurementTime_frac/60.))+"min.png").c_str());

      cmult->Print((dir+filename_absorb+"/Images/"+filename_absorb+"_mult_"+to_string((int)(measurementTime_frac/60.))+"min.pdf").c_str());
      cmult->Print((dir+filename_absorb+"/Images/"+filename_absorb+"_mult_"+to_string((int)(measurementTime_frac/60.))+"min.png").c_str());

      c1d->Print((dir+filename_absorb+"/Images/"+filename_absorb+"_mult_projection_"+to_string((int)(measurementTime_frac/60.))+"min.pdf").c_str());
      c1d->Print((dir+filename_absorb+"/Images/"+filename_absorb+"_mult_projection_"+to_string((int)(measurementTime_frac/60.))+"min.png").c_str());

      cCalibrationGamma->Print((dir+filename_absorb+"/Images/"+filename_absorb+"_Gamma_Calibration_"+to_string((int)(measurementTime_frac/60.))+"min.pdf").c_str());
      cCalibrationGamma->Print((dir+filename_absorb+"/Images/"+filename_absorb+"_Gamma_Calibration_"+to_string((int)(measurementTime_frac/60.))+"min.png").c_str());

      cAbsorbtionGamma->Print((dir+filename_absorb+"/Images/"+filename_absorb+"_Gamma_Absorption_"+to_string((int)(measurementTime_frac/60.))+"min.pdf").c_str());
      cAbsorbtionGamma->Print((dir+filename_absorb+"/Images/"+filename_absorb+"_Gamma_Absorption_"+to_string((int)(measurementTime_frac/60.))+"min.png").c_str());

      cMultGamma->Print((dir+filename_absorb+"/Images/"+filename_absorb+"_Gamma_Mult_"+to_string((int)(measurementTime_frac/60.))+"min.pdf").c_str());
      cMultGamma->Print((dir+filename_absorb+"/Images/"+filename_absorb+"_Gamma_Mult_"+to_string((int)(measurementTime_frac/60.))+"min.png").c_str());

      c1dGamma->Print((dir+filename_absorb+"/Images/"+filename_absorb+"_Gamma_projection_"+to_string((int)(measurementTime_frac/60.))+"min.pdf").c_str());
      c1dGamma->Print((dir+filename_absorb+"/Images/"+filename_absorb+"_Gamma_projection_"+to_string((int)(measurementTime_frac/60.))+"min.png").c_str());

      TFile *fOut = new TFile((dir+filename_absorb+"/Images/"+filename_absorb+"_mult_"+to_string((int)(measurementTime_frac/60.))+"min.root").c_str(), "RECREATE");
      h2dNCalib->Write();
      h2dNAbsorb->Write();
      h2dNCalib->Write();
      h2dGCalib->Write();
      h2dGAbsorb->Write();
      h2dGMult->Write();
      c1d->Write();
      fOut->Write();
      fOut->Close();
    }

    neutronRates_abs.clear();
    gammaRates_abs.clear();
    neutronRates_mult.clear();
    gammaRates_mult.clear();
    h2dNCalib->Clear();
    h2dNAbsorb->Clear();
    h2dNmult->Clear();
    // cImages->Clear();
    // delete cImages, h2dNCalib, h2dNAbsorb, h2dNCorrected;
    cout << endl;
  }

  theApp.Run();

  return 0;
}
