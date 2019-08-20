#include <math.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <sstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <functional>
#include <unistd.h>
#include <dirent.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TFile.h>
#include <TChain.h>
#include <TVector.h>
#include <TLine.h>
#include <TCanvas.h>
#include <stdint.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TF1.h>
#include <fstream>
#include "utilities.h"


using namespace std;

struct output{
	vector<double> vPeak;
	vector<double> vPSD;
	vector<double> vQS;
	vector<double> vQL;
	vector<int> vValid;
} outputData;

//////////////////////////////////////////////////////////////////////////////
// function that reads the ROOT file
output filereader(vector<string> file, double DeltaS=25, double DeltaL=270, bool draft=0)
{

	bool showPulse = 0;
	bool do_baseline_hist = 0;

	TCanvas *cwave = new TCanvas("cwave", "cwave", 900, 600);
  TCanvas *cBase = new TCanvas("cBase", "cBase", 900, 600);
  cBase->SetLogy();

	TGraph *grS = new TGraph(); // waveform
	TGraph *grPU = new TGraph(); // Pile Up
	TGraph *grP = new TGraph(); // peak
	TH1F *hB = new TH1F("hB", "hB", 200, 1510, 1522); // baseline histogram

	TChain *tree = new TChain("Events");
	for (int iFile=0; iFile<file.size(); iFile++) {
		tree->Add(file[iFile].c_str());
		cout << "Opening file: " << file[iFile] << endl;
	}

	if(tree == nullptr)
  {
		cout << "Tree could not be loaded. Closing program." << endl;
		exit(EXIT_FAILURE);
	}

	//////////////////////////////////////////////////////////////////////////////
	// time separations of the readout boards
	float dt_DRS = 0.5; //ns
	float dt_CAEN = 2.; //ns
	float dt = dt_CAEN;

	// voltage scale
	double VScale = 2./16.384; // 2V over # ADC counts

	//////////////////////////////////////////////////////////////////////////////
	// data storage
	// vector<uint16_t> *vSamples = nullptr;
	// vector<double> vTime;
	// vector<double> vVolt;
	// int N = 768; // samples in one pulse waveform
	// int N_samples = 256;
	// vector<bool> *dp1 = nullptr;
	// vector<bool> *dp2 = nullptr;
  // bool isPUR = false;
  // uint16_t qlong = 0, qshort = 0;
  // uint32_t extras = 0;
  // uint32_t format;
  // uint8_t channel;
  // uint8_t device;

	vector<uint16_t> *vSamples = nullptr;
	vector<double> vTime;
	vector<double> vVolt;
	int N = 768; // samples in one pulse waveform
	int N_samples = 256;
  uint16_t energy = 0, energyShort = 0;
  uint16_t channel;
	uint16_t device;
  uint16_t nSamples;


	// analysis variables
	int i0 = 0;
	double t0 = 0;

	double baselength = 70; // 40 // ns
	double VBase = 0.;

	double Vpeak = 0;
	int ipeak = 0;

	double a, b; // linear coefficients for trigger point
	double tTrg = 0;
	double VTrg = 0;

	int iStart = 0;
	double DeltaOff = -40; // ns
	// double DeltaS = 25; // ns
	// double DeltaL = 270; // ns
	double QS = 0;
	double QL = 0;
	int iQS=0;
	int iQL=0;

	int iPileUp = 0;
	bool is_pileUp = false;
	double PileMax = 400; // mV
	double VPU = 0; // height of pile up

	int N_tot = 0;
	int N_det = 0;
  int N_PU = 0;

	int N_special = 0;

	vector<double> vPSD;
	vector<double> vPeak;
	vector<double> vQS;
	vector<double> vQL;
  vector<vector<double>> result;
	vector<double> vtBase;
	vector<double> vVBase;
	vector<int> vValid;
	int isValid=0;

	//////////////////////////////////////////////////////////////////////////////
	// load data
	tree->SetBranchAddress("device", &device);
	tree->SetBranchAddress("samples", &vSamples);
	tree->SetBranchAddress("channel", &channel);
	tree->SetBranchAddress("nSamples", &nSamples);
	// tree->SetBranchAddress("dp1", &dp1);
	// tree->SetBranchAddress("dp2", &dp2);
  // tree->SetBranchAddress("PUR", &isPUR);
  // tree->SetBranchAddress("qlong", &qlong);
  // tree->SetBranchAddress("qshort", &qshort);
  // tree->SetBranchAddress("extras", &extras);
	// tree->SetBranchAddress("format", &format);


	int treeSize;
	if (draft == false)
	{
		treeSize = tree->GetEntries();
	}
	else
	{
		treeSize = 1e4;
	}
	// cout << "treeSize: " << treeSize << endl;

	//////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////
	// event loop
	//////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////
	for(int iEvent=0; iEvent < treeSize; ++iEvent)
	// for(int iEvent=0; iEvent < 1e5; ++iEvent)
	{
		N_tot++;
		tree->GetEntry(iEvent);
		cout << '\r' << "Reading event: " << iEvent << " / " << treeSize << flush;


		// cout << "i event: " << iEvent << endl;

		// cout << "size: " << vSamples->size() << endl;
		// cout << (int) vSamples->at(0) << endl;
		// if (iEvent==0) {
		// 	N = vSamples->size();
		// 	cout << "N: " << N << endl;
		// }

		////////////////////////////////////////////////////////////////////////////
		// reset variables
		i0=0;
		t0=0;
		ipeak=0;
		iStart=0;
		iQS=0;
		iQL=0;
		VBase=0;
		Vpeak=0;
		VTrg=0;
		iPileUp=0;
		is_pileUp = false;
		VPU = 0;
		isValid=0;
		N=nSamples-1;
		// N=768;
		// cout << "nSamples: " << nSamples << endl;

		////////////////////////////////////////////////////////////////////////////
		// convert to double arrays
		for (int i=0; i<N; i++)
		{
			// if (iEvent==0) {
				// cout << "sample " << i << " : " << vSamples->at(i) << endl;
			// }
			vTime.push_back(i*dt);
			vVolt.push_back(vSamples->at(i) * VScale);
		}

		////////////////////////////////////////////////////////////////////////////
		// compute baseline
		while (t0 < baselength)
		{
			VBase += vVolt[i0];
			vtBase.push_back(t0);
			i0++;
			t0=i0*dt;
		}
		VBase /= i0;
		// cout << VBase << endl;
		hB->Fill(VBase);

		// substract baseline and flip
		for (int i=0; i<N; i++)
		{
			vVolt[i] -= VBase;
			vVolt[i] = -vVolt[i];
		}

		////////////////////////////////////////////////////////////////////////////
		// find peakmaximum
		ipeak=i0;
		Vpeak=vVolt[i0];
		for (int i=i0; i<N/2; i++)
		{
			if (vVolt[i] > Vpeak) {
				Vpeak = vVolt[i];
				ipeak = i;
			}
		}

		// cout << ipeak << endl;
		// cout << vTime[ipeak] << endl;
		// cout << vVolt[ipeak] << endl;

		////////////////////////////////////////////////////////////////////////////
		// find trigger point
		for (int i=i0; i<ipeak; i++)
		{
			if (vVolt[i] < Vpeak/2 && vVolt[i+1] > Vpeak/2)
			{
				a = (vVolt[i+1] - vVolt[i]) / (vTime[i+1] - vTime[i]);
				b = vVolt[i] - a * vTime[i];
				// tTrg = (vTime[i] + vTime[i+1])/2.;
				tTrg = 1/a * (Vpeak/2. - b);
				VTrg = Vpeak/2.;
				break;
			}
		}


		////////////////////////////////////////////////////////////////////////////
		// find integral start point
		iStart=0;
		for (int i=0; i<ipeak; i++)
		{
			if (vTime[i] > tTrg + DeltaOff)
			{
				iStart = i;
				break;
			}
		}

		////////////////////////////////////////////////////////////////////////////
		// check for pile up
		i0=ipeak;
		// look inside QL window, follow falling edge
		while (vVolt[i0+2] <= vVolt[i0] && vTime[i0] < tTrg + DeltaL)
		{
			i0+=2;
			// iPileUp = i0;
		}
		// search for maximum in QL
		while (vTime[i0] < tTrg + DeltaL)
		{
			if (vVolt[i0] > VPU)
			{
				VPU = vVolt[i0];
				iPileUp = i0;
			}
			i0++;
		}
		// check if indeed a pile up
		if (VPU > Vpeak*0.9 && VPU > 20) //0.6
		{
			N_PU++;
			// cout << "detected Pile Up!!" << endl;
			is_pileUp = true;
		}

		////////////////////////////////////////////////////////////////////////////
		// compute integrals
		if (!is_pileUp)
		{
			QS=0;
			QL=0;
			i0=iStart;
			while (vTime[i0] < tTrg + DeltaL)
			{
				if (vTime[i0] < tTrg + DeltaS)
				{
					QS += (vVolt[i0] + vVolt[i0+1]) / 2;
					iQS=i0;
				}
				if (vTime[i0] < tTrg + DeltaL)
				{
					QL += (vVolt[i0] + vVolt[i0+1]) / 2;
					iQL=i0;
				}
				i0++;
			}

			////////////////////////////////////////////////////////////////////////////
			// store PSD variable
			vPeak.push_back(Vpeak);
			vPSD.push_back((QL - QS) / QL);
			vQS.push_back(QS);
			vQL.push_back(QL);
			isValid=1;

			// cout << "QS : " << QS << endl;
			// cout << "QL : " << QL << endl;
			// cout << "PSD : " << vPSD[iEvent];
			// cout << "\tPeak : " << vPeak[iEvent] << endl;
		}
		vValid.push_back(isValid);

		if (showPulse && (QL - QS) / QL > 0.6 && Vpeak > 100.)// && N_special<1e2)
		{
			for (int i=0; i<vVolt.size(); i++)
			{
				grS->SetPoint(grS->GetN(), vTime[i], vVolt[i]);
			}
			if (is_pileUp)
			grPU->SetPoint(grPU->GetN(), vTime[iPileUp], vVolt[iPileUp]);

			grP->SetPoint(grP->GetN(), vTime[ipeak], vVolt[ipeak]);

			N_special++;
		}

		////////////////////////////////////////////////////////////////////////////
		// clear vectors for next event samples
		vTime.clear();
		vVolt.clear();
		vtBase.clear();
		vVBase.clear();
	}

	cout << endl;
	cout << "Total:\t\t" << N_tot << endl;
	cout << "Detected:\t" << N_det << endl;
	cout << "Pile Ups:\t" << N_PU << endl;

	if (showPulse)
	{
		cwave->cd();
		grS->SetMarkerStyle(20);
		grS->SetMarkerSize(0.3);
		grS->SetTitle("Pulse Waveform CAEN");
		grS->GetXaxis()->SetTitle("t / ns");
		grS->GetYaxis()->SetTitle("U / mV");
		grS->Draw("APL");

		grPU->SetMarkerStyle(20);
		grPU->SetMarkerSize(1);
		grPU->SetMarkerColor(kBlue);
		grPU->Draw("same P");

		grP->SetMarkerStyle(23); // triangle down
		grP->SetMarkerColor(kRed);
		grP->SetMarkerSize(1);
		grP->Draw("same p");

	}
	if (!showPulse) delete cwave;

	if (do_baseline_hist)
  {
		cBase->cd();
		hB->SetTitle("Baseline CAEN");
		hB->GetXaxis()->SetTitle("U / mV");
    hB->Draw();
		// cBase->Print("Plots/Baseline/CAEN_baseline.pdf");
		// cBase->Print("Plots/Baseline/CAEN_baseline.svg");
  }
	if (!do_baseline_hist) delete cBase;
	delete tree;

	outputData.vPeak = vPeak;
	outputData.vPSD = vPSD;
	outputData.vQS = vQS;
	outputData.vQL = vQL;
	outputData.vValid = vValid;

	// return {vPeak, vPSD, vQS, vQL, vValid};

	cout << "Detected " << N_special << " special events." << endl;
  return outputData;
}


// check for pile up
// i0=ipeak;
// while (vTime[i0] < tTrg + DeltaL)
// {
// 	if (vVolt[i0+2] < vVolt[i0])
// 	{
// 		i0+=2;
// 	}
// 	else
// 	{
// 		while (vVolt[i0+1] > vVolt[i0])
// 		{
// 			i0++;
// 		}
// 		if (vVolt[i0] > PileMax)
// 		{
// 			iPileUp = i0;
// 			is_pileUp = true;
// 			cout << "detected Pile Up!!" << endl;
// 		}
// 	}
// 	i0++;
// }
