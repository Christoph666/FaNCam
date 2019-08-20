////////////////////////////////////////////////////////////////////////////////
/// CAEN CoMPASS binary file to ROOT file converter                          ///
/// Version 1.0                                                              ///
/// AUTHOR: Christoph Guenther                                               ///
////////////////////////////////////////////////////////////////////////////////

#include <fstream>
#include <iostream>
#include <arpa/inet.h>
#include <vector>
#include <math.h>
#include <stdio.h>

#include <TROOT.h>
#include <TFile.h>
#include <TChain.h>

using namespace std;

///////////////////////////////////////////////////////////////////////////////
// usage: ./CompassConverter measurement_name blocks_per_buffer
// e.g.: ./CompassConverter run_1 100
///////////////////////////////////////////////////////////////////////////////

int main (int argc, char** argv)
{
	string measurement = argv[argc-2];
	int blocks_per_buffer = atoi(argv[argc-1]);

	// cout << "measurement: " << measurement << endl;
	// cout << "blocks_per_buffer: " << blocks_per_buffer << endl;

	// vector<string> vMeasurement = {"Messung1", "Messung1_1"};
	string data_dir = "/home/christoph/Documents/MasterThesis/Data/CoMPASS/"; // change this if necessary
	// string data_dir = "/home/christoph/Documents/MasterThesis/Data/CoMPASS/Stolberg/"; // change this if necessary
	string dir = data_dir+measurement+"/UNFILTERED/";
	vector<string> board_ID {"50", "58"};

	uint16_t board;
	uint16_t channel;
	double timeStamp;
	uint16_t energy;
	uint16_t energyShort;
	uint16_t flags;
	uint32_t nSamples;
	uint16_t sample;
	vector<uint16_t> *samples = new vector<uint16_t>;
	int nEvents = 0;
	// int blocks_per_buffer = 100;
	int blockCounter = 0;
	int bufferCounter = 0;


	for (int iMeasurement=0; iMeasurement<1; iMeasurement++) {
		// string dir = data_dir+measurement+"/DAQ/"+measurement+"_pulse_"+to_string((int)iMeasurement)+"/UNFILTERED/";

		for (int iDEV=0; iDEV<2; iDEV++) {
			for (int iCHN=0; iCHN<8; iCHN++) {

				cout << "Output file: " << data_dir+measurement+(string)"/"+measurement+(string)"_"+to_string((int)iDEV)+(string)"_"+to_string((int)iCHN)+(string)".root" << endl;
				TFile *file = new TFile((data_dir+measurement+(string)"/"+measurement+(string)"_"+to_string((int)iDEV)+(string)"_"+to_string((int)iCHN)+(string)".root").c_str(), "RECREATE");

				// cout << endl << "Output file: " << data_dir+measurement+(string)"/"+measurement+"_pulse_"+to_string((int)iMeasurement)+(string)"_"+to_string((int)iDEV)+(string)"_"+to_string((int)iCHN)+(string)".root" << endl;
				// TFile *file = new TFile((data_dir+measurement+(string)"/"+measurement+"_pulse_"+to_string((int)iMeasurement)+(string)"_"+to_string((int)iDEV)+(string)"_"+to_string((int)iCHN)+(string)".root").c_str(), "RECREATE");
				TTree *tree = new TTree("Events", "Events");

				tree->Branch("device", &board, "device/s");
				tree->Branch("channel", &channel, "channel/s");
				tree->Branch("triggerTimeTag", &timeStamp, "triggerTimeTag/D");
				tree->Branch("energy", &energy, "energy/s");
				tree->Branch("energyShort", &energyShort, "energyShort/s");
				tree->Branch("flags", &flags, "flags/s");
				tree->Branch("nSamples", &nSamples, "nSamples/i");
				tree->Branch("samples", "std::vector<uint16_t>", &samples);

				string filename = (string)"CH_"+to_string(iCHN) + (string)"@N6730B_"+board_ID[iDEV]+(string)"_Data_"+measurement+(string)".bin";
				// string filename = (string)"CH_" + to_string(iCHN) + (string)"@N6730B_"+board_ID[iDEV]+(string)"_Data_"+measurement+(string)".bin";
				// string filename = to_string(iCHN) + (string)"@N6730B #3-11-"+board_ID[iDEV]+(string)"_Data_"+measurement+"_pulse_"+to_string((int)iMeasurement)+(string)".bin";
				cout << "Reading file: " << dir+filename << endl;

				ifstream is;
				is.open((dir+filename).c_str(), ios::binary);

				if ((is.good())) {

					blockCounter=0;
					bufferCounter = 0;

					// get length of file:
					int fileSize;
					is.seekg (0, ios::end);
					fileSize = is.tellg();
					is.seekg (0, ios::beg);
					cout << "File size: " << fileSize/1000 << " kB" << endl;


					// compute number of samples in one block
					char buffer0[24];
					is.read(buffer0, 24);
					// for (int j=0; j<22; j++) {
					// 	cout << "byte " << j << " :" << (unsigned short)(unsigned char) buffer0[j] << endl;
					// }
					uint32_t nSamples_0=0;
					int iPow=0;
					for (int i=20; i<24; i++) {
						nSamples_0 += (unsigned short)(unsigned char) buffer0[i] * pow(2,iPow*8);
						// cout << "byte " << i << " :" << (unsigned short)(unsigned char) buffer0[i]<< endl;
						// cout << "pow " << pow(2,iPow*8) << endl;
						iPow++;
					}
					is.seekg (0, ios::beg);

					cout << "Samples per waveform: " 	<< nSamples_0 << endl;
					// delete buffer0;

					int bytes_per_block = 24+nSamples_0*2;
					cout << "Bytes per block: " 	<< bytes_per_block << " B" << endl;

					int nBlocks = fileSize / bytes_per_block;
					cout << "Blocks in file: " 	<< nBlocks << endl;

					// allocate memory
					char buffer[blocks_per_buffer*bytes_per_block];

					// read the data in multiple buffers
					while(1) {
						if (fileSize==0) {
							cout << "File is empty." << endl << endl;
							break;
						}

						is.read(buffer, blocks_per_buffer*bytes_per_block);
						// cout << "Reading buffer number: " 	<< bufferCounter << endl;
						// cout << "file pos: " << is.tellg() << endl;

						int iStart=0;

						for (int iBlock=0; iBlock<blocks_per_buffer; iBlock++) {
							// for (int iBlock=0; iBlock<10; iBlock++) {

							cout << '\r' << "Reading block: " << blockCounter << " / " << nBlocks << flush;

							board=0;
							channel=0;
							timeStamp=0;
							energy=0;
							energyShort=0;
							flags=0;
							nSamples=0;
							if (nEvents>0){
								samples->clear();
							}

							iStart = iBlock*bytes_per_block;
							// cout << "iStart : " << iStart << endl;

							iPow=0;
							for (int i=0+iStart; i<2+iStart; i++) {
								board += (unsigned short)(unsigned char)buffer[i] * pow(2,iPow*8);
								iPow++;
							}
							iPow=0;
							for (int i=2+iStart; i<4+iStart; i++) {
								channel += (unsigned short)(unsigned char)buffer[i] * pow(2,iPow*8);
								iPow++;
							}
							iPow=0;
							for (int i=4+iStart; i<12+iStart; i++) {
								timeStamp += (unsigned short)(unsigned char)buffer[i] * pow(2,iPow*8);
								iPow++;
							}
							iPow=0;
							for (int i=12+iStart; i<14+iStart; i++) {
								energy += (unsigned short)(unsigned char)buffer[i] * pow(2,iPow*8);
								iPow++;
							}
							iPow=0;
							for (int i=14+iStart; i<16+iStart; i++) {
								energyShort += (unsigned short)(unsigned char)buffer[i] * pow(2,iPow*8);
								iPow++;
							}
							iPow=0;
							for (int i=16+iStart; i<20+iStart; i++) {
								flags += (unsigned short)(unsigned char)buffer[i] * pow(2,iPow*8);
								iPow++;
							}
							iPow=0;
							for (int i=20+iStart; i<24+iStart; i++) {
								nSamples += (unsigned short)(unsigned char)buffer[i] * pow(2,iPow*8);
								iPow++;
							}
							for (int j=1; j<nSamples_0; j++) {
								iPow=0;
								sample=0;
								for (int i=24+iStart+2*j; i<26+iStart+2*j; i++) {
									sample += (unsigned short)(unsigned char)buffer[i] * pow(2,iPow*8);
									iPow++;
								}
								// cout << "sample : " << sample << endl;
								samples->push_back(sample);
							}

							tree->Fill();

							if (nEvents < 0) {
								cout << endl << "\t-board: "			<< board << endl;
								cout << "\t-channel: "		<< channel << endl;
								cout << "\t-timeStamp: "	<< timeStamp << endl;
								cout << "\t-energy: "			<< energy << endl;
								cout << "\t-energyShort: "<< energyShort << endl;
								cout << "\t-flags: " 			<< flags << endl;
								cout << "\t-nSamples: " 	<< nSamples << endl;
								cout << "\t-sample0: " 		<< samples->at(0) << endl;
								cout << "\t-sample1: " 		<< samples->at(1) << endl;
								cout << endl;
							}

							nEvents++;
							blockCounter++;

							if (blockCounter == nBlocks) {
								cout << endl;
								cout << "End of file." << endl;
								break;
							}
						}
						if (blockCounter == nBlocks) {
							break;
						}
						bufferCounter++;
					}
					is.close();
					cout << endl;
				}
				else {
					cout << "File not found." << endl;
					file->Close();
					// cout << "Deleting file: " << data_dir+measurement+(string)"/"+measurement+"_pulse_"+to_string((int)iMeasurement)+(string)"_"+to_string((int)iDEV)+(string)"_"+to_string((int)iCHN)+(string)".root" << endl << endl;
					// remove((data_dir+measurement+(string)"/"+measurement+"_pulse_"+to_string((int)iMeasurement)+(string)"_"+to_string((int)iDEV)+(string)"_"+to_string((int)iCHN)+(string)".root").c_str());
					cout << "Deleting file: " << data_dir+measurement+(string)"/"+measurement+"_"+to_string((int)iDEV)+(string)"_"+to_string((int)iCHN)+(string)".root" << endl << endl;
					remove((data_dir+measurement+(string)"/"+measurement+"_"+to_string((int)iDEV)+(string)"_"+to_string((int)iCHN)+(string)".root").c_str());
				}
				if (nEvents>0){
					file->Write();
					file->Close();
					delete file, tree;
				}
			}
		}
	}

	cout << endl << "Found in total " 	<< nEvents << " events." << endl;
	if (nEvents>0) delete samples;

	return 0;
}
