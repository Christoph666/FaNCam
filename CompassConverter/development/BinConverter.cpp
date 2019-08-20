#include <fstream>
#include <iostream>
#include <arpa/inet.h>
#include <vector>
#include <math.h>

using namespace std;

// string dir = "/home/christoph/Documents/MasterThesis/Data/CoMPASS/test/DAQ/run_bin/UNFILTERED/";
// string filename = "0@N6730B #3-11-50_Data_run_bin.bin";

string dir = "/home/christoph/Documents/MasterThesis/Data/CoMPASS/run_AmBe/UNFILTERED/";
string filename = "1@N6730B #3-11-50_Data_run_AmBe.bin";

// struct nBuffer {
// 	uint16_t board;
// 	uint16_t channel;
// 	uint64_t timeStamp;
// 	uint16_t energy;
// 	uint16_t energyShort;
// 	uint16_t flags;
// 	uint32_t nSamples;
// 	uint16_t samples[496];
// } ;
//
// union {
// 	char*		cBuffer;
// 	nBuffer*	myBuffer;
// } myUnion;

// struct nBuffer {
// 	uint16_t board;
// 	uint16_t channel;
// 	uint64_t timeStamp;
// 	uint16_t energy;
// 	uint16_t energyShort;
// 	uint16_t flags;
// 	uint32_t nSamples;
// 	uint16_t samples[496];
// } ;
//
// union {
// 	unsigned char*		cBuffer;
// 	nBuffer*	myBuffer;
// } myUnion;

struct dataBlock {
	uint16_t board=0;
	uint16_t channel=0;
	uint64_t timeStamp=0;
	uint16_t energy=0;
	uint16_t energyShort=0;
	uint16_t flags=0;
	uint32_t nSamples=0;
	vector<uint16_t> samples;
};

dataBlock readBlock(char* buffer, dataBlock myDataBlock, int num_of_samples) {
	int iPow=0;
	for (int i=0; i<2; i++) {
		myDataBlock.board += (unsigned short)(unsigned char)buffer[i] * pow(2,iPow*8);
		iPow++;
	}
	iPow=0;
	for (int i=2; i<4; i++) {
		myDataBlock.channel += (unsigned short)(unsigned char)buffer[i] * pow(2,iPow*8);
		iPow++;
	}
	iPow=0;
	for (int i=4; i<12; i++) {
		myDataBlock.timeStamp += (unsigned short)(unsigned char)buffer[i] * pow(2,iPow*8);
		iPow++;
	}
	iPow=0;
	for (int i=12; i<14; i++) {
		myDataBlock.energy += (unsigned short)(unsigned char)buffer[i] * pow(2,iPow*8);
		iPow++;
	}
	iPow=0;
	for (int i=14; i<16; i++) {
		myDataBlock.energyShort += (unsigned short)(unsigned char)buffer[i] * pow(2,iPow*8);
		iPow++;
	}
	iPow=0;
	for (int i=16; i<18; i++) {
		myDataBlock.flags += (unsigned short)(unsigned char)buffer[i] * pow(2,iPow*8);
		iPow++;
	}
	iPow=0;
	for (int i=18; i<22; i++) {
		myDataBlock.nSamples += (unsigned short)(unsigned char)buffer[i] * pow(2,iPow*8);
		iPow++;
	}
	return myDataBlock;
}



int main ()
{

	cout << "Reading file: " << dir+filename << endl;

	// streampos size;
  // char * memblock;
	//
  // ifstream file ((dir+filename).c_str(), ios::in|ios::binary|ios::ate);
  // if (file.is_open())
  // {
  //   size = file.tellg();
  //   memblock = new char [size];
  //   file.seekg (0, ios::beg);
  //   // file.read (memblock, size);
	//
	// 	int x;
	// 	for (int i=0; i<10; i++) {
	// 		file.seekg(i);
	// 		file.read((&x), 1);
	// 		// cout << i << " : " << static_cast<uint16_t>(x) << endl;
	// 		cout << i << " : " << x << endl;
	// 	}
	//
	// 	file.close();
	// 	cout << "size: " << size << endl;
  //   cout << "the entire file content is in memory" << endl;
	//
  //   delete[] memblock;
  // }
  // else cout << "Unable to open file" << endl;

	// fstream file((dir+filename).c_str(), ios::in | ios::out | ios::binary);
	// file.seekg(3);
	// char x;
	// file.read((&x), 1);
	// cout<< "x: " << x << endl;
	// file.close();

	int length;

	ifstream is;
	// ofstream os;

	is.open((dir+filename).c_str(), ios::binary );

	// get length of file:
	is.seekg (0, ios::end);
	length = is.tellg();
	is.seekg (0, ios::beg);
	cout << "length: " << length << endl;

	// streampos begin,end;
	// begin = is.tellg();
  // is.seekg (0, ios::end);
  // end = is.tellg();
  // cout << "size is: " << (end-begin) << " bytes.\n";


	// allocate memory:
	// myUnion.cBuffer = new unsigned char [1024];

	uint32_t nSamples=0;
	int iPow=0;

	is.read(buffer, 22);
	for (int j=0; j<22; j++) {
			cout << "byte " << j << " :" << (unsigned short)(unsigned char)buffer[j] << endl;
	}

	for (int i=18; i<22; i++) {
		nSamples += (unsigned short)(unsigned char)buffer[i] * pow(2,iPow*8);
		iPow++;
	}
	cout << "nSamples: " 	<< nSamples << endl;


	char buffer[1024];


	// cout << "\t-board: "			<< board << endl;
	// cout << "\t-channel: "		<< channel << endl;
	// cout << "\t-timeStamp: "	<< timeStamp << endl;
	// cout << "\t-energy: "			<< energy << endl;
	// cout << "\t-energyShort: "<< energyShort << endl;
	// cout << "\t-flags: " 			<< flags << endl;
	// cout << "\t-nSamples: " 	<< nSamples << endl;
	// cout << "\t-sample0: " 		<< samples << endl;
	// cout << endl;


	// for (int i=0; i<1; i++) {
	//
	// 	// read data as a block:
	// 	// is.read (myUnion.cBuffer, 1024);
	// 	is.read(buffer, 1024);
	//
	//
	// 	cout << "reading block " 	<< i << endl;
	// 	for (int j=0; j<10; j++) {
	// 		cout << "byte " << j << " :" << (unsigned short)(unsigned char)buffer[j] << endl;
	// 	}
	// 	// cout << "\t-board: "			<< myUnion.myBuffer->board << endl;
	// 	// cout << "\t-channel: "		<< myUnion.myBuffer->channel << endl;
	// 	// cout << "\t-timeStamp: "	<< myUnion.myBuffer->timeStamp << endl;
	// 	// cout << "\t-energy: "			<< myUnion.myBuffer->energy << endl;
	// 	// cout << "\t-energyShort: "<< myUnion.myBuffer->energyShort << endl;
	// 	// cout << "\t-flags: " 			<< myUnion.myBuffer->flags << endl;
	// 	// cout << "\t-nSamples: " 	<< myUnion.myBuffer->nSamples << endl;
	// 	// cout << "\t-sample0: " 		<< myUnion.myBuffer->samples[0] << endl;
	// 	// cout << "\t-sample1: " 		<< myUnion.myBuffer->samples[1] << endl;
	// 	// cout << endl;
	// }

	is.close();

	// delete[] myUnion.cBuffer;
	return 0;
}
