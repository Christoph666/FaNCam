#include <iostream>
#include <vector>
#include <fstream>
#include <bitset>
#include <sys/stat.h>
// #include <unistd.h>
// #include <chrono>
// #include <ctime>
// #include <sstream>
// #include <algorithm>
// #include <typeinfo>

using namespace std;

// #define buffer_length 2

string dir = "/home/christoph/Documents/MasterThesis/Data/CoMPASS/test/DAQ/run_bin/UNFILTERED/";
string filename = "0@N6730B #3-11-50_Data_run_bin.bin";

const int buffer_length = 128;
// char *buffer;
// buffer = new char[buffer_length];
// char buffer[buffer_length];
char buffer[128];
struct stat results;

int main(){
  uint16_t board;
  uint16_t channel;
  uint64_t timeStamp;
  uint16_t energy;
  uint16_t energyShort;
  uint16_t flags;
  uint32_t nSamples;
  uint16_t oneSample;
  vector<uint16_t> samples;

  // if (stat((dir+filename).c_str(), &results) == 0)
  //   cout << "The size of the file in bytes is: "
  //   << results.st_size << endl;
  // else {
  //   cout << "An error occured." << endl;
  // }

  ifstream myFile((dir+filename).c_str(), ios::in | ios::binary);

  // get length of file
	myFile.seekg(0, myFile.end);
	int length = myFile.tellg();
	myFile.seekg(0, myFile.beg);
  cout << "length: " << length << endl;

  // read from file
  myFile.read(buffer, 2);
  // board = atoi(buffer);
  board = buffer[0] | buffer[1] << 8;
  cout << endl;
  cout << "board:" << board << endl;
  bitset<16> x1(board);
  cout << x1 << endl;

  cout << myFile.tellg() << endl;

  myFile.read(buffer, 2);
  // channel = atoi(buffer);
  channel = buffer[0] | buffer[1] << 8;
  cout << endl;
  cout << "channel:" << channel << endl;
  bitset<16> x2(channel);
  cout << x2 << endl;

  cout << myFile.tellg() << endl;

  myFile.read(buffer, 8);
  // timeStamp = atol(buffer);
  timeStamp = buffer[0] | buffer[1] | buffer[2]| buffer[3]| buffer[4]| buffer[5]| buffer[6]| buffer[7] << 8;
  cout << endl;
  cout << "timeStamp:" << timeStamp << endl;
  // cout << "timeStamp:" << hex(buffer) << endl;
  bitset<64> x3(timeStamp);
  cout << x3 << endl;

  cout << myFile.tellg() << endl;


  myFile.close();
  // delete buffer;

}


// // read binary file
// vector<vector<double>> read_binary (string name, const double vdiv=20E-3, const double voffset=0, const double samplingrate=2E9, const double triggerdelay=0) {
//
// 	ifstream io;
// 	io.open(name, ios::in | ios::binary);
//
// 	vector<vector<double>> data;
// 	vector<double> tvec;
// 	vector<double> Uvec;
// 	const int buffer_length = 1024;
// 	char *buffer;
// 	buffer = new char[buffer_length];
// 	int count = 0;
//
// 	if (!io.good()) {
// 		cout << "Datei nicht gefunden !" << endl;
// 		exit(-1);
// 	}
//
// 	// get length of file
// 	io.seekg(0, io.end);
// 	int length = io.tellg();
// 	io.seekg(0, io.beg);
//
// 	for ( int i = 0; i< (int) length/buffer_length; i++ ) {
// 		io.read(buffer, buffer_length);
// 		for (int k = 0; k<buffer_length; k++) {
// 			double time = (double) count / samplingrate + triggerdelay;
// 			double volt = (double) buffer[k] / 255.0 * vdiv * 7.92 - voffset;
// 			tvec.push_back(time);
// 			Uvec.push_back(volt);
// 			count++;
// 		}
// 	}
//
// 	io.close();
//
// 	// delete buffer ...
// 	delete buffer;
// 	data = {tvec, Uvec};
// 	return data;
// }
