#ifndef _INC_GUARD
#define _INC_GUARD

#include <vector>
#include <iostream>
#include <fstream>
#include "TMath.h"
#include "TComplex.h"
#include "TGaxis.h"
#include "TPaveText.h"
#include "TF1.h"


using namespace std;

/* Linspace Vektoren erstellen */
vector<Double_t> linspace (const Double_t &start, const Double_t &end, const unsigned &n) {
	Double_t step = (end - start)/(n-1);
	vector<Double_t> x;

	for (unsigned i=0; i<n; i++) {
		x.push_back(start + i * step);
	}
	return x;
}

/* Konst Vector erstellen */
vector<Double_t> konst_vector (const Double_t &eintrag, const unsigned &n){
	vector<Double_t> x;

	for (unsigned i=0;i<n;i++){
		x.push_back(eintrag);
	}
	return x;
}

/* liest eine Textdatei mit n spalten gefüllt mit zahlen aus */
vector<vector<Double_t> > read(string name, const unsigned spalten)
{
	ifstream data(name, ios_base::in);
	if(!data.good())
	{
		cout << "Datei nicht gefunden !" << endl;
		exit(-1);
	}
	vector<vector<Double_t> > v(spalten);
	vector<Double_t> tmp(spalten);
	bool leer = false;
	while (!leer)
	{
		for(unsigned j=0;j<spalten;j++)
		{
			if(!(data >> tmp[j]))
			{
				leer = true;
			}
		}
		if(!leer)
		{
			for(unsigned j=0;j<spalten;j++)
			{
				v[j].push_back(tmp[j]);
			}
		}
	}
	return v;
}

/* gibt einen vector mit einträgen Double_t in der Konsole aus */
void print_vector(const vector<Double_t> &v){
	int size = v.size();
	if(size==0){
		cout << "vector hat die Laenge 0 " << endl;
	}
	else{
		cout << "[ ";
		for(int i=0;i<size-1;i++){
			cout << v[i] << ", ";
		}
		cout << v[size-1] << " ]" << endl;
	}
}

/* gewichtetes Mittel */
vector<Double_t> gew_mit(vector<Double_t> x, vector<Double_t> ex){
	int len = x.size();       	//Laenge des Eingangsvektors
	Double_t s = 0;               	//Hilfsvektor
	Double_t wsum = 0;            	//Summe der w-Eintraege
	Double_t sx;              	//Fehler auf den Mittelwert
	Double_t xm;              	//Mittelwert
	vector<Double_t> w;       	//Fehler**2
	vector<Double_t> result;  	//Ergebnisvektor

	for (int i = 0; i<len; i++){
		w.push_back(1./(ex[i]*ex[i]));
		s += (w[i] * x[i]);
		wsum += w[i];
	}

	xm = s/wsum;
	sx = sqrt(1./wsum);
	result.push_back(xm);
	result.push_back(sx);

	return result;
}

/* Mittel : gibt  Mittelwert, Fehler auf die Einzelmessung, Fehler auf den Mittelwert  wieder*/
vector<Double_t> mit(const vector<Double_t> &x)
{
	const unsigned len = x.size();
	Double_t x_mid = 0;
	Double_t sig_i_2 = 0;
	vector<Double_t> result;
	for (unsigned i=0;i<len;i++)
	{
		x_mid += x[i]/len;
	}
	for (unsigned i=0;i<len;i++)
	{
		sig_i_2 += pow(x[i]-x_mid,2)/(len-1);
	}
	Double_t sig_x = sqrt(sig_i_2 / len);
	result.push_back(x_mid);
	result.push_back(sqrt(sig_i_2));
	result.push_back(sig_x);
	return result;
}



/* Überladung der operatoren *,+,- */
vector<Double_t> operator* (const vector<Double_t> &v1, const vector<Double_t> &v2){
	if(!(v1.size()==v2.size())){
		cout << "vectoren haben nicht die selbe groesse !!" << endl;
		exit(-1);
	}
	int size = v1.size();
	vector<Double_t> v3;
	for(int i=0;i<size;i++){
		v3.push_back(v1[i] * v2[i]);
	}
	return v3;
}

vector<Double_t> operator* (Double_t zahl, vector<Double_t> v){
	int size = v.size();

	for(int i=0;i<size;i++){
		v[i] = zahl * v[i];
	}
	return v;
}

vector<Double_t> operator* (vector<Double_t> v, Double_t zahl){
	int size = v.size();

	for(int i=0;i<size;i++){
		v[i] = v[i] * zahl;
	}
	return v;
}

vector<Double_t> operator+ (const vector<Double_t> &v1,const vector<Double_t> &v2){
	if(!(v1.size()==v2.size())){
		cout << "vectoren haben nicht die selbe groesse !!" << endl;
		exit(-1);
	}
	int size = v1.size();
	vector<Double_t> v3;
	for(int i=0;i<size;i++){
		v3.push_back(v1[i] + v2[i]);
	}
	return v3;
}

vector<Double_t> operator+ (Double_t zahl, vector<Double_t> v){
	int size = v.size();
	for(int i=0;i<size;i++){
		v[i] = zahl + v[i];
	}
	return v;
}

vector<Double_t> operator+ (vector<Double_t> v, Double_t zahl){
	int size = v.size();

	for(int i=0;i<size;i++){
		v[i] = v[i] + zahl;
	}
	return v;
}

vector<Double_t> operator- (const vector<Double_t> &v1, const vector<Double_t> &v2){
	if(!(v1.size()==v2.size())){
		cout << "vectoren haben nicht die selbe groesse !!" << endl;
		exit(-1);
	}
	int size = v1.size();
	vector<Double_t> v3;
	for(int i=0;i<size;i++){
		v3.push_back(v1[i] - v2[i]);
	}
	return v3;
}

vector<Double_t> operator- (Double_t zahl, vector<Double_t> v){
	int size = v.size();

	for(int i=0;i<size;i++){
		v[i] = zahl - v[i];
	}
	return v;
}

vector<Double_t> operator- (vector<Double_t> v, Double_t zahl){
	int size = v.size();

	for(int i=0;i<size;i++){
		v[i] = v[i] - zahl;
	}
	return v;
}

vector<Double_t> operator/ (Double_t zahl, vector<Double_t> v){
	int size = v.size();

	for(int i=0;i<size;i++){
		v[i] = zahl / v[i];
	}
	return v;
}
/* ende der überladung */



/* funktionen von vektoren */
vector<Double_t> exp(vector<Double_t> v)
{
	int size = v.size();
	for(int i=0;i<size;i++)
	{
		v[i] = exp(v[i]);
	}
	return v;
}

vector<Double_t> log(vector<Double_t> v)
{
	int size = v.size();
	for(int i=0;i<size;i++)
	{
		v[i] = log(v[i]);
	}
	return v;
}

vector<Double_t> sqrt(vector<Double_t> v)
{
	int size = v.size();
	for(int i=0;i<size;i++)
	{
		v[i] = sqrt(v[i]);
	}
	return v;
}

/* slice funktionen */
vector<Double_t> slice1 (const vector<Double_t> &x, const Double_t &start, const Double_t &end)
{
	const Int_t len = x.size();
	vector<Double_t> resvec;
	for (Int_t i = 0; i<len; i++)
	{
		if(x[i]>start && x[i]<end)
		{
			resvec.push_back(x[i]);
		}
	}
	return resvec;
}

vector<vector<Double_t> > slice2 (const vector<Double_t> &x, const vector<Double_t> &y, const Double_t &start, const Double_t &end)
{
	if (end > start)
	{
		vector<vector<Double_t> > resvecvec(2);
		const Int_t len = x.size();
		for (Int_t i = 0; i<len; i++)
		{
			if(x[i]>start && x[i]<end)
			{
				resvecvec[0].push_back(x[i]);
				resvecvec[1].push_back(y[i]);
			}
		}
		return resvecvec;
	}
	else
	{
		cout << "da ist was beim slicen schief gegangen :(" << endl;
		exit(-1);
	}
}

vector<vector<Double_t> > slice4 (const vector<Double_t> &x, const vector<Double_t> &y, const vector<Double_t> &ex, const vector<Double_t> &ey, const Double_t &start, const Double_t &end)
{
	if (end > start)
	{
		vector<vector<Double_t> > resvecvec(4);
		const unsigned len = x.size();
		for (unsigned i = 0; i<len; i++)
		{
			if(x[i]>start && x[i]<end)
			{
				resvecvec[0].push_back(x[i]);
				resvecvec[1].push_back(y[i]);
				resvecvec[2].push_back(ex[i]);
				resvecvec[3].push_back(ey[i]);
			}
		}
		return resvecvec;
	}
	else
	{
		cout << "da ist was beim slicen schief gegangen :(" << endl;
		exit(-1);
	}
}

// Abstandsfunktionen von vectoren
vector<Double_t> abs(const vector<Double_t> &x){
	const unsigned size = x.size();
	vector<Double_t> abs;
	for(unsigned i=0;i<size;i++){
		abs.push_back(x[i]*x[i]);
	}
	return sqrt(abs);
}

vector<Double_t> abs(const vector<TComplex> &x){
	const unsigned size = x.size();
	vector<Double_t> abs;
	for(unsigned i=0;i<size;i++){
		abs.push_back(TComplex::Abs(x[i]));
	}
	return abs;
}

// diskrete fouriertransformation
vector<TComplex> DFT(const vector<TComplex> &x){
	const Double_t Pi = 3.141592654;
	const unsigned N = x.size();
	TComplex I = TComplex(0, 1);
	vector<TComplex> X;
	for(unsigned k=0;k<N;k++){
		TComplex X_k = TComplex(0, 0);
		for(unsigned n=0;n<N;n++){
			X_k +=  x[n] * TComplex::Exp(-I * ((2.*Pi/N) * k * n) );
		}
		X.push_back(X_k);
	}
	return X;
}

vector<TComplex> DFT_b(const vector<TComplex> &x){
	const Double_t Pi = 3.141592654;
	const unsigned N = x.size();
	TComplex I = TComplex(0, 1);
	vector<TComplex> X;
	for(unsigned k=0;k<N;k++){
		TComplex X_k = TComplex(0, 0);
		for(unsigned n=0;n<N;n++){
			X_k += x[n] * TComplex::Exp(I * ((2.*Pi/N) * k * n) );
		}
		X.push_back(X_k);
	}
	for(unsigned i=0;i<N;i++){
		X[i] = X[i]/Double_t(N);
	}
	return X;
}

vector<TComplex> DFT(const vector<Double_t> &x){
	const Double_t Pi = 3.141592654;
	const unsigned N = x.size();
	TComplex I = TComplex(0, 1);
	vector<TComplex> X;
	for(unsigned k=0;k<N;k++){
		TComplex X_k = TComplex(0, 0);
		for(unsigned n=0;n<N;n++){
			X_k +=  x[n] * TComplex::Exp(-I * ((2.*Pi/N) * k * n) );
		}
		X.push_back(X_k);
	}
	return X;
}

vector<TComplex> DFT_b(const vector<Double_t> &x){
	const Double_t Pi = 3.141592654;
	const unsigned N = x.size();
	TComplex I = TComplex(0, 1);
	vector<TComplex> X;
	for(unsigned k=0;k<N;k++){
		TComplex X_k = TComplex(0, 0);
		for(unsigned n=0;n<N;n++){
			X_k += x[n] * TComplex::Exp(I * ((2.*Pi/N) * k * n) );
		}
		X.push_back(X_k);
	}
	for(unsigned i=0;i<N;i++){
		X[i] = X[i]/Double_t(N);
	}
	return X;
}

// liest dateien für versuch T6 ein
vector<Double_t> read_T6(string pfd){
	ifstream data(pfd);
	if(!data.good()){
		cout << "Datei nicht gefunden !" << endl;
		exit(-1);
	}
	vector<Double_t> v;
	string line;
	string end2 = "<<DATA>>" ;
	if(data.is_open()){
		while(data){
			getline(data,line);
			line.pop_back();
			if(line==end2){
				while(data){
					Double_t tmp;
					data >> tmp;
					v.push_back(tmp);
				}
				v.pop_back();
				return v;
			}
		}
	}
	return v;
}

// setzt alle einträge eines vectors zwischen start und end auf 0, Double_t
vector<Double_t> set_0(vector<Double_t> v, const unsigned start, const unsigned end){
	for(unsigned i=start;i<end;i++){
		v[i] = 0;
	}
	return v;
}
vector<TComplex> set_0(vector<TComplex> v, const unsigned start, const unsigned end){
	for(unsigned i=start;i<end;i++){
		v[i] = 0;
	}
	return v;
}





/*
vector<Double_t> Gauss(const Double_t &A,const Double_t &sigma,const Double_t &mid,const vector<Double_t> &x){
	unsigned size = x.size();
	vector<Double_t> erg;
	for(unsigned i=0;i<size;i++){
		erg.push_back( A * TComplex::Exp( - (x[i] - mid) * (x[i] - mid) / ( 2. * sigma * sigma ) ) );
	}
	return erg;
}
*/

//fancy Paramter Box
void parbox2 (Double_t chiq=0, Double_t *par=0, Double_t *parerr=0){
	TPaveText *b1 = new TPaveText(0.67,0.73,0.89,0.90,"NDC");
	char buffer[250];
	sprintf(buffer,"#chi^{2}/ndf = %.2f", chiq);
	b1->AddText(0., 5., buffer);
	sprintf(buffer,"A = %.2f #pm %.2f", par[1],parerr[1]);
	b1->AddText(5., 10., buffer);
	sprintf(buffer,"B = %.2f #pm %.2f", par[0],parerr[0]);
	b1->AddText(5., 10., buffer);
	b1->Draw();
}

void parbox3 (Double_t chiq=0, Double_t *par=0, Double_t *parerr=0){
	TPaveText *b1 = new TPaveText(0.67,0.73,0.89,0.90,"NDC");
	char buffer[250];
	sprintf(buffer,"#chi^{2}/ndf = %.2f", chiq);
	b1->AddText(0., 5., buffer);
	sprintf(buffer,"A = %.2f #pm %.2f", par[0],parerr[0]);
	b1->AddText(5., 10., buffer);
	sprintf(buffer,"#mu = %.2f #pm %.2f", par[1],parerr[1]);
	b1->AddText(5., 10., buffer);
	sprintf(buffer,"#sigma = %.2f #pm %.2f", par[2],parerr[2]);
	b1->AddText(5., 10., buffer);
	b1->Draw();
}

void parbox6 (Double_t chiq=0, Double_t *par=0, Double_t *parerr=0){
	TPaveText *b1 = new TPaveText(0.67,0.65,0.90,0.85,"NDC");
	char buffer[250];
	sprintf(buffer,"#chi^{2}_{1}/ndf = %.2f", chiq); //3.14->GetChisquare()
	b1->AddText(0., 5., buffer);
	sprintf(buffer,"A_{1} = %.2f #pm %.2f", par[0],parerr[0]);
	b1->AddText(5., 10., buffer);
	sprintf(buffer,"#mu_{1} = %.2f #pm %.2f", par[1],parerr[1]);
	b1->AddText(5., 10., buffer);
	sprintf(buffer,"#sigma_{1} = %.2f #pm %.2f", par[2],parerr[2]);
	b1->AddText(5., 10., buffer);
	sprintf(buffer,"A_{2} = %.2f #pm %.2f", par[3],parerr[3]);
	b1->AddText(5., 10., buffer);
	sprintf(buffer,"#mu_{2} = %.2f #pm %.2f", par[4],parerr[4]);
	b1->AddText(5., 10., buffer);
	sprintf(buffer,"#sigma_{2} = %.2f #pm %.2f", par[5],parerr[5]);
	b1->AddText(5., 10., buffer);
	b1->Draw();
}

Double_t max(const vector<Double_t> &v){
	Double_t tmp = v[0];
	const Int_t size = v.size();
	for(Int_t i=0;i<size;i++){
		if(tmp <= v[i]){
			tmp = v[i];
		}
	}
	return tmp;
}


vector<Double_t> real(const vector<TComplex> &x){
	vector<Double_t> re;
	const unsigned size = x.size();
	for(unsigned i=0;i<size;i++){
		re.push_back(x[i].Re());
	}
	return re;
}

Double_t Gau(Double_t *x, Double_t *par){
	return par[0] * exp( - pow(x[0] - par[1],2) / (2. * pow(par[2], 2)) );
}

Double_t double_Gau(Double_t *x, Double_t *par){
	return Gau(x, par) + Gau(x, &par[3]);
}

vector<Double_t> double_Gau(const vector<Double_t> &x, Double_t *par){
	vector<Double_t> d_gau;
	const Int_t size = x.size();
	for(Int_t i=0;i<size;i++){
		d_gau.push_back( par[0] * exp( - pow(x[i] - par[1],2) / (2. * pow(par[2], 2)) )
						+par[3] * exp( - pow(x[i] - par[4],2) / (2. * pow(par[5], 2)) ) );
	}
	return d_gau;
}

Double_t dblGaussian(Double_t *x, Double_t *par){
	return par[0] * exp(-pow(x[0] - par[1], 2) / (2.*pow(par[2], 2)) )
	+ par[3] * exp(-pow(x[0] - par[4], 2) / (2.*pow(par[5], 2)) );
}

Double_t tripleGaussian(Double_t *x, Double_t *par){
	return par[0] * exp(-pow(x[0] - par[1], 2) / (2.*pow(par[2], 2)) )
	+ par[3] * exp(-pow(x[0] - par[4], 2) / (2.*pow(par[5], 2)) )
	+ par[6] * exp(-pow(x[0] - par[7], 2) / (2.*pow(par[8], 2)) );
}

Double_t fLin(Double_t *x, Double_t *par){
	return par[1] * x[0] + par[0];
}

vector<Double_t> Gau(const vector<Double_t> &x, Double_t *par)
{
	const Int_t size = x.size();
	vector<Double_t> aus;
	for(Int_t i=0;i<size;i++){
		aus.push_back( par[0] * exp(- 0.5 * pow((x[i]-par[1])/par[2], 2.)) );
	}
	return aus;
}


Double_t linear(Double_t *x, Double_t *par)
{
	return par[1] * x[0] + par[0];
}

Double_t linear(Double_t x, Double_t par0, Double_t par1)
{
	return par1 * x + par0;
}

double interpolate1d(double x0, vector<double> x, vector<double> y)
{
	if (x0<x[0] | x0>x[x.size()-1])
	{
		cout << "Interpolate error: Out of range.\tx0 = " << x0 << endl;
		return 0;
	}
	double a, b, y0;
	for (int i=0; i<x.size()-1; i++)
	{
		if (x[i]<x0 && x[i+1]>=x0)
		{
			a = (y[i+1]-y[i]) / (x[i+1]-x[i]);
			b = y[i] - a*x[i];
			y0 = a*x0 + b;
		}
	}
	return y0;
}

double min2v (vector<double> x1, vector<double> x2) {
	vector<double> x;
	x.reserve( x1.size() + x2.size() ); // preallocate memory
	x.insert( x.end(), x1.begin(), x1.end() );
	x.insert( x.end(), x2.begin(), x2.end() );
	int n = x.size();
	double result=x[0];
	for (int i=0; i<n; i++) {
		if (x[i] < result) {
			result = x[i];
		}
	}
	return result*0.9;
}

double max2v (vector<double> x1, vector<double> x2) {
	vector<double> x;
	x.reserve( x1.size() + x2.size() ); // preallocate memory
	x.insert( x.end(), x1.begin(), x1.end() );
	x.insert( x.end(), x2.begin(), x2.end() );
	int n = x.size();
	double result=x[0];
	for (int i=0; i<n; i++) {
		if (x[i] > result) {
			result = x[i];
		}
	}
	return result*1.1;
}

#endif
