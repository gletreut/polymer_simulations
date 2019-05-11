//----------------------------------------------------------------------------
//  2014
//  Sampling of polymer Boltzmann ensemble with Monte-Carlo.
//
//  G. Le Treut - IPhT, CEA
//
//
//  Compilation with:
//  g++ -std=c++11 main.cpp
//----------------------------------------------------------------------------

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <limits>
#include <map>

#include "model.cpp"
#include "utils.cpp"
#include "geometry.cpp"

using namespace std;
using namespace string_utils;

//** GLOBAL
double macheps=std::numeric_limits<double>::epsilon();

//** MAIN
int main(int argc, char *argv[]){
	string pathtoinput, pathtooutput, name;
	stringstream convert;
	vector<string> parsechain;
	ifstream fin;
	ofstream fpol,fen;
	map<string,double> params;
	double dt, T, lp;
	long long int iter;
	int dump,N;
	vector<vector<double> > kmap;

//-----------------------------------------------
//** INITIALIZATION
//-----------------------------------------------
	cout << "Enter parameters values in the form: <key> <value> for the following keys: (iter,N,T,lp)" << endl;
	load_params<double>(cin,params);
	cout << "You entered the following parameters:" << endl;
	print_map(cout,params);

//	parsechain=split(pathtoinput,'.'); // remove extension
//	if (parsechain[0] == "/")
//		name="./";
//	else
//		name=parsechain[0];
//	for (size_t i=1; i<parsechain.size()-1; i++)
//		name+="." + parsechain[i];

	convert.clear();
	convert.str("");
	convert << "pol.xyz";
  pathtooutput = convert.str();
	fpol.open(pathtooutput.c_str());
	fpol << left;
  cout << "output file: " << pathtooutput << endl;

	fin.open(pathtoinput.c_str());
	load_map<double>(fin,kmap);
	fin.close();

	iter=params["iter"];
	dump=params["dump"];
	T=params["T"];
	lp=params["lp"];
	N=params["N"];
  print_map(cout, params);

//-----------------------------------------------
//	METROPOLIS MONTE-CARLO
//-----------------------------------------------
	Polymer pol(N,lp);
	for (int i=0; i<=iter; i++){
		pol.gen();
		pol.accept(T);

		if (i%dump == 0){
			cout << setw(5) << "i=" << setw(10) << i;
			cout << setw(5) << "E=" << setw(20) << pol.E << endl;

			fpol << setw(10) << pol.N << endl;
			fpol << "Atoms. Timestep:" << setw(10) << i << endl;
			pol.print_current(fpol);
		}
	}
	//load_map<double>(fin,cmap_raw);


//-----------------------------------------------
//	EXIT
//-----------------------------------------------
	cout << "normal exit" << endl;
	return 0;
}
