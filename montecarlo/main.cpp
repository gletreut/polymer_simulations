//----------------------------------------------------------------------------
//  2014
//  Sampling of polymer Boltzmann ensemble with Monte-Carlo.
//
//  G. Le Treut - IPhT, CEA
//
//
//  Build with `make` (requires GSL); run with:
//  ../bin/simu_montecarlo < ../configs/config_montecarlo
//----------------------------------------------------------------------------

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <limits>
#include <map>

#include "model.h"
#include "utils.h"

using namespace std;
using namespace string_utils;

//** GLOBAL
double macheps=std::numeric_limits<double>::epsilon();

//** MAIN
int main(int argc, char *argv[]){
	string pathtooutput;
	stringstream convert;
	ofstream fpol;
	map<string,double> params;
	double T, lp;
	long long int iter;
	int dump,N;
	unsigned long seed;

//-----------------------------------------------
//** INITIALIZATION
//-----------------------------------------------
	cout << "Enter parameters values in the form: <key> <value> for the following keys: (iter,N,T,lp,dump[,seed])" << endl;
	load_params<double>(cin,params);
	cout << "You entered the following parameters:" << endl;
	print_map(cout,params);

	convert.clear();
	convert.str("");
	convert << "pol.xyz";
  pathtooutput = convert.str();
	fpol.open(pathtooutput.c_str());
	fpol << left;
  cout << "output file: " << pathtooutput << endl;

	iter=params["iter"];
	dump=params["dump"];
	T=params["T"];
	lp=params["lp"];
	N=params["N"];
	seed=params.count("seed") ? (unsigned long) params["seed"] : 12345;
  print_map(cout, params);

//-----------------------------------------------
//	METROPOLIS MONTE-CARLO
//-----------------------------------------------
	Polymer pol(N,lp,seed);
	for (long long int i=0; i<=iter; i++){
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


//-----------------------------------------------
//	EXIT
//-----------------------------------------------
	cout << "normal exit" << endl;
	return 0;
}
