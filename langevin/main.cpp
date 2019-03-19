//*******************************************************************************
//*
//* Langevin Dynamics
//*	main.cpp
//*
//*
//* Author: Guillaume Le Treut
//*	CEA - 2014
//*
//* Compilation with:
//* g++ -std=c++11 main.cpp
//*
//*******************************************************************************
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
	string pathtoinput, pathtooutput, pathtoconf="", name;
	stringstream convert;
	vector<string> parsechain;
	ifstream fin;
	ofstream fpol,fen;
	map<string,double> params;
	double dt, T, Ke, lp;
	long long int iter;
	int dump, N;

//-----------------------------------------------
//** INITIALIZATION
//-----------------------------------------------
	if ( (argc != 1) && (argc != 2) ){
		cout << "SYNTAX ERROR: arguments are [<initial conf>]" << endl;
		return 1;
	}
		if (argc == 2){
			convert.clear();
			convert.str(argv[1]);
			convert >> pathtoconf;
		}

	cout << "Enter parameters values in the form: <key> <value> for the following keys: (iter,dt,T,N,Ke,lp,dump)" << endl;
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
	dt=params["dt"];
	T=params["T"];
	Ke=params["Ke"];
	lp=params["lp"];
	N=params["N"];

//-----------------------------------------------
//	INTEGRATION
//-----------------------------------------------
	Polymer pol(N,Ke,lp);
	if (pathtoconf != ""){
		cout << "pathtoconf:" << pathtoconf << endl;
		fin.open(pathtoconf.c_str());
		pol.init(fin);
		pol.recenter();
		fin.close();
	}
	for (int i=0; i<=iter; i++){
		pol.step(T,dt);
		pol.update();

		if (i%dump == 0){
			cout << setw(5) << "i=" << setw(10) << i;
			cout << setw(5) << "E=" << setw(20) << pol.energy() << endl;

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
