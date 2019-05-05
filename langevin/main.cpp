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

// external library
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>

// local library
#include "global.h"
#include "utils.h"
#include "linalg.h"
#include "model.h"
#include "stepper.h"

using namespace std;
using namespace utils;

//** GLOBAL
double macheps=std::numeric_limits<double>::epsilon();
vector<string> required_parameters = vector<string>({"lx", "ly", "lz", "eps_LJ","npart_max", "dt", "itermax", "idump_thermo", "idump_thermo"});
/* random number generator type */
const gsl_rng_type* RDT = gsl_rng_ranlxs1;
//const gsl_rng_type* RDT = gsl_rng_mt19937;
const size_t seed = 123;


//** FUNCTION
void check_params_keys(map<string,double> myparams, vector<string> list);

//** MAIN
int main(int argc, char *argv[]){
	string pathtoinput, pathtooutput, pathtoconf="", name;
	stringstream convert;
	vector<string> parsechain;
	ifstream fin;
	ofstream fpol,fen;
	map<string,double> params;
  gsl_rng *rng(0);

  MDWorld *world(0);
	double dt, T, b, lp;
	size_t iter, idump, N;

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

	cout << "Enter parameters values in the form: <key> <value> for the following keys: (iter,dt,T,N,b,lp,dump)" << endl;
	load_map<double>(cin,params);
	cout << "You entered the following parameters:" << endl;
	print_map<double>(cout,params);

  check_params_keys(params, required_parameters);

	convert.clear();
	convert.str("");
	convert << "pol.xyz";
  pathtooutput = convert.str();
	fpol.open(pathtooutput.c_str());
	fpol << left;
  cout << "output file: " << pathtooutput << endl;

	iter=params["iter"];      // number of iterations
	idump=params["idump"];    // dump interval
	dt=params["dt"];          // integration time step

	T=params["T"];            // temperature
	b=params["b"];            // bond length
	lp=params["lp"];          // persistence length
	N=params["N"];            // number of monomers

  /* initialize random number generator */
  rng = gsl_rng_alloc(RDT);
  gsl_rng_set(rng, seed);

  /* initialize MDWorld */
  world = new MDWorld(params["npart"], params["lx"], params["ly"], params["lz"], params["sig_hard_core"]);

  /* initialize MDStepper */

  /* initialize force fields */

  /* initialize simulation */
  //init_positions();
  //init_velocities(rng);

  /* integration */
//-----------------------------------------------
//	INTEGRATION
//-----------------------------------------------
//	Polymer pol(N,Ke,lp);
//	if (pathtoconf != ""){
//		cout << "pathtoconf:" << pathtoconf << endl;
//		fin.open(pathtoconf.c_str());
//		pol.init(fin);
//		pol.recenter();
//		fin.close();
//	}
//	for (int i=0; i<=iter; i++){
//		pol.step(T,dt);
//		pol.update();
//
//		if (i%dump == 0){
//			cout << setw(5) << "i=" << setw(10) << i;
//			cout << setw(5) << "E=" << setw(20) << pol.energy() << endl;
//
//			fpol << setw(10) << pol.N << endl;
//			fpol << "Atoms. Timestep:" << setw(10) << i << endl;
//			pol.print_current(fpol);
//		}
//	}
//	//load_map<double>(fin,cmap_raw);


//-----------------------------------------------
//	EXIT
//-----------------------------------------------
	cout << "normal exit" << endl;
  delete world;
  gsl_rng_free(rng);
	return 0;
}

void check_params_keys(map<string,double> myparams, vector<string> list)
{
  for (vector<string>::iterator it = list.begin(); it != list.end(); ++it){
    map<string,double>::iterator iter = myparams.find(*it);
    if (iter == myparams.end()){
      cout << "key: " << *it << endl;
      throw invalid_argument("Required key not found in parameter map!");
    }
  }
  return;
}
