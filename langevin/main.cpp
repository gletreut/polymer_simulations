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
vector<string> required_parameters = vector<string>({"lx", "ly", "lz", "eps_LJ","npart_max", "dt", "itermax", "idump_thermo", "idump_pos"});
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
	ofstream fpos_dat, fpos_xyz,fthermo;
	map<string,double> params;
  gsl_rng *rng(0);

  MDWorld *world(0);
  MDStepper *stepper(0);
  ForceField *ffield(0);
	double dt, T, b, lp;
	size_t iter, itermax, idump_pos, idump_thermo, N;

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
	convert << "thermo.dat";
  pathtooutput = convert.str();
	fthermo.open(pathtooutput.c_str());
	fthermo << left;
  cout << "output file: " << pathtooutput << endl;

	convert.clear();
	convert.str("");
	convert << "positions.dat";
  pathtooutput = convert.str();
	fpos_dat.open(pathtooutput.c_str());
	fpos_dat << left;
  cout << "output file: " << pathtooutput << endl;

	convert.clear();
	convert.str("");
	convert << "positions.xyz";
  pathtooutput = convert.str();
	fpos_xyz.open(pathtooutput.c_str());
	fpos_xyz << left;
  cout << "output file: " << pathtooutput << endl;

	itermax=params["itermax"];                // number of iterations
	idump_thermo=params["idump_thermo"];      // dump interval
	idump_pos=params["idump_pos"];            // dump interval
	dt=params["dt"];                          // integration time step

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
  stepper = new MDStepper_VVerlet(world, dt, 0);

  /* initialize force fields */
  ffield = new ConfinmentBox(-0.5*world->m_lx, +0.5*world->m_lx, -0.5*world->m_ly,+0.5*world->m_ly,-0.5*world->m_lz,+0.5*world->m_lz,1.0,1.0);
  world->m_ffields.push_back(ffield);

  /* initialize simulation */
  world->init_positions();
  world->init_velocities(rng);
  world->update_energy_kinetics();
  stepper->call_forces(); // compute forces at current position

  /* integration */
  for (iter = 0; iter < itermax; ++iter){
    /* dump */
    if (iter % idump_thermo == 0){
      // thermo dump
      fthermo << setw(10) << fixed << setprecision(0) << noshowpos << iter;
      world->dump_thermo(fthermo);
      fthermo << endl;
    }

    if (iter % idump_pos == 0){
      // conf dump -- standard
      world->dump_pos(fpos_dat, true, true, true);
      fpos_dat << endl;
      // conf dump -- xyz format
			fpos_xyz << setw(10) << world->m_npart << endl;
			fpos_xyz << "Atoms. Timestep:" << setw(10) << iter << endl;
      world->dump_pos(fpos_xyz);
      // velocity dump
      //world->dump_vel(cout);
      //cout << endl;
    }

    /* step */
    stepper->step();

    /* update energy */
    // potential energy already updated when calling the forces in the stepper
    world->update_energy_kinetics();

  }

  /* final dump */
  // thermo
  ofstream fthermo_final("thermo_final.dat");
  fthermo_final << left;
  world->dump_thermo(fthermo_final);
  fthermo_final << endl;
  // positions
  ofstream fpos_final("positions_final.dat");
  fpos_final << left << dec << fixed;
  world->dump_pos(fpos_final);

//-----------------------------------------------
//	EXIT
//-----------------------------------------------
	cout << "normal exit" << endl;
  delete stepper;
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
