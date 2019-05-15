//*******************************************************************************
//*
//* Langevin Dynamics
//*	main.cpp
//*
//*
//* Author: Guillaume Le Treut
//*	CEA - 2014
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
// 1) GSL
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>

// 2) BOOST
#include <boost/filesystem.hpp>

// local library
#include "global.h"
#include "utils.h"
#include "linalg.h"
#include "model.h"
#include "stepper.h"

using namespace std;
using namespace utils;
namespace fs = boost::filesystem;

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
	string pathtoinput, pathtooutput, pathtoconf, path, trajname;
	stringstream convert;
  fs::path state_path;
	vector<string> parsechain;
	ifstream fin;
	ofstream fpos_dat, fpos_xyz,fthermo;
	map<string,double> params;
  gsl_rng *rng(0);

  MDWorld *world(0);
  MDStepper *stepper(0);
  ForceField *ffield(0);
	double dt, T, b, lp;
	size_t iter, itermax, idump_pos, idump_thermo, N, iterwidth;

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

  // I/O
  // directories
  // // root directory for trajectories
  path = "trajectory";
  fs::path traj_dir(path.c_str());
  if (fs::create_directory(traj_dir))
  {
    cout << "Directory created: " << traj_dir.string() << endl;
  }
  // // directory for dat format
  path = "dat";
  fs::path traj_dat = traj_dir;
  traj_dat /=  fs::path(path.c_str());
  cout << traj_dat.string() << endl;
  if (fs::create_directory(traj_dat))
  {
    cout << "Directory created: " << traj_dat.string() << endl;
  }
  // // directory for xyz format
  path = "xyz";
  fs::path traj_xyz = traj_dir;
  traj_xyz /=  fs::path(path.c_str());
  cout << traj_xyz.string() << endl;
  if (fs::create_directory(traj_xyz))
  {
    cout << "Directory created: " << traj_xyz.string() << endl;
  }
  // // trajname
  trajname = "state";

  // // thermo
	convert.clear();
	convert.str("");
	convert << "thermo.dat";
  pathtooutput = convert.str();
	fthermo.open(pathtooutput.c_str());
	fthermo << left;
  cout << "output file: " << pathtooutput << endl;

  // parameters
	itermax=params["itermax"];                // number of iterations
	idump_thermo=params["idump_thermo"];      // dump interval
	idump_pos=params["idump_pos"];            // dump interval
	dt=params["dt"];                          // integration time step
  iterwidth = size_t(log10(itermax)+0.5);

	T=params["T"];                            // temperature
	b=params["b"];                            // bond length
	lp=params["lp"];                          // persistence length
	N=params["N"];                            // number of monomers

  /* initialize random number generator */
  rng = gsl_rng_alloc(RDT);
  gsl_rng_set(rng, seed);

  /* initialize MDWorld */
  world = new MDWorld(params["npart"], params["lx"], params["ly"], params["lz"], params["sig_hard_core"]);

  /* initialize MDStepper */
  //stepper = new MDStepper_VVerlet(world, dt, world->m_gamma);
  stepper = new LangStepper_VVerlet(world, rng, dt, world->m_gamma);

  /* initialize force fields */
  /** confinement **/
  //ffield = new ConfinmentBox(-0.5*world->m_lx, +0.5*world->m_lx, -0.5*world->m_ly,+0.5*world->m_ly,-0.5*world->m_lz,+0.5*world->m_lz,1.0,1.0);
  //ffield = new ConfinmentSphere(params["radius_conf"],1.0,1.0);
  //world->m_ffields.push_back(ffield);
  /** polymer **/
//  ffield = new PolymerGaussian(0,world->m_npart,1.0);
//  world->m_ffields.push_back(ffield);
  ffield = new PolymerHarmonic(0,world->m_npart,1000.0,1.);
  world->m_ffields.push_back(ffield);
  ffield = new PolymerKratkyPorod(0,world->m_npart,3.0);
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
      // generate file name
      convert.clear();
      convert.str("");
      convert << trajname << setw(iterwidth) << setfill('0') << setprecision(0) << fixed << dec << iter;
      // conf dump -- standard
      state_path = traj_dat;
      state_path /=  convert.str();
      state_path += ".dat";
      world->dump_dat(state_path.string());
      // conf dump -- xyz format
      state_path = traj_xyz;
      state_path /=  convert.str();
      state_path += ".xyz";
      world->dump_xyz(state_path.string(), iter);
    }

    /* step */
    stepper->step();

    /* update energy */
    // potential energy already updated when calling the forces in the stepper
    world->update_energy_kinetics();

  }

  /* final dump */
  // thermo dump
  fthermo << setw(10) << fixed << setprecision(0) << noshowpos << iter;
  world->dump_thermo(fthermo);
  fthermo << endl;
  // configuration dump
  convert.clear();
  convert.str("");
  convert << trajname << setw(iterwidth) << setfill('0') << setprecision(0) << fixed << dec << iter;
  // conf dump -- standard
  state_path = traj_dat;
  state_path /=  convert.str();
  state_path += ".dat";
  world->dump_dat(state_path.string());
  // conf dump -- xyz format
  state_path = traj_xyz;
  state_path /=  convert.str();
  state_path += ".xyz";
  world->dump_xyz(state_path.string(), iter);

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
