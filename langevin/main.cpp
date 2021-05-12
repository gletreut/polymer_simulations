//*******************************************************************************
//*
//* Langevin Dynamics
//*	main.cpp
//*
//*
//* Author: Guillaume Le Treut
//*	UCSD - 2019
//*
//*******************************************************************************
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <limits>
#include <vector>
#include <map>

// external library
// 1) GSL
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>

// 2) BOOST
#include <boost/filesystem.hpp>

// 3) yaml-cpp
#include <yaml-cpp/yaml.h>

// local library
#include "global.h"
#include "utils.h"
#include "linalg.h"
#include "model.h"
#include "stepper.h"
#include "integration.h"
#include "yaml_config.h"

using namespace std;
using namespace utils;
namespace fs = boost::filesystem;
namespace yc = yaml_config;

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
	string pathtoinput, pathtooutput, pathtoinitconf, pathname, trajname;
	stringstream convert;
  fs::path state_path;
	vector<string> parsechain;
	ifstream fin;
	ofstream fpos_dat, fpos_xyz,fthermo;
	map<string,double> params;
  YAML::Node config;
  gsl_rng *rng(0);

  MDWorld *world(0);
  MDStepper *stepper(0);
  IntegrationParams iparams;
	size_t iter;

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
			convert >> pathtoinitconf;
		}
    else {
      pathtoinitconf = "";
    }

  /* load parameters */
  cout << "Enter parameters:" << endl;
  config = YAML::Load(cin);  // load yaml from standart input
  cout << config << endl;

  /* initialize random number generator */
  rng = gsl_rng_alloc(RDT);
  yc::init_rng(config, rng);

  /* initialize MDWorld */
  world = new MDWorld();
  yc::init_world(config, rng, world);
  world->print_infos(cout);

  /* initialize MDStepper */
  yc::init_stepper(config, rng, world, stepper);
  stepper->print_infos(cout);

  /* load force fields */
  yc::init_forcefields(config, world);

  /* load constraints */
  yc::init_constraints(config, world);

  /* init integration parameters */
  yc::init_integration_params(config, iparams);
  iparams.print_infos(cout);

  /* I/O */
  // root directory for trajectories
  pathname = "trajectory";
  fs::path traj_dir(pathname.c_str());
  if (fs::create_directory(traj_dir))
  {
    cout << "Directory created: " << traj_dir.string() << endl;
  }
  // directory for dat format
  pathname = "dat";
  fs::path traj_dat = traj_dir;
  traj_dat /=  fs::path(pathname.c_str());
  cout << traj_dat.string() << endl;
  if (iparams.pos_dat) {
    if (fs::create_directory(traj_dat))
    {
      cout << "Directory created: " << traj_dat.string() << endl;
    }

    // remove all existing configurations
    for (fs::directory_iterator end_dir_it, it(traj_dat); it!=end_dir_it; ++it) {
      fs::remove_all(it->path());
    }
  }

  // directory for xyz format
  pathname = "xyz";
  fs::path traj_xyz = traj_dir;
  traj_xyz /=  fs::path(pathname.c_str());
  if (iparams.pos_xyz) {
    cout << traj_xyz.string() << endl;
    if (fs::create_directory(traj_xyz))
    {
      cout << "Directory created: " << traj_xyz.string() << endl;
    }
    // remove all existing configurations
    for (fs::directory_iterator end_dir_it, it(traj_xyz); it!=end_dir_it; ++it) {
      fs::remove_all(it->path());
    }
  }

  // trajname
  trajname = "state";

  // thermo
	convert.clear();
	convert.str("");
	convert << "thermo.dat";
  pathtooutput = convert.str();
	fthermo.open(pathtooutput.c_str());
	fthermo << left;
  cout << "output file: " << pathtooutput << endl;

  /* initialize simulation */
  if (pathtoinitconf != "") {
    fs::path path(pathtoinitconf.c_str());
    string ext;

    ext = path.extension().string();
    if (ext == ".dat") {
      cout << "Init from DAT format..." << endl;
      cout << "Initializing positions from input" << endl;
      cout << "Initializing velocities from input" << endl;
      world->load_dat(pathtoinitconf);
    }
    else if (ext == ".xyz") {
      cout << "Init from XYZ format..." << endl;
      cout << "Initializing positions from input" << endl;
      world->load_xyz(pathtoinitconf);
    }
    else {
      cout << path.string() << endl;
      throw invalid_argument("Extension not recognized.");
    }
  }

  world->set_constraints();
  world->update_energy_kinetics();
  stepper->call_forces(); // compute forces at current position

  /* integration */
  for (iter = 0; iter < iparams.itermax; ++iter){
    /* dump */
    if (iter % iparams.idump_thermo == 0){
      // thermo dump
      fthermo << setw(10) << fixed << setprecision(0) << noshowpos << iter;
      world->dump_thermo(fthermo);
      fthermo << endl;
    }
    if (iter % iparams.idump_pos == 0){
      // generate file name
      convert.clear();
      convert.str("");
      convert << trajname << setw(iparams.iterwidth) << setfill('0') << setprecision(0) << fixed << dec << iter;
      // conf dump -- standard
      if (iparams.pos_dat) {
        state_path = traj_dat;
        state_path /=  convert.str();
        state_path += ".dat";
        world->dump_dat(state_path.string());
      }
      // conf dump -- xyz format
      if (iparams.pos_xyz) {
        state_path = traj_xyz;
        state_path /=  convert.str();
        state_path += ".xyz";
        world->dump_xyz(state_path.string(), iter);
      }
    }

    /* neighbor list */
    if (iter % iparams.ineighbor_build == 0) {
      world->build_neighbors();
      // world->dump_neighbor(cout);
    }

    /* step */
    // cout << "iter = " << iter << endl;
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
  convert << trajname << setw(iparams.iterwidth) << setfill('0') << setprecision(0) << fixed << dec << iter;
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
  world->dump_pos(fpos_final, true, true, true, false);

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

