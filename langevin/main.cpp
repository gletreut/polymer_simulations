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
#include "neighborlist.h"
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

//** MAIN
int main(int argc, char *argv[]){
	string pathtooutput, pathtoinitconf, pathname, trajname;
	stringstream convert;
  fs::path traj_dat, traj_xyz, traj_neighbor, traj_polarity, state_path;
	ofstream fthermo;
  YAML::Node config;
  gsl_rng *rng(0);

  MDWorld *world(0);
  MDStepper *stepper(0);
  IntegrationParams iparams;
  bool nlist_updated=false;
  size_t iter=0;

  /* test: start
  ofstream ftest;
  ftest.open("test.dat", ofstream::out);
  ftest << "";
  ftest.close();
  // test: end */

//-----------------------------------------------
// INITIALIZATION
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
  traj_dat = traj_dir;
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
  traj_xyz = traj_dir;
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

  // directory for neighbor lists
  pathname = "neighbor";
  traj_neighbor = traj_dir;
  traj_neighbor /=  fs::path(pathname.c_str());
  if (iparams.neighbor) {
    cout << traj_neighbor.string() << endl;
    if (fs::create_directory(traj_neighbor))
    {
      cout << "Directory created: " << traj_neighbor.string() << endl;
    }
    // remove all existing configurations
    for (fs::directory_iterator end_dir_it, it(traj_neighbor); it!=end_dir_it; ++it) {
      fs::remove_all(it->path());
    }
  }

  // directory for polarity vector format
  pathname = "polarity";
  traj_polarity = traj_dir;
  traj_polarity /=  fs::path(pathname.c_str());
  cout << traj_polarity.string() << endl;
  if (fs::create_directory(traj_polarity))
  {
    cout << "Directory created: " << traj_polarity.string() << endl;
  }
  // remove all existing configurations
  for (fs::directory_iterator end_dir_it, it(traj_polarity); it!=end_dir_it; ++it) {
    fs::remove_all(it->path());
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

  /* initialize neighbor lists */
  for (vector<NeighborList*>::iterator it=world->m_neighbors.begin(); it!=world->m_neighbors.end(); ++it){
    (*it)->build(world->m_x, world->m_bonds);
  }
  nlist_updated=true;

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
      // polarity vector dump
      state_path = traj_polarity;
      state_path /=  convert.str();
      state_path += ".dat";
      world->dump_polarity_vectors(state_path.string());
    }

    if (nlist_updated) {
      convert.clear();
      convert.str("");
      convert << trajname << setw(iparams.iterwidth) << setfill('0') << setprecision(0) << fixed << dec << iter;
      state_path = traj_neighbor;
      state_path /=  convert.str();
      state_path += ".dat";
      world->dump_neighbors(state_path.string());
    }

    /* step */
    // cout << "iter = " << iter << endl;
    stepper->step();

    /* update energy */
    // potential energy already updated when calling the forces in the stepper
    world->update_energy_kinetics();

    /* update neighbor list */
    // if ( world->check_neighbor_lists() || (iter % iparams.ineighbor_build == 0) ) {
    // the check method doesn't work well so far. To make it work I would need to distinguish
    // rskin from rcut. rskin could be a fraction of rcut. Then yes, re-build when one atom
    // has moved more than the skin distance.
    if ( iter % iparams.ineighbor_build == 0 ) {
      world->update_neighbor_lists(iter);
      nlist_updated = true;
    }
    else
      nlist_updated = false;
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
  // polarity dump -- xyz format
  state_path = traj_polarity;
  state_path /=  convert.str();
  state_path += ".dat";
  world->dump_polarity_vectors(state_path.string());

  // thermo
  ofstream fthermo_final("thermo_final.dat");
  fthermo_final << left;
  world->dump_thermo(fthermo_final);
  fthermo_final << endl;

  // positions
  ofstream fpos_final("positions_final.dat");
  fpos_final << left << dec << fixed;
  world->dump_pos(fpos_final, true, true, true, false);

  /* print statistics */
  cout << "STATISTICS" << endl;
  for (vector<NeighborList*>::iterator it=world->m_neighbors.begin(); it!=world->m_neighbors.end(); ++it){
    cout << setw(40) << "number of builds:";
    cout << setw(20) << setprecision(0) << noshowpos << (*it)->m_ibuild_count;
    cout << endl;
    cout << setw(40) << "average build interval:";
    cout << setw(20) << setprecision(1) << noshowpos << double((*it)->m_ibuild_sum)/((*it)->m_ibuild_count);
    cout << endl;
    cout << endl;
  }

//-----------------------------------------------
//	EXIT
//-----------------------------------------------
  delete stepper;
  delete world;
  gsl_rng_free(rng);

	cout << "normal exit" << endl;
	return 0;
}
