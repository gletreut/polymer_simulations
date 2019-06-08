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

// 3) yaml-cpp
#include <yaml-cpp/yaml.h>

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
std::string yamltype_str(YAML::NodeType::value type);
void yaml_init_world(YAML::Node config, gsl_rng *rng, MDWorld* &world);
void yaml_init_rng(YAML::Node config, gsl_rng *rng);
void yaml_init_stepper(YAML::Node config, gsl_rng *rng, MDWorld *world, MDStepper* &stepper);
void yaml_init_forcefields(YAML::Node config, MDWorld *world);

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
  ForceField *ffield(0);
  Constraint *constraint(0);
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
			convert >> pathtoinitconf;
		}
    else {
      pathtoinitconf = "";
    }

  /* load parameters */
  config = YAML::Load(cin);  // load yaml from standart input
  cout << config << endl;

  /* initialize random number generator */
  rng = gsl_rng_alloc(RDT);
  yaml_init_rng(config, rng);

  /* initialize MDWorld */
  world = new MDWorld();
  yaml_init_world(config, rng, world);
  world->print_infos(cout);

  /* initialize MDStepper */
  yaml_init_stepper(config, rng, world, stepper);
  stepper->print_infos(cout);

  /* load force fields */
  yaml_init_forcefields(config, world);
  return 0;

	cout << "Enter parameters values in the form: <key> <value> for the following keys: (iter,dt,T,N,b,lp,dump)" << endl;
	load_map<double>(cin,params);
	cout << "You entered the following parameters:" << endl;
	print_map<double>(cout,params);

  check_params_keys(params, required_parameters);

  // I/O
  // directories
  // // root directory for trajectories
  pathname = "trajectory";
  fs::path traj_dir(pathname.c_str());
  if (fs::create_directory(traj_dir))
  {
    cout << "Directory created: " << traj_dir.string() << endl;
  }
  // // directory for dat format
  pathname = "dat";
  fs::path traj_dat = traj_dir;
  traj_dat /=  fs::path(pathname.c_str());
  cout << traj_dat.string() << endl;
  if (fs::create_directory(traj_dat))
  {
    cout << "Directory created: " << traj_dat.string() << endl;
  }
  // // remove all existing configurations
  for (fs::directory_iterator end_dir_it, it(traj_dat); it!=end_dir_it; ++it) {
    fs::remove_all(it->path());
  }
  // // directory for xyz format
  pathname = "xyz";
  fs::path traj_xyz = traj_dir;
  traj_xyz /=  fs::path(pathname.c_str());
  cout << traj_xyz.string() << endl;
  if (fs::create_directory(traj_xyz))
  {
    cout << "Directory created: " << traj_xyz.string() << endl;
  }
  // // remove all existing configurations
  for (fs::directory_iterator end_dir_it, it(traj_xyz); it!=end_dir_it; ++it) {
    fs::remove_all(it->path());
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
  iterwidth = size_t(log10(itermax)) + 1;

	T=params["T"];                            // temperature
	b=params["b"];                            // bond length
	lp=params["lp"];                          // persistence length
	N=params["N"];                            // number of monomers

  /* initialize force fields */
  /** confinement **/
  //ffield = new ConfinmentBox(-0.5*world->m_lx, +0.5*world->m_lx, -0.5*world->m_ly,+0.5*world->m_ly,-0.5*world->m_lz,+0.5*world->m_lz,1.0,1.0);
  //ffield = new ConfinmentSphere(params["radius_conf"],1.0,1.0);
  //world->m_ffields.push_back(ffield);
  /** polymer **/
  ffield = new PolymerGaussian(0,world->m_npart,1.0);
  world->m_ffields.push_back(ffield);
  //ffield = new PolymerHarmonic(0,world->m_npart,200.0,1.);
  //world->m_ffields.push_back(ffield);
  //ffield = new PolymerFENE(0,world->m_npart, 30.0, 1.5, 1.0, 1.0);
  //world->m_ffields.push_back(ffield);
  //ffield = new PolymerKratkyPorod(0,world->m_npart,3.0);
  //world->m_ffields.push_back(ffield);
  /** GEM **/
  ffield = new GEMField(0, world->m_npart, 1., "kmat.in");
  world->m_ffields.push_back(ffield);

  /* initialize constraints */
  /** attach monomer 0 to origin **/
  //constraint = new PointConstraint(0,0.0,0.0,0.0);
  //world->m_constraints.push_back(constraint);
  /** constrain axis of monomers (0,1) **/
  //constraint = new AxisConstraint(1,0.0,0.0,0.0,1.0,0.0,0.0);
  //world->m_constraints.push_back(constraint);
  /** constrain plane of monomers (0,1,2) **/
  //constraint = new PlaneConstraint(2,0.0,0.0,0.0,0.0,0.0,1.0);
  //world->m_constraints.push_back(constraint);

  /* initialize simulation */
  if (pathtoinitconf == "") {
    cout << "Initializing positions on a lattice" << endl;
    world->init_positions_lattice(1.123);
    cout << "Initializing velocities randomly" << endl;
    world->init_velocities(rng);
  }
  else {
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
      cout << "Initializing velocities randomly" << endl;
      world->init_velocities(rng);
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

std::string yamltype_str(YAML::NodeType::value type){
  if (type == YAML::NodeType::Undefined) {
    return "Undefined";
  }
  else if (type == YAML::NodeType::Null) {
    return "Null";
  }
  else if (type == YAML::NodeType::Scalar) {
    return "Scalar";
  }
  else if (type == YAML::NodeType::Sequence) {
    return "Sequence";
  }
  else if (type == YAML::NodeType::Map) {
    return "Map";
  }
  else {
    return "Error";
  }
}

void yaml_init_world(YAML::Node config, gsl_rng *rng, MDWorld* &world) {
  /*
   * Re-initialize the MDWorld object based on the input parameters.
   */

  string rootkey;
  YAML::Node lineup;

  rootkey = "MDWorld";

  // check that rootkey is in the config
  if (! config[rootkey]) {
    cout << config;
    cout << "root key: " << rootkey << endl;
    throw invalid_argument("Missing the root key!");
  }
  lineup = config[rootkey];

  // initialize new world based on the number of particles
  world->clear();
  world->m_npart = lineup["npart"].as<size_t>();
  world->init();

  // box dimensions
  YAML::Node box = lineup["box"];
  if (! box ) {
    throw invalid_argument("Missing key: box");
  }
  if (! (box.Type() == YAML::NodeType::Sequence) ){
    cout << box;
    throw invalid_argument("<box> must be a sequence.");
  }
  world->m_lx = box[0].as<double>();
  world->m_ly = box[1].as<double>();
  world->m_lz = box[2].as<double>();

  // thermodynamics parameters
  if (lineup["gamma"]){
    world->m_gamma = lineup["gamma"].as<double>();
  }
  if (lineup["temp"]){
    world->m_temp = lineup["temp"].as<double>();
  }
  if (lineup["mass"]){
    world->m_mass = lineup["mass"].as<double>();
  }

  // positions initialization
  YAML::Node init_pos = lineup["init_pos"];
  if (! init_pos){
    throw invalid_argument("Missing key: init_pos");
  }
  if (! (init_pos.Type() == YAML::NodeType::Map) ){
    cout << init_pos;
    throw invalid_argument("<init_pos> must be dictionary.");
  }
  YAML::Node mode, dist;
  mode = init_pos["mode"];
  if (! mode){
    throw invalid_argument("<mode> key missing.");
  }
  dist = init_pos["dist"];
  if (! dist){
    throw invalid_argument("<dist> key missing.");
  }
  string modestr = mode.as<string>();

  if (modestr == "lattice") {
        cout << "Initializing positions on a lattice" << endl;
        world->init_positions_lattice(dist.as<double>());
  }
  else {
    cout << "mode: " << mode.as<string>() << endl;
    throw invalid_argument("Unrecognized <mode>");
  }

  // velocities initialization
  cout << "Initializing velocities randomly" << endl;
  world->init_velocities(rng);

  return;
}

void yaml_init_rng(YAML::Node config, gsl_rng *rng) {
  /*
   * initialize random number generator
   */

  double seed;

  cout << left << dec << fixed;
  if (! config["seed"]) {
    cout << config;
    cout << "<seed> not found --> default value." << endl;
    seed = 123;
  }
  else {
    seed = config["seed"].as<size_t>();
  }
  cout << setw(20) << "seed" << setw(20) << setprecision(0) << noshowpos << seed << endl;
  gsl_rng_set(rng, seed);

  return;
}

void yaml_init_stepper(YAML::Node config, gsl_rng *rng, MDWorld *world, MDStepper* &stepper) {
  /*
   * initialize a stepper
   */
  double dt;
  string rootkey;
  string method;
  YAML::Node lineup;

  rootkey = "stepper";
  // check that rootkey is in the config
  if (! config[rootkey]) {
    cout << config;
    cout << "root key: " << rootkey << endl;
    throw invalid_argument("Missing the root key!");
  }
  lineup = config[rootkey];

  // integration timestep
  YAML::Node dtnode = lineup["dt"];
  if (! dtnode) {
    throw invalid_argument("<dt> key missing.");
  }
  dt = dtnode.as<double>();

  //method
  YAML::Node mnode = lineup["method"];
  if (! mnode) {
    throw invalid_argument("<method> key missing.");
  }
  method = mnode.as<string>();

  // instantiate the stepper
  if (method == "MD_VVerlet") {
    cout << "Initializing method: MD_VVerlet" << endl;
    stepper = new MDStepper_VVerlet(world, dt, world->m_gamma);
  }
  else if (method == "Lang_VVerlet") {
    cout << "Initializing method: Lang_VVerlet" << endl;
    stepper = new LangStepper_VVerlet(world, rng, dt, world->m_gamma);
  }
  else {
    throw invalid_argument("<method> must be one of the following: MD_VVerlet, Lang_VVerlet.");
  }

  return;
}

void yaml_init_forcefields(YAML::Node config, MDWorld *world) {
  /*
   * Initialize forcefields
   */

  string rootkey;
  YAML::Node lineup;
  ForceField *ffield(0);

  rootkey = "forcefields";

  // check that rootkey is in the config
  if (! config[rootkey]) {
    cout << config;
    cout << "root key: " << rootkey << endl;
    throw invalid_argument("Missing the root key!");
  }
  lineup = config[rootkey];

  // check that we have a map
  if ( lineup.Type() == YAML::NodeType::Null){ // if null just return (no forcefields)
    cout << "No force field detected" << endl;
    return;
  }
  else if ( lineup.Type() != YAML::NodeType::Map){ // error
    throw invalid_argument("<forcefield> must be a map");
  }

  // iterate on nodes and add forcefields
  for (YAML::const_iterator it=lineup.begin(); it!=lineup.end(); ++it){
    string key;
    YAML::Node FNode;
    key = it->first.as<string>();
    FNode = it->second;
    if (key == "ConfinmentBox") {
      double sigma, epsilon;
      cout << "Adding force field of type: " << key << endl;
      cout << "Parameters:" << endl;
      cout << FNode << endl;
      sigma = FNode["sigma"].as<double>();
      epsilon = FNode["epsilon"].as<double>();
      ffield = new ConfinmentBox(-0.5*world->m_lx, +0.5*world->m_lx, -0.5*world->m_ly,+0.5*world->m_ly,-0.5*world->m_lz,+0.5*world->m_lz,sigma,epsilon);
      world->m_ffields.push_back(ffield);
    }
    else if (key == "ConfinmentSphere") {
      double radius, sigma, epsilon;
      cout << "Adding force field of type: " << key << endl;
      cout << "Parameters:" << endl;
      cout << FNode << endl;
      radius = FNode["radius"].as<double>();
      sigma = FNode["sigma"].as<double>();
      epsilon = FNode["epsilon"].as<double>();
      ffield = new ConfinmentSphere(radius, sigma, epsilon);
      world->m_ffields.push_back(ffield);
    }
    else if (key == "PolymerGaussian") {
      size_t offset, N;
      double b;
      cout << "Adding force field of type: " << key << endl;
      cout << "Parameters:" << endl;
      cout << FNode << endl;
      offset = FNode["offset"].as<size_t>();
      N = FNode["N"].as<size_t>();
      b = FNode["b"].as<double>();
      ffield = new PolymerGaussian(offset, N, b);
      world->m_ffields.push_back(ffield);
    }
    else if (key == "PolymerHarmonic") {
      size_t offset, N;
      double ke, r0;
      cout << "Adding force field of type: " << key << endl;
      cout << "Parameters:" << endl;
      cout << FNode << endl;
      offset = FNode["offset"].as<size_t>();
      N = FNode["N"].as<size_t>();
      ke = FNode["ke"].as<double>();
      r0 = FNode["r0"].as<double>();
      ffield = new PolymerHarmonic(offset, N, ke, r0);
      world->m_ffields.push_back(ffield);
    }
  else {
    cout << "Unrecognized force field type: " << key << endl;
  }

  }
  /* initialize force fields */
  /** confinement **/
  //ffield = new ConfinmentBox(-0.5*world->m_lx, +0.5*world->m_lx, -0.5*world->m_ly,+0.5*world->m_ly,-0.5*world->m_lz,+0.5*world->m_lz,1.0,1.0);
  //ffield = new ConfinmentSphere(params["radius_conf"],1.0,1.0);
  //world->m_ffields.push_back(ffield);
  /** polymer **/
  //ffield = new PolymerGaussian(0,world->m_npart,1.0);
  //world->m_ffields.push_back(ffield);
  //ffield = new PolymerHarmonic(0,world->m_npart,200.0,1.);
  //world->m_ffields.push_back(ffield);
  //ffield = new PolymerFENE(0,world->m_npart, 30.0, 1.5, 1.0, 1.0);
  //world->m_ffields.push_back(ffield);
  //ffield = new PolymerKratkyPorod(0,world->m_npart,3.0);
  //world->m_ffields.push_back(ffield);
  /** GEM **/
  //ffield = new GEMField(0, world->m_npart, 1., "kmat.in");
  //world->m_ffields.push_back(ffield);

  return;
}
