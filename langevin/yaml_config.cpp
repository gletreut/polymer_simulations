//----------------------------------------------------------------------------
//  2019-03-29
//  G. Le Treut - Jun Lab, UC San Diego.
//
//  file: yaml_config.cpp
//----------------------------------------------------------------------------

#include "yaml_config.h"
using namespace std;

//****************************************************************************
//* yaml_config namespace
//****************************************************************************
string yamltype_str(YAML::NodeType::value type){
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

void yaml_config::init_world(YAML::Node config, gsl_rng *rng, MDWorld* &world) {
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
  YAML::Node npart = lineup["npart"];
  if (! npart ) {
    throw invalid_argument("Missing key: npart");
  }
  world->clear();
  world->m_npart = npart.as<size_t>();
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

void yaml_config::init_rng(YAML::Node config, gsl_rng *rng) {
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

void yaml_config::init_stepper(YAML::Node config, gsl_rng *rng, MDWorld *world, MDStepper* &stepper) {
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

void yaml_config::init_forcefields(YAML::Node config, MDWorld* &world) {
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
    else if (key == "PolymerFENE") {
      size_t offset, N;
      double ke, rc, sigma, epsilon;
      cout << "Adding force field of type: " << key << endl;
      cout << "Parameters:" << endl;
      cout << FNode << endl;
      offset = FNode["offset"].as<size_t>();
      N = FNode["N"].as<size_t>();
      ke = FNode["ke"].as<double>();
      rc = FNode["rc"].as<double>();
      sigma = FNode["sigma"].as<double>();
      epsilon = FNode["epsilon"].as<double>();
      ffield = new PolymerFENE(offset, N, ke, rc, sigma, epsilon);
      world->m_ffields.push_back(ffield);
    }
    else if (key == "PolymerKratkyPorod") {
      size_t offset, N;
      double lp;
      cout << "Adding force field of type: " << key << endl;
      cout << "Parameters:" << endl;
      cout << FNode << endl;
      offset = FNode["offset"].as<size_t>();
      N = FNode["N"].as<size_t>();
      lp = FNode["lp"].as<double>();
      ffield = new PolymerKratkyPorod(offset, N, lp);
      world->m_ffields.push_back(ffield);
    }
    else if (key == "GEMField") {
      size_t offset, N;
      double b;
      string kmat_fpath;
      cout << "Adding force field of type: " << key << endl;
      cout << "Parameters:" << endl;
      cout << FNode << endl;
      offset = FNode["offset"].as<size_t>();
      N = FNode["N"].as<size_t>();
      b = FNode["b"].as<double>();
      kmat_fpath = FNode["kmat_fpath"].as<string>();
      ffield = new GEMField(offset, N, b, kmat_fpath);
      world->m_ffields.push_back(ffield);
    }
    else {
      cout << "Unrecognized force field type: " << key << endl;
    }

  }
  return;
}

void yaml_config::init_constraints(YAML::Node config, MDWorld* &world) {
  /*
   * Initialize constraints
   */

  string rootkey;
  YAML::Node lineup;
  Constraint *constraint(0);

  rootkey = "constraints";

  // check that rootkey is in the config
  if (! config[rootkey]) {
    cout << config;
    cout << "root key: " << rootkey << endl;
    throw invalid_argument("Missing the root key!");
  }
  lineup = config[rootkey];

  // check that we have a map
  if ( lineup.Type() == YAML::NodeType::Null){ // if null just return (no forcefields)
    cout << "No constraint detected" << endl;
    return;
  }
  else if ( lineup.Type() != YAML::NodeType::Map){ // error
    throw invalid_argument("<constraints> must be a map");
  }

  // iterate on nodes and add forcefields
  for (YAML::const_iterator it=lineup.begin(); it!=lineup.end(); ++it){
    string key;
    YAML::Node FNode;
    key = it->first.as<string>();
    FNode = it->second;
    if (key == "PointConstraint") {
      size_t index;
      double rx, ry, rz;
      cout << "Adding constraint of type: " << key << endl;
      cout << "Parameters:" << endl;
      cout << FNode << endl;
      index = FNode["index"].as<size_t>();
      rx = FNode["point"][0].as<double>();
      ry = FNode["point"][1].as<double>();
      rz = FNode["point"][2].as<double>();
      constraint = new PointConstraint(index, rx, ry, rz);
      world->m_constraints.push_back(constraint);
    }
    else if (key == "AxisConstraint") {
      size_t index;
      double rx,ry,rz,ux,uy,uz;
      cout << "Adding constraint of type: " << key << endl;
      cout << "Parameters:" << endl;
      cout << FNode << endl;
      index = FNode["index"].as<size_t>();
      rx = FNode["point"][0].as<double>();
      ry = FNode["point"][1].as<double>();
      rz = FNode["point"][2].as<double>();
      ux = FNode["axis"][0].as<double>();
      uy = FNode["axis"][1].as<double>();
      uz = FNode["axis"][2].as<double>();
      constraint = new AxisConstraint(index, rx, ry, rz, ux, uy, uz);
      world->m_constraints.push_back(constraint);
    }
    else if (key == "PlaneConstraint") {
      size_t index;
      double rx,ry,rz,ux,uy,uz;
      cout << "Adding constraint of type: " << key << endl;
      cout << "Parameters:" << endl;
      cout << FNode << endl;
      index = FNode["index"].as<size_t>();
      rx = FNode["point"][0].as<double>();
      ry = FNode["point"][1].as<double>();
      rz = FNode["point"][2].as<double>();
      ux = FNode["direction"][0].as<double>();
      uy = FNode["direction"][1].as<double>();
      uz = FNode["direction"][2].as<double>();
      constraint = new PlaneConstraint(index, rx, ry, rz, ux, uy, uz);
      world->m_constraints.push_back(constraint);
    }
    else {
      cout << "Unrecognized constraint type: " << key << endl;
    }

  }
  return;
}

void yaml_config::init_integration_params(YAML::Node config, IntegrationParams &iparams) {
  /*
   * Initialize integration parameters.
   */

  string rootkey;
  YAML::Node lineup;
  Constraint *constraint(0);

  rootkey = "integration";

  // check that rootkey is in the config
  if (! config[rootkey]) {
    cout << config;
    cout << "root key: " << rootkey << endl;
    throw invalid_argument("Missing the root key!");
  }
  lineup = config[rootkey];

  // check that we have a map
  if ( lineup.Type() == YAML::NodeType::Null){ // if null just return
    cout << "Must provide the integration parameters" << endl;
    return;
  }
  else if ( lineup.Type() != YAML::NodeType::Map){ // error
    throw invalid_argument("<integration> must be a map");
  }

  // itermax
  YAML::Node itermax = lineup["itermax"];
  if (! itermax ) {
    throw invalid_argument("Missing key: itermax");
  }
  iparams.itermax = size_t(itermax.as<double>());

  // idump_thermo
  YAML::Node idump_thermo = lineup["thermo"]["idump"];
  if (! idump_thermo ) {
    throw invalid_argument("Missing key: <thermo><idump>");
  }
  iparams.idump_thermo = size_t(idump_thermo.as<double>());

  // idump_pos
  YAML::Node idump_pos = lineup["positions"]["idump"];
  if (! idump_pos ) {
    throw invalid_argument("Missing key: <positions><idump>");
  }
  iparams.idump_pos = size_t(idump_pos.as<double>());

  // idump xyz
  if (lineup["positions"]["xyz"]){
    iparams.pos_xyz = lineup["positions"]["xyz"].as<bool>();
  }

  // idump dat
  if (lineup["positions"]["dat"]){
    iparams.pos_dat = lineup["positions"]["dat"].as<bool>();
  }

  // update computed parameters
  iparams.init();

  return;
}
