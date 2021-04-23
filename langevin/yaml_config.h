#ifndef _YAML_CONFIG_H
#define _YAML_CONFIG_H
//----------------------------------------------------------------------------
//  2019-06-22
//  G. Le Treut - Jun Lab, UC San Diego.
//
//  file: utils.h
//----------------------------------------------------------------------------

#include <cstdlib>
#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <map>

// external library
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <yaml-cpp/yaml.h>

// local library
#include "model.h"
#include "stepper.h"
#include "integration.h"

// yaml_config
namespace yaml_config {
  std::string yamltype_str(YAML::NodeType::value type);
  void init_world(YAML::Node config, gsl_rng *rng, MDWorld* &world);
  void init_rng(YAML::Node config, gsl_rng *rng);
  void init_stepper(YAML::Node config, gsl_rng *rng, MDWorld *world, MDStepper* &stepper);
  void init_forcefields(YAML::Node config, MDWorld* &world);
  void init_constraints(YAML::Node config, MDWorld* &world);
  void init_integration_params(YAML::Node config, IntegrationParams &iparams);
};

#endif
