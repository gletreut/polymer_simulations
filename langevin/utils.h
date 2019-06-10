#ifndef _UTILS_H
#define _UTILS_H
//----------------------------------------------------------------------------
//  2019-03-29
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

// utils
namespace utils {
	std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems);
	std::vector<std::string> split(const std::string &s, char delim);

  template<typename T>
  void load_matrix(std::istream &mystream, std::vector<std::vector<T> > &cmap);

  template<typename T>
  void print_matrix(std::ostream &mystream, std::vector<std::vector<T> > &cmap);

  template<typename T>
  void load_map(std::istream &mystream, std::map<std::string,T> &params);

  template<typename T>
  void print_map(std::ostream &mystream, std::map<std::string,T> &params);

  void print_vector(std::ostream &mystream, gsl_vector *v);
  void print_matrix(std::ostream &mystream, gsl_matrix *m);
  void load_matrix(std::istream &mystream, gsl_matrix *m);
};

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

// include template definitions
#include "utils.tpp"

#endif
