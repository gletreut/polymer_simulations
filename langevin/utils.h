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

// other libraries
#include "linalg.h"

// utils
namespace utils {
  void print_vector(std::ostream &mystream, gsl_vector *v);
};

#endif
