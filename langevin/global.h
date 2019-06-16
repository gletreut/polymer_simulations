#ifndef GLOBAL_H
#define GLOBAL_H

#include <cstdlib>
#include <fstream>

#include <gsl/gsl_rng.h>

/* random number generator type */
extern const gsl_rng_type*  RDT;
extern const size_t seed;

/* test file streams */
extern std::ofstream fforce;
#endif
