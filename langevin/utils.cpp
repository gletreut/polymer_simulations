//----------------------------------------------------------------------------
//  2019-03-29
//  G. Le Treut - Jun Lab, UC San Diego.
//
//  file: utils.cpp
//----------------------------------------------------------------------------
#include <cmath>
#include "utils.h"
using namespace std;

//****************************************************************************
//* utils namespace
//****************************************************************************
void utils::print_vector(ostream &mystream, gsl_vector *v){
  mystream << left << dec << fixed;

  for (size_t i=0; i<v->size; ++i){
    mystream << noshowpos << setprecision(0) << setw(10) << i;
    mystream << showpos << setprecision(8) << setw(18) << gsl_vector_get(v,i);
    mystream << endl;
  }

  return;
}
