//*******************************************************************************
//*
//* Langevin Dynamics
//*	forcefield.h
//*
//* Author: Guillaume Le Treut
//*	UCSD - 2019
//*
//*******************************************************************************

#ifndef _FORCEFIELD_H
#define _FORCEFIELD_H

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <vector>
#include <limits>
#include <stdexcept>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include "linalg.h"

class ForceField {
  /*
   * Class defining a force field.
   */
  virtual void energy_force(gsl_matrix *x, double *u, gsl_matrix *force) = 0;
};

class PolymerFENE : ForceField {
  /*
   * Class defining a polymer force field.
   */

  public:
    /* attributes */
    size_t m_offset;    // index of first monomer
    size_t m_N;         // number of monomers
    double m_ke;        // energy / b^2 for FENE bond
    double m_rc;        // cutoff for FENE bond
    double m_sigma;     // hard-core distance for LJ
    double m_eps;       // energy scale for LJ

    /* constructor and destructor */
    PolymerFENE(size_t offset, size_t N, double ke=30.0, double rc=1.5, double sigma=1.0, double eps=1.0);

    /* methods */
    void energy_force(gsl_matrix *x, double *u, gsl_matrix *force);

  private:
    double m_pref;
    double m_4eps;
    double m_48eps;
    double m_rc_LJ;
};


#endif
