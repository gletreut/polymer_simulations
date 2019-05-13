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
#include <gsl/gsl_rng.h>

#include "linalg.h"

class ForceField {
  /*
   * Virtual class defining a force field.
   */
  public:
    virtual ~ForceField();
    virtual void energy_force(gsl_matrix *x, double *u, gsl_matrix *forces) = 0;
};

class ConfinmentBox : public ForceField {
  /*
   * Class defining a force field representing a box confinment
   */
  public:
    double m_xlo, m_xhi, m_ylo, m_yhi, m_zlo, m_zhi;
    double m_sigma;     // hard-core distance for LJ
    double m_eps;       // energy scale for LJ

    /* constructor and destructor */
    ConfinmentBox(double xlo, double xhi, double ylo, double yhi, double zlo, double zhi, double sigma=1.0, double eps=1.0);
    ~ConfinmentBox();

    /* methods */
    double energy_LJ_scal(double r);
    double force_LJ_scal(double r);
    void energy_force(gsl_matrix *x, double *u, gsl_matrix *forces);

  private:
    double m_4eps;
    double m_fpref;
    double m_rc_LJ;
};

class ConfinmentSphere : public ForceField {
  /*
   * Class defining a force field representing a spherical confinment
   */
  public:
    double m_radius;
    double m_sigma;     // hard-core distance for LJ
    double m_eps;       // energy scale for LJ

    /* constructor and destructor */
    ConfinmentSphere(double radius, double sigma=1.0, double eps=1.0);
    ~ConfinmentSphere();

    /* methods */
    double energy_LJ_scal(double r);
    double force_LJ_scal(double r);
    void energy_force(gsl_matrix *x, double *u, gsl_matrix *forces);

  private:
    double m_4eps;
    double m_fpref;
    double m_rc_LJ;
};

class PolymerGaussian : public ForceField {
  /*
   * Class defining a polymer force field.
   * Gaussian backbone
   */

  public:
    /* attributes */
    size_t m_offset;    // index of first monomer
    size_t m_N;         // number of monomers
    double m_b;         // bond length b
    double m_ke;        // energy / b^2 for FENE bond

    /* constructor and destructor */
    PolymerGaussian(size_t offset, size_t N, double b);
    ~PolymerGaussian();

    /* methods */
    void energy_force(gsl_matrix *x, double *u, gsl_matrix *forces);

  private:
    double m_fpref;
};

class PolymerFENE : public ForceField {
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
    ~PolymerFENE();

    /* methods */
    void energy_force(gsl_matrix *x, double *u, gsl_matrix *forces);

  private:
    double m_pref;
    double m_4eps;
    double m_48eps;
    double m_rc_LJ;
};

class LangThermostat : public ForceField {
  /*
   * Class defining a random force field reproducing a a thermal bath.
   * <f> = 0
   * <f^2> = 2 m \gamma k_B T
   */

  public:
    /* attributes */
    double m_mass;         // mass
    double m_gamma;        // gamma
    double m_temp;         // temperature
    gsl_rng *m_rng;        // random number generator

    /* constructor and destructor */
    LangThermostat(double mass, double gamma, double temp, gsl_rng *rng);
    ~LangThermostat();

    /* methods */
    void energy_force(gsl_matrix *x, double *u, gsl_matrix *forces);

  private:
    double m_s12;
    double m_fpref;
};


#endif
