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
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <limits>
#include <stdexcept>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>

#include "neighborlist.h"
#include "linalg.h"
#include "utils.h"

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

class PolymerHarmonic : public ForceField {
  /*
   * Class defining a polymer force field.
   * Harmonic bonds backbone
   */

  public:
    /* attributes */
    size_t m_offset;    // index of first monomer
    size_t m_N;         // number of monomers
    double m_ke;        // energy / b^2 for FENE bond
    double m_r0;

    /* constructor and destructor */
    PolymerHarmonic(size_t offset, size_t N, double ke=100, double r0=1.);
    ~PolymerHarmonic();

    /* methods */
    void energy_force(gsl_matrix *x, double *u, gsl_matrix *forces);

  private:
    double m_fpref;
};

class PolymerFENE : public ForceField {
  /*
   * Class defining a polymer force field.
   * See LAMMPS documentation: https://lammps.sandia.gov/doc/bond_fene.html
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
    double energy_LJ_scal(double r);
    double force_LJ_scal(double r);
    void energy_force(gsl_matrix *x, double *u, gsl_matrix *forces);

  private:
    double m_pref;
    double m_4eps;
    double m_fpref_LJ;
    double m_rc_LJ;
};

class PolymerKratkyPorod : public ForceField {
  /*
   * Class defining a force field representing a polymer semi-flexibility
   * Kratky-Porod model
   */

  public:
    /* attributes */
    size_t m_offset;    // index of first monomer
    size_t m_N;         // number of monomers
    double m_lp;        // persistence length (in units of bond length)

    /* constructor and destructor */
    PolymerKratkyPorod(size_t offset, size_t N, double lp);
    ~PolymerKratkyPorod();

    /* methods */
    void energy_force(gsl_matrix *x, double *u, gsl_matrix *forces);
};

class GEMField : public ForceField {
  /*
   * Class defining a Gaussian Effective Model force field.
   * It consists of elastic interactions to be overlaid on a Polymer force
   * field. It is parametrized by a matrix W such that the potential is:
   *   \beta U(X, Y, Z) = 3/(2b^2) [ X^T W X + Y^T W Y + Z^T W Z],
   *   where X, Y and Z are vectors of size N.
   */

  public:
    /* attributes */
    size_t m_offset;    // index of first monomer
    size_t m_N;         // number of monomers
    double m_b;
    gsl_matrix *m_W;    // interaction potential matrix.

    /* constructor and destructor */
    GEMField(size_t offset, size_t N, double b, std::string filepath);
    ~GEMField();

    /* methods */
    void K2W(const gsl_matrix *K, gsl_matrix *W);
    void energy_force(gsl_matrix *x, double *u, gsl_matrix *forces);

  private:
    double m_fpref;
    gsl_vector *m_fa;
};

class SoftCore : public ForceField {
  /*
   * Class defining a force field for soft-core excluded volume.
   *
   */

  public:
    /* attributes */
    double m_A;                     // scale
    double m_rc;                    // cutoff of interaction
    double m_y;                     // hold the ratio pi / rc
    NeighborList *m_neighbors;      // neighbor list

    /* constructor and destructor */
    SoftCore(double A, double sigma, NeighborList* neighbors);
    ~SoftCore();

    /* methods */
    double energy_scal(double r);
    double force_scal(double r);
    void energy_force(gsl_matrix *x, double *u, gsl_matrix *forces);
};

class PairLJ : public ForceField {
  /*
   * Class defining a force field for Lennard-Jones pair potential
   *
   */

  public:
    /* attributes */
    double m_eps;                 // scale
    double m_sigma;               // hard-core distance
    double m_rc_LJ;               // cutoff
    NeighborList *m_neighbors;    // neighbor list

    /* constructor and destructor */
    PairLJ(double eps, double sigma, double rc_LJ, NeighborList* neighbors);
    ~PairLJ();

    /* methods */
    double energy_LJ_scal(double r);
    double force_LJ_scal(double r);
    void energy_force(gsl_matrix *x, double *u, gsl_matrix *forces);

  private:
    double m_4eps;
    double m_fpref;
    double m_u0;
};

#endif
