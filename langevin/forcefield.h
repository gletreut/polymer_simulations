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
#include <gsl/gsl_spmatrix.h>
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
    gsl_spmatrix_uint *m_bonds; // bond selector

    /* constructor and destructor */
    PolymerGaussian(size_t offset, size_t N, double b, gsl_spmatrix_uint *bonds);
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
    gsl_spmatrix_uint *m_bonds; // bond selector

    /* constructor and destructor */
    // PolymerHarmonic(size_t offset, size_t N, double ke=100, double r0=1.);
    PolymerHarmonic(size_t offset, size_t N, double ke, double r0, gsl_spmatrix_uint *bonds);
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
    gsl_spmatrix_uint *m_bonds; // bond selector

    /* constructor and destructor */
    // PolymerFENE(size_t offset, size_t N, double ke=30.0, double rc=1.5, double sigma=1.0, double eps=1.0);
    PolymerFENE(size_t offset, size_t N, double ke, double rc, double sigma, double eps, gsl_spmatrix_uint *bonds);
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
    gsl_spmatrix_uint *m_bonds; // bond selector

    /* constructor and destructor */
    PolymerKratkyPorod(size_t offset, size_t N, double lp, gsl_spmatrix_uint *bonds);
    ~PolymerKratkyPorod();

    /* methods */
    void energy_force(gsl_matrix *x, double *u, gsl_matrix *forces);
};



// GEMField Class ----------------------------------------------------------------------------------------------------------------------------------------------------------
// would need to deal with bonded interactions for the following potential (not sparse)
// class GEMField : public ForceField {
//   /*
//    * Class defining a Gaussian Effective Model force field.
//    * It consists of elastic interactions to be overlaid on a Polymer force
//    * field. It is parametrized by a matrix W such that the potential is:
//    *   \beta U(X, Y, Z) = 3/(2b^2) [ X^T W X + Y^T W Y + Z^T W Z],
//    *   where X, Y and Z are vectors of size N.
//    */
//
//   public:
//     /* attributes */
//     size_t m_offset;    // index of first monomer
//     size_t m_N;         // number of monomers
//     double m_b;
//     gsl_matrix *m_W;    // interaction potential matrix.
//
//     /* constructor and destructor */
//     GEMField(size_t offset, size_t N, double b, std::string filepath);
//     ~GEMField();
//
//     /* methods */
//     void K2W(const gsl_matrix *K, gsl_matrix *W);
//     void energy_force(gsl_matrix *x, double *u, gsl_matrix *forces);
//
//   private:
//     double m_fpref;
//     gsl_vector *m_fa;
// };
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------


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
    double m_4eps;
    double m_48eps;
    double m_fpref;
    double m_u0;

    /* constructor and destructor */
    PairLJ(double eps, double sigma, double rc_LJ, NeighborList* neighbors);
    ~PairLJ();

    /* methods */
    double energy_LJ_scal(double r);
    double force_LJ_scal(double r);
    void energy_force(gsl_matrix *x, double *u, gsl_matrix *forces);

};

class PolarPairLJ : public PairLJ {
  /*
   * Class defining a force field for Lennard-Jones pair potential (Polarity Vector Included)
   *
   */

  public:
    /* attributes */
    std::vector<std::pair<size_t, size_t> > m_chain_ends;
    double m_alpha;
    gsl_matrix *m_pol_vec;

    /* constructor and destructor */
    PolarPairLJ(double eps, double sigma, double rc_LJ, double alpha, NeighborList* neighbors, std::vector<std::pair<size_t, size_t> > chain_ends);
    ~PolarPairLJ();

    /* methods */
    double V_LJ(double rhat);
    double V_LJ_prime(double rhat);
    double get_sigma(double phi, double ti, double tj);
    double get_sigma_dphi(double phi, double ti, double tj);
    double energy_LJ_scal(const gsl_vector *xi, const gsl_vector *xj, const gsl_vector *ni, const gsl_vector *nj);
    void force_LJ_scal(const gsl_vector *xi, const gsl_vector *xj, const gsl_vector *ni, const gsl_vector *nj, gsl_vector *force);
    void energy_force(gsl_matrix *x, double *u, gsl_matrix *forces);

};

class PolarPair48 : public ForceField {
  /*
   * Class defining a force field as in H. Li and G. Lykotrafitis. “Two-component coarse-grained molecular-dynamics model
   * for the human erythrocyte membrane”. In: Biophysical journal 102.1 (2012).
   */

  public:
    /* attributes */
    double m_eps;
    double m_sigma;
    double m_r0;
    double m_rc;
    double m_alpha;
    NeighborList *m_neighbors;
    std::vector<std::pair<size_t, size_t> > m_chain_ends;
    gsl_matrix *m_pol_vec;

    /* constructor and destructor */
    PolarPair48(double eps, double sigma, double rc, double alpha, NeighborList* neighbors, std::vector<std::pair<size_t, size_t> > chain_ends);
    ~PolarPair48();

    /* methods */
    double get_A(double phi, double ti, double tj);
    double get_A_dphi(double phi, double ti, double tj);
    double get_energy(const gsl_vector *xi, const gsl_vector *xj, const gsl_vector *ni, const gsl_vector *nj);
    void get_force(const gsl_vector *xi, const gsl_vector *xj, const gsl_vector *ni, const gsl_vector *nj, gsl_vector *force);
    void energy_force(gsl_matrix *x, double *u, gsl_matrix *forces);

};

class PolarPair48_2site : public ForceField {
  /*
   * Class defining a force field as in H. Li and G. Lykotrafitis. “Two-component coarse-grained molecular-dynamics model
   * for the human erythrocyte membrane”. In: Biophysical journal 102.1 (2012).
   * Object with 2 sites aligned along the polarity vectors.
   */

  public:
    /* attributes */
    double m_eps_nn, m_eps_ww, m_eps_nw;
    double m_sigma_nn, m_sigma_ww, m_sigma_nw;
    double m_r0_nn, m_r0_ww, m_r0_nw;
    double m_rc_nn, m_rc_ww, m_rc_nw;
    double m_rho_n, m_rho_w;
    NeighborList *m_neighbors;
    std::vector<std::pair<size_t, size_t> > m_chain_ends;
    gsl_matrix *m_pol_vec;

    /* constructor and destructor */
    PolarPair48_2site(double eps_nn, double eps_ww, double eps_nw,
                double sigma_nn, double sigma_ww, double sigma_nw,
                double rc_nn, double rc_ww, double rc_nw,
                double rho_n, double rho_w,
                NeighborList* neighbors, std::vector<std::pair<size_t, size_t> > chain_ends);
    ~PolarPair48_2site();

    /* methods */
    double potential(double r, double eps, double rc, double r0);
    double potential_deriv(double r, double eps, double rc, double r0);
    double get_energy(const gsl_vector *xi, const gsl_vector *xj, const gsl_vector *ni, const gsl_vector *nj);
    double get_energy_force(const gsl_vector *xi, const gsl_vector *xj, const gsl_vector *ni, const gsl_vector *nj, gsl_vector *force);
    // void get_force(const gsl_vector *xi, const gsl_vector *xj, const gsl_vector *ni, const gsl_vector *nj, gsl_vector *force);
    void energy_force(gsl_matrix *x, double *u, gsl_matrix *forces);

};
#endif
