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
#include <algorithm>
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

class ConfinementBox : public ForceField {
  /*
   * Class defining a force field representing a box confinement
   */
  public:
    double m_xlo, m_xhi, m_ylo, m_yhi, m_zlo, m_zhi;
    double m_sigma;     // hard-core distance for LJ
    double m_eps;       // energy scale for LJ

    /* constructor and destructor */
    ConfinementBox(double xlo, double xhi, double ylo, double yhi, double zlo, double zhi, double sigma=1.0, double eps=1.0);
    ~ConfinementBox();

    /* methods */
    double energy_LJ_scal(double r);
    double force_LJ_scal(double r);
    void energy_force(gsl_matrix *x, double *u, gsl_matrix *forces);

  private:
    double m_4eps;
    double m_fpref;
    double m_rc_LJ;
};

class ConfinementSphere : public ForceField {
  /*
   * Class defining a force field representing a spherical confinement
   */
  public:
    double m_radius;
    double m_sigma;     // hard-core distance for LJ
    double m_eps;       // energy scale for LJ

    /* constructor and destructor */
    ConfinementSphere(double radius, double sigma=1.0, double eps=1.0);
    ~ConfinementSphere();

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
    double m_ke;        // spring constant (energy / length^2)
    double m_r0;        // rest length
    gsl_spmatrix_uint *m_bonds; // bond selector

    /* constructor and destructor */
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

class PolarChains : public ForceField {
  /*
   * Base class for pair force fields between polar beads, in 2D.
   *
   * The polarity of bead k is the in-plane normal to its chain:
   *   n_k = ez x t_k / |t_k|,  t_k = x_{k+} - x_{k-},
   * with (k+, k-) = (k+1, k-1) for interior beads, (k+1, k) at the chain start
   * and (k, k-1) at the chain end.
   *
   * The pair energy depends on n_k, hence on x_{k+} and x_{k-}. With
   * g_k = dU/dn_k and J_k = dn_k/dt_k, the force therefore includes
   *   F_{k+} -= J_k^T g_k,  F_{k-} += J_k^T g_k,
   * where, in 2D, J_k^T g_k = -(t_k/|t_k| . g_k) n_k / |t_k|.
   *
   * Every particle must belong to exactly one chain of at least 2 beads.
   */

  public:
    /* attributes */
    NeighborList *m_neighbors;
    std::vector<std::pair<size_t, size_t> > m_chain_ends;
    gsl_matrix *m_pol_vec;        // n_k
    gsl_matrix *m_pol_grad;       // g_k = dU/dn_k
    std::vector<double> m_tnorm;  // |t_k|

    /* constructor and destructor */
    PolarChains(NeighborList *neighbors, std::vector<std::pair<size_t, size_t> > chain_ends);
    virtual ~PolarChains();

    /* methods */
    // Pair energy; sets fj = force on j at fixed polarities (the force on i is -fj),
    // gi = dU/dn_i and gj = dU/dn_j.
    virtual double pair_energy_force(const gsl_vector *xi, const gsl_vector *xj,
                                     const gsl_vector *ni, const gsl_vector *nj,
                                     gsl_vector *fj, gsl_vector *gi, gsl_vector *gj) = 0;
    void energy_force(gsl_matrix *x, double *u, gsl_matrix *forces);
    void update_polarity(const gsl_matrix *x);

  protected:
    void add_polarity_forces(gsl_matrix *forces);
};

class PolarPairLJ : public PolarChains {
  /*
   * Lennard-Jones pair potential with a polarity-dependent diameter:
   *   U = V(r/s) - V(rc),  V(u) = 4 eps (u^-12 - u^-6),
   *   s = sigma (1 + alpha (1 + a)/2),  a = (n_i - n_j).e_ij / 2  in [-1, 1],
   * with e_ij the unit vector from i to j. The cutoff rc is in units of s.
   */

  public:
    /* attributes */
    double m_eps;
    double m_sigma;
    double m_rc_LJ;
    double m_alpha;

    /* constructor and destructor */
    PolarPairLJ(double eps, double sigma, double rc_LJ, double alpha, NeighborList* neighbors, std::vector<std::pair<size_t, size_t> > chain_ends);
    ~PolarPairLJ();

    /* methods */
    double V_LJ(double u);
    double V_LJ_prime(double u);
    double pair_energy_force(const gsl_vector *xi, const gsl_vector *xj,
                             const gsl_vector *ni, const gsl_vector *nj,
                             gsl_vector *fj, gsl_vector *gi, gsl_vector *gj);
};

class PolarPair48 : public PolarChains {
  /*
   * Pair potential as in H. Li and G. Lykotrafitis. "Two-component coarse-grained molecular-dynamics model
   * for the human erythrocyte membrane". In: Biophysical journal 102.1 (2012):
   *   U = eps x^4 (x^4 - 2A) for r < rc,  x = (rc - r)/(rc - r0),  r0 = 2^(1/6) sigma,
   *   A = 1 - alpha (1 + a),  a = (n_i - n_j).e_ij / 2.
   * The cutoff rc is absolute.
   */

  public:
    /* attributes */
    double m_eps;
    double m_sigma;
    double m_r0;
    double m_rc;
    double m_alpha;

    /* constructor and destructor */
    PolarPair48(double eps, double sigma, double rc, double alpha, NeighborList* neighbors, std::vector<std::pair<size_t, size_t> > chain_ends);
    ~PolarPair48();

    /* methods */
    double pair_energy_force(const gsl_vector *xi, const gsl_vector *xj,
                             const gsl_vector *ni, const gsl_vector *nj,
                             gsl_vector *fj, gsl_vector *gi, gsl_vector *gj);
};

class PolarPair48_2site : public PolarChains {
  /*
   * Each bead carries a narrow site at x - rho_n n and a wide site at x + rho_w n.
   * Sites interact by U = eps x^4 (x^4 - 2) for r < rc, x = (rc - r)/(rc - r0),
   * r0 = 2^(1/6) sigma, with (eps, sigma, rc) depending on the pair of site types.
   * The cutoffs rc are absolute.
   */

  public:
    /* attributes */
    double m_eps_nn, m_eps_ww, m_eps_nw;
    double m_sigma_nn, m_sigma_ww, m_sigma_nw;
    double m_r0_nn, m_r0_ww, m_r0_nw;
    double m_rc_nn, m_rc_ww, m_rc_nw;
    double m_rho_n, m_rho_w;

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
    double pair_energy_force(const gsl_vector *xi, const gsl_vector *xj,
                             const gsl_vector *ni, const gsl_vector *nj,
                             gsl_vector *fj, gsl_vector *gi, gsl_vector *gj);
};
#endif
