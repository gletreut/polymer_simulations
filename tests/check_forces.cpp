//----------------------------------------------------------------------------
//  check_forces.cpp
//
//  Force/energy consistency for every force field: compare the force F
//  returned by energy_force() with the central difference
//      F_fd = -[U(x + eps e) - U(x - eps e)] / (2 eps),
//  where U is the full energy returned by the same energy_force() call.
//  Polar force fields are checked on one chain and on several unsorted chains,
//  and must reject chain_ends that do not partition the particles.
//  Also check that Point/Axis/Plane constraints hold after one Langevin step.
//
//  Build and run: make -C langevin check
//----------------------------------------------------------------------------
#include <cmath>
#include <cstdio>
#include <functional>
#include <stdexcept>
#include <string>
#include <vector>
#include <utility>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_spmatrix.h>

#include "forcefield.h"
#include "neighborlist.h"
#include "constraint.h"
#include "model.h"
#include "stepper.h"

using namespace std;

static const double FD_EPS = 1.0e-6;
static const double TOL_FORCE = 1.0e-5;
static const double TOL_CONSTRAINT = 1.0e-12;

typedef vector<pair<size_t, size_t> > ChainEnds;
typedef function<double(gsl_matrix *x, gsl_matrix *f)> EnergyForce;

//----------------------------------------------------------------------------
// configurations
//----------------------------------------------------------------------------
static double uniform(gsl_rng *rng, double a, double b) {
  return a + (b - a)*gsl_rng_uniform(rng);
}

static double dist(const gsl_matrix *x, size_t i, size_t j) {
  double s = 0.;
  for (size_t a=0; a<3; ++a) {
    double d = gsl_matrix_get(x, i, a) - gsl_matrix_get(x, j, a);
    s += d*d;
  }
  return sqrt(s);
}

static bool pairs_ok(const gsl_matrix *x, double rmin, const vector<double> &cutoffs) {
  /* all pairs farther than rmin, and not within 1e-3 of any cutoff */
  for (size_t i=0; i<x->size1; ++i)
    for (size_t j=i+1; j<x->size1; ++j) {
      double r = dist(x, i, j);
      if (r < rmin) return false;
      for (size_t k=0; k<cutoffs.size(); ++k)
        if (fabs(r - cutoffs[k]) < 1.0e-3) return false;
    }
  return true;
}

static gsl_matrix* random_chain(gsl_rng *rng, size_t N, size_t dim,
                                double rmin=0.85, vector<double> cutoffs=vector<double>()) {
  /* random walk with bond lengths in [0.95, 1.05]; z = 0 when dim = 2 */
  gsl_matrix *x = gsl_matrix_calloc(N, 3);
  do {
    gsl_matrix_set_all(x, 0.);
    for (size_t i=1; i<N; ++i) {
      double b = uniform(rng, 0.95, 1.05);
      double u[3] = {0., 0., 0.};
      if (dim == 2) {
        double t = uniform(rng, 0., 2.*M_PI);
        u[0] = cos(t); u[1] = sin(t);
      } else {
        double cz = uniform(rng, -1., 1.), t = uniform(rng, 0., 2.*M_PI);
        double s = sqrt(1. - cz*cz);
        u[0] = s*cos(t); u[1] = s*sin(t); u[2] = cz;
      }
      for (size_t a=0; a<3; ++a)
        gsl_matrix_set(x, i, a, gsl_matrix_get(x, i-1, a) + b*u[a]);
    }
  } while (!pairs_ok(x, rmin, cutoffs));
  return x;
}

//----------------------------------------------------------------------------
// finite-difference check
//----------------------------------------------------------------------------
static double g_fmax = 0.;  // max|F| of the last fd_error call

static double fd_error(const EnergyForce &ef, gsl_matrix *x, size_t dim) {
  /* max|F - F_fd| / max(max|F|, 1) */
  size_t N = x->size1;
  gsl_matrix *f = gsl_matrix_calloc(N, 3);
  gsl_matrix *ftp = gsl_matrix_calloc(N, 3);
  ef(x, f);

  double fmax = 0., emax = 0.;
  for (size_t i=0; i<N; ++i)
    for (size_t a=0; a<dim; ++a) {
      double x0 = gsl_matrix_get(x, i, a);
      gsl_matrix_set(x, i, a, x0 + FD_EPS);
      double up = ef(x, ftp);
      gsl_matrix_set(x, i, a, x0 - FD_EPS);
      double um = ef(x, ftp);
      gsl_matrix_set(x, i, a, x0);

      double fa = gsl_matrix_get(f, i, a);
      double ffd = -(up - um)/(2.*FD_EPS);
      fmax = fmax > fabs(fa) ? fmax : fabs(fa);
      emax = emax > fabs(fa - ffd) ? emax : fabs(fa - ffd);
    }

  gsl_matrix_free(f);
  gsl_matrix_free(ftp);
  g_fmax = fmax;
  return emax / (fmax > 1. ? fmax : 1.);
}

static EnergyForce wrap(ForceField *ff) {
  /* full energy and force as the code computes them */
  return [ff](gsl_matrix *x, gsl_matrix *f) {
    double u = 0.;
    gsl_matrix_set_all(f, 0.);
    ff->energy_force(x, &u, f);
    return u;
  };
}

static EnergyForce wrap_frozen_polarity(PolarChains *ff) {
  /*
   * Diagnostic: energy and force with the polarity vectors frozen at their
   * current values (computed once by energy_force at the reference x).
   * Checks the pair force alone, without the d(polarity)/dx terms.
   */
  return [ff](gsl_matrix *x, gsl_matrix *f) {
    double u = 0.;
    gsl_vector *ftp = gsl_vector_calloc(3);
    gsl_vector *gi = gsl_vector_calloc(3), *gj = gsl_vector_calloc(3);
    gsl_matrix_set_all(f, 0.);
    NeighborList *nl = ff->m_neighbors;
    for (size_t k=0; k<nl->m_npair; ++k) {
      size_t n = gsl_matrix_uint_get(nl->m_pairs, k, 0);
      size_t m = gsl_matrix_uint_get(nl->m_pairs, k, 1);
      gsl_vector_view xn = gsl_matrix_row(x, n), xm = gsl_matrix_row(x, m);
      gsl_vector_view pn = gsl_matrix_row(ff->m_pol_vec, n), pm = gsl_matrix_row(ff->m_pol_vec, m);
      gsl_vector_view fn = gsl_matrix_row(f, n), fm = gsl_matrix_row(f, m);
      u += ff->pair_energy_force(&xn.vector, &xm.vector, &pn.vector, &pm.vector, ftp, gi, gj);
      gsl_vector_sub(&fn.vector, ftp);
      gsl_vector_add(&fm.vector, ftp);
    }
    gsl_vector_free(ftp);
    gsl_vector_free(gi);
    gsl_vector_free(gj);
    return u;
  };
}

static void set_polarity(ForceField *ff, gsl_matrix *x) {
  /* one energy_force call sets the polarity vectors at x */
  gsl_matrix *f = gsl_matrix_calloc(x->size1, 3);
  wrap(ff)(x, f);
  gsl_matrix_free(f);
}

static int nfail = 0;

static void report(const string &name, double err, bool counted=true) {
  bool pass = err < TOL_FORCE;
  if (counted && !pass) nfail += 1;
  printf("%-34s %12.3e %12.3e  %s\n", name.c_str(), g_fmax, err, counted ? (pass ? "PASS" : "FAIL") : (pass ? "(pass)" : "(fail)"));
}

template <class Make>
static void expect_invalid_argument(const string &name, Make make) {
  /* PASS if make() throws std::invalid_argument */
  bool pass = false;
  try { make(); }
  catch (const invalid_argument &) { pass = true; }
  catch (...) {}
  if (!pass) nfail += 1;
  printf("%-34s %12s %12s  %s\n", name.c_str(), "-", "-", pass ? "PASS" : "FAIL");
}

static NeighborList* all_pairs(gsl_matrix *x, gsl_spmatrix_uint *bonds) {
  /* neighbor list holding every non-bonded pair */
  size_t N = x->size1;
  NeighborList *nl = new NeighborList(N*N, 1.0e3, N);
  nl->build(x, bonds);
  return nl;
}

//----------------------------------------------------------------------------
// constraint check
//----------------------------------------------------------------------------
static void check_constraints(gsl_rng *rng) {
  const size_t N = 10;
  const double p[3] = {0.3, -0.2, 0.1};
  double u[3] = {1., 2., 3.};
  double un = sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
  for (size_t a=0; a<3; ++a) u[a] /= un;

  MDWorld *world = new MDWorld(N, 20., 20., 20., 1.0, 1.0, 1.0, 3);
  gsl_matrix *x = random_chain(rng, N, 3);
  gsl_matrix_memcpy(world->m_x, x);
  gsl_matrix_free(x);
  world->init_velocities(rng);
  world->m_ffields.push_back(new PolymerHarmonic(0, N, 200., 1., world->m_bonds));
  world->m_ffields.push_back(new PolymerKratkyPorod(0, N, 3., world->m_bonds));
  world->m_constraints.push_back(new PointConstraint(0, p[0], p[1], p[2]));
  world->m_constraints.push_back(new AxisConstraint(1, p[0], p[1], p[2], u[0], u[1], u[2]));
  world->m_constraints.push_back(new PlaneConstraint(2, p[0], p[1], p[2], u[0], u[1], u[2]));

  LangStepper_VVerlet stepper(world, rng, 0.01, world->m_gamma);
  world->set_constraints();
  stepper.call_forces();
  stepper.step();

  // d = x - p; residuals: |d| (point), |d - (d.u)u| (axis), |d.u| (plane)
  double d[3][3], du[3] = {0., 0., 0.};
  for (size_t i=0; i<3; ++i)
    for (size_t a=0; a<3; ++a) {
      d[i][a] = gsl_matrix_get(world->m_x, i, a) - p[a];
      du[i] += d[i][a]*u[a];
    }
  double res_point = 0., res_axis = 0.;
  for (size_t a=0; a<3; ++a) {
    res_point += d[0][a]*d[0][a];
    res_axis += (d[1][a] - du[1]*u[a])*(d[1][a] - du[1]*u[a]);
  }
  double res[3] = {sqrt(res_point), sqrt(res_axis), fabs(du[2])};
  const char *names[3] = {"PointConstraint", "AxisConstraint", "PlaneConstraint"};
  for (size_t i=0; i<3; ++i) {
    bool pass = res[i] < TOL_CONSTRAINT;
    if (!pass) nfail += 1;
    printf("%-34s %12s %12.3e  %s\n", names[i], "-", res[i], pass ? "PASS" : "FAIL");
  }
  delete world;
}

//----------------------------------------------------------------------------
// main
//----------------------------------------------------------------------------
int main() {
  const size_t N = 10;
  gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(rng, 12345);
  gsl_spmatrix_uint *bonds = gsl_spmatrix_uint_alloc(N, N);
  ChainEnds chain_ends(1, make_pair(size_t(0), N-1));
  // several chains, unsorted, including a 2-bead chain
  const ChainEnds chains3 = {{6, 9}, {0, 1}, {2, 5}};

  printf("%-34s %12s %12s  %s\n", "check", "max|F|", "error", "status");

  // --- confinement: random points near the walls, away from them
  {
    const double L = 1.8;  // half box size
    gsl_matrix *x = gsl_matrix_calloc(N, 3);
    for (size_t i=0; i<N; ++i)
      for (size_t a=0; a<3; ++a)
        gsl_matrix_set(x, i, a, uniform(rng, -L + 0.88, L - 0.88));
    ConfinementBox ff(-L, L, -L, L, -L, L, 1.0, 1.0);
    report("ConfinementBox", fd_error(wrap(&ff), x, 3));
    gsl_matrix_free(x);
  }
  {
    const double R = 2.0;
    gsl_matrix *x = gsl_matrix_calloc(N, 3);
    for (size_t i=0; i<N; ++i) {
      double r, v[3];
      do {
        r = 0.;
        for (size_t a=0; a<3; ++a) { v[a] = uniform(rng, -R, R); r += v[a]*v[a]; }
        r = sqrt(r);
      } while (r > R - 0.88 || r < 0.5);
      for (size_t a=0; a<3; ++a) gsl_matrix_set(x, i, a, v[a]);
    }
    ConfinementSphere ff(R, 1.0, 1.0);
    report("ConfinementSphere", fd_error(wrap(&ff), x, 3));
    gsl_matrix_free(x);
  }

  // --- bonded: random 3D chains
  {
    gsl_matrix *x = random_chain(rng, N, 3);
    PolymerGaussian ff(0, N, 1.0, bonds);
    report("PolymerGaussian", fd_error(wrap(&ff), x, 3));
    gsl_matrix_free(x);
  }
  {
    gsl_matrix *x = random_chain(rng, N, 3);
    PolymerHarmonic ff(0, N, 200., 1.0, bonds);
    report("PolymerHarmonic", fd_error(wrap(&ff), x, 3));
    gsl_matrix_free(x);
  }
  {
    gsl_matrix *x = random_chain(rng, N, 3);
    PolymerFENE ff(0, N, 30., 1.5, 1.0, 1.0, bonds);
    report("PolymerFENE", fd_error(wrap(&ff), x, 3));
    gsl_matrix_free(x);
  }
  {
    gsl_matrix *x = random_chain(rng, N, 3);
    PolymerKratkyPorod ff(0, N, 3.0, bonds);
    report("PolymerKratkyPorod", fd_error(wrap(&ff), x, 3));
    gsl_matrix_free(x);
  }

  // --- pair potentials: all pairs in the neighbor list (no bonds excluded)
  gsl_spmatrix_uint *nobonds = gsl_spmatrix_uint_alloc(N, N);
  {
    const double sigma = 1.0;
    SoftCore *ff = 0;
    gsl_matrix *x = random_chain(rng, N, 3, 0.85, vector<double>(1, pow(2., 1./6)*sigma));
    NeighborList *nl = all_pairs(x, nobonds);
    ff = new SoftCore(10., sigma, nl);
    report("SoftCore", fd_error(wrap(ff), x, 3));
    delete ff; delete nl; gsl_matrix_free(x);
  }
  {
    const double rc = 2.5;
    gsl_matrix *x = random_chain(rng, N, 3, 0.85, vector<double>(1, rc));
    NeighborList *nl = all_pairs(x, nobonds);
    PairLJ ff(1.0, 1.0, rc, nl);
    report("PairLJ", fd_error(wrap(&ff), x, 3));
    delete nl; gsl_matrix_free(x);
  }

  // --- polar pair potentials: 2D chains (z = 0), chain_ends cover all beads
  {
    gsl_matrix *x = random_chain(rng, N, 2, 0.95);
    NeighborList *nl = all_pairs(x, nobonds);
    PolarPairLJ ff(1.0, 1.0, 2.5, 0.2, nl, chain_ends);
    report("PolarPairLJ", fd_error(wrap(&ff), x, 2));
    set_polarity(&ff, x);
    report("  PolarPairLJ (frozen polarity)", fd_error(wrap_frozen_polarity(&ff), x, 2), false);
    PolarPairLJ ff3(1.0, 1.0, 2.5, 0.2, nl, chains3);
    report("PolarPairLJ (3 chains)", fd_error(wrap(&ff3), x, 2));
    delete nl; gsl_matrix_free(x);
  }
  {
    gsl_matrix *x = random_chain(rng, N, 2, 0.95, vector<double>(1, 2.0));
    NeighborList *nl = all_pairs(x, nobonds);
    const double w48 = 2.0 - pow(2., 1./6);   // rc = 2
    PolarPair48 ff(1.0, 1.0, w48, 0.3, nl, chain_ends);
    report("PolarPair48", fd_error(wrap(&ff), x, 2));
    set_polarity(&ff, x);
    report("  PolarPair48 (frozen polarity)", fd_error(wrap_frozen_polarity(&ff), x, 2), false);
    PolarPair48 ff3(1.0, 1.0, w48, 0.3, nl, chains3);
    report("PolarPair48 (3 chains)", fd_error(wrap(&ff3), x, 2));
    PolarPair48 ffa(1.0, 1.0, w48, 0.8, nl, chain_ends);
    report("PolarPair48 (alpha = 0.8)", fd_error(wrap(&ffa), x, 2));
    delete nl; gsl_matrix_free(x);
  }
  {
    gsl_matrix *x = random_chain(rng, N, 2, 0.95);
    NeighborList *nl = all_pairs(x, nobonds);
    const double w2 = 1.5/0.6 - pow(2., 1./6);   // rc = 1.5
    PolarPair48_2site ff(1.0, 0.5, 0.8,     // eps nn, ww, nw
                         0.6, 0.6, 0.6,     // sigma
                         w2, w2, w2,        // width
                         0.3, 0.3,          // rho n, w
                         nl, chain_ends);
    report("PolarPair48_2site", fd_error(wrap(&ff), x, 2));
    set_polarity(&ff, x);
    report("  PolarPair48_2site (frozen pol.)", fd_error(wrap_frozen_polarity(&ff), x, 2), false);
    PolarPair48_2site ff3(1.0, 0.5, 0.8, 0.6, 0.6, 0.6, w2, w2, w2, 0.3, 0.3, nl, chains3);
    report("PolarPair48_2site (3 chains)", fd_error(wrap(&ff3), x, 2));
    delete nl; gsl_matrix_free(x);
  }

  // --- polar chain_ends that do not partition 0..N-1 are rejected
  {
    const pair<const char*, ChainEnds> bad[3] = {
      {"gap", {{0, 4}, {6, 9}}},
      {"overlap", {{0, 5}, {4, 9}}},
      {"singleton", {{0, 0}, {1, 9}}},
    };
    NeighborList nl(100, 1.0e3, 10);
    for (const auto &b : bad)
      expect_invalid_argument(string("PolarPair48 bad ends: ") + b.first,
                              [&] { PolarPair48 ff(1.0, 1.0, 0.9, 0.3, &nl, b.second); });
    const ChainEnds ends = {{0, 9}};
    expect_invalid_argument("PolarPair48 width = 0",
                            [&] { PolarPair48 ff(1.0, 1.0, 0.0, 0.3, &nl, ends); });
    expect_invalid_argument("PolarPair48 sigma = 0",
                            [&] { PolarPair48 ff(1.0, 0.0, 0.9, 0.3, &nl, ends); });
    expect_invalid_argument("PolarPair48_2site width < 0",
                            [&] { PolarPair48_2site ff(1.0, 1.0, 1.0, 0.6, 0.6, 0.6, 1.0, -1.0, 1.0,
                                                       0.3, 0.3, &nl, ends); });
  }

  // --- constraints after one Lang_VVerlet step
  check_constraints(rng);

  printf("%d check(s) failed\n", nfail);
  gsl_spmatrix_uint_free(bonds);
  gsl_spmatrix_uint_free(nobonds);
  gsl_rng_free(rng);
  return nfail == 0 ? 0 : 1;
}
