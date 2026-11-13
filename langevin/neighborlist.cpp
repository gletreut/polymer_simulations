#include "neighborlist.h"

using namespace std;

//****************************************************************************
// Neighborlist
//****************************************************************************
NeighborList::NeighborList(size_t nmax, double rskin, size_t npart) :
  m_nmax(nmax), m_rskin(rskin), m_rcut(0.)
{
  m_npair = 0;
  m_iter = 0;
  m_iterref = 0;
  m_ibuild_count = 0;
  m_ibuild_sum = 0;
  m_rmax = 10.*m_rskin;

  // initialize holder for pairs
  m_pairs = gsl_matrix_uint_calloc(m_nmax, 2);

  // initialize holder for positions
  m_xref = gsl_matrix_calloc(npart,3);
}

NeighborList::~NeighborList()
{
  /* delete gsl objects */
  if (m_pairs != nullptr){
    gsl_matrix_uint_free(m_pairs);
  }

  if (m_xref != nullptr){
    gsl_matrix_free(m_xref);
  }
}

void
NeighborList::set_range(double rcut){
  /*
   * Register the interaction range of a force field using this list.
   * The skin rskin - rcut must be positive.
   */
  if (!(m_rskin > rcut)) {
    stringstream msg;
    msg << "Neighbor list: rskin = " << m_rskin << " must exceed the interaction range " << rcut << ".";
    throw invalid_argument(msg.str());
  }
  if (rcut > m_rcut)
    m_rcut = rcut;
}

void
NeighborList::update_counter(size_t iter){
  /*
   * Update the neighbor list counters
   */

    // update build counter
    m_ibuild_sum += (iter - m_iterref);
    m_ibuild_count += 1;
    m_iterref = iter;

  return;
}

void
NeighborList::build(gsl_matrix *x, gsl_spmatrix_uint *bonds){
  /*
   * Build the neighbor list by collecting all pairs separated by less
   * than the skin distance. Discard pairs which are bonded.
   */
  double r;
  gsl_vector *xtp(0);

  // initialize
  m_npair = 0; // reset the number of pairs of neighbors.
  xtp = gsl_vector_calloc(3);

  // iterate over particles
  for (size_t n=0; n<m_xref->size1; ++n) {
    gsl_vector_view xn = gsl_matrix_row(x, n);
    gsl_vector_view xn_ref = gsl_matrix_row(m_xref, n);
    gsl_vector_memcpy(&xn_ref.vector, &xn.vector);

    for (size_t m=n+1; m<m_xref->size1; ++m) {
      size_t bonded = gsl_spmatrix_uint_get(bonds, n, m);
      if (bonded == 1)
        continue;

      gsl_vector_view xm = gsl_matrix_row(x, m);
      gsl_vector_memcpy(xtp, &xm.vector);
      linalg_daxpy(-1, &xn.vector, xtp);
      r = linalg_dnrm2(xtp);
      if (r < m_rskin) {
        if (m_npair == m_nmax)
          throw length_error("Maximum number of neighbors reached! Increase npair_max.");
        gsl_matrix_uint_set(m_pairs, m_npair, 0, n);
        gsl_matrix_uint_set(m_pairs, m_npair, 1, m);
        m_npair += 1;
      }
    }
  }

  /* exit */
  gsl_vector_free(xtp);

  return;
}

bool
NeighborList::check(gsl_matrix *x){
  /*
   * Check whether the neighbor list needs to be rebuilt.
   * A pair can come closer by at most twice the largest displacement since
   * the last build, so the list stays complete while
   *   max_i |x_i - x_i^ref| <= (rskin - rcut)/2.
   */

  double dmax = 0.;

  for (size_t n=0; n<m_xref->size1; ++n) {
    double d2 = 0.;
    for (size_t a=0; a<3; ++a) {
      double dx = gsl_matrix_get(x, n, a) - gsl_matrix_get(m_xref, n, a);
      d2 += dx*dx;
    }
    if (d2 > dmax*dmax)
      dmax = sqrt(d2);
  }

  if (dmax > m_rmax)
    throw invalid_argument("Displacement is too large! Decrease the time step or increase rskin!");

  return dmax > 0.5*(m_rskin - m_rcut);
}

void
NeighborList::dump(std::ostream &mystream, bool positions) {
  /*
   * Dump neighbor list
   */

  mystream << left << dec << fixed;
	for (size_t i=0; i<m_npair; i++){
    size_t n = gsl_matrix_uint_get(m_pairs,i,0);
    size_t m = gsl_matrix_uint_get(m_pairs,i,1);
    mystream << setw(10) << setprecision(0) << noshowpos << i;
    mystream << setw(10) << setprecision(0) << noshowpos << n;
    mystream << setw(10) << setprecision(0) << noshowpos << m;
    if (positions){
      for (size_t a=0; a<3; ++a){
        mystream << setw(18) << setprecision(8) << showpos << gsl_matrix_get(m_xref,n,a);
      }
      for (size_t a=0; a<3; ++a){
        mystream << setw(18) << setprecision(8) << showpos << gsl_matrix_get(m_xref,m,a);
      }
    }
		mystream << endl;
	}

  return;
}
