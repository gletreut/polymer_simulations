#include "forcefield.h"

using namespace std;

//****************************************************************************
// PolymerFENE
//****************************************************************************
PolymerFENE::PolymerFENE(size_t offset, size_t N, double ke, double rc, double sigma, double eps) :
  m_offset(offset),
  m_N(N),
  m_ke(ke),
  m_rc(rc),
  m_sigma(sigma),
  m_eps(eps)
{
  m_pref = -0.5*m_ke*m_rc*m_rc;
  m_4eps = 4.0*m_eps;
  m_48eps = 48.0*m_eps;
  m_rc_LJ = m_sigma*pow(2.,1./6);

}

  /* methods */
void PolymerFENE::energy_force(gsl_matrix *x, double *u, gsl_matrix *force){
  /*
   * Compute the potential energy and force (minus gradient) of a FENE chain.
   */

  size_t n;
  double r,fr,y,y2, y6,y12,fnorm;
  gsl_vector *bond(0);

  // initialize
  bond = gsl_vector_calloc(3);

  // set energy and force to zero
  //*u = 0;
  //gsl_matrix_set_all(force, 0.0);

  // iterate over the bonds (r_\{i+1\} - r_i)
  for (size_t i=m_offset; i<m_N-1; i++){
    n = m_offset + i;
    gsl_vector_view vn = gsl_matrix_row(x, n);
    gsl_vector_view vnp = gsl_matrix_row(x, n+1);

    gsl_vector_memcpy(bond, &vnp.vector);

    linalg_daxpy(-1., &vn.vector, bond);
    r = linalg_dnrm2(bond);
    fr = r/m_rc;
    y = m_sigma/r;
    y2 = y*y;
    y6 = y2*y2*y2;
    y12 = y6*y6;

    // check length of bond
    if ( !(fr < 1.) ){
      stringstream convert;
      convert.clear();
      convert << "FENE bond too long: " << "monomers " << n << "," << n+1 << "  r = " << r;
      throw runtime_error(convert.str());
    }

    // energy
    *u += m_pref*log(1.0-fr*fr);          // fene contribution
    if ( fabs(r) < m_rc_LJ)
      *u += m_4eps*(y12-y6) + m_eps;      // LJ contribution

    // forces
    //  FENE
    gsl_vector_view fn = gsl_matrix_row(force, n);
    gsl_vector_view fnp = gsl_matrix_row(force, n+1);
    fnorm = m_ke / (1.0 - fr*fr);
    linalg_daxpy(-fnorm, bond, &fn.vector);
    linalg_daxpy(fnorm, bond, &fnp.vector);
    //  LJ
    if ( fabs(r) < m_rc_LJ) {
      fnorm = m_48eps / (r*r) * (y12-0.5*y6);
      linalg_daxpy(-fnorm, bond, &fn.vector);
      linalg_daxpy(fnorm, bond, &fnp.vector);
    }
  }

  /* exit */
  gsl_vector_free(bond);
  return;
}
