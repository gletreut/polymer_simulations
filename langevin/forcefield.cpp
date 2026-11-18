#include "forcefield.h"
using namespace std;

//****************************************************************************
// ForceField
//****************************************************************************
ForceField::~ForceField(){
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------

//****************************************************************************
// ConfinementBox
//****************************************************************************
ConfinementBox::ConfinementBox(double xlo, double xhi, double ylo, double yhi, double zlo, double zhi, double sigma, double eps) :
  m_xlo(xlo),
  m_xhi(xhi),
  m_ylo(ylo),
  m_yhi(yhi),
  m_zlo(zlo),
  m_zhi(zhi),
  m_sigma(sigma),
  m_eps(eps)

{
  m_4eps = 4.0*m_eps;
  m_fpref = 48.0*m_eps / m_sigma;
  m_rc_LJ = m_sigma*pow(2.,1./6);
}

ConfinementBox::~ConfinementBox(){
}

double ConfinementBox::energy_LJ_scal(double r){
  /*
   * compute the energy of the  LJ interaction when the algebraic
   * distance is r.
   */

  double x, x6, x12;

  if (fabs(r) < m_rc_LJ) {
    x = m_sigma/r;
    x6 = x*x*x*x*x*x;
    x12 = x6*x6;
    return m_4eps*(x12 - x6) + m_eps;
  }
  else {
    return 0.0;
  }
}

double ConfinementBox::force_LJ_scal(double r){
  /*
   * compute the algebraic norm of the LJ force applied when the algebraic
   * distance is r.
   */
  double x, x6, x12;

  if (fabs(r) < m_rc_LJ) {
    x = m_sigma/r;
    x6 = x*x*x*x*x*x;
    x12 = x6*x6;
    return m_fpref*x*(x12 - 0.5*x6);
  }
  else {
    return 0.0;
  }
}

void ConfinementBox::energy_force(gsl_matrix *x, double *u, gsl_matrix *forces){
  /*
   * Compute the potential energy and force (minus gradient) produced by a
   * confinement box
   */

  size_t N;
  double rx,ry,rz,r,fnorm;
  gsl_vector *force(0);

  // initialize
  N = x->size1;
  force = gsl_vector_calloc(3);

  // iterate over the positions (r_\{i\})
  for (size_t i=0; i<N; ++i){
    gsl_vector_view xi = gsl_matrix_row(x, i);
    gsl_vector_view fi = gsl_matrix_row(forces, i);
    rx = gsl_vector_get(&xi.vector,0);
    ry = gsl_vector_get(&xi.vector,1);
    rz = gsl_vector_get(&xi.vector,2);

    /* xlo */
    gsl_vector_set_all(force,0.0);
    gsl_vector_set(force,0,1.0);  // ex
    r = rx - m_xlo;
    *u += energy_LJ_scal(r);
    fnorm = force_LJ_scal(r);  // algebraic norm positive
    linalg_daxpy(fnorm, force, &fi.vector);

    /* xhi */
    gsl_vector_set_all(force,0.0);
    gsl_vector_set(force,0,1.0);  // ex
    r = rx - m_xhi;
    *u += energy_LJ_scal(r);
    fnorm = force_LJ_scal(r);  // algebraic norm positive
    linalg_daxpy(fnorm, force, &fi.vector);

    /* ylo */
    gsl_vector_set_all(force,0.0);
    gsl_vector_set(force,1,1.0);  // ey
    r = ry - m_ylo;
    *u += energy_LJ_scal(r);
    fnorm = force_LJ_scal(r);  // algebraic norm positive
    linalg_daxpy(fnorm, force, &fi.vector);

    /* yhi */
    gsl_vector_set_all(force,0.0);
    gsl_vector_set(force,1,1.0);  // ey
    r = ry - m_yhi;
    *u += energy_LJ_scal(r);
    fnorm = force_LJ_scal(r);  // algebraic norm positive
    linalg_daxpy(fnorm, force, &fi.vector);

    /* zlo */
    gsl_vector_set_all(force,0.0);
    gsl_vector_set(force,2,1.0);  // ez
    r = rz - m_zlo;
    *u += energy_LJ_scal(r);
    fnorm = force_LJ_scal(r);  // algebraic norm positive
    linalg_daxpy(fnorm, force, &fi.vector);

    /* zhi */
    gsl_vector_set_all(force,0.0);
    gsl_vector_set(force,2,1.0);  // ez
    r = rz - m_zhi;
    *u += energy_LJ_scal(r);
    fnorm = force_LJ_scal(r);  // algebraic norm positive
    linalg_daxpy(fnorm, force, &fi.vector);
  }

  /* exit */
  gsl_vector_free(force);
  return;
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------

//****************************************************************************
// ConfinementSphere
//****************************************************************************
ConfinementSphere::ConfinementSphere(double radius, double sigma, double eps) :
  m_radius(radius),
  m_sigma(sigma),
  m_eps(eps)

{
  m_4eps = 4.0*m_eps;
  m_fpref = 48.0*m_eps / m_sigma;
  m_rc_LJ = m_sigma*pow(2.,1./6);
}

ConfinementSphere::~ConfinementSphere(){
}

double ConfinementSphere::energy_LJ_scal(double r){
  /*
   * compute the energy of the  LJ interaction when the algebraic
   * distance is r.
   */

  double x, x6, x12;

  if (fabs(r) < m_rc_LJ) {
    x = m_sigma/r;
    x6 = x*x*x*x*x*x;
    x12 = x6*x6;
    return m_4eps*(x12 - x6) + m_eps;
  }
  else {
    return 0.0;
  }
}

double ConfinementSphere::force_LJ_scal(double r){
  /*
   * compute the algebraic norm of the LJ force applied when the algebraic
   * distance is r.
   */
  double x, x6, x12;

  if (fabs(r) < m_rc_LJ) {
    x = m_sigma/r;
    x6 = x*x*x*x*x*x;
    x12 = x6*x6;
    return m_fpref*x*(x12 - 0.5*x6);
  }
  else {
    return 0.0;
  }
}

void ConfinementSphere::energy_force(gsl_matrix *x, double *u, gsl_matrix *forces){
  /*
   * Compute the potential energy and force (minus gradient) produced by a
   * confinement sphere
   */

  size_t N;
  double rx,ry,rz,r,fnorm;
  gsl_vector *force(0);

  // initialize
  N = x->size1;
  force = gsl_vector_calloc(3);

  // iterate over the positions (r_\{i\})
  for (size_t i=0; i<N; ++i){
    gsl_vector_view xi = gsl_matrix_row(x, i);
    gsl_vector_view fi = gsl_matrix_row(forces, i);
    rx = gsl_vector_get(&xi.vector,0);
    ry = gsl_vector_get(&xi.vector,1);
    rz = gsl_vector_get(&xi.vector,2);
    r = sqrt(rx*rx+ry*ry+rz*rz);

    gsl_vector_set_all(force,0.0);
    // er
    gsl_vector_set(force,0,rx/r);
    gsl_vector_set(force,1,ry/r);
    gsl_vector_set(force,2,rz/r);
    // distance
    r = r - m_radius;
    *u += energy_LJ_scal(r);
    fnorm = force_LJ_scal(r);  // algebraic norm positive
    linalg_daxpy(fnorm, force, &fi.vector);
  }

  /* exit */
  gsl_vector_free(force);
  return;
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------

//****************************************************************************
// PolymerGaussian
//****************************************************************************
PolymerGaussian::PolymerGaussian(size_t offset, size_t N, double b, gsl_spmatrix_uint *bonds) :
  m_offset(offset),
  m_N(N),
  m_b(b),
  m_bonds(bonds)
{
  /* initialize some parameters */
  m_ke = 1.5/(m_b*m_b);
  m_fpref = 2.*m_ke;

  for (size_t i=0; i<m_N-1; ++i){
    size_t n = m_offset + i;
    gsl_spmatrix_uint_set(m_bonds, n, n+1, 1);
  }
}

PolymerGaussian::~PolymerGaussian(){
}

void PolymerGaussian::energy_force(gsl_matrix *x, double *u, gsl_matrix *forces){
  /*
   * Compute the potential energy and force (minus gradient) of a Gaussian chain.
   */

  size_t n;
  double r;
  gsl_vector *bond(0);

  // initialize
  bond = gsl_vector_calloc(3);

  // iterate over the bonds (r_\{i+1\} - r_i)
  for (size_t i=0; i<m_N-1; ++i){
    n = m_offset + i;
    gsl_vector_view vn = gsl_matrix_row(x, n);
    gsl_vector_view vnp = gsl_matrix_row(x, n+1);

    // compute bond
    gsl_vector_memcpy(bond, &vnp.vector);
    linalg_daxpy(-1., &vn.vector, bond);
    r = linalg_dnrm2(bond);

    // energy
    *u += m_ke*r*r;

    // forces
    gsl_vector_view fn = gsl_matrix_row(forces, n);
    gsl_vector_view fnp = gsl_matrix_row(forces, n+1);
    linalg_daxpy(m_fpref, bond, &fn.vector);
    linalg_daxpy(-m_fpref, bond, &fnp.vector);
  }

  /* exit */
  gsl_vector_free(bond);
  return;
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------

//****************************************************************************
// PolymerHarmonic
//****************************************************************************
PolymerHarmonic::PolymerHarmonic(size_t offset, size_t N, double ke, double r0, gsl_spmatrix_uint *bonds) :
  m_offset(offset),
  m_N(N),
  m_ke(ke),
  m_r0(r0),
  m_bonds(bonds)
{
  /* initialize some parameters */
  m_fpref = 2.*m_ke;

  for (size_t i=0; i<m_N-1; ++i){
    size_t n = m_offset + i;
    gsl_spmatrix_uint_set(m_bonds, n, n+1, 1);
  }
}

PolymerHarmonic::~PolymerHarmonic(){
}

void PolymerHarmonic::energy_force(gsl_matrix *x, double *u, gsl_matrix *forces){
  /*
   * Compute the potential energy and force (minus gradient) of a Gaussian chain.
   */

  size_t n;
  double r,dr,fr;
  gsl_vector *bond(0);

  // initialize
  bond = gsl_vector_calloc(3);

  // iterate over the bonds (r_\{i+1\} - r_i)
  for (size_t i=0; i<m_N-1; ++i){
    n = m_offset + i;
    gsl_vector_view vn = gsl_matrix_row(x, n);
    gsl_vector_view vnp = gsl_matrix_row(x, n+1);

    // compute bond
    gsl_vector_memcpy(bond, &vnp.vector);
    linalg_daxpy(-1., &vn.vector, bond);
    r = linalg_dnrm2(bond);
    dr = r-m_r0;

    // energy
    *u += m_ke*dr*dr;

    // forces
    fr = (1.-m_r0/r);
    gsl_vector_view fn = gsl_matrix_row(forces, n);
    gsl_vector_view fnp = gsl_matrix_row(forces, n+1);
    linalg_daxpy(m_fpref*fr, bond, &fn.vector);
    linalg_daxpy(-m_fpref*fr, bond, &fnp.vector);
  }

  /* exit */
  gsl_vector_free(bond);
  return;
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------

//****************************************************************************
// PolymerFENE
//****************************************************************************
PolymerFENE::PolymerFENE(size_t offset, size_t N, double ke, double rc, double sigma, double eps, gsl_spmatrix_uint *bonds) :
  m_offset(offset),
  m_N(N),
  m_ke(ke),
  m_rc(rc),
  m_sigma(sigma),
  m_eps(eps),
  m_bonds(bonds)
{
  m_pref = -0.5*m_ke*m_rc*m_rc;
  m_4eps = 4.0*m_eps;
  m_fpref_LJ = 48.0*m_eps / m_sigma;
  m_rc_LJ = m_sigma*pow(2.,1./6);

  for (size_t i=0; i<m_N-1; ++i){
    size_t n = m_offset + i;
    gsl_spmatrix_uint_set(m_bonds, n, n+1, 1);
  }
}

PolymerFENE::~PolymerFENE(){
}

double PolymerFENE::energy_LJ_scal(double r){
  /* methods */
  /*
   * compute the energy of the  LJ interaction when the algebraic
   * distance is r.
   */

  double x, x6, x12;

  if (fabs(r) < m_rc_LJ) {
    x = m_sigma/r;
    x6 = x*x*x*x*x*x;
    x12 = x6*x6;
    return m_4eps*(x12 - x6) + m_eps;
  }
  else {
    return 0.0;
  }
}

double PolymerFENE::force_LJ_scal(double r){
  /*
   * compute the algebraic norm of the LJ force applied when the algebraic
   * distance is r.
   */
  double x, x6, x12;

  if (fabs(r) < m_rc_LJ) {
    x = m_sigma/r;
    x6 = x*x*x*x*x*x;
    x12 = x6*x6;
    return m_fpref_LJ*x*(x12 - 0.5*x6);
  }
  else {
    return 0.0;
  }
}

void PolymerFENE::energy_force(gsl_matrix *x, double *u, gsl_matrix *forces){
  /*
   * Compute the potential energy and force (minus gradient) of a FENE chain.
   */

  size_t n;
  double r,fr,fnorm;
  gsl_vector *bond(0);
  stringstream convert;

  // initialize
  bond = gsl_vector_calloc(3);

  // iterate over the bonds (r_\{i+1\} - r_i)
  for (size_t i=0; i<m_N-1; ++i){
    n = m_offset + i;
    gsl_vector_view vn = gsl_matrix_row(x, n);
    gsl_vector_view vnp = gsl_matrix_row(x, n+1);

    gsl_vector_memcpy(bond, &vnp.vector);
    linalg_daxpy(-1., &vn.vector, bond);

    r = linalg_dnrm2(bond);
    fr = r/m_rc;

    // check length of bond
    if ( !(fr < 1.) ){
      convert.clear();
      convert << "FENE bond too long: " << "monomers " << n << "," << n+1 << "  r = " << r << endl;
      convert << "vn: " << endl;
      utils::print_vector(convert, &vn.vector);
      convert << "vnp: " << endl;
      utils::print_vector(convert, &vnp.vector);
      convert << "bond: " << endl;
      utils::print_vector(convert, bond);
      throw runtime_error(convert.str());
    }

    // energy
    *u += m_pref*log(1.0-fr*fr);          // fene contribution
    *u += energy_LJ_scal(r);

    // forces
    //  FENE
    gsl_vector_view fn = gsl_matrix_row(forces, n);
    gsl_vector_view fnp = gsl_matrix_row(forces, n+1);
    fnorm = -m_ke / (1.0 - fr*fr);
    linalg_daxpy(-fnorm, bond, &fn.vector);
    linalg_daxpy(+fnorm, bond, &fnp.vector);
    //  LJ
    fnorm = force_LJ_scal(r);  // algebraic norm
    fnorm /= r;                // bond is not normalized to 1
    linalg_daxpy(-fnorm, bond, &fn.vector);
    linalg_daxpy(+fnorm, bond, &fnp.vector);
  }

  /* exit */
  gsl_vector_free(bond);
  return;
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------

//****************************************************************************
// PolymerKratkyPorod
//****************************************************************************
PolymerKratkyPorod::PolymerKratkyPorod(size_t offset, size_t N, double lp, gsl_spmatrix_uint *bonds) :
  m_offset(offset),
  m_N(N),
  m_lp(lp),
  m_bonds(bonds)
{
  for (size_t i=1; i<m_N-1; ++i){
    size_t n = m_offset + i;
    gsl_spmatrix_uint_set(m_bonds, n-1, n, 1);
    gsl_spmatrix_uint_set(m_bonds, n, n+1, 1);
    gsl_spmatrix_uint_set(m_bonds, n-1, n+1, 1);
  }
}

PolymerKratkyPorod::~PolymerKratkyPorod() {
}

void PolymerKratkyPorod::energy_force(gsl_matrix *x, double *u, gsl_matrix *forces){
  /* methods */
  /*
   * Compute the potential energy and force (minus gradient) of a Kratky-Porod
   * chain potential.
   * \beta U = l_p \sum_{i=1)^{N-1} (1 - \cos(\theta_i)),
   *   \cos(\theta_i)) = (u_i.u_{i+1})/(||u_i|| ||u_{i+1}||)
   *
   * Each joint i generate a force on the 3 surrounding monomers such that
   * (showing only x component):
   * f_{i-1} / lp = -(x_{i+1} - x_{i})/(||u_i|| ||u_{i+1}||) + (x_{i} - x_{i-1}) \cos(\theta_i)/||u_i||^2
   * f_{i+1} / lp = -(x_{i+1} - x_{i})\cos(\theta_i)/||u_{i+1}||^2 + (x_{i} - x_{i-1})/(||u_i|| ||u_{i+1}||)
   * f_{i} = - f_{i-1} - f{i+1}
   */

  size_t n;
  double ci,b2,bn2,bb;
  gsl_vector *bond(0), *bond_n(0), *fm(0), *fp(0);

  // initialize
  bond = gsl_vector_calloc(3);
  bond_n = gsl_vector_calloc(3);
  fm = gsl_vector_calloc(3);
  fp = gsl_vector_calloc(3);

  // initialize bond
  n = m_offset;
  gsl_vector_view v = gsl_matrix_row(x, n);
  gsl_vector_view vn = gsl_matrix_row(x, n+1);
  gsl_vector_memcpy(bond_n, &vn.vector);
  linalg_daxpy(-1., &v.vector, bond_n);

  // iterate over the bonds joints
  for (size_t i=1; i<m_N-1; ++i){
    n = m_offset + i;
    gsl_vector_memcpy(bond,bond_n);
    gsl_vector_view v = gsl_matrix_row(x, n);
    gsl_vector_view vn = gsl_matrix_row(x, n+1);
    gsl_vector_memcpy(bond_n, &vn.vector);
    linalg_daxpy(-1., &v.vector, bond_n);

    // compute quantities
    b2 = linalg_ddot(bond,bond);
    bn2 = linalg_ddot(bond_n,bond_n);
    bb = sqrt(b2)*sqrt(bn2);
    ci = linalg_ddot(bond,bond_n)/bb;  // cos(\theta_i)

    // energy
    *u += m_lp*(1.0-ci);

    // forces
    gsl_vector_view fnm = gsl_matrix_row(forces, n-1);
    gsl_vector_view fn = gsl_matrix_row(forces, n);
    gsl_vector_view fnp = gsl_matrix_row(forces, n+1);
    // force on n-1
    gsl_vector_set_all(fm,0.0);
    linalg_daxpy(-1./bb,bond_n,fm);
    linalg_daxpy(ci/b2,bond,fm);
    // force on n+1
    gsl_vector_set_all(fp,0.0);
    linalg_daxpy(-ci/bn2,bond_n,fp);
    linalg_daxpy(1./bb,bond,fp);

    // updates
    linalg_daxpy(m_lp, fm, &fnm.vector);
    linalg_daxpy(m_lp, fp, &fnp.vector);
    linalg_daxpy(-m_lp, fm, &fn.vector);
    linalg_daxpy(-m_lp, fp, &fn.vector);
  }

  /* exit */
  gsl_vector_free(bond);
  gsl_vector_free(bond_n);
  gsl_vector_free(fm);
  gsl_vector_free(fp);
  return;
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------

//****************************************************************************
// SoftCore
//****************************************************************************
SoftCore::SoftCore(double A, double sigma, NeighborList *neighbors) :
  m_A(A),
  m_neighbors(neighbors)
{
  double pi = atan(1.)*4;

  m_rc = pow(2., 1./6) * sigma;
  m_y = pi / m_rc;
  m_neighbors->set_range(m_rc);
}

SoftCore::~SoftCore() {
}

double SoftCore::energy_scal(double r){
  /* methods */
  /*
   * compute the energy of the  interaction when the algebraic
   * distance is r.
   */

  if (fabs(r) < m_rc) {
    return m_A*(1. + cos(r*m_y) );
  }
  else {
    return 0.0;
  }
}

double SoftCore::force_scal(double r){
  /*
   * compute the algebraic norm of the force applied when the algebraic
   * distance is r.
   */
  if (fabs(r) < m_rc) {
    return m_A*m_y*sin(r*m_y);
  }
  else {
    return 0.0;
  }
}

void SoftCore::energy_force(gsl_matrix *x, double *u, gsl_matrix *forces){
  /*
   * Compute the potential energy and force (minus gradient) of a soft-core potential.
   * \beta U(r) = A (1 + \cos(\pi r / \xi)),
   * where r is the separation between 2 particles and \xi is the cutoff.
   *
   * The force generated is in the direction of r and has the algebraic norm:
   * f = A \pi / \xi \sin(\pi r / \xi)
   */

  size_t n,m;
  double fnorm, r;
  gsl_vector *xtp(0);

  // initialize
  xtp = gsl_vector_calloc(3);

  for (size_t i=0; i<m_neighbors->m_npair; ++i){
    n = gsl_matrix_uint_get(m_neighbors->m_pairs, i, 0);
    m = gsl_matrix_uint_get(m_neighbors->m_pairs, i, 1);

    gsl_vector_view xn = gsl_matrix_row(x, n);
    gsl_vector_view xm = gsl_matrix_row(x, m);
    gsl_vector_view fn = gsl_matrix_row(forces, n);
    gsl_vector_view fm = gsl_matrix_row(forces, m);

    gsl_vector_memcpy(xtp, &xn.vector);
    linalg_daxpy(-1., &xm.vector, xtp);               // xtp = xn-xm
    r = linalg_dnrm2(xtp);

    // energy
    *u += energy_scal(r);

    // force
    fnorm = force_scal(r);
    linalg_daxpy(fnorm/r, xtp, &fn.vector);
    linalg_daxpy(-fnorm/r, xtp, &fm.vector);
  }

  /* exit */
  gsl_vector_free(xtp);
  return;
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------

//****************************************************************************
// PairLJ
//****************************************************************************
PairLJ::PairLJ(double eps, double sigma, double rc_LJ, NeighborList *neighbors) :
  m_eps(eps),
  m_sigma(sigma),
  m_rc_LJ(rc_LJ),
  m_neighbors(neighbors)

{
  m_4eps = 4.0*m_eps;
  m_48eps = 48.0*m_eps;
  m_fpref = m_48eps / m_sigma;

  double x = 1./m_rc_LJ;
  double x6 = x*x*x*x*x*x;
  double x12 = x6*x6;
  m_u0 = m_4eps*(x12 - x6);

  m_neighbors->set_range(m_rc_LJ*m_sigma);
}

PairLJ::~PairLJ(){
}

double PairLJ::energy_LJ_scal(double r){
  /*
   * compute the energy of the  LJ interaction when the algebraic
   * distance is r.
   */

  double x, x6, x12;

  if (fabs(r) < m_rc_LJ*m_sigma) {
    x = m_sigma/r;
    x6 = x*x*x*x*x*x;
    x12 = x6*x6;
    return m_4eps*(x12 - x6) - m_u0;
  }
  else {
    return 0.0;
  }
}

double PairLJ::force_LJ_scal(double r){
  /*
   * compute the algebraic norm of the LJ force applied when the algebraic
   * distance is r.
   */
  double x, x6, x12;

  if (fabs(r) < m_rc_LJ*m_sigma) {
    x = m_sigma/r;
    x6 = x*x*x*x*x*x;
    x12 = x6*x6;
    return m_fpref*x*(x12 - 0.5*x6);
  }
  else {
    return 0.0;
  }
}

void PairLJ::energy_force(gsl_matrix *x, double *u, gsl_matrix *forces){
  /*
   * Compute the potential energy and force (minus gradient) produced by pair Lennard-Jones interactions.
   */

  size_t n, m;
  double fnorm, r;
  gsl_vector *xtp(0);

  // initialize
  n = 0;
  m = 0;
  r = 0.;
  fnorm = 0.;
  xtp = gsl_vector_calloc(3);

  for (size_t i=0; i<m_neighbors->m_npair; ++i){
    n = gsl_matrix_uint_get(m_neighbors->m_pairs, i, 0);
    m = gsl_matrix_uint_get(m_neighbors->m_pairs, i, 1);

    gsl_vector_view xn = gsl_matrix_row(x, n);
    gsl_vector_view xm = gsl_matrix_row(x, m);
    gsl_vector_view fn = gsl_matrix_row(forces, n);
    gsl_vector_view fm = gsl_matrix_row(forces, m);

    gsl_vector_memcpy(xtp, &xn.vector);
    linalg_daxpy(-1., &xm.vector, xtp);               // xtp = xn-xm
    r = linalg_dnrm2(xtp);

    // energy
    *u += energy_LJ_scal(r);

    // force
    fnorm = force_LJ_scal(r);
    linalg_daxpy(fnorm/r, xtp, &fn.vector);
    linalg_daxpy(-fnorm/r, xtp, &fm.vector);
  }

  /* exit */
  gsl_vector_free(xtp);
  return;
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------

//****************************************************************************
// PolarChains
//****************************************************************************
PolarChains::PolarChains(NeighborList *neighbors, std::vector<std::pair<size_t, size_t> > chain_ends) :
  m_neighbors(neighbors),
  m_chain_ends(chain_ends)
{
  /* every particle must belong to exactly one chain of at least 2 beads */
  vector<pair<size_t, size_t> > chains(m_chain_ends);
  sort(chains.begin(), chains.end());
  size_t npart = 0;
  for (size_t c=0; c<chains.size(); ++c){
    if (chains[c].first != npart)
      throw invalid_argument("Polar force field: chain_ends must cover particles 0..npart-1 without gaps or overlaps.");
    if (chains[c].second <= chains[c].first)
      throw invalid_argument("Polar force field: each chain needs at least 2 beads.");
    npart = chains[c].second + 1;
  }
  if (npart == 0)
    throw invalid_argument("Polar force field: chain_ends is empty.");

  m_pol_vec = gsl_matrix_calloc(npart, 3);
  m_pol_grad = gsl_matrix_calloc(npart, 3);
  m_tnorm.assign(npart, 0.0);
}

PolarChains::~PolarChains()
{
  gsl_matrix_free(m_pol_vec);
  gsl_matrix_free(m_pol_grad);
}

static inline void chain_neighbors(size_t k, size_t istart, size_t iend, size_t &kp, size_t &km){
  /* beads defining the tangent t_k = x_{kp} - x_{km} */
  kp = (k == iend) ? k : k+1;
  km = (k == istart) ? k : k-1;
}

void PolarChains::update_polarity(const gsl_matrix *x){
  /*
   * Compute n_k = ez x t_k / |t_k| for all beads, and reset g_k.
   */
  if (x->size1 != m_pol_vec->size1)
    throw invalid_argument("Polar force field: chain_ends must cover all particles.");

  for (size_t c=0; c<m_chain_ends.size(); ++c){
    size_t istart = m_chain_ends[c].first;
    size_t iend = m_chain_ends[c].second;
    for (size_t k=istart; k<=iend; ++k){
      size_t kp, km;
      chain_neighbors(k, istart, iend, kp, km);
      double tx = gsl_matrix_get(x, kp, 0) - gsl_matrix_get(x, km, 0);
      double ty = gsl_matrix_get(x, kp, 1) - gsl_matrix_get(x, km, 1);
      double tn = sqrt(tx*tx + ty*ty);
      m_tnorm[k] = tn;
      gsl_matrix_set(m_pol_vec, k, 0, -ty/tn);
      gsl_matrix_set(m_pol_vec, k, 1, tx/tn);
      gsl_matrix_set(m_pol_vec, k, 2, 0.0);
    }
  }
  gsl_matrix_set_zero(m_pol_grad);
}

void PolarChains::add_polarity_forces(gsl_matrix *forces){
  /*
   * Add the forces due to the dependence of n_k on x_{k+} and x_{k-}:
   *   F_{k+} += c n_k,  F_{k-} -= c n_k,  c = (t_k/|t_k| . g_k) / |t_k|.
   */
  for (size_t c=0; c<m_chain_ends.size(); ++c){
    size_t istart = m_chain_ends[c].first;
    size_t iend = m_chain_ends[c].second;
    for (size_t k=istart; k<=iend; ++k){
      size_t kp, km;
      chain_neighbors(k, istart, iend, kp, km);
      double nx = gsl_matrix_get(m_pol_vec, k, 0);
      double ny = gsl_matrix_get(m_pol_vec, k, 1);
      // unit tangent: t_k/|t_k| = (ny, -nx)
      double coef = (ny*gsl_matrix_get(m_pol_grad, k, 0) - nx*gsl_matrix_get(m_pol_grad, k, 1)) / m_tnorm[k];
      gsl_matrix_set(forces, kp, 0, gsl_matrix_get(forces, kp, 0) + coef*nx);
      gsl_matrix_set(forces, kp, 1, gsl_matrix_get(forces, kp, 1) + coef*ny);
      gsl_matrix_set(forces, km, 0, gsl_matrix_get(forces, km, 0) - coef*nx);
      gsl_matrix_set(forces, km, 1, gsl_matrix_get(forces, km, 1) - coef*ny);
    }
  }
}

void PolarChains::energy_force(gsl_matrix *x, double *u, gsl_matrix *forces){
  /*
   * Compute the potential energy and force (minus gradient) of the pair interactions.
   */
  gsl_vector *fj = gsl_vector_calloc(3);
  gsl_vector *gi = gsl_vector_calloc(3);
  gsl_vector *gj = gsl_vector_calloc(3);

  update_polarity(x);

  for (size_t p=0; p<m_neighbors->m_npair; ++p)
  {
    size_t n = gsl_matrix_uint_get(m_neighbors->m_pairs, p, 0);
    size_t m = gsl_matrix_uint_get(m_neighbors->m_pairs, p, 1);

    gsl_vector_const_view xn = gsl_matrix_const_row(x, n);
    gsl_vector_const_view xm = gsl_matrix_const_row(x, m);
    gsl_vector_const_view pn = gsl_matrix_const_row(m_pol_vec, n);
    gsl_vector_const_view pm = gsl_matrix_const_row(m_pol_vec, m);
    gsl_vector_view fn = gsl_matrix_row(forces, n);
    gsl_vector_view fm = gsl_matrix_row(forces, m);
    gsl_vector_view gn = gsl_matrix_row(m_pol_grad, n);
    gsl_vector_view gm = gsl_matrix_row(m_pol_grad, m);

    *u += pair_energy_force(&xn.vector, &xm.vector, &pn.vector, &pm.vector, fj, gi, gj);
    gsl_vector_sub(&fn.vector, fj);
    gsl_vector_add(&fm.vector, fj);
    gsl_vector_add(&gn.vector, gi);
    gsl_vector_add(&gm.vector, gj);
  }

  add_polarity_forces(forces);

  gsl_vector_free(fj);
  gsl_vector_free(gi);
  gsl_vector_free(gj);
  return;
}

namespace {
  struct PolarPairGeometry {
    /*
     * Separation and orientation variables of a polar pair:
     *   r = |x_j - x_i|,  e = (x_j - x_i)/r,  d = n_i - n_j,  a = d.e / 2.
     */
    double r, a;
    double e[3], d[3];

    PolarPairGeometry(const gsl_vector *xi, const gsl_vector *xj, const gsl_vector *ni, const gsl_vector *nj){
      r = 0.0;
      a = 0.0;
      for (size_t k=0; k<3; ++k){
        e[k] = gsl_vector_get(xj, k) - gsl_vector_get(xi, k);
        d[k] = gsl_vector_get(ni, k) - gsl_vector_get(nj, k);
        r += e[k]*e[k];
      }
      r = sqrt(r);
      for (size_t k=0; k<3; ++k){
        e[k] /= r;
        a += 0.5*d[k]*e[k];
      }
    }

    void set_gradients(double U_r, double U_a, gsl_vector *fj, gsl_vector *gi, gsl_vector *gj) const {
      /*
       * Given U_r = dU/dr and U_a = dU/da:
       *   fj = -dU/dx_j = -(U_r e + U_a da/dx_j),  da/dx_j = (d - 2a e)/(2r),
       *   gi = dU/dn_i = U_a e/2,  gj = dU/dn_j = -U_a e/2.
       */
      for (size_t k=0; k<3; ++k){
        double da_dxj = 0.5*(d[k] - 2.0*a*e[k])/r;
        gsl_vector_set(fj, k, -(U_r*e[k] + U_a*da_dxj));
        gsl_vector_set(gi, k, 0.5*U_a*e[k]);
        gsl_vector_set(gj, k, -0.5*U_a*e[k]);
      }
    }
  };
}

//****************************************************************************
// PolarPairLJ
//****************************************************************************
PolarPairLJ::PolarPairLJ(double eps, double sigma, double rc_LJ, double alpha, NeighborList *neighbors, std::vector<std::pair<size_t, size_t> > chain_ends) :
  PolarChains(neighbors, chain_ends),
  m_eps(eps),
  m_sigma(sigma),
  m_rc_LJ(rc_LJ),
  m_alpha(alpha)
{
  // largest diameter s over orientations a in [-1, 1]
  m_neighbors->set_range(m_rc_LJ*m_sigma*std::max(1.0, 1.0 + m_alpha));
}

PolarPairLJ::~PolarPairLJ()
{
}

double PolarPairLJ::V_LJ(double u){
  /* u = r / s is the dimensionless separation */
  double u6 = pow(u, -6);
  return 4.0*m_eps*(u6*u6 - u6);
}

double PolarPairLJ::V_LJ_prime(double u){
  /* derivative of V_LJ with respect to u */
  double u6 = pow(u, -6);
  return -48.0*m_eps/u*(u6*u6 - 0.5*u6);
}

double PolarPairLJ::pair_energy_force(const gsl_vector *xi, const gsl_vector *xj,
                                      const gsl_vector *ni, const gsl_vector *nj,
                                      gsl_vector *fj, gsl_vector *gi, gsl_vector *gj){
  /*
   * U = V(r/s) - V(rc),  s = sigma (1 + alpha (1 + a)/2).
   *   dU/dr = V'(u)/s,  dU/da = -V'(u) (u/s) ds/da,  ds/da = sigma alpha/2.
   */
  PolarPairGeometry p(xi, xj, ni, nj);
  double s = m_sigma*(1.0 + 0.5*m_alpha*(1.0 + p.a));
  double u = p.r/s;

  gsl_vector_set_zero(fj);
  gsl_vector_set_zero(gi);
  gsl_vector_set_zero(gj);
  if (u >= m_rc_LJ)
    return 0.0;

  double dV = V_LJ_prime(u);
  p.set_gradients(dV/s, -dV*(u/s)*0.5*m_sigma*m_alpha, fj, gi, gj);
  return V_LJ(u) - V_LJ(m_rc_LJ);
}

//****************************************************************************
// PolarPair48
//****************************************************************************
PolarPair48::PolarPair48(double eps, double sigma, double width, double alpha, NeighborList* neighbors, std::vector<std::pair<size_t, size_t> > chain_ends) :
  PolarChains(neighbors, chain_ends),
  m_eps(eps),
  m_sigma(sigma),
  m_width(width),
  m_alpha(alpha)
{
  if (!(m_sigma > 0.0 && m_width > 0.0))
    throw std::invalid_argument("PolarPair48: sigma and width must be positive.");
  m_r0 = pow(2., 1./6)*m_sigma;
  m_rc = m_r0 + m_width*m_sigma;
  m_neighbors->set_range(m_rc);
}

PolarPair48::~PolarPair48()
{
}

double PolarPair48::pair_energy_force(const gsl_vector *xi, const gsl_vector *xj,
                                      const gsl_vector *ni, const gsl_vector *nj,
                                      gsl_vector *fj, gsl_vector *gi, gsl_vector *gj){
  /*
   * U = eps x^4 (x^4 - 2A),  x = (rc - r)/(rc - r0),  A = 1 - alpha (1 + a).
   *   dU/dr = -8 eps x^3 (x^4 - A)/(rc - r0),  dU/da = 2 alpha eps x^4.
   */
  PolarPairGeometry p(xi, xj, ni, nj);

  gsl_vector_set_zero(fj);
  gsl_vector_set_zero(gi);
  gsl_vector_set_zero(gj);
  if (p.r >= m_rc)
    return 0.0;

  double A = 1.0 - m_alpha*(1.0 + p.a);
  double x = (m_rc - p.r)/(m_rc - m_r0);
  double x3 = x*x*x;
  double x4 = x3*x;
  p.set_gradients(-8.0*m_eps*x3*(x4 - A)/(m_rc - m_r0), 2.0*m_alpha*m_eps*x4, fj, gi, gj);
  return m_eps*x4*(x4 - 2.0*A);
}

//****************************************************************************
// PolarPair48_2site
//****************************************************************************
PolarPair48_2site::PolarPair48_2site(double eps_nn, double eps_ww, double eps_nw,
    double sigma_nn, double sigma_ww, double sigma_nw,
    double width_nn, double width_ww, double width_nw,
    double rho_n, double rho_w,
    NeighborList* neighbors, std::vector<std::pair<size_t, size_t> > chain_ends) :
  PolarChains(neighbors, chain_ends),
  m_eps_nn(eps_nn),
  m_eps_ww(eps_ww),
  m_eps_nw(eps_nw),
  m_sigma_nn(sigma_nn),
  m_sigma_ww(sigma_ww),
  m_sigma_nw(sigma_nw),
  m_width_nn(width_nn),
  m_width_ww(width_ww),
  m_width_nw(width_nw),
  m_rho_n(rho_n),
  m_rho_w(rho_w)
{
  if (!(m_sigma_nn > 0.0 && m_sigma_ww > 0.0 && m_sigma_nw > 0.0
        && m_width_nn > 0.0 && m_width_ww > 0.0 && m_width_nw > 0.0))
    throw std::invalid_argument("PolarPair48_2site: sigma and width must be positive.");
  m_r0_nn = pow(2., 1./6)*m_sigma_nn;
  m_r0_ww = pow(2., 1./6)*m_sigma_ww;
  m_r0_nw = pow(2., 1./6)*m_sigma_nw;
  m_rc_nn = m_r0_nn + m_width_nn*m_sigma_nn;
  m_rc_ww = m_r0_ww + m_width_ww*m_sigma_ww;
  m_rc_nw = m_r0_nw + m_width_nw*m_sigma_nw;

  // sites sit at most |rho| from their bead center
  double rn = fabs(m_rho_n), rw = fabs(m_rho_w);
  m_neighbors->set_range(std::max(std::max(m_rc_nn + 2.0*rn, m_rc_ww + 2.0*rw), m_rc_nw + rn + rw));
}

PolarPair48_2site::~PolarPair48_2site()
{
}

double PolarPair48_2site::potential(double r, double eps, double rc, double r0){
  /*
   * Compute the potential as a function of the distance and the value of the parameters.
   */
  if (r < rc) {
    double x = (rc - r) / (rc - r0);
    double x4 = x*x*x*x;
    return eps * x4 * (x4 - 2.0);
  }
  else {
    return 0.0;
  }
}

double PolarPair48_2site::potential_deriv(double r, double eps, double rc, double r0){
  /*
   * Compute the derivative of the potential as a function of the distance and the value of the parameters.
   */
  if (r < rc) {
    double x = (rc - r) / (rc - r0);
    double x4 = x*x*x*x;
    return - 8.0 * eps * x4 / (rc-r) * (x4 - 1.0);
  }
  else {
    return 0.0;
  }
}

double PolarPair48_2site::pair_energy_force(const gsl_vector *xi, const gsl_vector *xj,
                                            const gsl_vector *ni, const gsl_vector *nj,
                                            gsl_vector *fj, gsl_vector *gi, gsl_vector *gj){
  /*
   * Sum over the 4 site pairs. A site of bead k sits at x_k + c n_k, with
   * c = -rho_n (narrow) or +rho_w (wide). For a site pair with separation
   *   y = (x_j + c_j n_j) - (x_i + c_i n_i)
   * and force f = -U'(|y|) y/|y| on the site of j:
   *   fj += f,  dU/dn_j += -c_j f,  dU/dn_i += c_i f.
   */
  struct SitePair { double eps, rc, r0, ci, cj; };
  const SitePair site_pairs[4] = {
    {m_eps_nn, m_rc_nn, m_r0_nn, -m_rho_n, -m_rho_n},   // narrow-narrow
    {m_eps_ww, m_rc_ww, m_r0_ww, +m_rho_w, +m_rho_w},   // wide-wide
    {m_eps_nw, m_rc_nw, m_r0_nw, -m_rho_n, +m_rho_w},   // narrow (i) - wide (j)
    {m_eps_nw, m_rc_nw, m_r0_nw, +m_rho_w, -m_rho_n},   // wide (i) - narrow (j)
  };

  double energy = 0.0;
  gsl_vector_set_zero(fj);
  gsl_vector_set_zero(gi);
  gsl_vector_set_zero(gj);

  for (size_t s=0; s<4; ++s){
    const SitePair &sp = site_pairs[s];
    double y[3], r = 0.0;
    for (size_t k=0; k<3; ++k){
      y[k] = gsl_vector_get(xj, k) + sp.cj*gsl_vector_get(nj, k)
           - gsl_vector_get(xi, k) - sp.ci*gsl_vector_get(ni, k);
      r += y[k]*y[k];
    }
    r = sqrt(r);
    if (r >= sp.rc)
      continue;

    energy += potential(r, sp.eps, sp.rc, sp.r0);
    double fnorm = -potential_deriv(r, sp.eps, sp.rc, sp.r0)/r;
    for (size_t k=0; k<3; ++k){
      double f = fnorm*y[k];
      gsl_vector_set(fj, k, gsl_vector_get(fj, k) + f);
      gsl_vector_set(gj, k, gsl_vector_get(gj, k) - sp.cj*f);
      gsl_vector_set(gi, k, gsl_vector_get(gi, k) + sp.ci*f);
    }
  }

  return energy;
}
