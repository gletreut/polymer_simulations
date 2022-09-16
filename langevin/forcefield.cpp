#include "forcefield.h"
using namespace std;

//****************************************************************************
// ForceField
//****************************************************************************
ForceField::~ForceField(){
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------



//****************************************************************************
// ConfinmentBox
//****************************************************************************
ConfinmentBox::ConfinmentBox(double xlo, double xhi, double ylo, double yhi, double zlo, double zhi, double sigma, double eps) :
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

ConfinmentBox::~ConfinmentBox(){
}

double ConfinmentBox::energy_LJ_scal(double r){
  /*
   * compute the energy of the  LJ interaction when the algebric
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

double ConfinmentBox::force_LJ_scal(double r){
  /*
   * compute the algebric norm of the LJ force applied when the algebric
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

void ConfinmentBox::energy_force(gsl_matrix *x, double *u, gsl_matrix *forces){
  /*
   * Compute the potential energy and force (minus gradient) produced by a
   * confinment box
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
    fnorm = force_LJ_scal(r);  // algebric norm postive
    linalg_daxpy(fnorm, force, &fi.vector);

    // // test
    // if (fabs(fnorm) > 0) {
    //   cout << setw(10) << setprecision(0) << i;
    //   cout << setw(20) << setprecision(6) << fnorm;
    //   cout << setw(20) << setprecision(6) << rx;
    //   cout << setw(20) << setprecision(6) << m_xlo;
    //   cout << setw(20) << setprecision(6) << r;
    //   cout << endl;
    //   throw invalid_argument("Test");
    // }
    // // test

    /* xhi */
    gsl_vector_set_all(force,0.0);
    gsl_vector_set(force,0,1.0);  // ex
    r = rx - m_xhi;
    *u += energy_LJ_scal(r);
    fnorm = force_LJ_scal(r);  // algebric norm postive
    linalg_daxpy(fnorm, force, &fi.vector);

    /* ylo */
    gsl_vector_set_all(force,0.0);
    gsl_vector_set(force,1,1.0);  // ey
    r = ry - m_ylo;
    *u += energy_LJ_scal(r);
    fnorm = force_LJ_scal(r);  // algebric norm postive
    linalg_daxpy(fnorm, force, &fi.vector);

    /* yhi */
    gsl_vector_set_all(force,0.0);
    gsl_vector_set(force,1,1.0);  // ey
    r = ry - m_yhi;
    *u += energy_LJ_scal(r);
    fnorm = force_LJ_scal(r);  // algebric norm postive
    linalg_daxpy(fnorm, force, &fi.vector);

    /* zlo */
    gsl_vector_set_all(force,0.0);
    gsl_vector_set(force,2,1.0);  // ez
    r = rz - m_zlo;
    *u += energy_LJ_scal(r);
    fnorm = force_LJ_scal(r);  // algebric norm postive
    linalg_daxpy(fnorm, force, &fi.vector);

    /* zhi */
    gsl_vector_set_all(force,0.0);
    gsl_vector_set(force,2,1.0);  // ez
    r = rz - m_zhi;
    *u += energy_LJ_scal(r);
    fnorm = force_LJ_scal(r);  // algebric norm postive
    linalg_daxpy(fnorm, force, &fi.vector);
  }

  /* exit */
  gsl_vector_free(force);
  return;
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------



//****************************************************************************
// ConfinmentSphere
//****************************************************************************
ConfinmentSphere::ConfinmentSphere(double radius, double sigma, double eps) :
  m_radius(radius),
  m_sigma(sigma),
  m_eps(eps)

{
  m_4eps = 4.0*m_eps;
  m_fpref = 48.0*m_eps / m_sigma;
  m_rc_LJ = m_sigma*pow(2.,1./6);
}

ConfinmentSphere::~ConfinmentSphere(){
}

double ConfinmentSphere::energy_LJ_scal(double r){
  /*
   * compute the energy of the  LJ interaction when the algebric
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

double ConfinmentSphere::force_LJ_scal(double r){
  /*
   * compute the algebric norm of the LJ force applied when the algebric
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

void ConfinmentSphere::energy_force(gsl_matrix *x, double *u, gsl_matrix *forces){
  /*
   * Compute the potential energy and force (minus gradient) produced by a
   * confinment box
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
    fnorm = force_LJ_scal(r);  // algebric norm postive
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
   * compute the energy of the  LJ interaction when the algebric
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
   * compute the algebric norm of the LJ force applied when the algebric
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

  // set energy and force to zero
  //*u = 0;
  //gsl_matrix_set_all(forces, 0.0);

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
      //cout << convert.str();
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
    fnorm = force_LJ_scal(r);  // algebric norm
    fnorm /= (r*r);            // because bond is not normalized to 1
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

// GEMField Class -----------------------------------------------------------------------------------------------------------------------------------------------------------
//****************************************************************************
// GEMField
//****************************************************************************
// GEMField::GEMField(size_t offset, size_t N, double b, string filepath) :
//   m_offset(offset),
//   m_N(N),
//   m_b(b)
// {
//   /* declarations */
//   gsl_matrix *K(0);
//   ifstream fin;
//
//   /* initialize the matrix */
//   m_W = gsl_matrix_calloc(m_N,m_N);
//   K = gsl_matrix_calloc(m_N,m_N);
//   gsl_matrix_set_all(m_W,0.0);
//   gsl_matrix_set_all(K,0.0);
//
//   /* initialize the temporary vector */
//   m_fa = gsl_vector_calloc(m_N);
//
//   /* load couplings */
//   fin.open(filepath.c_str());
//   utils::load_matrix(fin, K);
//   fin.close();
//
//   /* compute quadratic matrix of interactions */
//   K2W(K,m_W);
//
//   /* initialize prefactor for force */
//   m_fpref = 3./(m_b*m_b);
//
//   /* exit */
//   gsl_matrix_free(K);
// }
//
// GEMField::~GEMField() {
//   gsl_matrix_free(m_W);
//   gsl_vector_free(m_fa);
// }
//
// /* methods */
// void GEMField::K2W(const gsl_matrix *K, gsl_matrix *W){
//   /*
//    * Compute quadratic matrix of interactions from coupling matrix.
//    */
//   double kij, w;
//
//   gsl_matrix_set_all(W,0.0);
//   for (size_t i=0; i<m_N; ++i){
//     for (size_t j=i; j<m_N; ++j){
//       // coupling for i,j
//       kij = gsl_matrix_get(K,i,j);
//
//       // ri^2 term
//       w = gsl_matrix_get(W,i,i);
//       w += kij;
//       gsl_matrix_set(W,i,i,w);
//
//       // rj^2 term
//       w = gsl_matrix_get(W,j,j);
//       w += kij;
//       gsl_matrix_set(W,j,j,w);
//
//       // ri.rj term
//       w = gsl_matrix_get(W,i,j);
//       w += -kij;
//       gsl_matrix_set(W,i,j,w);
//
//       // rj.ri term
//       w = gsl_matrix_get(W,j,i);
//       w += -kij;
//       gsl_matrix_set(W,j,i,w);
//     }
//   }
//
//   return;
// }
// void GEMField::energy_force(gsl_matrix *x, double *u, gsl_matrix *forces){
//   /*
//    * Compute the potential energy and force (minus gradient) of a GEM potential.
//    *   \beta U(X, Y, Z) = 3/(2b^2) [ X^T W X + Y^T W Y + Z^T W Z],
//    *   where X, Y and Z are vectors of size N.
//    */
//
//   // make sub-matrix corresponding to atoms belonging to the GEM
//   gsl_matrix_view xsub = gsl_matrix_submatrix(x, m_offset, 0, m_N, x->size2);
//   gsl_matrix_view fsub = gsl_matrix_submatrix(forces, m_offset, 0, m_N, x->size2);
//
//   // iterate on the x,y and z coordinates and compute the force and energy
//   // contributions
//   for (size_t d=0; d<3; ++d){
//     gsl_vector_view xa = gsl_matrix_column(&xsub.matrix,d);
//     gsl_vector_view fa = gsl_matrix_column(&fsub.matrix,d);
//
//     // update force
//     linalg_dgemv(0, -m_fpref, m_W, &xa.vector, 0.0, m_fa);  // fa = -3/b^2 W.X
//     linalg_daxpy(1.0, m_fa, &fa.vector);
//
//     // update energy
//     *u += -0.5*linalg_ddot(&xa.vector, m_fa);               // U = 3/(2b^2) X^T W X
//   }
//
//   return;
// }
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
}

SoftCore::~SoftCore() {
}

double SoftCore::energy_scal(double r){
  /* methods */
  /*
   * compute the energy of the  interaction when the algebric
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
   * compute the algebric norm of the force applied when the algebric
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
}

PairLJ::~PairLJ(){
}

double PairLJ::energy_LJ_scal(double r){
  /*
   * compute the energy of the  LJ interaction when the algebric
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
   * compute the algebric norm of the LJ force applied when the algebric
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
// PolarPairLJ
//****************************************************************************
PolarPairLJ::PolarPairLJ(double eps, double sigma, double rc_LJ, double alpha, NeighborList *neighbors, std::vector<std::pair<size_t, size_t> > chain_ends) :
  PairLJ(eps, sigma, rc_LJ, neighbors),
  m_chain_ends(chain_ends),
  m_alpha(alpha)

{
  // m_4eps = 4.0*m_eps;
  // m_fpref = 48.0*m_eps / m_sigma;
  //
  // double x = 1./m_rc_LJ;
  // double x6 = x*x*x*x*x*x;
  // double x12 = x6*x6;
  // m_u0 = m_4eps*(x12 - x6);

  int imax = 0;
  for (vector<pair<size_t, size_t> >::iterator it=m_chain_ends.begin(); it!=m_chain_ends.end(); ++it)
  {
    int iend = it->second;                                   // iend assigned to iterator it that points to second (end of chain)
    if (iend > imax) imax = iend;
  }

  m_pol_vec = gsl_matrix_calloc(imax+1, 3);
}

PolarPairLJ::~PolarPairLJ()
{
  gsl_matrix_free(m_pol_vec);
}

double PolarPairLJ::V_LJ(double rhat){
  /*
   * rhat = r / sigma is the dimensionless separation
   */
  double u, u6, u12;

  u = 1./rhat;
  u6=u*u*u*u*u*u;
  u12=u6*u6;

  return m_4eps*(u12 - u6);
}

double PolarPairLJ::V_LJ_prime(double rhat){
  /*
   * rhat = r / sigma is the dimensionless separation
   */
  double u, u6, u12;

  u = 1./rhat;
  u6=u*u*u*u*u*u;
  u12=u6*u6;

  return -m_48eps*u*(u12 - 0.5*u6);
}

double PolarPairLJ::get_sigma(double phi, double ti, double tj){
  double aij;

  aij = -sin(phi - 0.5*(ti+tj))*sin(0.5*(tj-ti));  // in [-1,1]

  return m_sigma*( 1. + m_alpha * 0.5*(1. + aij) );
}

double PolarPairLJ::get_sigma_dphi(double phi, double ti, double tj){
  double aij_dphi;

  aij_dphi = 0.5*(sin(phi-tj) - sin(phi-ti));

  return m_sigma * m_alpha * 0.5 * aij_dphi;
}

double PolarPairLJ::energy_LJ_scal(const gsl_vector *xi, const gsl_vector *xj, const gsl_vector *ni, const gsl_vector *nj){
  /*
   * compute the energy of the interaction.
   * INPUT:
   *   xij = xj - xi.
   *   ni: polarity vector i
   *   nj: polarity vector j
   */

  double r, phi, ti, tj, sigma, rhat, rhatc, energy;
  gsl_vector *xij(0);

  xij = gsl_vector_calloc(3);

  gsl_vector_memcpy(xij, xj);
  linalg_daxpy(-1.0, xi, xij);
  r = linalg_dnrm2(xij);
  phi = utils::get_angle_2d(xij);
  ti = utils::get_angle_2d(ni);
  tj = utils::get_angle_2d(nj);
  sigma = get_sigma(phi, ti, tj);

  rhat = r / sigma;
  rhatc = m_rc_LJ;

  if (rhat < rhatc) {
    energy =  V_LJ(rhat) - V_LJ(rhatc);
  }
  else {
    energy = 0.0;
  }

  gsl_vector_free(xij);
  return energy;
}

void PolarPairLJ::force_LJ_scal(const gsl_vector *xi, const gsl_vector *xj, const gsl_vector *ni, const gsl_vector *nj, gsl_vector *force){
  /*
   * compute the force applied by i on j.
   * INPUT:
   *   xi: coordinates i
   *   xj: coordinates j
   *   ni: polarity vector i
   *   nj: polarity vector j
   *   force: vector.
   */

  double r, phi, ti, tj, sigma, sigma_dphi, rhat, rhatc, vx, vy, fpref;
  gsl_vector *er(0), *ephi(0);

  er = gsl_vector_calloc(3);
  ephi = gsl_vector_calloc(3);

  // compute separation, vectors
  gsl_vector_memcpy(er, xj);
  linalg_daxpy(-1.0, xi, er);
  r = linalg_dnrm2(er);
  linalg_dscal(1./r, er);
  vx = gsl_vector_get(er, 0);
  vy = gsl_vector_get(er, 1);
  gsl_vector_set(ephi, 0, -vy);
  gsl_vector_set(ephi, 1, vx);

  // compute angles
  phi = utils::get_angle_2d(er);
  ti = utils::get_angle_2d(ni);
  tj = utils::get_angle_2d(nj);

  // construct er, ephi vectors
  // compute sigma, sigma derivative to phi
  sigma = get_sigma(phi, ti, tj);
  sigma_dphi = get_sigma_dphi(phi, ti, tj);

  // other variables
  rhat = r / sigma;
  rhatc = m_rc_LJ; // * sigma / sigma
  fpref = -V_LJ_prime(rhat) / sigma;

  gsl_vector_set_all(force, 0.0);
  if (rhat < rhatc) {
    gsl_vector_memcpy(force, er);
    linalg_daxpy(-sigma_dphi/sigma, ephi, force);
    linalg_dscal(fpref, force);
  }

  gsl_vector_free(er);
  gsl_vector_free(ephi);
  return;
}

void PolarPairLJ::energy_force(gsl_matrix *x, double *u, gsl_matrix *forces){
  /*
   * Compute the potential energy and force (minus gradient) produced by pair Lennard-Jones interactions.
   */

  size_t n, m;
  gsl_vector *xtp(0), *ftp(0);
  gsl_vector_view pn, pm;

  xtp = gsl_vector_calloc(3);
  ftp = gsl_vector_calloc(3);


  //* fixed orientation
  double freq = 2.*3.14158/9.;
  for (size_t n=0; n<m_pol_vec->size1; n++){
    double ux = cos(n*freq);
    double uy = sin(n*freq);
    gsl_matrix_set(m_pol_vec, n, 0, ux);
    gsl_matrix_set(m_pol_vec, n, 1, uy);
  }
  // */

  /* random orientation
  const gsl_rng_type* RDT = gsl_rng_ranlxs1;
  gsl_rng *rng = gsl_rng_alloc(RDT);
  for (size_t n=0; n<m_pol_vec->size1; n++){
    double uvar = gsl_rng_uniform(rng)*2.0*3.14158;
    // double vvar = gsl_rng_uniform(rng);
    // double rvar = sqrt(-2.*log(vvar));
    double rvar = 1.0;
    double ux = rvar * cos(uvar);
    double uy = rvar * sin(uvar);
    gsl_matrix_set(m_pol_vec, n, 0, ux);
    gsl_matrix_set(m_pol_vec, n, 1, uy);
  }
  gsl_rng_free(rng);
  // */

  /* chain polarity vector
  for (vector<pair<size_t, size_t> >::iterator it=m_chain_ends.begin(); it!=m_chain_ends.end(); ++it)
  {
    size_t istart = it->first;                                  // istart assigned to iterator it that points to first (beginning of chain)
    size_t iend = it->second;                                   // iend assigned to iterator it that points to second (end of chain)
    double pn_norm;
    gsl_matrix *rot(0);
    gsl_vector_view x0, x1, x2;

    // rotation matrix around ez axis
    rot = gsl_matrix_calloc(3,3);
    gsl_matrix_set(rot, 0, 0, 0.0);
    gsl_matrix_set(rot, 1, 0, 1.0);
    gsl_matrix_set(rot, 0, 1, -1.0);
    gsl_matrix_set(rot, 1, 1, 0.0);


    // INTERIOR BEADS ([istart + 1] to [iend - 1])
    for (size_t i = istart + 1; i < iend; i++)
    {

      x0 = gsl_matrix_row(x, i-1);
      x2 = gsl_matrix_row(x, i+1);
      pn = gsl_matrix_row(m_pol_vec,i);

      gsl_vector_set_all(xtp, 0.0);
      linalg_daxpy(1.0, &x2.vector, xtp);
      linalg_daxpy(-1.0, &x0.vector, xtp);
      pn_norm = linalg_dnrm2(xtp);
      linalg_dgemv(0, 1.0/pn_norm, rot, xtp, 0.0, &pn.vector);
    }

    // ENDPOINTS
    // CASE 1: START/FIRST BEAD [istart]
    pn = gsl_matrix_row(m_pol_vec, istart);
    pm = gsl_matrix_row(m_pol_vec, istart+1);

    // Row Vectors
    x0 = gsl_matrix_row(x, istart);
    x1 = gsl_matrix_row(x, istart+1);

    gsl_vector_set_all(xtp,0);
    linalg_daxpy(1, &x1.vector, xtp);
    linalg_daxpy(-1, &x0.vector, xtp);
    pn_norm = linalg_dnrm2(xtp);
    linalg_dgemv(0, 1.0/pn_norm, rot, xtp, 0.0, &pn.vector);

    // CASE 2: END/LAST BEAD [iend]
    pn = gsl_matrix_row(m_pol_vec, iend);
    pm = gsl_matrix_row(m_pol_vec, iend-1);

    x1 = gsl_matrix_row(x, iend);
    x0 = gsl_matrix_row(x, iend - 1);

    gsl_vector_set_all(xtp,0);
    linalg_daxpy(1, &x1.vector, xtp);
    linalg_daxpy(-1, &x0.vector, xtp);

    gsl_vector_set_all(xtp,0);
    linalg_daxpy(1, &x1.vector, xtp);
    linalg_daxpy(-1, &x0.vector, xtp);
    pn_norm = linalg_dnrm2(xtp);
    linalg_dgemv(0, 1.0/pn_norm, rot, xtp, 0.0, &pn.vector);

    gsl_matrix_free(rot);
  }
  // */

  /* test: start
  ofstream ftest;
  ftest.open("test.dat", ofstream::app);
  ftest << scientific << setprecision(8) << left;
  // test: end */

  // Loop neighbors
  for (size_t i=0; i<m_neighbors->m_npair; ++i)
  {
    n = gsl_matrix_uint_get(m_neighbors->m_pairs, i, 0);
    m = gsl_matrix_uint_get(m_neighbors->m_pairs, i, 1);
    pn = gsl_matrix_row(m_pol_vec, n);
    pm = gsl_matrix_row(m_pol_vec, m);

    gsl_vector_view xn = gsl_matrix_row(x, n);
    gsl_vector_view xm = gsl_matrix_row(x, m);
    gsl_vector_view fn = gsl_matrix_row(forces, n);
    gsl_vector_view fm = gsl_matrix_row(forces, m);

    // energy
    double utp = energy_LJ_scal(&xn.vector, &xm.vector, &pn.vector, &pm.vector);
    *u += utp;

    gsl_vector_set_all(ftp,0.0);
    force_LJ_scal(&xn.vector, &xm.vector, &pn.vector, &pm.vector, ftp);
    linalg_daxpy(-1.0, ftp, &fn.vector);
    linalg_daxpy(1.0, ftp, &fm.vector);

    /* test: start
    gsl_vector_memcpy(xtp, &xm.vector);
    linalg_daxpy(-1.0, &xn.vector, xtp);
    double r = linalg_dnrm2(xtp);
    double phi_tp = utils::get_angle_2d(xtp);
    double ti_tp = utils::get_angle_2d(&pn.vector);
    double tj_tp = utils::get_angle_2d(&pm.vector);
    double sigma = get_sigma(phi_tp, ti_tp, tj_tp);
    double fnorm = linalg_dnrm2(ftp);
    if (fabs(fnorm) > 1.0e-10){
      ftest << setw(18) << n;
      ftest << setw(18) << m;
      ftest << setw(18) << gsl_vector_get(&xn.vector, 0);
      ftest << setw(18) << gsl_vector_get(&xn.vector, 1);
      ftest << setw(18) << gsl_vector_get(&xm.vector, 0);
      ftest << setw(18) << gsl_vector_get(&xm.vector, 1);
      ftest << setw(18) << gsl_vector_get(&pn.vector, 0);
      ftest << setw(18) << gsl_vector_get(&pn.vector, 1);
      ftest << setw(18) << gsl_vector_get(&pm.vector, 0);
      ftest << setw(18) << gsl_vector_get(&pm.vector, 1);
      ftest << setw(18) << gsl_vector_get(ftp, 0);
      ftest << setw(18) << gsl_vector_get(ftp, 1);
      ftest << setw(18) << gsl_vector_get(ftp, 2);
      ftest << setw(18) << utp;
      ftest << setw(18) << sigma;
      ftest << setw(18) << r;
      ftest << setw(18) << fnorm;
      ftest << endl;
    }
    // test: end */
  }

  /* test: start
  ftest.close();
  // test: end */


  /* exit */
  // These commands will free the previously allocated vector (xtp) and matrix (pol_vec).
  gsl_vector_free(xtp);
  gsl_vector_free(ftp);
  return;
}

//****************************************************************************
// PolarPair48
//****************************************************************************
PolarPair48::PolarPair48(double eps, double sigma, double rc, double alpha, NeighborList* neighbors, std::vector<std::pair<size_t, size_t> > chain_ends) :
  m_eps(eps),
  m_sigma(sigma),
  m_rc(rc),
  m_alpha(alpha),
  m_neighbors(neighbors),
  m_chain_ends(chain_ends)

{
  m_r0 = pow(2., 1./6)*m_sigma;

  int imax = 0;
  for (vector<pair<size_t, size_t> >::iterator it=m_chain_ends.begin(); it!=m_chain_ends.end(); ++it)
  {
    int iend = it->second;                                   // iend assigned to iterator it that points to second (end of chain)
    if (iend > imax) imax = iend;
  }

  m_pol_vec = gsl_matrix_calloc(imax+1, 3);
}

PolarPair48::~PolarPair48()
{
  gsl_matrix_free(m_pol_vec);
}

double PolarPair48::get_A(double phi, double ti, double tj){
  double aij;

  aij = -sin(phi - 0.5*(ti+tj))*sin(0.5*(tj-ti));  // in [-1,1]
  return 1. - m_alpha*(aij + 1.);
}

double PolarPair48::get_A_dphi(double phi, double ti, double tj){
  double aij_dphi;

  aij_dphi = 0.5*(sin(phi-tj) - sin(phi-ti));

  return -m_alpha*aij_dphi;
}

double PolarPair48::energy_LJ_scal(const gsl_vector *xi, const gsl_vector *xj, const gsl_vector *ni, const gsl_vector *nj){
  double r, phi, ti, tj, A, energy, x, x4;
  gsl_vector *xij(0);

  xij = gsl_vector_calloc(3);

  gsl_vector_memcpy(xij, xj);
  linalg_daxpy(-1.0, xi, xij);
  r = linalg_dnrm2(xij);
  phi = utils::get_angle_2d(xij);
  ti = utils::get_angle_2d(ni);
  tj = utils::get_angle_2d(nj);
  A = get_A(phi, ti, tj);

  if (r < m_rc) {
    x = (m_rc - r) / (m_rc - m_r0);
    x4 = x*x*x*x;

    energy = m_eps * x4 * (x4 - 2*A);
  }

  else {
    energy = 0.0;
  }

  gsl_vector_free(xij);
  return energy;

}

void PolarPair48::force_LJ_scal(const gsl_vector *xi, const gsl_vector *xj, const gsl_vector *ni, const gsl_vector *nj, gsl_vector *force){
  /*
   * compute the force applied by i on j.
   * INPUT:
   *   xi: coordinates i
   *   xj: coordinates j
   *   ni: polarity vector i
   *   nj: polarity vector j
   *   force: vector.
   */

  double r, phi, ti, tj, A, A_dphi, vx, vy, x, x4, fpref;
  gsl_vector *er(0), *ephi(0);

  er = gsl_vector_calloc(3);
  ephi = gsl_vector_calloc(3);

  // compute separation, vectors
  gsl_vector_memcpy(er, xj);
  linalg_daxpy(-1.0, xi, er);
  r = linalg_dnrm2(er);
  linalg_dscal(1./r, er);
  vx = gsl_vector_get(er, 0);
  vy = gsl_vector_get(er, 1);
  gsl_vector_set(ephi, 0, -vy);
  gsl_vector_set(ephi, 1, vx);

  // compute angles
  phi = utils::get_angle_2d(er);
  ti = utils::get_angle_2d(ni);
  tj = utils::get_angle_2d(nj);

  // construct er, ephi vectors
  // compute sigma, sigma derivative to phi
  A = get_A(phi, ti, tj);
  A_dphi = get_A_dphi(phi, ti, tj);

  // other variables
  x = (m_rc - r) / (m_rc - m_r0);
  x4 = x*x*x*x;
  fpref = 8*m_eps/(m_rc - r)*x4;

  gsl_vector_set_all(force, 0.0);
  if (r < m_rc) {
    linalg_daxpy(fpref*(x4 - A), er, force);
    linalg_daxpy(A_dphi*2.*m_eps*x4/r, ephi, force);
  }

  gsl_vector_free(er);
  gsl_vector_free(ephi);
  return;
}

void PolarPair48::energy_force(gsl_matrix *x, double *u, gsl_matrix *forces){
  /*
   * Compute the potential energy and force (minus gradient) produced by pair Lennard-Jones interactions.
   */

  size_t n, m;
  gsl_vector *xtp(0), *ftp(0);
  gsl_vector_view pn, pm;

  xtp = gsl_vector_calloc(3);
  ftp = gsl_vector_calloc(3);


  /* fixed orientation
  double freq = 2.*3.14158/9.;
  for (size_t n=0; n<m_pol_vec->size1; n++){
    double ux = cos(n*freq);
    double uy = sin(n*freq);
    gsl_matrix_set(m_pol_vec, n, 0, ux);
    gsl_matrix_set(m_pol_vec, n, 1, uy);
  }
  // */

  /* random orientation
  const gsl_rng_type* RDT = gsl_rng_ranlxs1;
  gsl_rng *rng = gsl_rng_alloc(RDT);
  for (size_t n=0; n<m_pol_vec->size1; n++){
    double uvar = gsl_rng_uniform(rng)*2.0*3.14158;
    // double vvar = gsl_rng_uniform(rng);
    // double rvar = sqrt(-2.*log(vvar));
    double rvar = 1.0;
    double ux = rvar * cos(uvar);
    double uy = rvar * sin(uvar);
    gsl_matrix_set(m_pol_vec, n, 0, ux);
    gsl_matrix_set(m_pol_vec, n, 1, uy);
  }
  gsl_rng_free(rng);
  // */

  //* chain polarity vector
  for (vector<pair<size_t, size_t> >::iterator it=m_chain_ends.begin(); it!=m_chain_ends.end(); ++it)
  {
    size_t istart = it->first;                                  // istart assigned to iterator it that points to first (beginning of chain)
    size_t iend = it->second;                                   // iend assigned to iterator it that points to second (end of chain)
    double pn_norm;
    gsl_matrix *rot(0);
    gsl_vector_view x0, x1, x2;

    // rotation matrix around ez axis
    rot = gsl_matrix_calloc(3,3);
    gsl_matrix_set(rot, 0, 0, 0.0);
    gsl_matrix_set(rot, 1, 0, 1.0);
    gsl_matrix_set(rot, 0, 1, -1.0);
    gsl_matrix_set(rot, 1, 1, 0.0);


    // INTERIOR BEADS ([istart + 1] to [iend - 1])
    for (size_t i = istart + 1; i < iend; i++)
    {

      x0 = gsl_matrix_row(x, i-1);
      x2 = gsl_matrix_row(x, i+1);
      pn = gsl_matrix_row(m_pol_vec,i);

      gsl_vector_set_all(xtp, 0.0);
      linalg_daxpy(1.0, &x2.vector, xtp);
      linalg_daxpy(-1.0, &x0.vector, xtp);
      pn_norm = linalg_dnrm2(xtp);
      linalg_dgemv(0, 1.0/pn_norm, rot, xtp, 0.0, &pn.vector);
    }

    // ENDPOINTS
    // CASE 1: START/FIRST BEAD [istart]
    pn = gsl_matrix_row(m_pol_vec, istart);
    pm = gsl_matrix_row(m_pol_vec, istart+1);

    // Row Vectors
    x0 = gsl_matrix_row(x, istart);
    x1 = gsl_matrix_row(x, istart+1);

    gsl_vector_set_all(xtp,0);
    linalg_daxpy(1, &x1.vector, xtp);
    linalg_daxpy(-1, &x0.vector, xtp);
    pn_norm = linalg_dnrm2(xtp);
    linalg_dgemv(0, 1.0/pn_norm, rot, xtp, 0.0, &pn.vector);

    // CASE 2: END/LAST BEAD [iend]
    pn = gsl_matrix_row(m_pol_vec, iend);
    pm = gsl_matrix_row(m_pol_vec, iend-1);

    x1 = gsl_matrix_row(x, iend);
    x0 = gsl_matrix_row(x, iend - 1);

    gsl_vector_set_all(xtp,0);
    linalg_daxpy(1, &x1.vector, xtp);
    linalg_daxpy(-1, &x0.vector, xtp);

    gsl_vector_set_all(xtp,0);
    linalg_daxpy(1, &x1.vector, xtp);
    linalg_daxpy(-1, &x0.vector, xtp);
    pn_norm = linalg_dnrm2(xtp);
    linalg_dgemv(0, 1.0/pn_norm, rot, xtp, 0.0, &pn.vector);

    gsl_matrix_free(rot);
  }
  //*/

  /* test: start
  ofstream ftest;
  ftest.open("test.dat", ofstream::app);
  ftest << scientific << setprecision(8) << left;
  // test: end */

  // Loop neighbors
  for (size_t i=0; i<m_neighbors->m_npair; ++i)
  {
    n = gsl_matrix_uint_get(m_neighbors->m_pairs, i, 0);
    m = gsl_matrix_uint_get(m_neighbors->m_pairs, i, 1);
    pn = gsl_matrix_row(m_pol_vec, n);
    pm = gsl_matrix_row(m_pol_vec, m);

    gsl_vector_view xn = gsl_matrix_row(x, n);
    gsl_vector_view xm = gsl_matrix_row(x, m);
    gsl_vector_view fn = gsl_matrix_row(forces, n);
    gsl_vector_view fm = gsl_matrix_row(forces, m);

    // energy
    double utp = energy_LJ_scal(&xn.vector, &xm.vector, &pn.vector, &pm.vector);
    *u += utp;

    gsl_vector_set_all(ftp,0.0);
    force_LJ_scal(&xn.vector, &xm.vector, &pn.vector, &pm.vector, ftp);
    linalg_daxpy(-1.0, ftp, &fn.vector);
    linalg_daxpy(1.0, ftp, &fm.vector);

    /* test: start
    gsl_vector_memcpy(xtp, &xm.vector);
    linalg_daxpy(-1.0, &xn.vector, xtp);
    double r = linalg_dnrm2(xtp);
    double phi_tp = utils::get_angle_2d(xtp);
    double ti_tp = utils::get_angle_2d(&pn.vector);
    double tj_tp = utils::get_angle_2d(&pm.vector);
    double A = get_A(phi_tp, ti_tp, tj_tp);
    double fnorm = linalg_dnrm2(ftp);
    if (fabs(fnorm) > 1.0e-10){
      ftest << setw(18) << n;
      ftest << setw(18) << m;
      ftest << setw(18) << gsl_vector_get(&xn.vector, 0);
      ftest << setw(18) << gsl_vector_get(&xn.vector, 1);
      ftest << setw(18) << gsl_vector_get(&xm.vector, 0);
      ftest << setw(18) << gsl_vector_get(&xm.vector, 1);
      ftest << setw(18) << gsl_vector_get(&pn.vector, 0);
      ftest << setw(18) << gsl_vector_get(&pn.vector, 1);
      ftest << setw(18) << gsl_vector_get(&pm.vector, 0);
      ftest << setw(18) << gsl_vector_get(&pm.vector, 1);
      ftest << setw(18) << gsl_vector_get(ftp, 0);
      ftest << setw(18) << gsl_vector_get(ftp, 1);
      ftest << setw(18) << gsl_vector_get(ftp, 2);
      ftest << setw(18) << utp;
      ftest << setw(18) << A;
      ftest << setw(18) << r;
      ftest << setw(18) << fnorm;
      ftest << endl;
    }
    // test: end */
  }

  /* test: start
  ftest.close();
  // test: end */


  /* exit */
  // These commands will free the previously allocated vector (xtp) and matrix (pol_vec).
  gsl_vector_free(xtp);
  gsl_vector_free(ftp);
  return;
}
