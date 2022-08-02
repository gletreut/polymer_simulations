#include "forcefield.h"
using namespace std;

//****************************************************************************
// (1) ForceField
//****************************************************************************
ForceField::~ForceField(){
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------



//****************************************************************************
// (2) ConfinmentBox 
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
// (3) ConfinmentSphere
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
// (4) PolymerGaussian
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
// (5) PolymerHarmonic 
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
// (6) PolymerFENE
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
// (7) PolymerKratkyPorod 
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
// (8) SoftCore
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
// (9) PairLJ
//****************************************************************************
PairLJ::PairLJ(double eps, double sigma, double rc_LJ, NeighborList *neighbors) :
  m_eps(eps),
  m_sigma(sigma),
  m_rc_LJ(rc_LJ),
  m_neighbors(neighbors)

{
  m_4eps = 4.0*m_eps;
  m_fpref = 48.0*m_eps / m_sigma;

  double x = m_sigma/m_rc_LJ;
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

  if (fabs(r) < m_rc_LJ) {
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
// (10) PolarPairLJ
//****************************************************************************
PolarPairLJ::PolarPairLJ(double eps, double sigma, double rc_LJ, NeighborList *neighbors, std::vector<std::pair<size_t, size_t> > chain_ends) :
  PairLJ(eps, sigma, rc_LJ, neighbors),
  m_chain_ends(chain_ends)

{
}

PolarPairLJ::~PolarPairLJ()
{
}

double PolarPairLJ::energy_LJ_scal(double r, double sigma){
  /*
   * compute the energy of the  LJ interaction when the algebric
   * distance is r.
   */

  double x, x6, x12;

  if (fabs(r) < m_rc_LJ) {
    x = sigma/r;
    x6 = x*x*x*x*x*x;
    x12 = x6*x6;
    return m_4eps*(x12 - x6) - m_u0;
  }
  else {
    return 0.0;
  }
}

double PolarPairLJ::force_LJ_scal(double r, double sigma){
  /*
   * compute the algebric norm of the LJ force applied when the algebric
   * distance is r.
   */
  double x, x6, x12;

  if (fabs(r) < m_rc_LJ) {
    x = sigma/r;
    x6 = x*x*x*x*x*x;
    x12 = x6*x6;
    return m_fpref*x*(x12 - 0.5*x6);
  }
  else {
    return 0.0;
  }
}

void PolarPairLJ::energy_force(gsl_matrix *x, double *u, gsl_matrix *forces){
  /*
   * Compute the potential energy and force (minus gradient) produced by pair Lennard-Jones interactions.
   */
  
  size_t n, m;
  double fnorm, r, sigma;
  // MAYBE Add alpha variable?

  // pol_vec and xtp INITIALIZED - Done through a pointer variable that holds the memory address
  gsl_vector *xtp(0);
  gsl_matrix *pol_vec(0);

  // initialize
  n = 0;
  m = 0;
  r = 0.;
  fnorm = 0.;
  xtp = gsl_vector_calloc(3);

  // pol_vec - represented by m sub i
  pol_vec = gsl_matrix_calloc(x->size1, 3);

//____________________________________________________________________________
  // -------------- NESTED for loop -------------- //
  // Establishes and calculates [istart] and [iend]
  for (vector<pair<size_t, size_t> >::iterator it=m_chain_ends.begin(); it!=m_chain_ends.end(); ++it) 
  {
    size_t istart = it->first;
    size_t iend = it->first;
    // istart vectors
    gsl_vector_view pn;      // pn = m sub (0) = m
    gsl_vector_view pnn;     // pnn = m sub (1) = m'
    // iend vectors 
    gsl_vector_view pN;      // pN = m sub (N)
    gsl_vector_view pNN;     // pNN = m sub (N-1) 

    // SECOND for loop Calculate the Interior Beads ONLY ([istart + 1] to [iend - 1])
    for (size_t i = istart + 1; i < iend; i++)
    {
      gsl_vector_view x0 = gsl_matrix_row(x, i-1);  // Sets x0 = r sub (i-1)
      gsl_vector_view x1 = gsl_matrix_row(x, i);    // Sets x1 = r sub (i)
      gsl_vector_view x2 = gsl_matrix_row(x, i+1);  // Sets x2 = r sub (i+1)

      gsl_vector_set_all(xtp,0);          // xtp = (0,0,0) - Initiates xtp to 0 in 3D
                                          // NOW --> xtp = 0
      linalg_daxpy(0.5, &x0.vector, xtp); // xtp <-- (1/2)x0 + xtp = (1/2)(r sub (i-1))
                                          // NOW --> xtp = (1/2)(r sub (i-1))
      linalg_daxpy(0.5, &x2.vector, xtp); // xtp <-- (1/2)x2 + xtp = (1/2)(r sub (i+1)) + (1/2)(r sub (i-1))
                                          // NOW --> xtp = (1/2)(r sub (i+1)) + (1/2)(r sub (i-1))

      pn = gsl_matrix_row(pol_vec,i);                      // Sets pn = m sub(i) from pol_vec
      gsl_vector_memcpy(&pn.vector, &x1.vector);          //  m sub (i) = r sub (i)
      linalg_daxpy(-1,xtp,&pn.vector);                   //  m sub (i) = r sub (i) - xtp =  r sub (i) - ((r sub (i+1) + r sub (i-1))/2)
      double pn_norm = linalg_dnrm2(&pn.vector);        // || m sub (i) || = || r sub (i) - ((r sub (i+1) + r sub (i-1))/2) ||
      linalg_dscal((1/pn_norm), &pn.vector);           // (1 / || m sub (i) ||) * (m sub (i))
    }

    // ENDPOINTS (ENDBEADS) CALCULATION
    // ________________CASE 1: START/FIRST BEAD [istart]________________
    pn = gsl_matrix_row(pol_vec, istart);                       // pn = Vector m sub (0)
    pnn = gsl_matrix_row(pol_vec, istart+1);                    // pnn = Vector m sub (1)

    gsl_vector_view xstart = gsl_matrix_row(x, istart);         // Sets xstart = r sub (0)
    gsl_vector_view x1 = gsl_matrix_row(x, istart+1);           // Sets x1 = r sub (1), where 1 = i = istart + 1

    /* THE LINES in (1047-1048) are now in (1040-1041)
    //pn = gsl_matrix_row(pol_vec, istart);                     // pn = Vector m sub (0)
    //pnn = gsl_matrix_row(pol_vec, istart+1);                  // pnn = Vector m sub (1)
    */

    gsl_vector_set_all(xtp,0);                                  // xtp = (0,0,0)
    linalg_daxpy(1, &x1.vector, xtp);                           // xtp = r sub (1)
    linalg_daxpy(-1, &xstart.vector, xtp);                      // xtp = r sub (1) - r sub (0)
    double xtp_normS = linalg_dnrm2(xtp);                       // (xtp_normS) = || u sub (1) || = || r sub (1) - r sub (0) ||
    linalg_dscal((1/xtp_normS), xtp);                           // (1 / || u sub (1) ||) * (u sub (1)) --> = m sub (0) HAT


    // We have vector ( a ) with positions ( 0 )
    //                ( b )                ( 1 )
    double ux = gsl_vector_get(xtp,0);                          // a = (ux) for coordinate position (0)  
    double uy = gsl_vector_get(xtp,1);                          // b = (uy) for coordinate position (1)
    gsl_vector_set(&pn.vector,0,-uy);                           // For [pn.vector (m sub 0)], position (0) becomes (-uy) 
    gsl_vector_set(&pn.vector,1,ux);                            //   "   "   "                position (1) becomes (ux)  
                                                                // s Number   
    double s = linalg_ddot(&pn.vector, &pnn.vector);            // Number (s) = Dot product of (pn.vector) and (pnn.vector)
    int sgn = s/fabs(s);                                        // Number (sgn) = Number (s) divided by the Absolute Value of (s)
    double pn_norm = linalg_dnrm2(&pn.vector);                  // Number (pn_norm) = Magnitude of (pn.vector) = ||pn.vector|| 
    gsl_vector_memcpy(xtp,&pn.vector);                          // Copies the vector (xtp) into the vector (pn.vector)
    gsl_vector_set_all(&pn.vector,0);                           // Sets all elements of pn.vector to 0
    linalg_daxpy(sgn * 1/pn_norm, xtp, &pn.vector);             // Normailized Vector/ Unit Vector of [[m sub (0) HAT]]


    // ________________CASE 2: END/LAST BEAD [iend]________________
    pN = gsl_matrix_row(pol_vec, iend);                         // pN = Vector m sub (N)
    pNN = gsl_matrix_row(pol_vec, iend-1);                      // pNN = Vector m sub (N-1)

    gsl_vector_view xend = gsl_matrix_row(x, iend);             // Sets xend = r sub (N)
    gsl_vector_view xpenultimate = gsl_matrix_row(x, iend - 1); // Sets xpenultimate = r sub (N-1), where 1 = iend - 1 

    gsl_vector_set_all(xtp,0);                                  // xtp = (0,0,0)
    linalg_daxpy(1, &xend.vector, xtp);                         // xtp = r sub (N)
    linalg_daxpy(-1, &xpenultimate.vector, xtp);                // xtp = r sub (N) - r sub (N-1)
    double xtp_normE = linalg_dnrm2(xtp);                       // || u sub (N) || = || r sub (N) - r sub (N-1) ||
    linalg_dscal((1/xtp_normE), xtp);                           // (1 / || u sub (N) ||) * (m sub (N)) --> = m sub (N) HAT

    // We have vector ( a ) with positions ( 0 )
    //                ( b )                ( 1 )
    // N added to represent variables associated with END/LAST BEAD
    double uNx = gsl_vector_get(xtp,0);                         // a = (ux) for coordinate position (0)  
    double uNy = gsl_vector_get(xtp,1);                         // b = (uy) for coordinate position (1)
    gsl_vector_set(&pN.vector,0,-uNy);                          // For [pN.vector (m sub (N) HAT)], position (0) becomes (-uNy) 
    gsl_vector_set(&pN.vector,1,uNx);                           //   "   "   "                  position (1) becomes (uNx)     
                                                                // t Number
    double t = linalg_ddot(&pN.vector, &pNN.vector);            // Number (t) = Dot product of (pN.vector) and (pNN.vector)
    int tgn = t/fabs(t);                                        // Number (tgn) = Number (t) divided by the Absolute Value of (t)
    double pN_norm = linalg_dnrm2(&pN.vector);                  // Number (pN_norm) = Magnitude of (pN.vector) = ||pN.vector|| 
    gsl_vector_memcpy(xtp,&pN.vector);                          // Copies the vector (xtp) into the vector (pN.vector)
    gsl_vector_set_all(&pN.vector,0);                           // Sets all elements of pN.vector to 0
    linalg_daxpy(tgn * 1/pN_norm, xtp, &pN.vector);             // Normailized Vector/ Unit Vector of [[m sub (N) HAT]]

     //----------------------------------------------------------------------
    /* ASIDE FOR LATER
    Rotation Matrix 
    m sub (0)[0] = -u sub (1) [1]
    m sub (0)[1] = u sub (1) [1]
    m sub (0)[2] = 0
    or 
    rotate() function // the reverse order
    or 
    std::rotate(ObjectToRotate.begin(),
      ObjectToRotate.end()-1, // the reverse order
      ObjectToRotate.end());
    //*/
    //----------------------------------------------------------------------
  }
//____________________________________________________________________________

  // Loop neighbors
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

    // retrieve the polarity vector

    // Define pn, pm, pol
    
    /*
    pn = pol_vec[n]; 
    pm = pol_vec[m];
    pol = pm - pn;
    // xtp = theta sub i = (alpha) * (((1+cos(theta sub i))/2)+1 ) * (theta/2)
    // pol = theta sub j = (alpha) * (((1-cos(theta sub j))/2)+1 ) * (theta/2)
    sigma = m_sigma*(1 + (alpha * f(xtp, pol)));
    //*/

    // energy
    // *u += energy_LJ_scal(r);
    sigma = m_sigma; // sigma
    *u += energy_LJ_scal(r, sigma);

    // force
    // fnorm = force_LJ_scal(r);
    fnorm = force_LJ_scal(r, sigma);
    linalg_daxpy(fnorm/r, xtp, &fn.vector);
    linalg_daxpy(-fnorm/r, xtp, &fm.vector);
  }

  /* exit */
  gsl_vector_free(xtp);
  gsl_matrix_free(pol_vec);
  return;
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------



// END


