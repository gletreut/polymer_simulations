#include "forcefield.h"

using namespace std;

//****************************************************************************
// ForceField
//****************************************************************************
ForceField::~ForceField(){
}

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
  double rx,ry,rz,r,sint,cosp,sinp,fnorm;
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

//****************************************************************************
// PolymerGaussian
//****************************************************************************
PolymerGaussian::PolymerGaussian(size_t offset, size_t N, double b) :
  m_offset(offset),
  m_N(N),
  m_b(b)
{
  /* initialize some parameters */
  m_ke = 1.5/(m_b*m_b);
  m_fpref = 2.*m_ke;
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
    //  FENE
    gsl_vector_view fn = gsl_matrix_row(forces, n);
    gsl_vector_view fnp = gsl_matrix_row(forces, n+1);
    linalg_daxpy(m_fpref, bond, &fn.vector);
    linalg_daxpy(-m_fpref, bond, &fnp.vector);
  }

  /* exit */
  gsl_vector_free(bond);
  return;
}

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

PolymerFENE::~PolymerFENE(){
}

/* methods */
void PolymerFENE::energy_force(gsl_matrix *x, double *u, gsl_matrix *forces){
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
    gsl_vector_view fn = gsl_matrix_row(forces, n);
    gsl_vector_view fnp = gsl_matrix_row(forces, n+1);
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

//****************************************************************************
// LangThermostat
//****************************************************************************
LangThermostat::LangThermostat(double mass, double gamma, double temp, gsl_rng *rng) :
  m_mass(mass),
  m_gamma(gamma),
  m_temp(temp),
  m_rng(rng)
{
  /* initialize some parameters */
  m_fpref = sqrt(2.*m_mass*m_gamma*m_temp);
  m_s12 = sqrt(12.);
}

LangThermostat::~LangThermostat(){
}

void LangThermostat::energy_force(gsl_matrix *x, double *u, gsl_matrix *forces){
  /*
   * Compute the random force in a thermal bath.
   * Velocity Langevin.
   */

  double r;
  gsl_vector *rforce(0);

  // initialize
  rforce = gsl_vector_calloc(3);

  // iterate over atomes
  for (size_t i=0; i<x->size1; ++i){
    // compute random force
    for (size_t a=0; a<3; ++a){
      r = m_s12*(gsl_rng_uniform(m_rng) - 0.5); // uniform with mean zero and variance 1.
      gsl_vector_set(rforce,a,r);
    }

    cout << "fpref = " << m_fpref << endl;
    // update forces
    gsl_vector_view fi = gsl_matrix_row(forces, i);
    linalg_daxpy(m_fpref, rforce, &fi.vector);
  }

  /* exit */
  gsl_vector_free(rforce);
  return;
}
