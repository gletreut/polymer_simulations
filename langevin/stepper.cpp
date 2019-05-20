//----------------------------------------------------------------------------
//  2019-03-29
//  G. Le Treut - Jun Lab, UC San Diego.
//
//  file: stepper.cpp
//----------------------------------------------------------------------------
#include "stepper.h"

using namespace std;

//****************************************************************************
// MDStepper
//****************************************************************************
MDStepper::MDStepper(MDWorld *world, double dt, double gam) :
  m_world(world), m_dt(dt), m_gam(gam) {

  // copy the reference to world positions, velocities and forces vector
  m_x = m_world->m_x;
  m_v = m_world->m_v;
  m_f = m_world->m_forces;

  // other variables to be used
  m_mass_i = 1./m_world->m_mass;
}

MDStepper::~MDStepper(){
}

void MDStepper::call_forces(){
  /*
   * Force computation
   */
  m_world->update_energy_forces();    // update the forces
  return;
}

void MDStepper::call_constraints_x(){
  /*
   * Apply constraint projection for positions
   */
  m_world->constrain_x();        // apply constraints
  return;
}

void MDStepper::call_constraints_v(){
  /*
   * Apply constraint projection for velocities
   */
  m_world->constrain_v();        // apply constraints
  return;
}
//****************************************************************************
// MDStepper_VVerlet
//****************************************************************************
MDStepper_VVerlet::MDStepper_VVerlet(MDWorld *world, double dt, double gam) :
  MDStepper::MDStepper(world, dt, gam) {

  m_dt_half = 0.5*m_dt;
  m_gam_m = (1. - m_gam*m_dt_half);           // 1.-\gamma x dt / 2
  m_gam_i = 1./(1. + m_gam*m_dt_half);        // 1./(1.+\gamma x dt / 2)
  m_pref_force = m_mass_i*m_dt_half;    // 1/M . dt/2
}

MDStepper_VVerlet::~MDStepper_VVerlet(){
}

void MDStepper_VVerlet::step() {
  /*
   * This method uses the Velocity-Verlet algorithm to step the dynamics.
   * Reference: Tamar Schlick. Molecular Modeling and Simulation. Springer.
   * p437.
   */

  // the algorithm assume that the forces are updated with respect to the
  // current positions
  // update velocities -- first half-step
  linalg_dscal(m_gam_m, m_v);             // V <- (1 - \gamma dt/2) V
  linalg_daxpy(m_pref_force, m_f, m_v);   // V <- V + 1/M dt/2 * F
  call_constraints_v();
  // update positions
  linalg_daxpy(m_dt, m_v, m_x);           // X <- X + dt * V
  call_constraints_x();
  call_forces();                          // F <- forces(X)
  // update velocities -- second half-step
  linalg_daxpy(m_pref_force, m_f, m_v);   // V <- V + 1/M dt/2 * F
  linalg_dscal(m_gam_i, m_v);             // V <- 1/(1 + \gamma dt/2) V
  call_constraints_v();

  return;
}

//****************************************************************************
// LangStepper_VVerlet
//****************************************************************************
LangStepper_VVerlet::LangStepper_VVerlet(MDWorld *world, gsl_rng *rng, double dt, double gam) :
  MDStepper_VVerlet::MDStepper_VVerlet(world, dt, gam), m_rng(rng) {

  /* initialize some parameters */

  m_s12 = sqrt(12.);
  m_fbd_norm = m_mass_i*sqrt(2.*world->m_mass*m_gam*world->m_temp) * sqrt(m_dt_half);
}

LangStepper_VVerlet::~LangStepper_VVerlet(){
}

void LangStepper_VVerlet::add_brownian_forces(){
  /*
   * Compute and brownian forces to velocities
   */

  // declaration
  double r;
  gsl_vector *rforce(0);

  // initialize
  rforce = gsl_vector_calloc(3);

  // iterate over atomes
  for (size_t i=0; i<m_v->size1; ++i){
    // compute random force
    for (size_t a=0; a<3; ++a){
      r = m_s12*(gsl_rng_uniform(m_rng) - 0.5); // uniform with mean zero and variance 1.
      gsl_vector_set(rforce,a,r);
    }

    // update forces
    gsl_vector_view vi = gsl_matrix_row(m_v, i);
    linalg_daxpy(m_fbd_norm, rforce, &vi.vector);
  }

  /* exit */
  gsl_vector_free(rforce);

  return;
}

void LangStepper_VVerlet::step() {
  /*
   * This method uses the Velocity-Verlet algorithm to step the dynamics.
   * Reference: Tamar Schlick. Molecular Modeling and Simulation. Springer.
   * p437.
   */

  // the algorithm assume that the forces are updated with respect to the
  // current positions
  // update velocities -- first half-step
  linalg_dscal(m_gam_m, m_v);             // V <- (1 - \gamma dt/2) V
  linalg_daxpy(m_pref_force, m_f, m_v);   // V <- V + 1/M dt/2 * F
  add_brownian_forces();
  call_constraints_v();
  // update positions
  linalg_daxpy(m_dt, m_v, m_x);           // X <- X + dt * V
  call_constraints_x();
  call_forces();                          // F <- forces(X)
  // update velocities -- second half-step
  linalg_daxpy(m_pref_force, m_f, m_v);   // V <- V + 1/M dt/2 * F
  add_brownian_forces();
  linalg_dscal(m_gam_i, m_v);             // V <- 1/(1 + \gamma dt/2) V
  call_constraints_v();

  return;
}

