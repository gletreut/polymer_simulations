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
   * Initiate the stepper
   */
  m_world->update_energy_forces();    // update the forces
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
  // update positions
  linalg_daxpy(m_dt, m_v, m_x);           // X <- X + dt * V
  call_forces();                          // F <- forces(X)
  // update velocities -- second half-step
  linalg_daxpy(m_pref_force, m_f, m_v);   // V <- V + 1/M dt/2 * F
  linalg_dscal(m_gam_i, m_v);             // V <- 1/(1 + \gamma dt/2) V

  return;
}

MDStepper_VVerlet::~MDStepper_VVerlet(){
}

