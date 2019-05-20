#ifndef _STEPPER_H
#define _STEPPER_H
//----------------------------------------------------------------------------
//  2019-03-29
//  G. Le Treut - Jun Lab, UC San Diego.
//
//  file: stepper.h
//----------------------------------------------------------------------------
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdexcept>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>

#include "linalg.h"
#include "utils.h"
#include "model.h"

class MDStepper {
  public:
    /* attributes */
    double m_dt;
    double m_gam;
    double m_mass_i;
    MDWorld *m_world;
    gsl_matrix *m_x, *m_v, *m_f;

    /* constructor / destructor */
    MDStepper(MDWorld *world, double dt, double gam=1.0);
    virtual ~MDStepper();

    /* methods */
    void call_forces();
    void call_constraints_x();
    void call_constraints_v();
    virtual void step()=0;
};

class MDStepper_VVerlet: public MDStepper {
  public:
    /* attributes */
    double m_dt_half;
    double m_gam_m;
    double m_gam_i;
    double m_pref_force;

    /* constructor / destructor */
    MDStepper_VVerlet(MDWorld *world, double dt, double gam=0.0);
    ~MDStepper_VVerlet();

    /* methods */
    void step();
};

class LangStepper_VVerlet: public MDStepper_VVerlet {
  public:
    /* attributes */
    double m_fbd_norm;
    double m_s12;
    gsl_rng *m_rng;

    /* constructor / destructor */
    LangStepper_VVerlet(MDWorld *world, gsl_rng *rng, double dt, double gam=1.0);
    ~LangStepper_VVerlet();

    /* methods */
    void add_brownian_forces();
    void step();
};

#endif
