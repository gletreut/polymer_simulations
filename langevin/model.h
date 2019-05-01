//*******************************************************************************
//*
//* Langevin Dynamics
//*	model.h
//*	Definition of model.h
//*
//* Author: Guillaume Le Treut
//*	UCSD - 2019
//*
//*******************************************************************************

#ifndef _MODEL_H
#define _MODEL_H

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <limits>
#include <stdexcept>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

class MDWorld {
	public:
    /* attributes */
    size_t m_npart, m_npartmax;
    double m_mass;
    double m_energy_pot, m_energy_kin;
    gsl_matrix *m_x;        // positions
    gsl_matrix *m_v;        // velocities
    gsl_matrix *m_forces;   // forces

    /* constructor / destructor */
    LangevinObject(size_t npart);
    virtual ~LangevinObject();

    /* methods */
    virtual void compute_energy()=0; // energy
    virtual void deriv(const gsl_matrix *x, const gsl_matrix *v, gsl_matrix *dxdt,  gsl_matrix *dvdt)=0;   // compute derivatives
    virtual void update()=0; // update object
    void dump_conf(std::ostream &mystream);
    void dump_thermo(std::ostream &mystream);
};

//class Polymer : public LangevinObject {
//	public:
//    /* attributes */
//    double m_b; // bond length
//    double m_lp; // persistence length
//
//    /* constructor / destructor */
//    Polymer(size_t npart, double b, double lp);
//    ~Polymer();
//
//    /* methods */
//    virtual void force_bonds(std::vector<double[3]> &forces);
//    void force_bending(std::vector<double[3]> &forces);
//    void init_line();
//    void init_conf(std::istream &mystream);
//};

class LangevinSimulation {
  public:
    /* attributes */
    size_t m_itermax;
    double m_temp;
    double m_dt, m_s2dt;
    // rng: random number generator
    std::vector<LangevinObject> m_objs;
    // Stepper
    // StepperDiffusion

    /* constructor / destructor */
    LangevinSimulation(size_t itermax, double temp=1.0, double dt=0.001);
    ~LangevinSimulation();

    /* methods */
    void update();
};

#endif

/*
double distance(std::vector<double> &ri, std::vector<double> &rj);
void recenter();
//*/
