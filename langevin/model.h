//*******************************************************************************
//*
//* Langevin Dynamics
//*	model.h
//*	Definition of model.h
//*
//* Author: Guillaume Le Treut
//*	UCSD - 2019
//*
//
//  Model:
//    * the damping term is taken to unity: \gamma = 1. Note that the
//    diffusion coefficient is given by: m \gamma D = T
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

#include "forcefield.h"

class MDWorld {
	public:
    /* attributes */
    size_t m_npart, m_npartmax;
    double m_temp;
    double m_mass;
    double m_energy_pot, m_energy_kin;
    gsl_vector_uint *m_mobile;        // vector storing the status of one particle (mobile:1 or fixed:0)
    gsl_matrix *m_x;                  // positions
    gsl_matrix *m_v;                  // velocities
    gsl_matrix *m_forces;             // forces
    std::vector<ForceField*> ffields; // list of active force fields (minus gradient of positional potential)

    /* constructor / destructor */
    MDWorld(size_t npart);
    virtual ~MDWorld();

    /* methods */
    // initiation methods
    void init_positions();
    void init_velocities();

    // update methods
    //void update_positions();
    //void update_velocities();
    void update_energy_forces();

    // dump methods
    void dump_conf(std::ostream &mystream);
    void dump_thermo(std::ostream &mystream);
};

//class Polymer : public MDWorld {
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

//class LangevinSimulation {
//  public:
//    /* attributes */
//    size_t m_itermax;
//    double m_temp;
//    double m_dt, m_s2dt;
//    // rng: random number generator
//    std::vector<MDWorld> m_objs;
//    // Stepper
//    // StepperDiffusion
//
//    /* constructor / destructor */
//    LangevinSimulation(size_t itermax, double temp=1.0, double dt=0.001);
//    ~LangevinSimulation();
//
//    /* methods */
//    void update();
//};

/*
double distance(std::vector<double> &ri, std::vector<double> &rj);
void recenter();
//*/

#endif

