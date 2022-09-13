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
#include <fstream>
#include <cmath>
#include <vector>
#include <limits>
#include <stdexcept>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_rng.h>

#include "forcefield.h"
#include "constraint.h"
#include "neighborlist.h"
#include "utils.h"

class MDWorld {
	public:
    /* attributes */
    size_t m_npart, m_npartmax;
    double m_lx, m_ly, m_lz;                  // dimensions of the box
    double m_gamma;
    double m_temp;
    double m_mass;
    size_t m_dim;                             // dimension for the dynamics.
    double m_energy_pot, m_energy_kin;
    // gsl_vector_uint *m_mobile;                // vector storing the status of one particle (mobile:1 or fixed:0)
    gsl_matrix *m_x;                          // positions
    gsl_matrix *m_v;                          // velocities
    gsl_matrix *m_forces, *m_forces_tp;       // forces
    std::vector<ForceField*> m_ffields;       // list of active force fields (minus gradient of positional potential)
    std::vector<Constraint*> m_constraints;   // list of active constraints
    std::vector<NeighborList*> m_neighbors;   // list of neighbor lists
    gsl_spmatrix_uint *m_bonds;               // is there bonds (0 or 1) between particles. If 1 then non-bonded interactions not computed.

    /* constructor / destructor */
    MDWorld(size_t npart=1, double lx=1., double ly=1., double lz=1.,
            double gamma=1.0, double temp=1.0, double mass=1.0, size_t dim=3);
    virtual ~MDWorld();

    /* methods */
    // clearing/init
    void init();
    void clear();

    // initiation methods
    void init_positions_lattice(double delta);
    void init_velocities(gsl_rng *rng);
    void set_constraints();

    // update methods
    //void update_positions();
    //void update_velocities();
    void update_energy_forces();
    void update_energy_kinetics();
    void constrain_x();
    void constrain_v();
    bool check_neighbor_lists();
    void update_neighbor_lists(size_t iter);

    // dump methods
    void dump_pos(std::ostream &mystream, bool index=true, bool positions=true, bool velocities=false, bool forces=false);
    void dump_thermo(std::ostream &mystream);
    void dump_dat(std::string fileout);
    void dump_xyz(std::string fileout, size_t iter);
    void dump_neighbors(std::ostream &mystream);
    void dump_neighbors(std::string fileout);
    void dump_polarity_vectors(std::ostream &mystream);
    void dump_polarity_vectors(std::string fileout);

    // loading methods
    void load_dat(std::string filein);
    void load_xyz(std::string filein);

    // miscellaneous
    void print_infos(std::ostream &mystream);
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

