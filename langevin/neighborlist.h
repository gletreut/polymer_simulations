//*******************************************************************************
//*
//* Langevin Dynamics
//*	neighborlist.h
//*
//* Author: Guillaume Le Treut
//*	CZ Biohub - 2021
//*
//*******************************************************************************

#ifndef _NEIGHBORLIST_H
#define _NEIGHBORLIST_H

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <limits>
#include <stdexcept>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_rng.h>

#include "linalg.h"
#include "utils.h"

class NeighborList {
  /*
   * Class defining a neighbor list.
   */
  public:
    /* attributes */
    size_t m_nmax;
    size_t m_npair;
    size_t m_iter, m_iterref;
    size_t m_ibuild_count, m_ibuild_sum;
    double m_rskin;
    double m_rmax;
    gsl_matrix_uint *m_pairs;
    gsl_matrix *m_xref;

    /* constructor and destructor */
    NeighborList(size_t nmax, double rskin, size_t npart);
    ~NeighborList();

    /* methods */
    void update_counter(size_t iter);
    void build(gsl_matrix *x, gsl_spmatrix_uint *bonds);
    bool check(gsl_matrix *x);
    void dump(std::ostream &mystream, bool positions=true);
};


#endif
