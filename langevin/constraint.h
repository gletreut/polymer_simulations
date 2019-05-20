//*******************************************************************************
//*
//* Langevin Dynamics
//*	constraint.h
//*
//* Author: Guillaume Le Treut
//*	UCSD - 2019
//*
//*******************************************************************************

#ifndef _CONSTRAINT_H
#define _CONSTRAINT_H

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <vector>
#include <limits>
#include <stdexcept>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>

#include "linalg.h"
#include "utils.h"

class Constraint {
  /*
   * Virtual class defining a constraint.
   */
  public:
    virtual ~Constraint();
    virtual void constrain_x(gsl_matrix *x) = 0;
    virtual void constrain_v(gsl_matrix *v) = 0;
    virtual void constrain_force(gsl_matrix *forces) = 0;
};

class PointConstraint : public Constraint {
  /*
   * Class defining a point constraint.
   * An atom is attached at a point location.
   */
  public:
    size_t m_index;       // index of the attached atom.
    gsl_vector *m_x0;     // coordinates of the attachment.

    /* constructor and destructor */
    PointConstraint(size_t index, double rx, double ry, double rz);
    ~PointConstraint();

    /* methods */
    void constrain_x(gsl_matrix *x);
    void constrain_v(gsl_matrix *v);
    void constrain_force(gsl_matrix *forces);

  private:
    gsl_vector *m_v0;
    gsl_vector *m_f0;
};

class AxisConstraint : public Constraint {
  /*
   * Class defining an axis constraint.
   * An atom is constrained to move along one axis.
   */
  public:
    size_t m_index;       // index of the constrained atom.
    gsl_vector *m_x0;     // point belonging to the axis
    gsl_vector *m_u;      // unit vector with the direction of the axis

    /* constructor and destructor */
    AxisConstraint(size_t index, double rx, double ry, double rz, double ux, double uy, double uz);
    ~AxisConstraint();

    /* methods */
    void constrain_x(gsl_matrix *x);
    void constrain_v(gsl_matrix *v);
    void constrain_force(gsl_matrix *forces);
};

class PlaneConstraint : public Constraint {
  /*
   * Class defining a plan constraint.
   * An atom is constrained to move within a plane
   */
  public:
    size_t m_index;       // index of the constrained atom.
    gsl_vector *m_x0;     // point belonging to the plane
    gsl_vector *m_u;      // unit vector perpendicular to the plane

    /* constructor and destructor */
    PlaneConstraint(size_t index, double rx, double ry, double rz, double ux, double uy, double uz);
    ~PlaneConstraint();

    /* methods */
    void constrain_x(gsl_matrix *x);
    void constrain_v(gsl_matrix *v);
    void constrain_force(gsl_matrix *forces);

  private:
    gsl_vector *m_vtp;
};

#endif
