#define _USE_MATH_DEFINES

#ifndef _GEOMETRY_H
#define _GEOMETRY_H

//*******************************************************************************
//*
//* Montecarlo
//*	geometry.h
//*	Definition of the geometry namespace
//*
//* 	Author: Guillaume Le Treut
//*	CEA - 2014
//*
//*******************************************************************************

//* Standard libraries to include
#include <cstdlib>					// Common
#include <iostream>					// I/O library C++
#include <cmath>					// maths library C++
#include <vector>					// Vector class C++
#include <exception>					// Exception handler
#include <stdexcept>


namespace geometry
{
	/* Algebra Methods */

	double algebra_dot_vv(const std::vector<double> &v1, const std::vector<double> &v2);
	std::vector<double> algebra_cross_vv(const std::vector<double> &v1, const std::vector<double> &v2);
	std::vector<double> algebra_rotation_to_quaternion(const double &theta, const std::vector<double> &axis);
	std::vector<double> algebra_quaternion_conjugate(const double &quaternion);
	std::vector<double> algebra_prod_qq(const std::vector<double> &q1, const std::vector<double> &q2);
	std::vector<double> algebra_to_quaternion(const double &s, const std::vector<double> &v);
	std::vector<double> algebra_q_getVector(const std::vector<double> &q);
	double algebra_q_getScalar(const std::vector<double> &q);
	std::vector<double> algebra_transformation_rotation(const std::vector<double> &v_in, const std::vector<double> &axis, const double &angle);

        void algebra_transformation_crankshaft(const std::vector<std::vector<double> > &r, const int i, const int j, const double &theta, std::vector<std::vector<double> > &rn);
	void algebra_transformation_pivot(const std::vector<std::vector<double> > &r, const int i, const std::vector<double> &axis, const double &theta, std::vector<std::vector<double> > &rn);


}
#endif
