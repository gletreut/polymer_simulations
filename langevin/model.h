//*******************************************************************************
//*
//* Langevin Dynamics
//*	model.h
//*	Definition of model.h
//*
//* Author: Guillaume Le Treut
//*	CEA - 2014
//*
//*******************************************************************************

#ifndef _MODEL_H
#define _MODEL_H

#include <cstdlib>
#include <iostream>
#include <cmath>
#include <vector>
#include <limits>
#include <stdexcept>
#include <time.h>

#include "numerical_recipes/nr3.h"
#include "numerical_recipes/gamma.h"
#include "numerical_recipes/ran.h"
#include "numerical_recipes/deviates.h"


class LangevinObject {
	public:
	double E;
	virtual double energy()=0; // energy
	virtual void step(const double T, const double dt)=0; // one-timestep integration
	virtual void update()=0; // update objects
};

class Polymer : public LangevinObject {
	public:
	int N;
	double Ke, lp;
	std::vector<std::vector<double> > xyz, xyzn;
	Normaldev *ran;

	Polymer(int N_, double Ke_, double lp_);
	~Polymer();

	virtual double energy();
	virtual void step(const double T, const double dt);
	virtual void update();

	std::vector<double> force_gauss(const int k);
	std::vector<double> force_bending(const int k);
	std::vector<double> force_total(const int k);
	double distance(std::vector<double> &ri, std::vector<double> &rj);
	void recenter();
	void init();
	void init(istream &mystream);
	void print_current(ostream &mystream);

};

#endif
