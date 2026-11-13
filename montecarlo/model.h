//*******************************************************************************
//*
//* Montecarlo
//*	model.h
//*	Definition of the model
//*
//* 	Author: Guillaume Le Treut
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
#include <ostream>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


class MonteCarloObject {
	public:
	double E;
	virtual double energy()=0; // energy
	virtual void gen()=0; // generate new configuration
	virtual void accept(const double t)=0; // accept or not the move
					     // with temperature t
};

class Polymer : public MonteCarloObject {
	public:
	int N;
	double lp;
	std::vector<std::vector<double> > xyz, xyzn;
	gsl_rng *rng;

	Polymer(int N_, double lp_, unsigned long seed);
	~Polymer();

	virtual double energy();
	virtual void gen();
	virtual void accept(const double t);

	double distance(std::vector<double> &ri, std::vector<double> &rj);
	void init();
	void print_current(std::ostream &mystream);

};

#endif
