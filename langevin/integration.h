#ifndef _INTEGRATION_H
#define _INTEGRATION_H
//*******************************************************************************
//*
//* Langevin Dynamics
//*	integration.h
//*	Definition of integration.h
//*
//* Author: Guillaume Le Treut
//*	UCSD - 2019
//*
//* structure with parameter for integration.
//*******************************************************************************

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>

struct IntegrationParams {

  size_t itermax;
  size_t iterwidth;
  size_t idump_thermo;
  size_t idump_pos;
  size_t ineighbor_build;
  bool pos_xyz;
  bool pos_dat;

  IntegrationParams();

  void init();
  void print_infos(std::ostream &mystream);
};
#endif
