//*******************************************************************************
//*
//* Langevin Dynamics
//*	integration.cpp
//*	Implementation of integration.h
//*
//* Author: Guillaume Le Treut
//*	UCSD - 2019
//*
//*******************************************************************************
#include "integration.h"

using namespace std;

//****************************************************************************
//* IntegrationParams
//****************************************************************************
IntegrationParams::IntegrationParams() :
    itermax(0), idump_thermo(1), idump_pos(1), pos_xyz(true), pos_dat(false), ineighbor_build(10) {

      init();
    }

void
IntegrationParams::init() {
  iterwidth = size_t(log10(itermax)) + 1;
  return;
}

void
IntegrationParams::print_infos(ostream &mystream) {
  /*
   * print some information on the world object
   */

  mystream << left << dec << fixed;
  mystream << setw(20) << "itermax" << setw(20) << setprecision(0) << noshowpos << itermax << endl;
  mystream << setw(20) << "iterwidth" << setw(20) << setprecision(0) << noshowpos << iterwidth << endl;
  mystream << setw(20) << "idump_thermo" << setw(20) << setprecision(0) << noshowpos << idump_thermo << endl;
  mystream << setw(20) << "idump_pos" << setw(20) << setprecision(0) << noshowpos << idump_pos << endl;
  mystream << setw(20) << "ineighbor_build" << setw(20) << setprecision(0) << noshowpos << ineighbor_build << endl;
  mystream << setw(20) << "pos_xyz" << setw(20) << setprecision(0) << noshowpos << pos_xyz << endl;
  mystream << setw(20) << "pos_dat" << setw(20) << setprecision(0) << noshowpos << pos_dat << endl;
  return;
}

