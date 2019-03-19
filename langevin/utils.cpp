//*******************************************************************************
//*
//* Langevin Dynamics
//*	utils.cpp
//*	Implementation of utils.h
//*
//* Author: Guillaume Le Treut
//*	CEA - 2014
//*
//*******************************************************************************
#include "utils.h"

std::vector<std::string> & string_utils::split(const std::string &s, char delim, std::vector<std::string> &elems) {
	    std::stringstream ss(s);
	        std::string item;
		    while (std::getline(ss, item, delim)) {
			            elems.push_back(item);
				        }
		        return elems;
}


std::vector<std::string> string_utils::split(const std::string &s, char delim) {
	    std::vector<std::string> elems;
	        split(s, delim, elems);
		    return elems;
}
