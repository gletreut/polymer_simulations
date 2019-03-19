//*******************************************************************************
//*
//* Langevin Dynamics
//*	utils.cpp
//*	Definition of utils.h
//*
//* Author: Guillaume Le Treut
//*	CEA - 2014
//*
//*******************************************************************************

#ifndef _HOUSEKEEPING_H
#define _HOUSEKEEPING_H

#include <cstdlib>
#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <map>

namespace string_utils{
	std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems);
	std::vector<std::string> split(const std::string &s, char delim);
};

template<typename T>
void load_map(std::istream &mystream, std::vector<std::vector<T> > &cmap){
	std::string line;
	std::stringstream convert;
	int i,j,n;
	T val;
	typename
	std::vector<std::vector<T> >::iterator it;

	while (getline(mystream,line)){
		convert.clear();
		convert.str(line);

		if ( (convert >> i) && (convert >> j) && (convert >> val) ){
			n=i+1;
			if (j+1 > n) n=j+1;

			// size issues
			if (cmap.size() < n)
				cmap.resize(n);
			for (it=cmap.begin(); it!=cmap.end(); it++){
				if (it->size() < n)
					it->resize(n);
			}

			// put value in table
			cmap[i][j]=val;
		}
	}

	std::cout << "MAP WITH SIZE " << n << " x " << n << " IMPORTED" << std::endl;

}

template<typename T>
void print_map(std::ostream &mystream, std::vector<std::vector<T> > &cmap){

	mystream << std::left;
	//mystream << std::fixed << std::setprecision(8);
	mystream << std::scientific << std::setprecision(8);

	for (size_t i=0; i<cmap.size(); i++){
		for (size_t j=0; j<cmap.size(); j++){
			mystream << std::setw(10) << i;
			mystream << std::setw(10) << j;
			mystream << std::setw(20) << cmap[i][j];
			mystream << std::endl;
		}
	}
}

template<typename T>
void load_params(std::istream &mystream, std::map<std::string,T> &params){
	std::string line, key;
	std::stringstream convert;
	T val;
	std::cout << std::left;

	while (getline(mystream, line)){
		convert.clear();
		convert.str(line);

		if (convert >> key) {// there is a string
			if (convert >> val) {// there is a value
				params[key]=val;
				std::cout << std::setw(30) << key;
				std::cout << val << std::endl;
			}
		}
	}

	return;
}

template<typename T>
void print_map(std::ostream &mystream, std::map<std::string,T> &params){
	int lmax,l;
	typename
	std::map<std::string,T>::iterator it;

	lmax=0;
	for (it=params.begin(); it!=params.end(); it++){
		l=it->first.length();
		if (l>lmax) lmax=l;
	}

	lmax+=2;

	mystream << std::left;
	for (it=params.begin(); it!=params.end(); it++){
		mystream << std::setw(lmax) << it->first.c_str();
		mystream << it->second;
		mystream << std::endl;
	}
	return;
}

//template <typename T>
//struct VectorSorter {
//	std::vector<T> &myvec;
//
//	VectorSorter(std::vector<T> &values): myvec(values) {}
//
//	bool operator() (size_t i1, size_t i2) {
//		return myvec[i1] < myvec[i2];
//	}
//};
//
//template <typename T>
//std::vector<int> argsort(std::vector<T> &v){
//	std::vector<int> idx(v.size());
//	for (size_t i=0; i != idx.size(); i++) idx[i] = i;
//
//	VectorSorter<T> mycomp(v);
//	std::sort(idx.begin(), idx.end(),mycomp);
//
//	return idx;
//}


#endif
