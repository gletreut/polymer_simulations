#ifndef _UTILS_H
#define _UTILS_H
//----------------------------------------------------------------------------
//  2019-03-29
//  G. Le Treut - Jun Lab, UC San Diego.
//
//  file: utils.h
//----------------------------------------------------------------------------

#include <cstdlib>
#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <map>

namespace utils {
	std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems);
	std::vector<std::string> split(const std::string &s, char delim);

  template<typename T>
  void load_matrix(std::istream &mystream, std::vector<std::vector<T> > &cmap);

  template<typename T>
  void print_matrix(std::ostream &mystream, std::vector<std::vector<T> > &cmap);

  template<typename T>
  void load_map(std::istream &mystream, std::map<std::string,T> &params);

  template<typename T>
  void print_map(std::ostream &mystream, std::map<std::string,T> &params);

};

// include template definitions
#include "utils.tpp"

#endif
