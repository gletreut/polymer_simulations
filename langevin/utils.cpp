//----------------------------------------------------------------------------
//  2019-03-29
//  G. Le Treut - Jun Lab, UC San Diego.
//
//  file: utils.cpp
//----------------------------------------------------------------------------
#include "utils.h"
using namespace std;

vector<string>& utils::split(const string &s, char delim, vector<string> &elems) {
	    stringstream ss(s);
	        string item;
		    while (getline(ss, item, delim)) {
			            elems.push_back(item);
				        }
		        return elems;
}

vector<string> utils::split(const string &s, char delim) {
	    vector<string> elems;
	        split(s, delim, elems);
		    return elems;
}

//template <typename T>
//struct VectorSorter {
//	vector<T> &myvec;
//
//	VectorSorter(vector<T> &values): myvec(values) {}
//
//	bool operator() (size_t i1, size_t i2) {
//		return myvec[i1] < myvec[i2];
//	}
//};
//
//template <typename T>
//vector<int> argsort(vector<T> &v){
//	vector<int> idx(v.size());
//	for (size_t i=0; i != idx.size(); i++) idx[i] = i;
//
//	VectorSorter<T> mycomp(v);
//	sort(idx.begin(), idx.end(),mycomp);
//
//	return idx;
//}

