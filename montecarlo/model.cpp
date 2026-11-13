//*******************************************************************************
//*
//* Montecarlo
//*	model.cpp
//*	Implementation of the model
//*
//* 	Author: Guillaume Le Treut
//*	CEA - 2014
//*
//*******************************************************************************

#include <iomanip>

#include "model.h"
#include "geometry.h"

using namespace std;
using namespace geometry;

//void algebra_transformation_crankshaft(vector<vector<double> > xyz, int r1, int r2, double theta, vector<vector<double> > xyzn) {return;}
//void algebra_transformation_pivot(vector<vector<double> > xyz, int r1, vector<double> v, double theta, vector<vector<double> > xyzn) {return;}

Polymer::Polymer(int N_, double lp_, unsigned long seed) : N(N_), lp(lp_) {
		cout << left << setprecision(8) << fixed;

		cout << "seed=" << seed << endl;
		rng=gsl_rng_alloc(gsl_rng_mt19937);
		gsl_rng_set(rng, seed);

		init();
		xyzn=xyz;
		E=energy();
}

Polymer::~Polymer(){
	gsl_rng_free(rng);
}

double Polymer::energy(){
	/*
	 * Compute the energy of the new configuration
	 */

	double e=0.0, d;
	vector<double> un(3),u(3);

	// 1) chain
	for (size_t n=2; n<N; n++){
		for (size_t l=0; l<3; l++){
			un[l]=xyzn[n][l]-xyzn[n-1][l];
			u[l]=xyzn[n-1][l]-xyzn[n-2][l];
		}

		// curvature
		d=distance(un,u);
		e+=lp/2.0*d*d;
	}

	return e;
}
void Polymer::gen(){
	/*
	 * generate a new configuration from the current one
	 */
	double PI=3.14159265359, theta, vnrm;
	int rc, r1, r2;
	vector<double> v(3);

	rc=int(gsl_rng_uniform(rng)*2);

	if (rc == 0) {// Crankshaft move
		r1=int(gsl_rng_uniform(rng)*N);
		r2=r1;
		while (r1==r2)
			r2=int(gsl_rng_uniform(rng)*N);
		theta=2.0*PI*gsl_rng_uniform(rng);
		algebra_transformation_crankshaft(xyz,r1,r2,theta,xyzn);
	}
	else { // Pivot move
		r1=int(gsl_rng_uniform(rng)*N);
		theta=2.0*PI*gsl_rng_uniform(rng);
		v[0]=gsl_ran_gaussian(rng, 1.0);
		v[1]=gsl_ran_gaussian(rng, 1.0);
		v[2]=gsl_ran_gaussian(rng, 1.0);
		vnrm=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
		v[0]/=vnrm;
		v[1]/=vnrm;
		v[2]/=vnrm;
		algebra_transformation_pivot(xyz,r1,v,theta,xyzn);
	}

	return;
}
void Polymer::accept(const double t){
	double r, Enew;

	Enew=energy();
	r=gsl_rng_uniform(rng);
	if (r < exp(-(Enew-E)/t) ) {// accept
		xyz=xyzn;
		E=Enew;
	}
	return;
}

void Polymer::init(){
	xyz.resize(N);
	for (int i=0; i<N; i++){
		xyz[i].resize(3);
		xyz[i][0]=-N/2.0 + i;
		xyz[i][1]=0.0;
		xyz[i][2]=0.0;
	}
	cout << "POLYMER CONFIGURATION INITIALIZED" << endl;
	return;
}

double Polymer::distance(vector<double> &ri, vector<double> &rj){
	double sum;
	sum=0.0;

	if (ri.size() != rj.size()){
		cerr << ("ri and rj must have same size!") << endl;
		return 0.0;
	}

	for (size_t i=0; i<ri.size(); i++)
		sum+=(ri[i]-rj[i])*(ri[i]-rj[i]);
	sum=sqrt(sum);
	return sum;
}

void Polymer::print_current(ostream &mystream){
	mystream << left << fixed << setprecision(6);

	if (xyz.size() == 0)
		throw invalid_argument("xyz is empty!");

	for (size_t i=0; i<N;i++){
		mystream << setw(10) << i;
		mystream << setw(20) << xyz[i][0];
		mystream << setw(20) << xyz[i][1];
		mystream << setw(20) << xyz[i][2];
		mystream << endl;
	}
	//mystream << endl;

	return;
}
