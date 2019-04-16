//*******************************************************************************
//*
//* Langevin Dynamics
//*	model.cpp
//*	Implementation of model.h
//*
//* Author: Guillaume Le Treut
//*	CEA - 2014
//*
//*******************************************************************************
#include "model.h"

using namespace std;

//****************************************************************************
//* LangevinSimulation
//****************************************************************************
LangevinSimulation::LangevinSimulation(size_t itermax, double temp, double dt) :
  m_itermax(itermax), m_temp(temp), m_dt(dt) {
  m_s2dt = sqrt(2.*m_dt);
}

LangevinSimulation::~LangevinSimulation(){
}

void
LangevinSimulation::update() {
  /*
   * Update all objects
   */
  for (vector<LangevinObject>::iterator it=m_objs.begin(); it != m_objs.end(); ++it){
    it->update();
  }
  return;
}

//****************************************************************************
//* LangevinObject
//****************************************************************************
LangevinObject::LangevinObject(size_t npart) : m_npart(npart) {
  /* initializations */
  m_mass = gsl_vector_calloc(m_npart);
  gsl_vector_set_all(m_mass, 1.0);
  m_x = gsl_matrix_calloc(m_npart, 3);
  gsl_matrix_set_all(m_x, 0.0);
  m_v = gsl_matrix_calloc(m_npart, 3);
  gsl_matrix_set_all(m_v, 0.0);
  m_energy_pot=0.;
  m_energy_kin=0.;
}

LangevinObject::~LangevinObject() {
  gsl_vector_free(m_mass);
  gsl_matrix_free(m_x);
  gsl_matrix_free(m_v);
}

//void
//LangevinObject::dump_thermo(ostream &mystream){
//  /*
//   * Dump thermodynamic quantities.
//   */
//
//  mystream << left << dec;
//  mystream << setw(10) << setprecision(0) << fixed << noshowpos << m_npart;
//  mystream << setw(18) << setprecision(8) << scientific << showpos << m_energy_pot;
//  mystream << setw(18) << setprecision(8) << scientific << showpos << m_energy_kin;
//  return;
//}
//
//void
//LangevinObject::dump_conf(ostream &mystream){
//  /*
//   * Dump configuration
//   */
//
//	if (m_npart == 0)
//		throw invalid_argument("LangevinObject is empty!");
//
//  mystream << left << dec << fixed << setprecision(8) << showpos;
//
//	for (size_t i=0; i<m_npart;i++){
//		mystream << setw(10) << i;
//		mystream << setw(18) << m_x[i][0];
//		mystream << setw(18) << m_x[i][1];
//		mystream << setw(18) << m_x[i][2];
//		mystream << endl;
//	}
//  return;
//}

//****************************************************************************
//* Polymer
//****************************************************************************
//Polymer::Polymer(int N_, double Ke_, double lp_) : Ke(Ke_), lp(lp_), N(N_) {
//		long long int seconds;
//		int ndigit=4, seed;
//		cout << left << setprecision(8) << fixed;
//
//		seconds=time(NULL); // seconds since 1970
//		seed=int( seconds % int(pow(10,ndigit)) );
//		cout << "seed=" << seed << endl;
//		ran=new Normaldev(0.0,1.0,seed);
//
//		init();
//		recenter();
//		xyzn=xyz;
//		E=energy();
//}
//
//Polymer::~Polymer(){
//	delete ran;
//}
//
//void Polymer::compute_energy(){
//	/*
//	 * Compute the energy of the new configuration
//	 */
//
//	double d;
//	vector<double> un(3),u(3);
//
//  energy = 0.;
//	// 1) chain
//	d=distance(xyz[1],xyz[0]);
//	energy+=Ke*d*d;
//	for (size_t n=2; n<N; n++){
//		for (size_t l=0; l<3; l++){
//			un[l]=xyz[n][l]-xyz[n-1][l];
//			u[l]=xyz[n-1][l]-xyz[n-2][l];
//		}
//
//		// spring
//		d=distance(xyz[n],xyz[n-1]);
//		energy+=Ke*d*d;
//
//		// curvature
//		d=distance(un,u);
//		energy+=lp/2.0*d*d;
//	}
//
//
//	return;
//}
//
//vector<double> Polymer::force_gauss(const int k){
//	/* Compute the force applied to monomer k
//	 */
//	vector<double> F(3);
//	int D;
//
//	if ( (k<0) || !(k<N) )
//		throw invalid_argument("k must range from 0 to N-1");
//
//	for (D=0; D<3; D++){
//		F[D]=0.0;
//		if ( (k+1<N) )
//			F[D]+=2.0*Ke*(xyz[k+1][D]-xyz[k][D]);
//		if ( !(k-1<0) )
//			F[D]+=2.0*Ke*(xyz[k-1][D]-xyz[k][D]);
//	}
//
//	return F;
//}
//
//vector<double> Polymer::force_bending(const int k){
//	/* Compute the force applied to monomer k
//	 */
//	vector<double> F(3);
//	int kk,D;
//
//	if ( (k<0) || !(k<N) )
//		throw invalid_argument("k must range from 0 to N-1");
//
//	for (D=0; D<3; D++){
//		F[D]=0.0;
//		if ( (k+2<N) )
//			F[D]+=-lp*( xyz[k+2][D]-2.0*xyz[k+1][D]+xyz[k][D]);
//		if ( (k+1<N) && !(k-1<0) )
//			F[D]+=-lp*(-2.0*xyz[k+1][D]+4.0*xyz[k][D]-2.0*xyz[k-1][D]);
//		if ( !(k-2<0) )
//			F[D]+=-lp*(xyz[k][D]-2.0*xyz[k-1][D]+xyz[k-2][D]);
//	}
//
//	return F;
//}
//
//vector<double> Polymer::force_total(const int k){
//	/*
//	 * Compute total force applied to monomer k
//	 */
//	vector<double> Ftot(3), F;
//	int D;
//
//	for (D=0; D<3; D++)
//		Ftot[D]=0.0;
//
//	// Gaussian chain contribution
//	F=force_gauss(k);
//	for (D=0; D<3; D++)
//		Ftot[D]+=F[D];
//
//	// Bending contribution
//	F=force_bending(k);
//	for (D=0; D<3; D++)
//		Ftot[D]+=F[D];
//
//	return Ftot;
//}
//
//void Polymer::step(const double T, const double dt){
//	/*
//	 * generate a new configuration from the current one
//	 */
//	vector<double> F(3);
//	double r, s=sqrt(2.0*T*dt);
//
//	for (int k=1; k<N; k++){
//		// compute the force applied to monomer k;
//		F=force_total(k);
//
//		for (int D=0; D<3; D++){
//
//			// draw random number with unit Gaussian distribution
//			r=ran->dev();
//
//			// integration step
//			xyzn[k][D]=xyz[k][D]+F[D]*dt+s*r;
//		}
//	}
//	return;
//}
//void Polymer::update(){
//	xyz=xyzn;
//	return;
//}
//
//void Polymer::recenter(){
//	double x0,y0,z0;
//
//	x0=xyz[0][0];
//	y0=xyz[0][1];
//	z0=xyz[0][2];
//
//	for (size_t i=0;i<N;i++){
//		xyz[i][0]-=x0;
//		xyz[i][1]-=y0;
//		xyz[i][2]-=z0;
//	}
//	return;
//}
//
//void Polymer::init_line(){
//	xyz.resize(N);
//	for (int i=0; i<N; i++){
//		xyz[i].resize(3);
//		xyz[i][0]=-N/2.0 + i;
//		xyz[i][1]=0.0;
//		xyz[i][2]=0.0;
//	}
//	cout << "POLYMER CONFIGURATION INITIALIZED" << endl;
//	return;
//}
//
//void Polymer::init_conf(istream &mystream){
//	string line;
//	stringstream convert;
//	double x,y,z;
//	int ind;
//	vector<double> coord;
//
//	xyz.clear();
//	while (getline(mystream,line)){
//		convert.clear();
//		convert.str(line);
//
//		if ( (convert >> ind) &&
//	             (convert >> x) &&
//		     (convert >> y) &&
//		     (convert >> z) ){
//			coord.clear();
//			coord.push_back(x);
//			coord.push_back(y);
//			coord.push_back(z);
//			xyz.push_back(coord);
//		}
//
//	}
//
//	if (xyz.size() != N)
//		throw invalid_argument("wrong size of the input configuration.");
//
//	return;
//
//}
//
//double Polymer::distance(vector<double> &ri, vector<double> &rj){
//	double sum;
//	sum=0.0;
//
//	if (ri.size() != rj.size()){
//		cerr << ("ri and rj must have same size!") << endl;
//		return 0.0;
//	}
//
//	for (size_t i=0; i<ri.size(); i++)
//		sum+=(ri[i]-rj[i])*(ri[i]-rj[i]);
//	sum=sqrt(sum);
//	return sum;
//}
//
//void Polymer::print_current(ostream &mystream){
//	mystream << left << fixed << setprecision(6);
//
//	if (xyz.size() == 0)
//		throw invalid_argument("xyz is empty!");
//
//	for (size_t i=0; i<N;i++){
//		mystream << setw(10) << i;
//		mystream << setw(20) << xyz[i][0];
//		mystream << setw(20) << xyz[i][1];
//		mystream << setw(20) << xyz[i][2];
//		mystream << endl;
//	}
//	//mystream << endl;
//
//	return;
//}


//vector<double> Polymer::force_couplings(const int k){
//	/*
//	 * Compute the force applied to monomer k
//	 */
//
//	vector<double> F(3);
//	int D;
//
//	if ( (k<0) || !(k<N) )
//		throw invalid_argument("k must range from 0 to N-1");
//
//	for (D=0; D<3; D++){
//		F[D]=0.0;
//		for (int j=0; j<N; j++){
//			if ( k != j ){
//			//	cout << "delta[" << D << "]=" <<  (xyz[j][D]-xyz[k][D]);
//			//	cout << "  K[" << j << "][" << k << "]=" << K[j][k];
//			//	cout << endl;
//				F[D]+=2.0*Ke*K[j][k]*(xyz[j][D]-xyz[k][D]);
//			}
//		}
//	}
//
//	return F;
//}
