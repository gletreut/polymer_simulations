//*******************************************************************************
//*
//* Montecarlo
//*	geometry.cpp
//*	Implementation of geometry.h namespace
//*
//* Author: Guillaume Le Treut
//*	CEA - 2014
//*
//*******************************************************************************

//* Specific libraries to include
#include "geometry.h"

using namespace std;

namespace geometry{
  /* Methods */

  double algebra_dot_vv(const vector<double> &v1, const vector<double> &v2){
  	/* Return the scalar product of vector v1 and v2 */

  	int m = min(v1.size(), v2.size());
  	double s = 0.0;

  	for (size_t i=0; i<m; i++){
  		s += v1[i]*v2[i];
  	}

  	return s;
  }

  vector<double> algebra_cross_vv(const vector<double> &v1, const vector<double> &v2){
  	/* Return the vectorial product of v1 with v2 : v1^v2 */

  	vector<double> v;

  	try
  	{
  		v.push_back(v1[1]*v2[2] - v1[2]*v2[1]);
  		v.push_back(v1[2]*v2[0] - v1[0]*v2[2]);
  		v.push_back(v1[0]*v2[1] - v1[1]*v2[0]);
  	}

  	catch(exception const& e)
  	{
  		cerr << "ERROR : " << e.what() << endl;
  	}

  	return v;
  }

  vector<double> algebra_to_quaternion(const double &s, const vector<double> &v){
  	/* Return the quaternion q = (s, v)
  	 */
  	vector<double> q(v.begin(), v.end());
  	q.insert(q.begin(), s);

  	return q;
  }

  vector<double> algebra_rotation_to_quaternion(const double &theta, const vector<double> &axis){
	  /* Return the quaternion representing the rotation R=(axis, theta)
	   */

	  double nrm;
	  vector<double> q, v;

	  v = axis;
	  nrm = sqrt(algebra_dot_vv(v, v));
	  for (size_t i=0; i<v.size(); i++)
		  v[i]*=sin(theta/2.)/nrm;

	  //cout << "norme q = " << algebra_dot_vv(algebra_to_quaternion(cos(theta/2.), v), algebra_to_quaternion(cos(theta/2.), v)) << endl;

	  return algebra_to_quaternion(cos(theta/2.), v);
  }

  vector<double> algebra_quaternion_conjugate(const vector<double> &q){
	  /* Return the conjugate quaternion
	   */
	  vector<double> qout;

	  qout = q;
	  qout[1] *= -1.;
	  qout[2] *= -1.;
	  qout[3] *= -1.;

	  return qout;
  }

  vector<double> algebra_q_getVector(const vector<double> &q){
  	/* Return the vectorial part of the quaternion q = (s, v)
  	 */

  	return vector<double>(q.begin()+1, q.end());
  }

  double algebra_q_getScalar(const vector<double> &q){
  	/* Return the scalar part of the quaternion q = (s, v)
  	 */

  	return q[0];
  }

  vector<double> algebra_prod_qq(const vector<double> &q1, const vector<double> &q2){
  	/* Return the product of quaternions q1 x q2
  	 * with the convention q(scalar, vector)
  	 */

  	vector<double> q;
  	double s, t, r;
  	vector<double> u, v, w;

	double a1, b1, c1, d1;
	double a2, b2, c2, d2;

	if (q1.size() != q2.size())
		throw invalid_argument("q1 and q2 must have same size!");

  	try
  	{
  		s = q1[0];
  		t = q2[0];
  		u = vector<double>(q1.begin()+1, q1.end());
  		v = vector<double>(q2.begin()+1, q2.end());
  		r = s*t - algebra_dot_vv(u, v);

  		w = algebra_cross_vv(u, v);
		for (size_t i=0; i<u.size(); i++){
			w[i]+=t*u[i]+s*v[i];
		}

  		q = algebra_to_quaternion(r, w);
  	}

  	catch (exception const &e)
  	{
  		cerr << "ERROR : " << e.what() << endl;
  	}

  	return q;
  }

  vector<double> algebra_transformation_rotation(const vector<double> &v_in, const vector<double> &axis, const double &angle){
	/* Return the input vector rotated from angle around axis
	 */

	  vector<double> v_out, q, qc;

  	  q = algebra_rotation_to_quaternion(angle, axis);
	  qc = algebra_quaternion_conjugate(q);

          v_out = algebra_q_getVector(algebra_prod_qq(q, algebra_prod_qq(algebra_to_quaternion(0., v_in), qc)));

	  return v_out;
  }

  void algebra_transformation_crankshaft(const vector<vector<double> > &r,const int i, const int j, const double &theta, vector<vector<double> > &rn){
        /* Perform a Crankshaft move on input vector
	 * using rj - ri as rotation axis, with angle theta
	 */

	  int m, n;

	  m = min(i, j);
	  n = max(i, j);
	  vector<double> u(3), v(3), axis(3);

	  rn = r;
	  axis[0] = r[n][0] - r[m][0];
	  axis[1] = r[n][1] - r[m][1];
	  axis[2] = r[n][2] - r[m][2];

	  try{
		  for (size_t k=m+1; k<n; k++){
			  u[0] = r[k][0] - r[m][0];
			  u[1] = r[k][1] - r[m][1];
			  u[2] = r[k][2] - r[m][2];

			  v = algebra_transformation_rotation(u, axis, theta);

			  rn[k][0] = rn[m][0] + v[0];
			  rn[k][1] = rn[m][1] + v[1];
			  rn[k][2] = rn[m][2] + v[2];
		  }
	  }
  	  catch (exception const &e)
  	  {
  		  cerr << "ERROR : " << e.what() << endl;
  	  }

	  return;

  }

  void algebra_transformation_pivot(const vector<vector<double> > &r, const int i, const vector<double> &axis, const double &theta, vector<vector<double> > &rn){
	  /* Perform a Pivot move on input vector
	   * using (i, axis) as the axis, with angle theta
	   */

	  int N;
	  vector<double> u(3), v(3);
	  rn = r;

	  N = r.size();

	  try{
		   for (size_t k=i+1; k<N; k++){
			   u[0] = r[k][0] - r[i][0];
			   u[1] = r[k][1] - r[i][1];
			   u[2] = r[k][2] - r[i][2];

			   v = algebra_transformation_rotation(u, axis, theta);

		  	   rn[k][0] = rn[i][0] + v[0];
		  	   rn[k][1] = rn[i][1] + v[1];
		  	   rn[k][2] = rn[i][2] + v[2];
		   }
	  }

  	  catch (exception const &e)
  	  {
  	       	  cerr << "ERROR : " << e.what() << endl;
  	  }

	return;
  }

}

