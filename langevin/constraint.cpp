#include "constraint.h"

using namespace std;

//****************************************************************************
// Constraint
//****************************************************************************
Constraint::~Constraint(){
}

//****************************************************************************
// PointConstraint
//****************************************************************************
PointConstraint::PointConstraint(size_t index, double rx, double ry, double rz) :
  m_index(index) {

  /* initialize vectors */
  m_x0 = gsl_vector_calloc(3);
  m_v0 = gsl_vector_calloc(3);
  m_f0 = gsl_vector_calloc(3);

  gsl_vector_set(m_x0, 0, rx);
  gsl_vector_set(m_x0, 1, ry);
  gsl_vector_set(m_x0, 2, rz);
  gsl_vector_set_all(m_v0, 0.0);
  gsl_vector_set_all(m_f0, 0.0);
}

PointConstraint::~PointConstraint() {
  gsl_vector_free(m_x0);
  gsl_vector_free(m_v0);
  gsl_vector_free(m_f0);
}

void PointConstraint::constrain_x(gsl_matrix *x) {
  /*
   * Set the coordinates of the attached atom to the attachment point.
   */

  gsl_vector_view xn = gsl_matrix_row(x, m_index);
  gsl_vector_memcpy(&xn.vector, m_x0);

  return;
}

void PointConstraint::constrain_v(gsl_matrix *v) {
  /*
   * Set the velocities to zero.
   */

  gsl_vector_view vn = gsl_matrix_row(v, m_index);
  gsl_vector_memcpy(&vn.vector, m_v0);

  return;
}

void PointConstraint::constrain_force(gsl_matrix *forces) {
  /*
   * Project the forces so that there is no acceleration.
   */

  gsl_vector_view fn = gsl_matrix_row(forces, m_index);
  gsl_vector_memcpy(&fn.vector, m_f0);

  return;
}

//****************************************************************************
// AxisConstraint
//****************************************************************************
AxisConstraint::AxisConstraint(size_t index, double rx, double ry, double rz, double ux, double uy, double uz) :
  m_index(index) {

  /* initialize vectors */
  m_x0 = gsl_vector_calloc(3);
  m_u = gsl_vector_calloc(3);

  gsl_vector_set(m_x0, 0, rx);
  gsl_vector_set(m_x0, 1, ry);
  gsl_vector_set(m_x0, 2, rz);
  gsl_vector_set(m_u, 0, ux);
  gsl_vector_set(m_u, 1, uy);
  gsl_vector_set(m_u, 2, uz);

  /* normalize u */
  double unorm = linalg_dnrm2(m_u);
  linalg_dscal(1./unorm, m_u);

}

AxisConstraint::~AxisConstraint() {
  gsl_vector_free(m_x0);
  gsl_vector_free(m_u);
}

void AxisConstraint::constrain_x(gsl_matrix *x) {
  /*
   * Project the coordinates of the attached atom on the axis.
   */

  double rproj;

  gsl_vector_view xn = gsl_matrix_row(x, m_index);

  // make OM vector
  linalg_daxpy(-1., m_x0, &xn.vector);

  // compute axis projection
  rproj = linalg_ddot(&xn.vector, m_u);
  gsl_vector_memcpy(&xn.vector, m_x0);
  linalg_daxpy(rproj, m_u, &xn.vector);

  return;
}

void AxisConstraint::constrain_v(gsl_matrix *v) {
  /*
   * Project the velocity to the axis.
   */

  double vproj;

  gsl_vector_view vn = gsl_matrix_row(v, m_index);

  // compute axis projection
  vproj = linalg_ddot(&vn.vector, m_u);
  gsl_vector_set_all(&vn.vector, 0.0);
  linalg_daxpy(vproj, m_u, &vn.vector);

  return;
}

void AxisConstraint::constrain_force(gsl_matrix *forces) {
  /*
   * Project the force to the axis.
   */

  double fproj;

  gsl_vector_view fn = gsl_matrix_row(forces, m_index);

  // compute axis projection
  fproj = linalg_ddot(&fn.vector, m_u);
  gsl_vector_set_all(&fn.vector, 0.0);
  linalg_daxpy(fproj, m_u, &fn.vector);

  return;
}

//****************************************************************************
// PlaneConstraint
//****************************************************************************
PlaneConstraint::PlaneConstraint(size_t index, double rx, double ry, double rz, double ux, double uy, double uz) :
  m_index(index) {

  /* initialize vectors */
  m_x0 = gsl_vector_calloc(3);
  m_u = gsl_vector_calloc(3);
  m_vtp = gsl_vector_calloc(3);

  gsl_vector_set(m_x0, 0, rx);
  gsl_vector_set(m_x0, 1, ry);
  gsl_vector_set(m_x0, 2, rz);
  gsl_vector_set(m_u, 0, ux);
  gsl_vector_set(m_u, 1, uy);
  gsl_vector_set(m_u, 2, uz);

  /* normalize u */
  double unorm = linalg_dnrm2(m_u);
  linalg_dscal(1./unorm, m_u);

}

PlaneConstraint::~PlaneConstraint() {
  gsl_vector_free(m_x0);
  gsl_vector_free(m_u);
  gsl_vector_free(m_vtp);
}

void PlaneConstraint::constrain_x(gsl_matrix *x) {
  /*
   * Project the coordinates of the attached atom on the plane.
   */

  double rproj;

  gsl_vector_view xn = gsl_matrix_row(x, m_index);

  // compute projection along perpendicular axis
  gsl_vector_memcpy(m_vtp, &xn.vector);
  linalg_daxpy(-1., m_x0, m_vtp);
  rproj = linalg_ddot(m_vtp, m_u);

  // project to plane
  linalg_daxpy(-rproj, m_u, &xn.vector);

  return;
}

void PlaneConstraint::constrain_v(gsl_matrix *v) {
  /*
   * Project the velocity on the plane.
   */

  double vproj;

  gsl_vector_view vn = gsl_matrix_row(v, m_index);

  // compute projection along perpendicular axis
  vproj = linalg_ddot(&vn.vector, m_u);

  // project to plane
  linalg_daxpy(-vproj, m_u, &vn.vector);

  return;
}

void PlaneConstraint::constrain_force(gsl_matrix *forces) {
  /*
   * Project the force on the plane.
   */

  double fproj;

  gsl_vector_view fn = gsl_matrix_row(forces, m_index);

  // compute projection along perpendicular axis
  fproj = linalg_ddot(&fn.vector, m_u);

  // project to plane
  linalg_daxpy(-fproj, m_u, &fn.vector);

  return;
}
