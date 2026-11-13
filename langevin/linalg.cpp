//****************************************************************************
// linalg.cpp - Implementation of the definitions in linalg.h.
// Date: 2017-05-03
// Created by: Guillaume Le Treut
//****************************************************************************

/* custom files */
#include "linalg.h"

using namespace std;
/* linear algebra */
double
linalg_dnrm2(const gsl_vector *v){
  double vnorm;
#ifdef GSL
  vnorm = gsl_blas_dnrm2(v);
#elif defined(LAPACK) || defined(MKL)
  vnorm = cblas_dnrm2(v->size, v->data, 1);
#endif

  return vnorm;
}

double
linalg_dnrm2(const gsl_matrix *m){
  gsl_vector_view mview = gsl_vector_view_array(m->block->data,m->block->size);

  return linalg_dnrm2(&mview.vector);
}

double
linalg_ddot(const gsl_vector *a, const gsl_vector *b){
  double ddot(0);

#ifdef GSL
  gsl_blas_ddot(a,b,&ddot);
#elif defined(LAPACK) || defined(MKL)
  ddot = cblas_ddot(a->size, a->data, 1, b->data, 1);
#endif
  return ddot;
}

double
linalg_ddot(const gsl_matrix *a, const gsl_matrix *b){
  gsl_vector_view aview = gsl_vector_view_array(a->block->data,a->block->size);
  gsl_vector_view bview = gsl_vector_view_array(b->block->data,b->block->size);

  return linalg_ddot(&aview.vector, &bview.vector);
}

void
linalg_dscal(const double s, gsl_vector *a){
#ifdef GSL
  gsl_blas_dscal(s,a);
#elif defined(LAPACK) || defined(MKL)
  cblas_dscal(a->size, s, a->data, 1);
#endif
  return;
}

void
linalg_dscal(const double s, gsl_matrix *a){

  gsl_vector_view aview = gsl_vector_view_array(a->block->data,a->block->size);

  linalg_dscal(s,&aview.vector);
  return;
}

void
linalg_daxpy(const double alpha, const gsl_vector *a, gsl_vector *b){

#ifdef GSL
  gsl_blas_daxpy(alpha,a,b);
#elif defined(LAPACK) || defined(MKL)
  cblas_daxpy(a->size, alpha, a->data, 1, b->data, 1);
#endif
  return;
}

void
linalg_daxpy(const double alpha, const gsl_matrix *a, gsl_matrix *b){

  gsl_vector_view aview = gsl_vector_view_array(a->block->data,a->block->size);
  gsl_vector_view bview = gsl_vector_view_array(b->block->data,b->block->size);

  linalg_daxpy(alpha,&aview.vector,&bview.vector);
  return;
}

void
linalg_dgemv  (int transA, const double alpha,
               const gsl_matrix *A, const gsl_vector *v,
               const double beta, gsl_vector *w){
#ifdef GSL
  CBLAS_TRANSPOSE_t transA_;
  if (transA == 0)
    transA_ = CblasNoTrans;
  else
    transA_=CblasTrans;

  gsl_blas_dgemv(transA_,alpha,A,v,beta,w);
#elif defined(LAPACK) || defined(MKL)
  CBLAS_TRANSPOSE transA_;
  int M,N,LDA;
  M = A->size1;
  N = A->size2;
  LDA = N;
  if (transA == 0)
    transA_ = CblasNoTrans;
  else
    transA_ = CblasTrans;

  cblas_dgemv(CblasRowMajor, transA_,
                 M, N,
                 alpha, A->data, LDA,
                 v->data, 1,
                 beta, w->data, 1);
#endif
  return;
}
