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
//LangevinSimulation::LangevinSimulation(size_t itermax, double temp, double dt) :
//  m_itermax(itermax), m_temp(temp), m_dt(dt) {
//  m_s2dt = sqrt(2.*m_dt);
//}
//
//LangevinSimulation::~LangevinSimulation(){
//}
//
//void
//LangevinSimulation::update() {
//  /*
//   * Update all objects
//   */
//  for (vector<MDWorld>::iterator it=m_objs.begin(); it != m_objs.end(); ++it){
//    it->update();
//  }
//  return;
//}

//****************************************************************************
//* MDWorld
//****************************************************************************
MDWorld::MDWorld(size_t npart, double lx, double ly, double lz,
                 double gamma, double temp, double mass, size_t dim) :
  m_npart(npart), m_lx(lx), m_ly(ly), m_lz(lz),
  m_gamma(gamma), m_temp(temp), m_mass(mass), m_dim(dim){
  /* check */
  if (m_npart == 0){
    throw invalid_argument("Number of particles must be larger than zero!");
  }

  init();
}

MDWorld::~MDWorld() {
  clear();
}

void MDWorld::init() {

  // energy
  m_energy_pot=0.;
  m_energy_kin=0.;

  // gsl vectors
  m_x = gsl_matrix_calloc(m_npart, 3);
  gsl_matrix_set_all(m_x, 0.0);
  m_v = gsl_matrix_calloc(m_npart, 3);
  gsl_matrix_set_all(m_v, 0.0);
  m_forces = gsl_matrix_calloc(m_npart, 3);
  gsl_matrix_set_all(m_forces, 0.0);
  m_forces_tp = gsl_matrix_calloc(m_npart, 3);
  gsl_matrix_set_all(m_forces_tp, 0.0);

  // vectors
  m_ffields.clear();
  m_constraints.clear();
  m_neighbors.clear();

  // bond matrix
  m_bonds = gsl_spmatrix_uint_alloc(m_npart, m_npart);

  return;
}

void MDWorld::clear() {
  /* clear vectors and matrices */

  /* delete gsl objects */
  if (m_x != nullptr){
    gsl_matrix_free(m_x);
  }
  if (m_v != nullptr){
    gsl_matrix_free(m_v);
  }
  if (m_forces != nullptr){
    gsl_matrix_free(m_forces);
  }
  if (m_forces_tp != nullptr){
    gsl_matrix_free(m_forces_tp);
  }

  /* delete forcefields */
  for (vector<ForceField*>::iterator it=m_ffields.begin(); it!=m_ffields.end(); ++it){
    if (*it != nullptr){
      delete (*it);
    }
  }

  /* delete constraints */
  for (vector<Constraint*>::iterator it=m_constraints.begin(); it!=m_constraints.end(); ++it){
    if (*it != nullptr){
      delete (*it);
    }
  }

  /* delete neighbor lists */
  for (vector<NeighborList*>::iterator it=m_neighbors.begin(); it!=m_neighbors.end(); ++it){
    if (*it != nullptr){
      delete (*it);
    }
  }

  /* delete weight matrix */
  gsl_spmatrix_uint_free(m_bonds);

  return;
}

void MDWorld::init_positions_lattice(double delta){
  /*
   * Initialize positions on a lattice.
   */

  size_t counter;
  // size_t nx, ny, nz;
  //double delta;
  double x,y,z,u;
  int dirx, diry;

  u =pow(2.,1./6);
  // nx = (m_lx-2*u)/delta;
  // ny = (m_ly-2*u)/delta;
  // nz = (m_lz-2*u)/delta;

  dirx = diry = 1;
  x = -0.5*m_lx + u;
  y = -0.5*m_ly + u;
  z = m_dim==3?-0.5*m_lz + u:0.;
  counter = 0 ;
  for (;;) {

    // place one particle
    gsl_matrix_set(m_x, counter, 0, x);
    gsl_matrix_set(m_x, counter, 1, y);
    gsl_matrix_set(m_x, counter, 2, z);
    counter += 1;

    // if all molecules placed exit
    if (counter == m_npart) break;

    // update coordinates
    double xnew = x + dirx*delta;
    if (fabs(xnew) > 0.5*m_lx - u) {
      dirx *= -1;
      double ynew = y + diry*delta;
      if (fabs(ynew) > 0.5*m_ly - u) {
        diry *= -1;
        if (m_dim == 3) {
          double znew = z + delta;
          if (fabs(znew) > 0.5*m_lz-u) {
            // break;
            throw runtime_error("Box is too small to position all particles on a lattice!");
          }
          else {
            z = znew;
          }
        }
        else {
          // break;
          throw runtime_error("Box (2d) is too small to position all particles on a lattice!");
        }
      }
      else {
        y = ynew;
      }
    }
    else {
      x = xnew;
    }
  }
  // for (size_t ix = 1; ix < nx; ++ix){
  //   for (size_t iy = 1; iy < ny; ++iy){
  //     gsl_matrix_set(m_x, counter, 0, x);
  //     gsl_matrix_set(m_x, counter, 1, y);
  //     if (m_dim == 3) {
  //       for (size_t iz = 1; iz < nz; ++iz){
  //         z = -0.5*m_lz + u + iz*delta;
  //         gsl_matrix_set(m_x, counter, 2, z);
  //         counter += 1;
  //
  //         if (counter == m_npart)
  //           break;
  //       }
  //     }
  //     else {
  //       counter += 1;
  //
  //       if (counter == m_npart)
  //         break;
  //     }
  //   }
  //   if (counter == m_npart)
  //     break;
  // }
  //
  // if (counter < m_npart) {
  //   throw runtime_error("Box is too small to position all particles on a lattice!");
  // }

  cout << "finished lattice init" << endl;
  return;
}

void MDWorld::init_velocities(gsl_rng *rng){
  /*
   * Initialize velocities
   */
  // declarations
  double vsqm, vvar;
  gsl_vector *vm(0);
  double fs;

  // initializations
  vsqm = 0.;
  vm = gsl_vector_calloc(3);

  // random velocities
  for (size_t n=0; n<m_npart; ++n){
    gsl_vector_view v = gsl_matrix_row(m_v, n);
    for (size_t i=0; i<m_dim; i++){
      gsl_vector_set(&v.vector, i, gsl_rng_uniform(rng));
    }

    // update average velocity
    linalg_daxpy(1., &v.vector, vm);

    // update average square velocity
    vsqm += linalg_ddot(&v.vector, &v.vector);
  }

  // compute mean velocity
  linalg_dscal(1./m_npart, vm);

  // compute variance
  vsqm /= m_npart;
  vvar = vsqm - linalg_ddot(vm, vm);  // var = <v^2> - <v>^2

  // shift and rescale velocities
  if (m_npart == 1){
    fs = sqrt(m_dim * m_temp / vsqm);
      gsl_vector_view v = gsl_matrix_row(m_v, 0);
      linalg_dscal(fs, &v.vector);
  }
  else {
    fs = sqrt(m_dim * m_temp / vvar);
    for (size_t n=0; n<m_npart; ++n){
      gsl_vector_view v = gsl_matrix_row(m_v, n);
      // shift
      linalg_daxpy(-1., vm, &v.vector);

      // rescale
      linalg_dscal(fs, &v.vector);
    }
  }

  // exit
  gsl_vector_free(vm);
}

void MDWorld::set_constraints() {
  /*
   * Set constraints.
   */

  constrain_x();
  constrain_v();

  return;
}

void MDWorld::update_energy_forces(){
  /*
   * compute all forces
   */
  m_energy_pot = 0;
  gsl_matrix_set_all(m_forces, 0.);

  /* iterate over force fields */
  for (vector<ForceField*>::iterator it=m_ffields.begin(); it!=m_ffields.end(); ++it){
      // cout << typeid((*it)).name() << endl;
      (*it)->energy_force(m_x, &m_energy_pot, m_forces);
  }

  return;
}

void MDWorld::constrain_x(){
  /*
   * Apply constraints on positions
   */

  /* iterate over constraints */
  for (vector<Constraint*>::iterator it=m_constraints.begin(); it!=m_constraints.end(); ++it){
      (*it)->constrain_x(m_x);
  }

  return;
}

void MDWorld::constrain_v(){
  /*
   * Apply constraints on velocities
   */
  /* iterate over constraints */
  for (vector<Constraint*>::iterator it=m_constraints.begin(); it!=m_constraints.end(); ++it){
      (*it)->constrain_v(m_v);
  }

  return;
}

void MDWorld::update_energy_kinetics(){
  /*
   * compute all forces
   */

  m_energy_kin = 0.5*m_mass*linalg_ddot(m_v, m_v);

  return;
}

bool MDWorld::check_neighbor_lists(){
  /*
   * Check neighbor lists
   */

  bool res=false;

  /* iterate over neighbor lists */
  for (vector<NeighborList*>::iterator it=m_neighbors.begin(); it!=m_neighbors.end(); ++it){
      bool res_tp = (*it)->check(m_x);
      res = res || res_tp;
  }

  return res;
}

void MDWorld::update_neighbor_lists(size_t iter){
  /*
   * Build neighbor lists
   */

  /* iterate over neighbor lists */
  for (vector<NeighborList*>::iterator it=m_neighbors.begin(); it!=m_neighbors.end(); ++it){
      (*it)->build(m_x, m_bonds);
      (*it)->update_counter(iter);
  }

  return;
}

void
MDWorld::dump_thermo(ostream &mystream){
  /*
   * Dump thermodynamic quantities.
   */

  mystream << left << dec;
  mystream << setw(10) << setprecision(0) << fixed << noshowpos << m_npart;
  mystream << setw(18) << setprecision(8) << scientific << showpos << m_energy_pot;
  mystream << setw(18) << setprecision(8) << scientific << showpos << m_energy_kin;
  return;
}

void
MDWorld::print_infos(ostream &mystream) {
  /*
   * print some information on the world object
   */

  mystream << left << dec << fixed;
  mystream << setw(20) << "npart" << setw(20) << setprecision(0) << noshowpos << m_npart << endl;
  mystream << setw(20) << "gamma" << setw(20) << setprecision(6) << showpos << m_gamma << endl;
  mystream << setw(20) << "temp" << setw(20) << setprecision(6) << showpos << m_temp << endl;
  mystream << setw(20) << "mass" << setw(20) << setprecision(6) << showpos << m_mass << endl;
  mystream << setw(20) << "box" << setw(20) << setprecision(6) << showpos << m_lx << setw(20) << m_ly << setw(20) << m_lz << endl;
  return;
}

void
MDWorld::dump_pos(ostream &mystream, bool index, bool positions, bool velocities, bool forces){
  /*
   * Dump configuration
   */

	if (m_npart == 0)
		throw invalid_argument("MDWorld is empty!");

  mystream << left << dec << fixed;

	for (size_t i=0; i<m_npart;i++){
    if (index) {
      mystream << setw(10) << setprecision(0) << noshowpos << i;
    }
    if (positions) {
      mystream << setw(18) << setprecision(8) << showpos << gsl_matrix_get(m_x,i,0);
      mystream << setw(18) << setprecision(8) << showpos << gsl_matrix_get(m_x,i,1);
      mystream << setw(18) << setprecision(8) << showpos << gsl_matrix_get(m_x,i,2);
    }
    if (velocities) {
      mystream << setw(18) << setprecision(8) << showpos << gsl_matrix_get(m_v,i,0);
      mystream << setw(18) << setprecision(8) << showpos << gsl_matrix_get(m_v,i,1);
      mystream << setw(18) << setprecision(8) << showpos << gsl_matrix_get(m_v,i,2);
    }
    if (forces) {
      mystream << setw(18) << setprecision(8) << showpos << gsl_matrix_get(m_forces,i,0);
      mystream << setw(18) << setprecision(8) << showpos << gsl_matrix_get(m_forces,i,1);
      mystream << setw(18) << setprecision(8) << showpos << gsl_matrix_get(m_forces,i,2);
    }
		mystream << endl;
	}
  return;
}

void
MDWorld::dump_dat(std::string fileout) {
  /*
   * Dump configuration in .dat format.
   */
  ofstream fout;
  fout.open(fileout.c_str());

  dump_pos(fout, true, true, true, true);

  fout.close();
  return;
}

void
MDWorld::dump_xyz(std::string fileout, size_t iter) {
  /*
   * Dump configuration in .xyz format.
   */
  ofstream fout;
  fout.open(fileout.c_str());

  fout << left << dec << fixed << setw(10) << setprecision(0) << noshowpos << m_npart << endl;
  fout << "Atoms. Timestep:" << setw(10) << iter << endl;
  dump_pos(fout, true, true, false, false);

  fout.close();
  return;
}

void
MDWorld::load_dat(std::string filein) {
  /*
   * Load configuration in .dat format.
   */
  /* declarations */
  ifstream fin;
  stringstream convert;
  string line;

  size_t n;
  double rx, ry, rz, vx, vy, vz;

  /* initializations */
  fin.open(filein.c_str());

  /* reading file */
  while (getline(fin, line)) {
    convert.clear();
    convert.str(line);
    if (!(
          (convert >> n) &&
          (convert >> rx) &&
          (convert >> ry) &&
          (convert >> rz) &&
          (convert >> vx) &&
          (convert >> vy) &&
          (convert >> vz)
         )
       ) {
      throw runtime_error("Problem in configuration import");
    }

    /* checks */
    if ( !(n < m_npart) ){
      cout << "n = " << n << endl;
      throw runtime_error("Problem in configuration import: n too large!");
    }

    // set position
    gsl_matrix_set(m_x, n, 0, rx);
    gsl_matrix_set(m_x, n, 1, ry);
    gsl_matrix_set(m_x, n, 2, rz);

    // set velocity
    gsl_matrix_set(m_v, n, 0, vx);
    gsl_matrix_set(m_v, n, 1, vy);
    gsl_matrix_set(m_v, n, 2, vz);

    // verbose
    cout << "particle " << n << " loaded." << endl;
  }

  fin.close();
  return;
}

void
MDWorld::load_xyz(std::string filein) {
  /*
   * Load configuration in .xyz format.
   */
  /* declarations */
  ifstream fin;
  stringstream convert;
  string line;

  size_t n,skiprows;
  double rx, ry, rz;

  /* initializations */
  skiprows=2;
  fin.open(filein.c_str());

  /* reading file */
  for (size_t i=0; i<skiprows; ++i){
    getline(fin,line);
  }

  while (getline(fin, line)) {
    convert.clear();
    convert.str(line);
    if (!(
          (convert >> n) &&
          (convert >> rx) &&
          (convert >> ry) &&
          (convert >> rz)
         )
       ) {
      throw runtime_error("Problem in configuration import");
    }

    /* checks */
    if ( !(n < m_npart) ){
      cout << "n = " << n << endl;
      throw runtime_error("Problem in configuration import: n too large!");
    }

    // set position
    gsl_matrix_set(m_x, n, 0, rx);
    gsl_matrix_set(m_x, n, 1, ry);
    gsl_matrix_set(m_x, n, 2, rz);

    // verbose
    cout << "particle " << n << " loaded." << endl;
  }

  fin.close();
  return;
}

void
MDWorld::dump_neighbors(ostream &mystream){
  /*
   * Dump neighbor list
   */

  // initialize
  // size_t n, m;
  // double r;
  // gsl_vector *xtp(0);
  // xtp = gsl_vector_calloc(3);

	if (m_npart == 0)
		throw invalid_argument("MDWorld is empty!");


  mystream << left << dec << fixed;
  for (vector<NeighborList*>::iterator it=m_neighbors.begin(); it!=m_neighbors.end(); ++it){
    (*it)->dump(mystream);
    // for (size_t i=0; i<(*it)->m_npair; i++){
    //   n = gsl_matrix_uint_get((*it)->m_pairs,i,0);
    //   m = gsl_matrix_uint_get((*it)->m_pairs,i,1);
    //   gsl_vector_view xn = gsl_matrix_row(m_x, n);
    //   gsl_vector_view xn_ref = gsl_matrix_row((*it)->m_xref, n);
    //   gsl_vector_memcpy(xtp, &xn_ref.vector);
    //   linalg_daxpy(-1, &xn.vector, xtp);
    //   r = linalg_dnrm2(xtp);
    //
    //   mystream << setw(10) << setprecision(0) << noshowpos << i;
    //   mystream << setw(10) << setprecision(0) << noshowpos << n;
    //   mystream << setw(10) << setprecision(0) << noshowpos << m;
    //   mystream << setw(20) << setprecision(6) << noshowpos << r;
    //   for (size_t a=0; a<3; ++a){
    //     mystream << setw(20) << setprecision(6) << noshowpos << gsl_vector_get(&xn.vector,a);
    //   }
    //   for (size_t a=0; a<3; ++a){
    //     mystream << setw(20) << setprecision(6) << noshowpos << gsl_vector_get(&xn_ref.vector,a);
    //   }
    //   mystream << endl;
    // }

      mystream << endl;
  }
  // gsl_vector_free(xtp);
  return;
}

void
MDWorld::dump_neighbors(string fileout){
  /*
   * Dump neighbor list
   */
  ofstream fout;
  fout.open(fileout.c_str());

  dump_neighbors(fout);

  fout.close();

  return;
}

void
MDWorld::dump_polarity_vectors(ostream &mystream){
  /*
   * Dump polarity vectors
   */

  mystream << left << dec << fixed;
  for (vector<ForceField*>::iterator it=m_ffields.begin(); it!=m_ffields.end(); ++it){
    // PolarPairLJ
    {
      PolarPairLJ *ffield(0);
      ffield = dynamic_cast<PolarPairLJ*>(*it);
      if ( ffield != nullptr ){
        gsl_matrix *pol_vec = ffield->m_pol_vec;

        for (size_t i=0; i<pol_vec->size1;i++){
          mystream << setw(10) << setprecision(0) << noshowpos << i;
          mystream << setw(18) << setprecision(8) << showpos << gsl_matrix_get(pol_vec,i,0);
          mystream << setw(18) << setprecision(8) << showpos << gsl_matrix_get(pol_vec,i,1);
          mystream << setw(18) << setprecision(8) << showpos << gsl_matrix_get(pol_vec,i,2);
          mystream << endl;
        }
      }
    }

    // PolarPair48
    {
      PolarPair48 *ffield(0);
      ffield = dynamic_cast<PolarPair48*>(*it);
      if ( ffield != nullptr ){
        gsl_matrix *pol_vec = ffield->m_pol_vec;

        for (size_t i=0; i<pol_vec->size1;i++){
          mystream << setw(10) << setprecision(0) << noshowpos << i;
          mystream << setw(18) << setprecision(8) << showpos << gsl_matrix_get(pol_vec,i,0);
          mystream << setw(18) << setprecision(8) << showpos << gsl_matrix_get(pol_vec,i,1);
          mystream << setw(18) << setprecision(8) << showpos << gsl_matrix_get(pol_vec,i,2);
          mystream << endl;
        }
      }
    }
  }
  return;
}

void
MDWorld::dump_polarity_vectors(string fileout){
  /*
   * Dump polarity vectors
   */
  ofstream fout;
  fout.open(fileout.c_str());

  dump_polarity_vectors(fout);

  fout.close();

  return;
}
