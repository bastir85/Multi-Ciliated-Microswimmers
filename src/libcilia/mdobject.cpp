//#define DEBUG
#include "mdobject.h"
#include "MatrixIO.h"
#include "global.h"
#include <assert.h>
#include <fstream>

void MDObject::spring(int i, int j, double b, double k) {
  // Attention: No Atomic since we parallelize per object
  Row<double, 1> pos1 = positions.get_row_by_idx(i);
  Row<double, 1> pos0 = positions.get_row_by_idx(j);
  Row<double, 1> force = forces.get_row_by_idx(i);
  double fac = -2.0 * k * (1.0 - b / sqrt(pos1.sqdistance(pos0)));
  force[0] += (pos1[0] - pos0[0]) * fac;
  force[1] += (pos1[1] - pos0[1]) * fac;
  force[2] += (pos1[2] - pos0[2]) * fac;
}

Vector3d MDObject::spring(const Matrix<double, 1> &pos1,
                          const Matrix<double, 1> &pos0, double b, double k) {
  Vector3d force = {0., 0., 0.};
  double d = sqrt(pos1.sqdistance(pos0));
  if (d != 0) {
    double fac = -2.0 * k * (1.0 - b / d);
    force = (pos1 - pos0) * fac;
  }
  return force;
  // return (pos1 - pos0)*fac;
}

//-------------------------------------------------------
//--------------move particle-------------------------------
//-----------------------------------------------------
// noch nicht kontroliert !!!!
// und verlangt angabe der kraftpromasseneinheit!
void MDObject::md_loop_start() { // used to be partmove
// is new functions inner_loop
#pragma omp for schedule(static) collapse(2) nowait
  for (int i = 0; i < num_rods; ++i) {
    for (int j = startidx; j < rodlength; ++j) {
      velocities[i][j].add_mul3(forces[i][j], wc.ststep / 2.0);
      mpcd.move(positions[i][j].data(), velocities[i][j].data(), wc.ststep);
      forces[i][j] = 0.0;
    }

    // is new functions inner_loop
  }
}

void MDObject::md_loop(void) {
  // haubtmove, ststeps
#pragma omp for schedule(static) collapse(2) nowait
  for (int i = 0; i < num_rods; ++i) {
    for (int j = startidx; j < rodlength; ++j) {
      velocities[i][j].add_mul3(forces[i][j], wc.ststep);
      mpcd.move(positions[i][j].data(), velocities[i][j].data(), wc.ststep);
      forces[i][j] = 0.0;
    }
  }
}

void MDObject::md_loop_stop(void) {
// letzter v move
#pragma omp for schedule(static) collapse(2) nowait
  for (int i = 0; i < num_rods; ++i) {
    for (int j = startidx; j < rodlength; ++j) {
      velocities[i][j].add_mul3(forces[i][j], wc.ststep / 2.0);
      mpcd.add_impuls(positions[i][j].data(), velocities[i][j].data());
      forces[i][j] = 0.0;
    }
  }
}

Matrix<double, 2> MDObject::get_avarage_rod() {
  Matrix<double, 2> ava(rodlength, 3);
  for (int i = 0; i < num_rods; ++i) {
    ava += positions[i];
  }
  ava /= num_rods;
  return ava;
}
// geeneric IO Functions from cilia.cpp and sphere.cpp
void MDObject::writedata() {
  if (write_data == false)
    return;
  string datname = name + ".xyz";
  ofstream data(datname.c_str(), ios::out | ios::app);
  while (!data.good()) {
    data.clear();
    data.open(datname.c_str(), ios::out | ios::app);
    cerr << "writing to datname.c_str(),ios::out|ios::app faild" << endl;
    if (data.good())
      cerr << " writing works" << endl;
  }
  // iconversion to xyz data...
  data << rodlength * num_rods << endl;
  data << "rods: " << wc.globalt << endl;
  // Cilia xyz data...
  for (int i = 0; i < num_rods; ++i) {
    for (int j = 0; j < rodlength; ++j) {
      data << "C ";
      data << positions[i][j][0] << " " << positions[i][j][1] << " "
           << positions[i][j][2] << endl;
    }
  }
  data.close();
}

// checked against cilia, traj,sperm
void MDObject::load(istream &in) {
  string garb;
  int rodlength, numberofrods;
  REED2(numberofrods, in);
  REED2(rodlength, in);
  assert(this->rodlength == rodlength);
  assert(this->num_rods == numberofrods);
  in >> positions;
  in >> garb;
  cout << garb << endl;
  in >> velocities;
}

void MDObject::save() {
  string sname = name + ".dat";
  cout << "basic save " << name << endl;
  fstream out(sname, ios::out);
  out.precision(100);
  save(out);
}
void MDObject::save(fstream &out) {
  int numberofrods = num_rods;

  out.precision(100);
  WRITE2(numberofrods);
  WRITE2(rodlength);
  out << positions;
  out << "#velocities" << endl;
  out << velocities;
}

bool MDObject::open_datfile() {
  string filename = name + ".dat";
  fstream in(filename, ios::in);
  if (not in.good()) {
    return false;
  } else {
    load(in);
    loaded = true;
  }
  return true;
};

void MDObject::init() {
  for (int i = 0; i < num_rods; ++i) {
    for (int j = startidx; j < rodlength; ++j) {
      velocities[i][j][3] = partmass;
      mpcd.add_impuls(positions[i][j].data(), velocities[i][j].data());
    }
  }
}

Matrix<double, 2> MDObject::getBoundingBox() {
  Matrix<double, 2> minmax(2, 3);
  minmax[0] = positions[0][0]; // min
  minmax[1] = positions[0][0]; // max
  return minmax;
}

void MDObject::calc_total_force() {
  Vector3d total_force;
  Vector3d total_vel;
  total_force = 0;
  total_vel = 0;
  for (int j = 0; j < rodlength; ++j) {
    for (int i = 0; i < num_rods - 1; ++i) {
      total_force += forces[i][j];
    }
  }
  total_force += forces[num_rods - 1][0];

  if (dot_product(total_force, total_force) > 10e-20) {
    cout << "TF: " << total_force;
    cout << "  <<====" << endl;
  }
}
