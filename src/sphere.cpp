//#define DEBUG
#include "sphere.h"
#include "basicf.h"
#include "global.h"
#include <assert.h>
#include <fstream>
using namespace Numeric_lib;
using namespace std;

#include "MatrixIO.h"

Sphere::Sphere(double radius, double springconstant, double wx, double wy,
               double wz, int npart, MPCD &mpcd, GlobalForces *gf, string name)
    : MDObject(name, 1, npart, 0, mpcd, gf), springconstant(springconstant),
      radius(radius), tforce(3), ava_vel(3), max_flowF(2, 3), npart(npart) {
  cout << "construtct " << name << " " << npart << endl;
  int fill_factor = 0;
  if (npart == 163) {
    fill_factor = 2;
  } else if (npart == 643) {
    fill_factor = 3;
  } else {

    cerr << "Wrong nparts has to 163 or 643" << endl;
    throw;
  }
  neilist = new int *[npart];
  for (int i = 0; i < npart; i++) {
    neilist[i] = new int[7];
  }
  if (not open_datfile()) {
    create(wx, wy, wz, fill_factor);
  }

#ifdef SIMPLE_FLOWF
// mu0 =npart/(6*PI*eta*oseen_a);
#endif
  tforce = 0.0;
  cout << "np:" << npart << endl;
  cout << "hb:" << headbond << endl;
}

Sphere::~Sphere() {
  for (int i = 0; i < npart; i++) {
    delete neilist[i];
  }
  delete neilist;
};

void Sphere::save(fstream &out) {
  cout << "saving neiullist also" << endl;
  MDObject::save(out);

  out << "#neilist" << endl;
  for (int i = 0; i < npart; i++) {
    for (int j = 0; j < 7; j++) {
      out << neilist[i][j] << " ";
    }
    out << endl;
  }
  out << "#headbond" << endl;
  out << headbond << endl;
}

void Sphere::load(istream &in) {
  cout << "load neillist for " << rodlength << endl;
  string garb;
  MDObject::load(in);
  in >> garb;
  if (garb != "#neilist") {
    cout << "input error neilist";
    exit(1);
  }
  for (int i = 0; i < rodlength; i++) {
    for (int j = 0; j < 7; j++) {
      in >> neilist[i][j];
    }
  }
  in >> garb;
  if (garb != "#headbond") {
    cout << "input error headbond. got: " << garb << endl;
    exit(1);
  }
  in >> headbond;
}

void Sphere::recalcforce(void) {
  // radials
  //  #pragma omp for schedule (guided) nowait
  for (int i = 1; i < npart; i++) {
    spring(0, i, radius, 0.1 * springconstant);
    spring(i, 0, radius, 0.1 * springconstant);
  }
#ifdef REP_WALL
  assert(true);
  for (int i = 1; i < tnpart; i++) {
    double cutoff = REP_WALL;
    double sigma = cutoff * 0.890898718; // entspricht siqma=cutof/2^(1/6)
    if (zparticle[i] < cutoff) {
      zforces[i] += 24. * (2 * pow(sigma, 12) / pow(zparticle[i], 13) -
                           pow(sigma, 6) / pow(zparticle[i], 7));
    }
    double zdis = (sizez - 1) * anull - zparticle[i];
    if (zdis < cutoff) {
      zforces[i] += -24. * (2 * pow(sigma, 12) / pow(zdis, 13) -
                            pow(sigma, 6) / pow(zdis, 7));
    }
  }
#endif

  // kanten
  //  #pragma omp for schedule (guided) nowait
  for (int i = 1; i < npart; i++) {
    for (int j = 1; j < neilist[i][0] + 1; j++) {
      spring(i, neilist[i][j], headbond, springconstant);
    }
  }

} // end of recalcf

void Sphere::debug() {}

istream &operator>>(istream &in, Sphere &Sphere) {
  /*
   *  usage is >> Sphere;
   */
  Sphere.load(in);

  return in;
}

void Sphere::create(Matrix<double, 1> position, int fill_factor) {
  create(position[0], position[1], position[2], fill_factor);
}

void Sphere::create(double wx, double wy, double wz, int fill_factor) {
  positions = 0;
  velocities = 0.0;

  // c =========================c
  // c     temporary variables
  // c =========================c

  // c     basic icosahedron:
  // 0 at center...
  positions[0][1][0] = 0.;
  positions[0][1][1] = 0.;
  positions[0][1][2] = 1.13;
  positions[0][2][0] = cos(2.0 * PI / 5);
  positions[0][2][1] = sin(2.0 * PI / 5);
  positions[0][2][2] = 0.50230;
  positions[0][3][0] = -cos(PI - 2.0 * 2.0 * PI / 5);
  positions[0][3][1] = sin(PI - 2.0 * 2.0 * PI / 5);
  positions[0][3][2] = 0.50230;
  positions[0][4][0] = -cos(PI - 2 * 2.0 * PI / 5);
  positions[0][4][1] = -sin(PI - 2 * 2.0 * PI / 5);
  positions[0][4][2] = 0.5023;
  positions[0][5][0] = cos(2.0 * PI / 5);
  positions[0][5][1] = -sin(2.0 * PI / 5);
  positions[0][5][2] = 0.50230;
  positions[0][6][0] = 1.0;
  positions[0][6][1] = 0.0;
  positions[0][6][2] = 0.5023;
  positions[0][7][0] = cos(PI - 2 * 2.0 * PI / 5);
  positions[0][7][1] = sin(PI - 2 * 2.0 * PI / 5);
  positions[0][7][2] = -0.5023;
  positions[0][8][0] = -cos(2.0 * PI / 5);
  positions[0][8][1] = sin(2.0 * PI / 5);
  positions[0][8][2] = -0.5023;
  positions[0][9][0] = -1.0;
  positions[0][9][1] = 0.0;
  positions[0][9][2] = -0.5023;
  positions[0][10][0] = -cos(2.0 * PI / 5);
  positions[0][10][1] = -sin(2.0 * PI / 5);
  positions[0][10][2] = -0.5023;
  positions[0][11][0] = cos(PI - 4.0 * PI / 5);
  positions[0][11][1] = -sin(PI - 2 * 2.0 * PI / 5);
  positions[0][11][2] = -0.5023;
  positions[0][12][0] = 0.;
  positions[0][12][1] = 0.;
  positions[0][12][2] = -1.13;

  // c     scale such that bonds are unity
  double scale = 0.8513;
  for (int i = 1; i < 13; i++) {
    positions[0][i][0] *= scale;
    positions[0][i][1] *= scale;
    positions[0][i][2] *= scale;
  }
  npart = 12;
  int tel = 12;
  // fill to Icosahedron of order n
  for (int i = 0; i < fill_factor; i++) { // was 2
    create_neillist(npart);
    // Add more particles

    for (int j = 1; j < npart; j++) {
      for (int k = j + 1; k < npart + 1; k++) {
        for (int l = 1; l < neilist[j][0] + 1; l++) {
          if (neilist[j][l] == k) {
            tel += 1;
            positions[0][tel] = (positions[0][j] + positions[0][k]) / 2.0;
          }
        }
      }
    }
    npart = tel;

    scale = 2.0;
    for (int i = 1; i < npart + 1; i++) {
      positions[0][i] *= scale;
    }

  } // end of fill

  // is in init after create
  create_neillist(npart);

  npart++; // number of particles is also center particle

  // turning to point head1 in x direciton!!!
  // turning 90 around y Axis
  for (int i = 1; i < npart; i++) {
    double xi = positions[0][i][0];
    double zi = positions[0][i][2];
    positions[0][i][0] = zi;
    positions[0][i][2] = -xi;
  }

  // rotate and translate sperm
  // double alphac = cos(alpha);
  // double alphas = sin(alpha);

  // project on sphere with radius: and move to where()
  for (int i = 1; i < npart; i++) {
    scale = radius / sqrt(positions[0][i].sqdistance(positions[0][0]));
    positions[0][i] *= scale;

    double nx = positions[0][i][0]; //*alphac - positions[0][i][2]*alphas;
    double nz = positions[0][i][2]; //*alphas + positions[0][i][2]*alphac;

    positions[0][i][0] = nx + wx;
    positions[0][i][1] += wy;
    positions[0][i][2] = nz + wz;
  }
  positions[0][0][0] += wx;
  positions[0][0][1] += wy;
  positions[0][0][2] += wz;

  // calculate average bondlength
  headbond = 0;
  double normierung = 0;
  for (int i = 1; i < npart; i++) {
    for (int j = 1; j < neilist[i][0] + 1; j++) {
      double dis =
          sqrt(positions[0][i].sqdistance(positions[0][neilist[i][j]]));
      headbond += dis;
      normierung++;
    }
  }
  headbond = headbond / normierung;
}

void Sphere::create_neillist(double max_part) {
  double rcutsq = (1.1 * 1.1);
  // create list again:
  // 0 neigbourlist
  for (int j = 0; j < max_part + 1; j++) {
    for (int k = 0; k < 7; k++) {
      neilist[j][k] = 0;
    }
  }
  cout << "create with " << max_part << endl;
  //  make neighbourlist
  for (int j = 1; j < max_part; j++) {
    for (int k = j + 1; k < max_part + 1; k++) {
      if (positions[0][j].sqdistance(positions[0][k]) < rcutsq) {
        neilist[j][0] += 1;
        neilist[k][0] += 1;
        neilist[j][neilist[j][0]] = k;
        neilist[k][neilist[k][0]] = j;
      }
    }
  }
}

Matrix<double, 2> Sphere::getBoundingBox() {

  Matrix<double, 2> minmax(2, 3);
  minmax[0] = positions[0][0] - radius; // min
  minmax[1] = positions[0][0] + radius; // max

  return minmax;
}

double Sphere::getT() {
  assert(true);
  return 0.0;
}

double Sphere::getPhase() { return 0.0; }

void Sphere::calc_total_force() {
  tforce = 0.0;
  for (int i = 0; i < num_rods; ++i) {
    for (int j = startidx; j < rodlength; ++j) {
      tforce += forces[i][j];
    }
  }

  ava_vel = 0.0;
  for (int i = 0; i < num_rods; ++i) {
    for (int j = startidx; j < rodlength; ++j) {
      ava_vel += velocities[i][j];
    }
  }
  ava_vel /= num_rods * rodlength;
}

#ifdef SIMPLE_FLOWF
Vector3d Sphere::flowField(const Matrix<double, 1> &pos) {
  Vector3d velocity;
  velocity = tensor_velocity(positions[0][0], tforce, pos);

  if (dot_product(max_flowF[1], max_flowF[1]) < velocity[0] * velocity[0] +
                                                    velocity[1] * velocity[1] +
                                                    velocity[2] * velocity[2]) {
    max_flowF[0] = pos;
    Vector_lib::to_matrix(max_flowF[1], velocity);
  }
  return velocity;
}

void Sphere::update_velocities() {

  Vector3d flowF = gf->flowField(positions[0][0]) / mu0;

  for (int i = 0; i < num_rods; ++i) {
    for (int j = startidx; j < rodlength; ++j) {
      Vector3d tmp;

#ifdef NOISE
      double noise = 1.0 * sqrt(2 / mu0);
      newacc[0] += noise * normalrand() / partmass;
      newacc[1] += noise * normalrand() / partmass;
      newacc[2] += noise * normalrand() / partmass;
#endif
      tmp = forces[i][j] + flowF;
      Vector3d tmp2 =
          velocities[i][j] * 2.0 + (oldforces[i][j] + tmp) * wc.ststep;
      to_matrix(velocities[i][j], tmp2, 1.0 / (2 + wc.ststep / partmass / mu0));
      tmp -= velocities[i][j] / mu0 / partmass;
      to_matrix(oldforces[i][j], tmp);
    }
  }
}
#endif
