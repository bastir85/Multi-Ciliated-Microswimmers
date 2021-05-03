#include "mpcd.h"
#include "MatrixIO.h"
#include "basicf.h"
#include "global.h"
#include <assert.h>
#include <fstream>

#ifndef MPCD_CPP
#define MPCD_CPP

#define MPCD_Y_DIST 70.71

inline double rgamma(int Nc) {
  // Marsaglia-Tsang algorithm of gamma random number
  const double d =
      1.5 * Nc -
      1.833333333333333333333333333333333333333333333333333333333333333; // 3
                                                                         // (Nc
                                                                         // - 1)
                                                                         // / 2
                                                                         // - 1
                                                                         // / 3
  const double c = 0.33333333333333333333333333333333333333333333333333333333 /
                   sqrt(d); // (1.0 / 3.0) / sqrt(d);
  double x, v;
#ifdef PARALLEL
  int tid = omp_get_thread_num();
#else
  int tid = 0;
#endif
  for (;;) {
    do {
      x = normalrand();
      v = 1.0 + c * x;
    } while (v <= 0.0);
    v = v * v * v;
    double u = randi[tid]->randDblExc(), x2 = x * x, x4 = x2 * x2;
    if (u < 1.0 - 0.0331 * x4)
      return d * v;
    if (log(u) < 0.5 * x2 + d * (1.0 - v + log(v)))
      return d * v;
  }
} // gamma distributed relative kinetic energy

MPCD::MPCD(int sizex, int sizey, int sizez)
    : virtimpuls(sizex * sizey * sizez, 4),
      numberofboxes(sizex * sizey * sizez), size {
  sizex, sizey, sizez
}
#ifdef XY_WALL_MAP
, lowerwall(LW, (double)sizex / resolution),
    walls(WALLS, (double)sizex / resolution)
#endif

{
  cout << size[0] << size[1] << size[2] << endl;
#ifdef NOSLIP
  cout << "top with NOSLIP" << endl;
#ifdef SLIP
  assert(0);
#endif
#endif
#ifdef SLIP
  cout << "top with SLIP" << endl;
#ifdef NOSLIP
  assert(0);
#endif
#endif

#ifdef PARALLEL
  int max_threads = omp_get_max_threads();
  cout << max_threads << " " << numberofboxes << " ";
#else
  int max_threads = 1;
#endif
  boxes = std::vector<std::vector<std::vector<double *>>>(max_threads);
  for (int tid = 0; tid < max_threads; tid++) {
    boxes[tid] = std::vector<std::vector<double *>>(numberofboxes);
    for (int box_id = 0; box_id < numberofboxes; box_id++) {
      boxes[tid][box_id] = std::vector<double *>(16);
      boxes[tid][box_id].clear();
    }
  }

  l_energy = new double[numberofboxes];
  fstream in(filename, ios::in);
  if (in.good()) {
    in >> a[0];
    in >> a[1];
    in >> a[2];
  } else {
    cout << "init with random A vector" << endl;
    a[0] = rnd();
    a[1] = rnd();
    a[2] = rnd();
  }
  in.close();
}
#ifdef XY_WALL
void MPCD::set_cfg(const libconfig::Setting &cfg) {
  cfg.lookupValue("radius_of_curvature", radius_of_curvature);
  cout << "Radius of curvature = " << radius_of_curvature << endl;
}
#endif

#ifdef XY_WALL_MAP
inline bool MPCD::WALLS(double x, double y) {
  return LW(x) > y or LW(x) + MPCD_Y_DIST < y;
}
#endif

#ifdef XY_WALL
inline double MPCD::LW(double x) const {
  const double L = wc.sizex;
  const double a = radius_of_curvature;
  const double aS = a / sqrt(2);
  // cout << "L="<<L<< " a="<<radius_of_curvature<<" aS="<<aS<<endl;
  if ((x > aS) and (x <= L / 2 - aS))
    return x;
  if ((x > L / 2 - aS) and (x < L / 2 + aS))
    return sqrt(a * a - (x - L / 2) * (x - L / 2)) + L / 2 - 2 * aS;
  if ((x > L / 2 + aS) and (x <= L - aS))
    return L / 2 - (x - L / 2);
  if (x > L - aS)
    return 2 * aS - sqrt(a * a - (x - L) * (x - L));
  return 2 * aS - sqrt(a * a - x * x);
}
#endif

inline int MPCD::box(const double *position) {
  return int(floor((position[0] + a[0]) / anull) -
             floor((position[0] + a[0]) / anull / size[0]) * size[0] +
             size[0] * floor((position[1] + a[1]) / anull) -
             floor((position[1] + a[1]) / anull / size[1]) * size[0] * size[1] +
             size[0] * size[1] * floor((position[2] + a[2]) / anull) -
             floor((position[2] + a[2]) / anull / size[2]) * size[0] * size[1] *
                 size[2]);
}

void MPCD::next_tstep(void) {

#ifdef PARALLEL
  int max_threads = omp_get_max_threads();
#else
  int max_threads = 1;
#endif

#pragma omp for schedule(static) nowait
  for (int cid = 0; cid < numberofboxes; cid++) {
    double ux = 0.0;
    double uy = 0.0;
    double uz = 0.0;

    int num = 0;
    int masse = 0;
    l_energy[cid] = 0.0;

    for (int tid = 0; tid < max_threads; tid++) {
      for (unsigned int vid = 0; vid < boxes[tid][cid].size(); vid++) {
        double *velocity = boxes[tid][cid][vid];
        const double particle_mass = velocity[3];
        ux += velocity[0] * particle_mass;
        uy += velocity[1] * particle_mass;
        uz += velocity[2] * particle_mass;
        masse += particle_mass;
        l_energy[cid] += particle_mass * (velocity[0] * velocity[0] +
                                          velocity[1] * velocity[1] +
                                          velocity[2] * velocity[2]);
      }
      // count
      num += boxes[tid][cid].size();
    }

    if (num > 0) {
#if defined(SLIP) or defined(NOSLIP)
      // top or bottom boxes
      if ((cid < size[0] * size[1]) or
          (cid > size[0] * size[1] * (size[2] - 1))) {
        if (num < rho) {
          const int mass_diff = rho - num;
          virtimpuls[cid][1] = sqrt(mass_diff * temp) * normalrand();
          virtimpuls[cid][2] = sqrt(mass_diff * temp) * normalrand();
          virtimpuls[cid][3] = sqrt(mass_diff * temp) * normalrand();
          virtimpuls[cid][0] = mass_diff;

          ux += virtimpuls[cid][1];
          uy += virtimpuls[cid][2];
          uz += virtimpuls[cid][3];
          // add fluid parts mass =1
          masse += mass_diff;
          num += mass_diff;
          l_energy[cid] += (virtimpuls[cid][1] * virtimpuls[cid][1] +
                            virtimpuls[cid][2] * virtimpuls[cid][2] +
                            virtimpuls[cid][3] * virtimpuls[cid][3]) /
                           mass_diff;
        }
      }
#endif

      ux /= masse;
      uy /= masse;
      uz /= masse;
      double rel_energy = l_energy[cid] - (ux * ux + uy * uy + uz * uz) * masse;
      double sc;
      if (rel_energy > 0 and num > 1) {
        sc = sqrt(rgamma(num) / (rel_energy * 0.5));
      } else {
        assert(rel_energy * rel_energy < 1e-8);
        sc = 1.0;
      }

      Matrix<double, 2> RotMat = generate_collision_matrix();

      for (int tid = 0; tid < max_threads; tid++) {
        for (unsigned int vid = 0; vid < boxes[tid][cid].size(); vid++) {
          double *velocity = boxes[tid][cid][vid];
          double vrel[3];
          vrel[0] = velocity[0] - ux;
          vrel[1] = velocity[1] - uy;
          vrel[2] = velocity[2] - uz;

          velocity[0] = ux + (RotMat[0][0] * vrel[0] + RotMat[0][1] * vrel[1] +
                              RotMat[0][2] * vrel[2]) *
                                 sc;
          velocity[1] = uy + (RotMat[1][0] * vrel[0] + RotMat[1][1] * vrel[1] +
                              RotMat[1][2] * vrel[2]) *
                                 sc;
          velocity[2] = uz + (RotMat[2][0] * vrel[0] + RotMat[2][1] * vrel[1] +
                              RotMat[2][2] * vrel[2]) *
                                 sc;

          l_energy[cid] -= velocity[3] * (velocity[0] * velocity[0] +
                                          velocity[1] * velocity[1] +
                                          velocity[2] * velocity[2]);
        }
      }
    }

    for (int tid = 0; tid < max_threads; tid++) {
      boxes[tid][cid].clear();
    }
  }
#pragma omp for reduction(+ : diss_energy_total)
  for (int cid = 0; cid < numberofboxes; cid++) {
    diss_energy_total += l_energy[cid];
  }
}

void MPCD::shift_boxes() {
  a[0] = rnd();
  a[1] = rnd();
  a[2] = rnd();
}

// TODO: these locks use a lot of time maybe collect per thread and than sum
// over threads?!
void MPCD::add_impuls(const double *position, double *velocity) {
  int cell = box(position);
#ifdef PARALLEL
  int tid = omp_get_thread_num();
#else
  int tid = 0;
#endif

  boxes[tid][cell].push_back(velocity);

}

inline Matrix<double, 2> MPCD::generate_collision_matrix() {
  Matrix<double, 2> RotMat(3, 3);
#ifdef PARALLEL
  int tid = omp_get_thread_num();
#else
  int tid = 0;
#endif

  double phi = 2 * PI * randi[tid]->rand();
  double theta = 2.0 * randi[tid]->rand() - 1.0;

  double sqrt_1_min_theta2 = sqrt(1.0 - theta * theta);
  double R_x = sqrt_1_min_theta2 * cos(phi);
  double R_y = sqrt_1_min_theta2 * sin(phi);
  double R_z = theta;

  RotMat[0][0] = R_x * R_x + (1.0 - R_x * R_x) * calpha;
  RotMat[0][1] = R_x * R_y * (1.0 - calpha) - R_z * salpha;
  RotMat[0][2] = R_x * R_z * (1.0 - calpha) + R_y * salpha;
  RotMat[1][0] = R_x * R_y * (1.0 - calpha) + R_z * salpha;
  RotMat[1][1] = R_y * R_y + (1.0 - R_y * R_y) * calpha;
  RotMat[1][2] = R_y * R_z * (1.0 - calpha) - R_x * salpha;
  RotMat[2][0] = R_x * R_z * (1.0 - calpha) - R_y * salpha;
  RotMat[2][1] = R_y * R_z * (1.0 - calpha) + R_x * salpha;
  RotMat[2][2] = R_z * R_z + (1.0 - R_z * R_z) * calpha;
  return RotMat;
}

#ifdef BOUNDS
int MPCD::outofbounds(const Row<double, 1> &position) const {
#if XY_CIRCLE
  double lx =
      position[0] - floor(position[0] / size[0]) * size[0] - wc.sizex / 2.;
  double ly =
      position[1] - floor(position[1] / size[1]) * size[1] - wc.sizey / 2.;
  return position[2] < 0 or position[2] > (size[2] - 1) * anull or
         ly * ly < radius_of_curvature * radius_of_curvature - lx * lx;
#endif

#if defined(UPPERWALL) && defined(LOWERWALL)
  int sizex = size[0];
  int sizey = size[1];
  int sizez = size[2];
  if (position[2] < 0)
    return 1; // bottom wall
  if (position[2] > (size[2] - 1) * anull)
    return 2; // top wall
#ifdef XY_WALL
  double x = position[0]; // zizack function of x, periodic in sizex, von -1/4
                          // auf +1/4, minimum bei x%sizex==0
  if (LOWERWALL > position[1])
    return 3;
  if (UPPERWALL < position[1])
    return 4;
#endif
  return 0;

#else
#ifdef XY_WALL
  double y_wall = LW(position[0] - floor(position[0] / size[0]) * size[0]);
#endif
  return position[2] < 0 or position[2] > (size[2] - 1) * anull
#ifdef XY_WALL
         or y_wall > position[1] or position[1] > y_wall + MPCD_Y_DIST;
#else
      ;
#endif
#endif
}
#endif

void MPCD::move(double *position, double *velocity, const double ts) {
  position[0] += velocity[0] * ts;
  position[1] += velocity[1] * ts;
  position[2] += velocity[2] * ts;

  // Bounce Back:
#if defined(SLIP) || defined(NOSLIP)
  if (position[2] > (size[2] - 1) * anull) { // Watch order and sign!!!
#ifdef NOSLIP
    double t2 = (position[2] - (size[2] - 1) * anull) / velocity[2];
    velocity[0] = 2 * vwx - velocity[0];
    velocity[1] = 2 * vwy - velocity[1];
    velocity[2] *= (-1);

    position[0] += 2.0 * t2 * (velocity[0] - vwx);
    position[1] += 2.0 * t2 * (velocity[1] - vwy);
    position[2] += 2.0 * t2 * velocity[2];
    // impulsuebertrag auf top wall. ( = Kraft auf Topwall * dt)
    // TODO: wall implementaition!!!
    // wallximpuls+=2*vwx-2*xvel[i];
    // wallyimpuls+=2*vwy-2*yvel[i];
#endif
#ifdef SLIP
    // Slip boundary at the top:
    velocity[2] *= (-1.0);
    position[2] = -position[2] + (size[2] - 1) * anull * 2.0;
#endif
  }
  if (position[2] < 0) { // Watch order and sign!!!
    double t2 = (position[2] / velocity[2]);
    velocity[0] *= (-1.0);
    velocity[1] *= (-1.0);
    velocity[2] *= (-1.0);
    position[0] += velocity[0] * 2.0 * t2;
    position[1] += velocity[1] * 2.0 * t2;
    position[2] += velocity[2] * 2.0 * t2;
  }
#endif
// xy walls from sperm code
#ifdef BOUNDS
  if (outofbounds(position)) {
    velocity[0] *= (-1.0);
    velocity[1] *= (-1.0);
    velocity[2] *= (-1.0);
    // position += velocity*ts;
    position[0] += velocity[0] * ts;
    position[1] += velocity[1] * ts;
    position[2] += velocity[2] * ts;
  }
#if defined(SLIP) || defined(NOSLIP)
  cerr << "bounds and slip, noslip configured" << endl;
  throw "Config Error";
#endif
#endif
}

MPCD::~MPCD() {

  fstream out(filename, ios::out);
  out.precision(100);
  out << a;
  out.close();
  cout << "mpcd saved";
}
void MPCD::debug() {
  double ux = 0.0;
  double uy = 0.0;
  double uz = 0.0;
  double g_energy = 0.0;
  double dissp_energy = 0.0;
  int num = 0;
  int masse = 0;

  double virtimp[4] = {0.0};

#ifdef PARALLEL
  int max_threads = omp_get_max_threads();
  cout << max_threads << " ";
#else
  int max_threads = 1;
#endif
  for (int cid = 0; cid < numberofboxes; cid++) {
    for (int tid = 0; tid < max_threads; tid++) {
      for (unsigned int vid = 0; vid < boxes[tid][cid].size(); vid++) {
        double *velocity = boxes[tid][cid][vid];
        double m = velocity[3];
        ux += velocity[0] * m;
        uy += velocity[1] * m;
        uz += velocity[2] * m;

        masse += m;
        g_energy += m * (velocity[0] * velocity[0] + velocity[1] * velocity[1] +
                         velocity[2] * velocity[2]);
      }

      // count
      num += boxes[tid][cid].size();
      int cap = boxes[tid][cid].capacity();
      if (cap > 32) {
        cout << "cid=" << cid << "," << tid
             << " size=" << boxes[tid][cid].size()
             << "cap =" << boxes[tid][cid].capacity() << endl;
      }
    }

    virtimp[0] += virtimpuls[cid][0];
    virtimp[1] += virtimpuls[cid][1];
    virtimp[2] += virtimpuls[cid][2];
    virtimp[3] += virtimpuls[cid][3];
    dissp_energy += l_energy[cid];
  }

  ofstream data("mpcd.san", ios::out | ios::app);
  cout << "------------------------------------------------------------"
       << endl;
  cout << "Total Momentum \t= " << ux / masse << " " << uy / masse << " "
       << uz / masse << endl;
  cout << "Total Number \t= " << num << endl;
  cout << "Total Mass \t= " << masse << endl;
  cout << "kbT \t= " << g_energy / 3.0 / num << endl << endl;

  cout << "Energy diff = " << dissp_energy << endl << endl;
  cout << "Energy total = " << diss_energy_total << endl << endl;
  data << wc.globalt << " " << diss_energy_total << " " << num << " " << masse
       << " " << ux / masse << " " << uy / masse << " " << uz / masse << " "
       << g_energy / 3.0 / num << endl;
  data.close();
}
#endif
