#ifndef MPCD__H
#define MPCD__H

#include <cmath>
#include <vector>

#ifdef PARALLEL
#include <omp.h>
#endif

#include "Matrix.h"
#include "abbrev.cpp"
#include <libconfig.h++>


using namespace Numeric_lib;
using namespace std;

class MPCD {
  friend class Fluid;

public:
  MPCD(int sizex, int sizey, int sizez);
  ~MPCD();
  void debug();
  void move(double *position, double *velocity, const double ts);
  void add_impuls(const double *positions, double *velocity);

  void next_tstep(void);

  void shift_boxes();
#ifdef TEST_CASES
  Row<double, 1> test_getcellimp(int idx);
#endif

#ifdef XY_WALL
  void set_cfg(const libconfig::Setting &cfg);
#endif
  long double diss_energy_total = 0.0;

private:
#ifdef XY_WALL
  // static inline double UW(double x);
  inline double LW(double x) const;
  static inline bool WALLS(double x, double y);
#endif
#ifdef BOUNDS
  int outofbounds(const Row<double, 1> &position) const;
#endif
  inline int box(const double *position);
  inline Matrix<double, 2> generate_collision_matrix();
  Matrix<double, 2> virtimpuls;

  std::vector<std::vector<std::vector<double *>>> boxes;

  double a[3];
  const int numberofboxes;
  const int size[3];
  const string filename = "mpcd.dat";
  double *l_energy;

// OMP Variables
#ifdef PARALLEL
  omp_lock_t *locks;
#endif

  static constexpr double anull = 1.0;
  static constexpr int temp = 1;
  static constexpr int rho = 10;
  static constexpr double vwx = 0;
  static constexpr double vwy = 0;
  static constexpr double degree = 130;
  static constexpr double alpha = degree * PI / 180.0;
  static constexpr double salpha = 0.76604444;
  static constexpr double calpha = -0.6427876;
  static constexpr double talpha = -1.1917536;
  static constexpr int resolution = 400 * 10;
#ifdef XY_WALL_MAP
  LUT_function<0, resolution> lowerwall;
  LUT_function2d<0, resolution, 0, resolution> walls;
#endif
  //#ifdef XY_WALL
  double radius_of_curvature = 0.0;
  //#endif

  void virtpart(int i);
};
#endif
