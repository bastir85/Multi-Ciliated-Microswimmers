#ifndef FLUID__H
#define FLUID__H
#include "Matrix.h"
#include "mpcd.h"
#include <atomic>
#include <cassert>
#include <string.h>

using namespace Numeric_lib;
using namespace std;

extern std::atomic<int> KeyCount;
class Key;
class Iterator;

class Fluid {
public:
  Fluid(int numberofparticles, MPCD &mpcd, string filename);
  static Fluid *FluidFromFile(string filename, MPCD &mpcd);
  void move();
  ~Fluid();
  void create();

  /* IO Functions*/
  void load(fstream &in);

  void writedata(void);
  void save();
  void check();

private:
  MPCD &mpcd;
  Matrix<double, 2> positions;
  Matrix<double, 2> velocities;
  int numberofparticles;
  string filename;

  int *cells = NULL;
  const double partmass = 1.0;
};
#endif
