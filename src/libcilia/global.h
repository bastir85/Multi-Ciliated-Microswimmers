#ifndef GLOBAL__H
#define GLOBAL__H

#include "abbrev.cpp"
#include <libconfig.h++>
using namespace libconfig;
class World {
public:
  World();

  static const double gravity;
  static int sizex, sizey, sizez;
  double tstep;
  static double ststep;
  static const double anull;

  double globalt;
  double desiredt;
  double runtime;
  double omega;
  static const int measureinterval;

  const string parfile = "parameter";
  void loadparameter();
  void loadparameter(const Setting &root);
  void saveparameter();
  void saveparameter(const Setting &root);
};

extern World wc;

#define WRITE2(x) out << "#" << #x << " " << x << endl;
#define REED2(x, in)                                                           \
  in >> garb;                                                                  \
  if (strcmp(&garb[1], #x)) {                                                  \
    cout << "error with " << garb << " as " << #x;                             \
    exit(0);                                                                   \
  }                                                                            \
  in >> x;

#endif
