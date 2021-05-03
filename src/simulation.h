#ifndef SIMULATION__H
#define SIMULATION__H
#include "fluid.h"
#include "globalforces.h"
#include "mpcd.h"
#include <fstream>
#include <libconfig.h++>
#include <vector>

using namespace libconfig;

class Simulation {
public:
  Simulation(const Setting &cfg);
  Simulation(string filename);
  Simulation(int sizex, int sizey, int sizez);
  ~Simulation();
  void run();
  void equlibrate(double time);
  MDObject *add_md_object(const Setting &cfg, MPCD &mpcd);
  void add_md_object(MDObject *obj);

private:
  MPCD mpcd;
  int number_of_objs;
  std::vector<MDObject *> md_objects;
  Fluid *fluid = NULL;
#ifdef A_FLOWFIELD
  fstream phases_out;
  Matrix<double, 4> flowF;
#endif

  GlobalForces globalF;

  void load_config(const Setting &root);
  void load_config(string filename);

  void move_particles();

  void next_tstep();
  void saveparameter();
  void loadparameter();
  void recalc_forces();

  double max_force;
};
#endif
