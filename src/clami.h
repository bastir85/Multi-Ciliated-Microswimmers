#ifndef CLAMI__H
#define CLAMI__H
#include "cilia.h"
#include "mdobject.h"
#include "simulation.h"
#include "sphere.h"
#include "utils.h"

class Clami : public MDObject {
public:
  Clami(const Setting &cfg, MPCD &mpcd, Simulation &sim, GlobalForces *gf);
  ~Clami();
  void writedata();

  void recalcforce();
  void global_forces();
  Sphere *sphere;
  // Cilia *cilia2;
private:
  Simulation &sim;

  double springconstant = 20000.0;
  double bondlength = 0.5;
  int num_particles;

  std::vector<int> cilias_attach_idx;
  std::vector<std::vector<SpringCfg>> cilias_attachment_pnts;
  std::vector<Cilia *> cilias;

  double getPhase() { return 0.0; };
  double getT() { return 0.0; };
  Matrix<double, 2> getBoundingBox();
  void create_cilia(string name, Matrix<double, 1> &pos, double thetaX,
                    double thetaY, double thetaZ, double start_delay,
                    const Setting &cfg);
  void init_attachment_springs(Cilia *cilia, double depth);
  void add_springs_for_rods_idx(std::vector<SpringCfg> &springs, Cilia *cilia,
                                const int wo);

  void save(fstream &out);
  void load(istream &in);
};
#endif
