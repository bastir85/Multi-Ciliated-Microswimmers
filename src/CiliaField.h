#ifndef CILIA_FIELD_H
#define CILIA_FIELD_H
#include "cilia.h"
#include "mdobject.h"
#include "simulation.h"
#include "utils.h"

class CiliaField : public MDObject {
public:
  CiliaField(const Setting &cfg, MPCD &mpcd, Simulation &sim);
  ~CiliaField();
  void writedata();

  void recalcforce();
  void global_forces();

private:
  double springconstant = 20000.0;
  double bondlength = 0.5;
  int num_particles;

  std::vector<Cilia *> cilias;

  double getPhase() { return 0.0; };
  double getT() { return 0.0; };
  Matrix<double, 2> getBoundingBox();
  Cilia *create_cilia(string name, Matrix<double, 1> &pos, double thetaX,
                      double thetaY, double thetaZ, const Setting &cfg);

  Simulation &sim;
};
#endif
