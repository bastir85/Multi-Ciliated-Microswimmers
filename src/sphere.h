#ifndef SPHERE__H
#define SPHERE__H
#include "mdobject.h"

class Sphere : public MDObject {
  friend class GlobalForces;
  friend class TrajectorySphere;
  friend class Sperm;
  friend class Clami;
  friend class NewSperm;
  friend class Akashiwo;

public:
  Sphere(const Setting &cfg, MPCD &mpcd, GlobalForces *gf)
      : Sphere(cfg["radius"], cfg["springconstant"], cfg["position"][0],
               cfg["position"][1], cfg["position"][2], cfg["npart"], mpcd, gf,
               cfg["name"]){};
  Sphere(double radius, double springconstant, double wx, double wy, double wz,
         int npart, MPCD &mpcd, GlobalForces *gf, string name);
  ~Sphere();

  /* IO Functions */
  void load(istream &is);
  void save(fstream &out);
  // void writedata();

  void recalcforce();

  void debug();

  void create(double wx, double wy, double wz, int fill_factor);
  void create(Matrix<double, 1> position, int fill_factor);
  Matrix<double, 2> getBoundingBox();
  double getT();
  double getPhase();

private:
// Functions
#ifdef SIMPLE_FLOWF
  Vector_lib::Vector3d flowField(const Matrix<double, 1> &pos);
  void update_velocities();
#endif
  void calc_total_force();
  // Elliptic traj

  double springconstant = 20000.0;
  double radius; //=2
  double getForce(double time);
  Matrix<double, 1> tforce;
  Matrix<double, 1> ava_vel;

  Matrix<double, 2> max_flowF;

  int npart;
  int **neilist;
  double headbond; // needed (set in create but where else?)
  void create_neillist(double max_part);
};
#endif
