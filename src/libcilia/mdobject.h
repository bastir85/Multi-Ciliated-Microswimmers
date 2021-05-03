#ifndef MDOBJECT__H
#define MDOBJECT__H
#include "Vector.h"
#include "mpcd.h"
#include <libconfig.h++>
#include <string.h>
using namespace libconfig;
using namespace Numeric_lib;
using namespace std;

class GlobalForces;

class MDObject {
  friend class GlobalForces;

public:
  MDObject(string name, int num_rods, int rodlength, double beatstrength,
           MPCD &mpcd, GlobalForces *gf, int startidx = 0,
           bool write_data = true)
      : name(name), num_rods(num_rods), rodlength(rodlength),
        startidx(startidx), mpcd(mpcd), positions(num_rods, rodlength, 3),
        velocities(num_rods, rodlength, 4), forces(num_rods, rodlength, 3),
        oldforces(num_rods, rodlength, 3), beatstrength(beatstrength), gf(gf),
        write_data(write_data) {
    forces = 0.0;
    oldforces = 0.0;
  };

  virtual ~MDObject(){};
  /* IO Functions */
  virtual void load(istream &in);
  virtual void save(fstream &out);
  void save();
  bool is_loaded() { return loaded; };
  virtual void writedata();

  /* COLLIDE AND MOVE */
  virtual void md_loop_start(void);
  virtual void md_loop(void);
  virtual void md_loop_stop(void);

  virtual void recalcforce() = 0;
  virtual void global_forces(){};

  virtual void debug(){};
  virtual void switch_logic(){};

  double get_beatstrength() { return beatstrength; };
  void set_beatstrength(double strength) { beatstrength = strength; };

  virtual double getPhase() { return 0.0; };
  virtual Matrix<double, 2> getBoundingBox();
  virtual double getT() { return 0.0; };

  virtual void calc_total_force();
  virtual void nextT() { forces = 0.0; };
  virtual void init();
  void calculate_flow_field();

protected:
  string name;
  const int num_rods;
  const int rodlength;
  const int startidx;
  MPCD &mpcd;

  Matrix<double, 3> positions;
  Matrix<double, 3> velocities;
  Matrix<double, 3> forces;
  Matrix<double, 3> oldforces;

  double beatstrength;
  GlobalForces *gf;

  bool loaded = false;
  static constexpr int partmass = 5;

  void spring(int i, int j, double b, double k);
  Vector3d spring(const Matrix<double, 1> &pos1, const Matrix<double, 1> &pos0,
                  double b, double k);
  Matrix<double, 2> get_avarage_rod();
  bool open_datfile();
  double get_num_particles() { return rodlength * num_rods; }

  double mu0 = 1.0; //(6*PI*eta*0.5);
  const bool write_data;
};

#endif
