#ifndef CILIA__H
#define CILIA__H
#include "mdobject.h"

class Cilia : public MDObject {
  friend class GlobalForces;
  friend class Clami;
  friend class HexagonGrid;

public:
  Cilia(const Setting &cfg, MPCD &mpcd, GlobalForces *gf)
      : Cilia(cfg["rodlength"], cfg["beatstrength"], cfg["bending"],
              cfg["min_bend"], cfg["max_bend"], cfg["m_speed"],
              cfg["springconstant"], cfg["position"][0], cfg["position"][1],
              cfg["position"][2], cfg["thetaX"], cfg["thetaY"], cfg["thetaZ"],
              mpcd, gf, cfg["name"], 0, 2, cfg["write_data"]){};
  Cilia(const Setting &cfg, Matrix<double, 1> pos, double thetaX, double thetaY,
        double thetaZ, MPCD &mpcd, GlobalForces *gf, string name, int startidx,
        int num_anker)
      : Cilia(cfg["rodlength"], cfg["beatstrength"],

              cfg["bending"], cfg["min_bend"], cfg["max_bend"], cfg["m_speed"],
              cfg["springconstant"], pos[0], pos[1], pos[2], thetaX, thetaY,
              thetaZ, mpcd, gf, name, startidx, num_anker){};
  Cilia(int rodlength, double beatstrength, double bending, double min_curv,
        double max_curv, double m_speed, double springconstant, double wx,
        double wy, double wz, double thetaX, double thetaY, double thetaZ,
        MPCD &mpcd, GlobalForces *gf, string name, int startidx = 2,
        int num_anker = 2, bool write_data = true);
  ~Cilia();

  /* IO Functions */
  void load(istream &is);
  void save(fstream &out);

  void recalcforce();
  void switch_logic();

  void debug();

  void create(double wx, double wy, double wz, double thetaX, double thetaY,
              double thetaZ);

  double getPhase();
  Matrix<double, 2> getBoundingBox();
  double getT() { return 0; }

protected:
  double localbend(int a, int b);
  double length(int i);

  void addcross();

  Matrix<double, 1> newlength;
  Matrix<double, 1> oldlength;

  double midpnt;
  const double start;
  double bending = 0.0150;

  double min_bend = -1.1;
  double max_bend = 1.1;
  double m_speed = 10.0;

  double Tstroke = 0.0;

private:
  double deltae_md_mpcd;

#ifdef MEM_AVA
  double ***vel_ava;
#endif

  // local controll variable
  static constexpr double bondlength = 0.5;
  double crosslength = sqrt(2.0) * bondlength;
  double springconstant = 20000.0;
  int num_anker = 1;

  static constexpr double recover = 0.0;
  static constexpr double crossstrength = 1.0;
};
#endif
