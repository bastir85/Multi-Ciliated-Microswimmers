#include "global.h"

#ifdef FLUID_FLOW_G
const double World::gravity = FLUID_FLOW_G;
#else
const double World::gravity = 0.0;
#endif
int World::sizex = 100;
int World::sizey = 100;
int World::sizez = 50;
double World::ststep = 0.001;
const double World::anull = 1.0;

const int World::measureinterval = SIM_MEASUREINTERVAL;

World wc;
//= {0.0,10,10,10,0.05,1.0};
//
World::World() {}

void World::loadparameter(const Setting &root) {
  wc.globalt = root["globalt"];
  wc.desiredt = root["desiredt"];
  wc.omega = root["omega"];
  wc.sizex = root["sizex"];
  wc.sizey = root["sizey"];
  wc.sizez = root["sizez"];
  wc.runtime = root["runtime"];
  wc.tstep = root["tstep"];
};

void World::loadparameter() {
  // lodes Data
  ifstream in(parfile.c_str());
  string garb; // neccesarry for REED
  double globalt, omega, sizex, sizey, sizez;
  REED(globalt);
  REED(desiredt);
  REED(runtime);
  REED(sizex);
  REED(sizey);
  REED(sizez);
  REED(omega);
  in.close();
  wc.globalt = globalt;
  wc.omega = omega;
  wc.sizex = sizex;
  wc.sizey = sizey;
  wc.sizez = sizez;
};

void World::saveparameter() {
  fstream out(parfile.c_str(), ios::out);
  double globalt = wc.globalt;
  WRITE(globalt);
  WRITE(desiredt);
  WRITE(runtime);
  double omega = wc.omega;
  double sizex = wc.sizex;
  double sizey = wc.sizey;
  double sizez = wc.sizez;
  WRITE(sizex);
  WRITE(sizey);
  WRITE(sizez);
  WRITE(omega)
  out.close();
}

void World::saveparameter(const Setting &root) {
  root["globalt"] = wc.globalt;
  root["desiredt"] = wc.desiredt;
  root["runtime"] = wc.runtime;
  root["omega"] = wc.omega;
  root["sizex"] = wc.sizex;
  root["sizey"] = wc.sizey;
  root["sizez"] = wc.sizez;
}
