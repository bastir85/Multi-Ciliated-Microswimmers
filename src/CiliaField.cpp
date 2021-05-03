#define DEBUG
#include "CiliaField.h"
#include "basicf.h"
#include "global.h"
#include <assert.h>
#include <fstream>
using namespace Numeric_lib;
using namespace std;

#include "MatrixIO.h"

CiliaField::CiliaField(const Setting &cfg, MPCD &mpcd, Simulation &sim)
    : sim(sim), MDObject(cfg["name"], 0, 0, 0.0, mpcd, gf, 0) {
  double radius;
  Matrix<double, 1> pos(3);
  Matrix<double, 1> posT(3);
  cfg.lookupValue("radius", radius);
  string filename = name + ".dat";
  fstream in(filename, ios::in);
  loaded = in.good();
  try {
    pos[0] = cfg["position"][0];
    pos[1] = cfg["position"][1];
    pos[2] = cfg["position"][2];

    cfg.lookupValue("springconstant", springconstant);
    cfg["springconstant"] = 20000.0;

    int Nxl = cfg["num_equators"].getLength();
    int Nyl = cfg["y_thetas"].getLength();
    int num_offset = cfg["theta_offset"].getLength();
    if ((Nxl != Nyl) or (Nxl != num_offset)) {
      cerr << "Num of offsets mismatch" << endl;
      throw "config error";
    }
    cout << Nxl << ":" << Nyl << endl;
    int cil_idx = 0;
    for (int i = 0; i < Nyl; i++) {
      double thetaY = cfg["y_thetas"][i];
      thetaY *= PI / 180.;
      double theta_offset = cfg["theta_offset"][i];
      int Nx = cfg["num_equators"][i];
      for (int j = 0; j < Nx; j++) {
        double thetaX = theta_offset + 2 * PI * j / Nx;
        cout << "tx: " << thetaX << " ty: " << thetaY << endl;
        string cilia_name = name + "_cilia" + to_string(cil_idx);
        Cilia *cilia =
            create_cilia(cilia_name, pos, thetaX, thetaY, 0.0, cfg["cilia"]);
        cil_idx++;
        bool test = thetaX < 2 * PI + theta_offset;
        cout << i << "-" << thetaX << " " << test;
        cout << endl << endl;
      }
    }
  } catch (const SettingTypeException &exp) {
    cerr << "error with:" << exp.getPath() << endl;
    throw "config error";
  }

  num_particles = 0;
  for (unsigned int i = 0; i < cilias.size(); ++i) {
    num_particles += cilias[i]->positions.size() / 3;
  }
  if (loaded) {
    cout << "----------------" << loaded;
    load(in);
  }
  in.close();
}

Cilia *CiliaField::create_cilia(string name, Matrix<double, 1> &pos,
                                double thetaX, double thetaY, double thetaZ,
                                const Setting &cfg) {
  Matrix<double, 1> pos_cilia(3);
  int wo = 3; // add third layer to the surface of the sphere.

  double radius = wo * bondlength;
  pos_cilia[0] = radius * sin(thetaY);
  pos_cilia[1] = -radius * cos(thetaY) * sin(thetaX);
  pos_cilia[2] = radius * cos(thetaX) * cos(thetaY);
  pos_cilia += pos;

  cfg["name"] = name;
  cfg["position"][0] = pos_cilia[0];
  cfg["position"][1] = pos_cilia[1];
  cfg["position"][2] = pos_cilia[2];
  cfg["thetaX"] = thetaX;
  cfg["thetaY"] = thetaY;
  cfg["thetaZ"] = thetaZ;

  MDObject *obj = sim.add_md_object(cfg, mpcd);

  Cilia *cilia = dynamic_cast<Cilia *>(obj);
  if (not loaded) {
    cout << "create attachment springs" << endl;
    init_attachment_springs(cilia, wo * bondlength);
  }

  cilias.push_back(cilia);
  return cilia;
}

void CiliaField::recalcforce(void) {}

void CiliaField::global_forces(void) {
  const double l0 = bondlength * sin(PI / 3.0) * 2. / 3.;

  for (unsigned int i = 0; i < cilias.size(); ++i) {
    Cilia *cilia = cilias[i];
    for (unsigned int j = 0; j < cilias_attachment_pnts[i].size(); ++j) {
      SpringCfg &spring_cfg = cilias_attachment_pnts[i][j];
      Vector3d force = spring(spring_cfg.pos0, spring_cfg.pos1, spring_cfg.l0,
                              springconstant);
      spring_cfg.force0 += force;
      spring_cfg.force1 -= force;
    }
  }
} // end of recalcf

Matrix<double, 2> CiliaField::getBoundingBox() {
  Matrix<double, 2> minmax(2, 3);
  // negative so we don't have an overlap
  minmax[0] = 1.0;  // min
  minmax[1] = -1.0; // max
  return minmax;
}

void CiliaField::writedata() {
  string datname = name + ".xyz";
  ofstream data(datname.c_str(), ios::out | ios::app);
  while (!data.good()) {
    data.clear();
    data.open(datname.c_str(), ios::out | ios::app);
    cerr << "writing to datname.c_str(),ios::out|ios::app faild" << endl;
    if (data.good())
      cerr << " writing works" << endl;
  }
  // conversion to xyz data...
  data << num_particles << endl;

  data << "rods: " << wc.globalt << endl;

  for (unsigned int i = 0; i < cilias.size(); ++i) {
    Cilia *cilia = cilias[i];
    for (int j = 0; j < cilia->positions.dim2(); ++j) {
      data << "H ";
      data << cilia->positions[0][j][0] << " " << cilia->positions[0][j][1]
           << " " << cilia->positions[0][j][2] << endl;
    }

    for (int i = 1; i < cilia->positions.dim1(); ++i) {
      for (int j = 0; j < cilia->positions.dim2(); ++j) {
        data << "O ";
        data << cilia->positions[i][j][0] << " " << cilia->positions[i][j][1]
             << " " << cilia->positions[i][j][2] << endl;
      }
    }
  }

  data.close();
}

CiliaField::~CiliaField(){};
