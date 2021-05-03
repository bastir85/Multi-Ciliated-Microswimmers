#define DEBUG
#include "clami.h"
#include "global.h"
#include <assert.h>
#include <stdexcept>

#include "basicf.h"
#include <fstream>
using namespace Numeric_lib;
using namespace std;

#include "MatrixIO.h"

class SphereNoWrite : public Sphere {
public:
  SphereNoWrite(double radius, double springconstant, double wx, double wy,
                double wz, MPCD &mpcd, GlobalForces *gf, string name)
      : Sphere(radius, springconstant, wx, wy, wz, 643, mpcd, gf, name){};
  void writedata(){};
};

class CiliaNoWrite : public Cilia {
public:
  CiliaNoWrite(const Setting &cfg, Matrix<double, 1> pos, double thetaX,
               double thetaY, double thetaZ, MPCD &mpcd, GlobalForces *gf,
               string name, int num_anker)
      : Cilia(cfg, pos, thetaX, thetaY, thetaZ, mpcd, gf, name, 0,
              num_anker){}; // was 16
  void writedata(){};
};

Clami::Clami(const Setting &cfg, MPCD &mpcd, Simulation &sim, GlobalForces *gf)
    : MDObject(cfg["name"], 0, 0, 0.0, mpcd, gf, 0), sim(sim) {
  double radius;
  Matrix<double, 1> pos(3);
  Matrix<double, 1> posT(3);
  cfg.lookupValue("radius", radius);
  if (radius <= 0.0) {
    throw std::invalid_argument("received negative value");
  }
  string filename = name + ".dat";
  fstream in(filename, ios::in);
  loaded = in.good();
  try {
    pos[0] = cfg["position"][0];
    pos[1] = cfg["position"][1];
    pos[2] = cfg["position"][2];

    cfg.lookupValue("springconstant", springconstant);
    cfg["springconstant"] = 20000.0;

    sphere = new SphereNoWrite(radius, springconstant, pos[0], pos[1], pos[2],
                               mpcd, gf, name + "_sphere");
    sim.add_md_object(sphere);

    int Nxl = cfg["num_equators"].getLength();
    int Nyl = cfg["y_thetas"].getLength();
    int num_offset = cfg["theta_offset"].getLength();
    if ((Nxl != Nyl) or (Nxl != num_offset)) {
      cerr << "Num of offsets mismatch" << endl;
      throw "config error";
    }
    int num_cilia = 0;
    for (int i = 0; i < Nyl; i++) {
      int num_cilia_eq = cfg["num_equators"][i];
      num_cilia += num_cilia_eq;
    }

    if (cfg.exists("start_delays")) {
      if (num_cilia != cfg["start_delays"].getLength()) {
        cerr << "Phase lacks" << num_offset << "mismatch num of cilia "
             << num_cilia << endl;
        throw "config error";
      }
    } else {
      cout << "!!! Using default phase lag of ZERO for all cilia " << endl;
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
        double start_delay = 0.0;
        if (cfg.exists("start_delays")) {
          start_delay = cfg["start_delays"][cil_idx];
        }
        cout << "tx: " << thetaX << " ty: " << thetaY
             << " start_delay:" << start_delay << endl;
        string cilia_name = name + "_cilia" + to_string(cil_idx);
        create_cilia(cilia_name, pos, thetaX, thetaY, cfg["cilia"]["thetaZ"],
                     start_delay, cfg["cilia"]);
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

  num_particles = sphere->positions.size() / 3;
  for (unsigned int i = 0; i < cilias.size(); ++i) {
    num_particles += cilias[i]->positions.size() / 3;
  }
  cout << "sP @" << sphere->positions[0][0].data() << endl;
  if (loaded) {
    cout << "-----------------" << loaded;
    load(in);
  }
  in.close();
}

void Clami::create_cilia(string name, Matrix<double, 1> &pos, double thetaX,
                         double thetaY, double thetaZ, double start_delay,
                         const Setting &cfg) {
  Matrix<double, 1> pos_cilia(3);
  int wo = 3; // add third layer to the surface of the sphere.

  double radius = sphere->radius - wo * bondlength;
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
  cfg["start_delay"] = start_delay;

  MDObject *obj = sim.add_md_object(cfg, mpcd);

  Cilia *cilia = dynamic_cast<Cilia *>(obj);
  if (not loaded) {
    cout << "create attachment springs" << endl;
    init_attachment_springs(cilia, wo * bondlength);
  }

  cilias.push_back(cilia);
}

void Clami::add_springs_for_rods_idx(std::vector<SpringCfg> &springs,
                                     Cilia *cilia, const int wo) {
  for (int j = 0; j < cilia->positions.dim1(); ++j) {
    // Attach  rod to the closest surface point!
    int idx = sphere->positions[0].get_closest_idx(cilia->positions[j][wo]);
    Row<double, 1> posS = sphere->positions[0][idx];
    double l1 = sqrt(posS.sqdistance(cilia->positions[j][wo]));
    // cout << "1-surface: " << wo <<"<-->" << idx << "(" << l1 << ")" <<endl;
    springs.emplace_back(posS, cilia->positions[j][wo], sphere->forces[0][idx],
                         cilia->forces[j][wo], l1);

    // Attach  rod to all its neighbors!
    for (int n = 1; n < sphere->neilist[idx][0]; ++n) {
      int attachSphereIdx = sphere->neilist[idx][n];
      Row<double, 1> posSphere = sphere->positions[0][attachSphereIdx];
      double sqdis = posSphere.sqdistance(cilia->positions[j][wo]);
      l1 = sqrt(sqdis);
      springs.emplace_back(
          sphere->positions[0][attachSphereIdx], cilia->positions[j][wo],
          sphere->forces[0][attachSphereIdx], cilia->forces[j][wo], l1);
    }
  }
}

void Clami::init_attachment_springs(Cilia *cilia, double depth) {
  std::vector<SpringCfg> springs;

  int wo = int(depth / bondlength);
  add_springs_for_rods_idx(springs, cilia, 0);
  add_springs_for_rods_idx(springs, cilia, wo);

  cilias_attachment_pnts.push_back(std::move(springs));
}

void Clami::recalcforce(void) {}

void Clami::global_forces(void) {
  for (unsigned int i = 0; i < cilias.size(); ++i) {
    for (unsigned int j = 0; j < cilias_attachment_pnts[i].size(); ++j) {
      SpringCfg &spring_cfg = cilias_attachment_pnts[i][j];
      Vector3d force = spring(spring_cfg.pos0, spring_cfg.pos1, spring_cfg.l0,
                              springconstant);
      spring_cfg.force0 += force;
      spring_cfg.force1 -= force;
    }
  }
} // end of recalcf

Matrix<double, 2> Clami::getBoundingBox() {
  Matrix<double, 2> minmax(2, 3);
  // negative so we don't have an overlap
  minmax[0] = 1.0;  // min
  minmax[1] = -1.0; // max
  return minmax;
}

void Clami::writedata() {
  string datname = name + ".xyz";
  ofstream data(datname.c_str(), ios::out | ios::app);
  while (!data.good()) {
    data.clear();
    data.open(datname.c_str(), ios::out | ios::app);
    cerr << "writing to datname.c_str(),ios::out|ios::app faild" << endl;
    if (data.good())
      cerr << " writing works" << endl;
  }

  Vector3d xcm, vcm, n, b, p;
  xcm = vcm = n = b = p = 0.0;
  const int even =
      cilias[0]->positions.dim2() % 2; // assume all cilia same size here!

  // conversion to xyz data...
  data << cilias.size() * cilias[0]->positions.dim1() *
              (cilias[0]->positions.dim2() / 2 + even)
       << endl;

  // data << (sphere->positions.size()+cilia1->positions.size())/3<<endl;
  data << "rods: " << wc.globalt << endl;
  // Sphere
  for (int i = 0; i < sphere->positions.dim1(); ++i) {
    for (int j = 0; j < sphere->positions.dim2(); ++j) {
      // data << "C ";
      // data << sphere->positions[i][j][0] << " "  <<
      // sphere->positions[i][j][1] << " " << sphere->positions[i][j][2] << endl;
      xcm += sphere->positions[i][j];
      vcm += sphere->velocities[i][j];
    }
  }
  for (unsigned int i = 0; i < cilias.size(); ++i) {
    Cilia *cilia = cilias[i];
    for (int j = 0; j < cilia->positions.dim2(); ++j) {
      if ((j + 1) % 2 == even) {
        data << "H ";
        data << cilia->positions[0][j][0] << " " << cilia->positions[0][j][1]
             << " " << cilia->positions[0][j][2] << endl;
      }
      xcm += cilia->positions[0][j];
      vcm += cilia->velocities[0][j];
    }

    for (int i = 1; i < cilia->positions.dim1(); ++i) {
      for (int j = 0; j < cilia->positions.dim2(); ++j) {
        if ((j + 1) % 2 == even) {
          data << "O ";
          data << cilia->positions[i][j][0] << " " << cilia->positions[i][j][1]
               << " " << cilia->positions[i][j][2] << endl;
        }
        xcm += cilia->positions[i][j];
        vcm += cilia->velocities[i][j];
      }
    }
  }

  data.close();

  // Writeout calcualted values
  xcm /= num_particles;
  vcm /= num_particles;
  n = sphere->positions[0][1] - sphere->positions[0][0];
  b = sphere->positions[0][32] - sphere->positions[0][0];
  p = sphere->positions[0][137] - sphere->positions[0][0];

  datname = name + ".san";
  data.open(datname.c_str(), ios::out | ios::app);
  data << wc.globalt << " " << xcm[0] << " " << xcm[1] << " " << xcm[2] << " "
       << vcm[0] << " " << vcm[1] << " " << vcm[2] << " " << n[0] << " " << n[1]
       << " " << n[2] << " " << b[0] << " " << b[1] << " " << b[2] << " "
       << p[0] << " " << p[1] << " " << p[2] << " "
       << mpcd.diss_energy_total / wc.measureinterval / wc.tstep << endl;
  mpcd.diss_energy_total = 0.0;
  data.close();
}

Clami::~Clami(){};

void Clami::save(fstream &out) {
  cout << "saving springs only" << endl;
  // MDObject::save(out);

  out << "#springs" << endl;
  out << cilias.size() << endl;
  for (unsigned int i = 0; i < cilias.size(); ++i) {
    for (unsigned int j = 0; j < cilias_attachment_pnts[i].size(); ++j) {
      SpringCfg &spring_cfg = cilias_attachment_pnts[i][j];
      // cil_idx pos_cilia_idx pos_sphere_idx l0
      out << i << " " << get_index(spring_cfg.pos0, sphere->positions) << " "
          << get_index(spring_cfg.pos1, cilias[i]->positions) << " "
          << spring_cfg.l0 << endl;
    }
  }
}

void Clami::load(istream &in) {
  cout << "load springs" << endl;
  string garb;
  // MDObject::load(in);
  unsigned int num_cilia;
  unsigned int cil_idx;
  in >> garb;
  if (garb != "#springs") {
    cout << "input error springs" << endl;
    exit(1);
  }
  in >> num_cilia;
  if (num_cilia != cilias.size()) {
    cout << "Mismatch between created cilia and dat file." << endl;
    exit(1);
  }
  in >> cil_idx;
  for (unsigned int i = 0; i < cilias.size(); ++i) {
    std::vector<SpringCfg> springs;
    Cilia *cilia = cilias[cil_idx];
    if (cil_idx != i) {
      cout << "Mismatch between created cilia and dat file index." << endl;
      exit(1);
    }
    while (cil_idx == i and not in.eof()) {
      int sphere_pidx, cil_pidx;
      double l0;
      in >> sphere_pidx;
      in >> cil_pidx;
      in >> l0;
      cout << "addS: " << sphere_pidx << "<->" << cil_pidx << "(" << l0 << ")"
           << endl;
      springs.emplace_back(sphere->positions.get_row_by_idx(sphere_pidx),
                           cilia->positions.get_row_by_idx(cil_pidx),
                           sphere->forces.get_row_by_idx(sphere_pidx),
                           cilia->forces.get_row_by_idx(cil_pidx), l0);
      in >> cil_idx;
    }
    cilias_attachment_pnts.push_back(std::move(springs));
  }
}
