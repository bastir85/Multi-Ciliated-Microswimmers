#include "fluid.h"
#include "basicf.h"
#include "global.h"

#include "MatrixIO.h"

#include <atomic>
#include <iterator>

using namespace std;

Fluid::Fluid(int numberofparticles, MPCD &mpcd, string filename)
    : mpcd(mpcd), positions(numberofparticles, 3),
      velocities(numberofparticles, 4), numberofparticles(numberofparticles),
      filename(filename)

{
  cells = new int[numberofparticles];
}

Fluid *Fluid::FluidFromFile(string filename, MPCD &mpcd) {
  fstream in(filename, ios::in);
  Fluid *fl;
  int num_parts = wc.sizex * wc.sizey * wc.sizez * mpcd.rho;
  if (not in.good()) {
    cout << "create fluid rh=10 " << wc.sizex << "x" << wc.sizey << "x"
         << wc.sizez << endl;
    cout << "#part: " << num_parts;
    fl = new Fluid(num_parts, mpcd, filename);
    fl->create();
  } else {
    string garb;
    int numberofparticles;
    REED2(numberofparticles, in)
    if (num_parts != numberofparticles) {
      cerr << "!!! Loaded fluid does not match box size !!!" << wc.sizex << "x"
           << wc.sizey << "x" << wc.sizez << endl;
      if (numberofparticles > num_parts) {
        throw "Fluidfile to big.";
      }
    }
    cout << "nFL:" << numberofparticles << endl;
    fl = new Fluid(numberofparticles, mpcd, filename);
    fl->load(in);
  }

  return fl;
}

Fluid::~Fluid() {
  cout << "deleting fluid";
  delete[] cells;
}

void Fluid::create() {
  positions = 0;
  velocities = 0;
  cout << "NpS:" << numberofparticles << endl;
  for (int i = 0; i < numberofparticles; ++i) {
    positions[i][0] = rnd() * wc.sizex * wc.anull;
    positions[i][1] = rnd() * wc.sizey * wc.anull;
    positions[i][2] = rnd() * wc.sizez * wc.anull;
    velocities[i][0] = normalrand() * sqrt(mpcd.temp);
    velocities[i][1] = normalrand() * sqrt(mpcd.temp);
    velocities[i][2] = normalrand() * sqrt(mpcd.temp);
    velocities[i][3] = partmass;
#ifdef BOUNDS
    if (mpcd.outofbounds(positions[i])) {
      i = i - 1;
      numberofparticles -= 1;
    }
#endif
  }
  cout << "Npf:" << numberofparticles << endl;

  Vector3d cm;

  cm = 0.0;
  for (int i = 0; i < numberofparticles; ++i) {
    cm += velocities[i];
  }
  cm /= numberofparticles;
  cout << cm << endl;
  for (int i = 0; i < numberofparticles; ++i) {
    velocities[i][0] -= cm[0];
    velocities[i][1] -= cm[1];
    velocities[i][2] -= cm[2];
  }
  cout << "cm_fluid " << cm;
  cm = 0.0;
  for (int i = 0; i < numberofparticles; ++i) {
    cm += velocities[i];
  }
  cout << " -> " << cm << endl;
  for (int i = 0; i < numberofparticles; ++i) {
    mpcd.add_impuls(positions[i].data(), velocities[i].data());
  }
}

void Fluid::save() {

  fstream out(filename, ios::out);
  out << std::scientific;
  WRITE2(numberofparticles);
  out << positions.slice(0, numberofparticles);
  out << "#velocities" << endl;
  out << velocities.slice(0, numberofparticles);
}

void Fluid::check() {
  for (int i = 0; i < numberofparticles; i++) {
#ifdef BOUNDS
    if (mpcd.outofbounds(positions[i])) {
      cout << i << " is out of bounds!" << endl;
    }
#endif
    cout << velocities[i];
  }
}

void Fluid::load(fstream &in) {
  string garb;
  in >> positions;
  in >> garb;
  in >> velocities;
#ifdef BOUNDS
  for (int i = 0; i < numberofparticles; i++) {
    while (mpcd.outofbounds(positions[i])) {
      cerr << "fluid particle " << i << " is out of bounds, trying to move..."
           << endl;
      cerr << positions[i];
      double dx = 0.1 * (rnd() - 0.5);
      double dy = 0.1 * (rnd() - 0.5);
      double dz = 0.1 * (rnd() - 0.5);
      positions[i][0] += dx;
      positions[i][1] += dy;
      positions[i][2] += dz;
      if (mpcd.outofbounds(positions[i])) {
        positions[i][0] -= dx;
        positions[i][1] -= dy;
        positions[i][2] -= dz;
      }
    }
  }
#endif

  for (int i = 0; i < numberofparticles; ++i) {
    mpcd.add_impuls(positions[i].data(), velocities[i].data());
  }
}

void Fluid::move() {
#pragma omp for schedule(static) nowait
  for (int i = 0; i < numberofparticles; ++i) {
    mpcd.move(positions[i].data(), velocities[i].data(), wc.tstep);
#ifdef FLUID_FLOW_G
    velocities[i][0] += wc.gravity * cos(wc.globalt * wc.omega) * wc.tstep;
#endif
#ifdef FLUID_FLOW_SEGMENT
    if (fmod(cells[i], wc.sizex) - wc.sizex / 2 == 0) {
      velocities[i][0] += wc.gravity * cos(wc.globalt * wc.omega) * wc.tstep;
      // cout << "applied fluid for part at" << positions[i][0] << endl;
    }
#endif
    mpcd.add_impuls(positions[i].data(), velocities[i].data());
  }
}

void Fluid::writedata() {
  /*    const string datname  ="fl_datas.xyz";
      ofstream data(datname.c_str(),ios::out|ios::app);
      while (!data.good()){
          data.clear();
          data.open(datname.c_str(),ios::out|ios::app);
          cerr << "writing to datname.c_str(),ios::out|ios::app faild" << endl;
          if (data.good())
              cerr << " writing works" << endl;
      }
      //iconversion to xyz data...
      for (int i=0;i< numberofparticles;++i){
        data << positions[i][0] << " "  << positions[i][1] << " " <<
     positions[i][2]
             << " " << cells[i] << endl;
      }
      data.close();
  */

  double ltemp = 0.0;
  for (int i = 0; i < numberofparticles; ++i) {
    for (int q = 0; q < 3; ++q) {
      ltemp += velocities[i][q] * velocities[i][q];
    }
  }
  ltemp = ltemp / numberofparticles / 3.0;

  string datname = "temp.san";
  ofstream data(datname.c_str(), ios::out | ios::app);
  while (!data.good()) {
    data.clear();
    data.open(datname.c_str(), ios::out | ios::app);
    cerr << "writing to datname.c_str(),ios::out|ios::app faild" << endl;
    if (data.good())
      cerr << " writing works" << endl;
  }
  data << ltemp << endl;
  data.close();
}
