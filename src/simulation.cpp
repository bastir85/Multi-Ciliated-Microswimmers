#include "simulation.h"
#include "MatrixIO.h"
#include "basicf.h"
#include "global.h"
#include <assert.h>

#include "ConstantBeatCilia.h"
#include "cilia.h"
#include "clami.h"
#include "sphere.h"

#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>

using namespace std::chrono;

Simulation::Simulation(const Setting &cfg)
    : Simulation(cfg["sizex"], cfg["sizey"], cfg["sizez"]) {
  load_config(cfg);
#ifdef XY_WALL
  mpcd.set_cfg(cfg);
#endif
  globalF.set_md_objects(md_objects.data(), number_of_objs);

  fluid = Fluid::FluidFromFile("fluid_n.dat", mpcd);
  cout << "init md objects" << endl;
  for (int i = 0; i < number_of_objs; ++i) {
    md_objects[i]->init();
  }

  cout << "inital config" << endl;
  mpcd.debug();
};

Simulation::Simulation(int sizex, int sizey, int sizez)
    : mpcd(sizex, sizey, sizez), number_of_objs(0), md_objects(0), max_force(0)
#ifdef A_FLOWFIELD
      ,
      phases_out("phases.log", ios::out), flowF(100 + 1, sizex, sizey, 4)
#endif
                                              {};

Simulation::Simulation(string filename)
    : mpcd(wc.sizex, wc.sizey, wc.sizez), number_of_objs(0), md_objects(0),
      max_force(0)
#ifdef A_FLOWFIELD
      ,
      phases_out("phases.log", ios::out), flowF(100 + 1, wc.sizex, wc.sizey, 4)
#endif
{
  load_config(filename);
  globalF.set_md_objects(md_objects.data(), number_of_objs);

  fluid = Fluid::FluidFromFile("fluid_n.dat", mpcd);
}

void Simulation::run() {
  int started = time(0);
  int lastsave = time(0);
  double anfangsgt = wc.globalt;

  cout << "start simulation:" << endl;
#ifdef A_FLOWFIELD
  fstream out("flowF.dat", ios::out);
#endif
#ifdef SIM_LOG_FPS
  fstream fps_out("fps.san", ios::out);
#endif
#ifdef DEBUG2
  for (int i = 0; i < number_of_objs; ++i) {
    md_objects[i]->writedata();
  }
#endif
  for (; wc.globalt < wc.desiredt;) {
    milliseconds deltaT =
        duration_cast<milliseconds>(system_clock::now().time_since_epoch());

    for (int smalstep = 0; smalstep < wc.measureinterval; smalstep++) {
      next_tstep();

#if defined(A_FLOWFIELD)
      int phase_bin = floor(md_objects[0]->getPhase() / 2 / PI * 100);
      mpcd.getFlowField(flowF[phase_bin]);
      smalstep += vsmalstep - 1; // since we already have one at the outer loop
                                 // (works w/o fine greadeing)
#endif
    }
    max_force = 0.0;

#ifdef WRITE_FLUID
    mpcd.writedata();
    fluid->writedata();
#endif

#if defined(DEBUG)
    for (int i = 0; i < number_of_objs; ++i) {
      md_objects[i]->debug();
    }
    mpcd.debug();
    fluid->check();
#endif

    for (int i = 0; i < number_of_objs; ++i) {
      md_objects[i]->writedata();
    }

    milliseconds current =
        duration_cast<milliseconds>(system_clock::now().time_since_epoch());
    double delta = (current - deltaT).count();
    cout << "T =" << wc.globalt;
    if (delta > 0) {
      cout << " ("
           << 1000.0 * wc.tstep * wc.measureinterval /
                  (current - deltaT).count()
           << "fps)" << endl;
    }
#ifdef SIM_LOG_FPS
    fps_out << 1.0 * wc.measureinterval / (current - deltaT).count() << endl;
#endif

    if (((time(0) - lastsave) > 7200) || ((time(0) - wc.runtime) > started)) {
      if ((time(0) - wc.runtime) > started) {
        break;
      }
      cout << "save stats" << endl;
      wc.saveparameter();
      for (int i = 0; i < number_of_objs; ++i) {
        md_objects[i]->save();
      }
      fluid->save();
#ifdef A_FLOWFIELD
      string sname = "flowF.dat";
      fstream out(sname, ios::out);
      out.precision(6);
      out << "#" << 100 << " " << wc.sizex << " " << wc.sizey << endl;
      out << flowF;
      fflush(stdout);
      out.close();
#endif
      lastsave = time(0);
    }
  }
#ifdef A_FLOWFIELD
  out.close();
#endif

  cout << endl << endl;
  cout << "end simulation:" << endl;
  cout << wc.globalt - anfangsgt << " r:" << (time(0) - started) << endl;
  cout << "fps: " << (wc.globalt - anfangsgt) / (time(0) - started) << endl;
#ifdef SIM_LOG_FPS
  fps_out << (wc.globalt - anfangsgt) / wc.tstep / (time(0) - started) << endl;
#endif
}

void Simulation::equlibrate(double time) {
  vector<double> lagtimes(number_of_objs);
  vector<double> beatstrengthes(number_of_objs);
  if (time == 0.0)
    time = md_objects[0]->getT() * 1.5;
  for (int i = 1; i < number_of_objs; ++i) {
    lagtimes[i] = md_objects[i]->getT();
    beatstrengthes[i] = md_objects[i]->get_beatstrength();
    md_objects[i]->set_beatstrength(0.0);
    cout << "Lagtime:" << lagtimes[i] << endl;
  }

  for (; wc.globalt < time;) { // runs until desired time is reached
    for (int smalstep = 0; smalstep < wc.measureinterval; smalstep++) {
      for (int i = 0; i < number_of_objs; ++i) {
        if ((lagtimes[i] - wc.globalt) * (lagtimes[i] - wc.globalt) <
            wc.tstep * wc.tstep) {
          cout << "activate " << i << wc.globalt << " " << lagtimes[i] << endl;
          md_objects[i]->set_beatstrength(beatstrengthes[i]);
        }
      }
      next_tstep();
    }
    for (int i = 0; i < number_of_objs; ++i) {
      md_objects[i]->writedata();
    }
    cout << wc.globalt << endl;
  }
}

#ifdef OMP_PARALLEL_FORCES
inline void Simulation::recalc_forces() {
  for (int i = 0; i < number_of_objs; ++i) {
    md_objects[i]->recalcforce();
  }
  globalF.recalc_forces();

#pragma omp for schedule(guided)
  for (int i = 0; i < number_of_objs; ++i) {
    md_objects[i]->global_forces();
  }
}
#endif

#ifdef OMP_PARALLEL_MDOBJ
inline void Simulation::recalc_forces() {
#pragma omp for schedule(guided)
  for (int i = 0; i < number_of_objs; ++i) {
    md_objects[i]->recalcforce();
  }
  globalF.recalc_forces();

#pragma omp for schedule(guided)
  for (int i = 0; i < number_of_objs; ++i) {
    md_objects[i]->global_forces();
  }
}
#endif

void Simulation::next_tstep() {
  const static int n = int((wc.tstep + 0.00000001) / wc.ststep);

#pragma omp parallel
  {
    mpcd.next_tstep();
    recalc_forces();
#pragma omp barrier
    for (int i = 0; i < number_of_objs; ++i) {
      md_objects[i]->md_loop_start();
    }
#pragma omp barrier
    recalc_forces();
  }

  // inner MD Loop.
  for (int k = 1; k < n; k++) {
    wc.globalt += wc.ststep;

#pragma omp parallel
    {
      for (int i = 0; i < number_of_objs; ++i) {
        md_objects[i]->md_loop();
      }
// Wait for all velocity verlet -> update v,x
#pragma omp barrier
      recalc_forces();
    }
  }

  wc.globalt += wc.ststep;
  mpcd.shift_boxes();

#pragma omp parallel
  {
#if !defined(MPCD_RANDOM)
    fluid->move();
#endif

    for (int i = 0; i < number_of_objs; ++i) {
      md_objects[i]->md_loop_stop();
    }

#pragma omp for schedule(guided)
    for (int i = 0; i < number_of_objs; ++i) {
      md_objects[i]->switch_logic();
    }
  }
}

void Simulation::load_config(string filename) {
  Config cfg;
  cout << "loading: " << filename;
  try {
    cfg.readFile(filename.c_str());
  } catch (const FileIOException &fioex) {
    std::cerr << "I/O error while reading file." << std::endl;
    throw "Config Error";
  } catch (const ParseException &pex) {
    std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
              << " - " << pex.getError() << std::endl;
    throw "Config Error";
  }

  try {
    string name = cfg.lookup("name");
    cout << "Simulation " << name << " MD Objects:" << endl << endl;
  } catch (const SettingNotFoundException &nfex) {
    cerr << "No 'name' setting in configuration file." << endl;
    throw "Config Error";
  }

  const Setting &root = cfg.getRoot();
  load_config(root);
}

void Simulation::load_config(const Setting &root) {
  try {
    const Setting &c_md_objects = root["md_objects"];
    int count = c_md_objects.getLength();

    cout << setw(15) << left << "type"
         << "  " << setw(10) << left << "name"
         << "  " << setw(5) << left << "new"
         << "  " << endl;

    // md_objects = new MDObject*[count];
    for (int i = 0; i < count; ++i) {
      const Setting &c_obj = c_md_objects[i];

      string type, name;

      if (!(c_obj.lookupValue("type", type) &&
            c_obj.lookupValue("name", name))) {
        cout << name << "skipped" << endl;
        continue;
      }

      MDObject *obj = add_md_object(c_obj, mpcd);
      cout << setw(15) << left << type << "  " << setw(10) << left << name
           << "  " << setw(7) << left << obj->is_loaded() << endl
           << endl;
    }

    number_of_objs = md_objects.size();
    cout << "# md objects:" << number_of_objs << endl << endl;
  } catch (const SettingNotFoundException &nfex) {
    cerr << nfex.getPath() << endl;
    throw;
  }
}

void Simulation::add_md_object(MDObject *obj) { md_objects.push_back(obj); }

MDObject *Simulation::add_md_object(const Setting &cfg, MPCD &mpcd) {
  const string type = cfg["type"];
  MDObject *obj = NULL;
  cout << "add:" << type << endl;

  if (type == "Clami") {
    obj = new Clami(cfg, mpcd, *this, &globalF);
  } else if (type == "Sphere") {
    obj = new Sphere(cfg, mpcd, &globalF);
  } else if (type == "Cilia") {
    obj = new Cilia(cfg, mpcd, &globalF);
  } else if (type == "ConstantBeatCilia") {
    obj = new ConstantBeatCilia(cfg, mpcd, &globalF);
  } else {
    cerr << "Unkown type:" << type << endl;
    throw "Syntax Error";
  }

  add_md_object(obj);

  return obj;
}

Simulation::~Simulation() {
#ifndef SIM_NOSAVE
  for (int i = 0; i < number_of_objs; ++i) {
    md_objects[i]->save();
  }

  if (fluid != NULL)
    fluid->save();
#endif
  for (int i = 0; i < number_of_objs; ++i) {
    delete md_objects[i];
  }

#ifdef A_FLOWFIELD
  phases_out.close();

  string sname = "flowF.dat";
  fstream out(sname, ios::out);
  out.precision(6);
  cout << "#" << 100 << " " << wc.sizex << " " << wc.sizey << endl;
  out << flowF;
  fflush(stdout);
#endif
#if !defined(MPCD_RANDOM)
  delete fluid;
#endif
}
