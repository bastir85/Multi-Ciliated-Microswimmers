//#define REPELL 0.66
// umkreisradius gl dreicek=a/sqrt(3)
// REPELL=0.66 enspr c.a. radius von 0.3 (umkreisradius = 0.29)
#define DEBUG

#include <cstdlib>
#include <fenv.h>
#include <fstream>
#include <iostream>
#include <math.h>
#include <string>
#include <time.h>

#include <signal.h>
#include <unistd.h>

using namespace std;

#ifdef PARALLEL
#include <omp.h>
#endif
#define CILIA

#include "abbrev.cpp"
#include "basicf.cpp"
#include "global.cpp"
#include "simulation.h"
void save_seed(string filename) {
#ifdef PARALLEL
  int max_threads = omp_get_max_threads();
  cout << max_threads << " ";
#else
  int max_threads = 1;
#endif
  MTRand::uint32 tmp[MTRand::SAVE];
  ofstream data(filename.c_str(), ios::out | ios::app);

  for (int i = 0; i < max_threads; ++i) {
    randi[i]->save(tmp);
    for (int j = 0; j < MTRand::SAVE; ++j) {
      data << tmp[j] << " ";
    }
    data << endl;
  }
  data.close();
}

void From_cfgFile(string filename) {
  // Loads md_objects from filenmae.ini
  // init md_objects and num_objects!

  Config cfg;

  cout << "loading: " << filename << endl;
  // Read the file. If there is an error, report it and exit.
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

  save_seed("inital_randi.dat");
  // Get the store name.
  const Setting &root = cfg.getRoot();
  try {
    string name = root["name"];
    cout << "Simulation " << name << endl
         << endl
         << " MD Objects:" << endl
         << endl;
    // init global object here for now
    wc.loadparameter(root);
    cout << "T=" << wc.globalt << "->" << wc.desiredt << endl;
    Simulation sim(root);
#ifdef EQUILIBRATE
    if (wc.globalt == 0.0) {
      cout << "equilibrate for T=" << EQUILIBRATE << endl;
      sim.equlibrate(EQUILIBRATE);
    }
#endif
    cout << wc.globalt;
    sim.run();

  } catch (const SettingNotFoundException &nfex) {
    cerr << "No " << nfex.getPath() << " setting in configuration file."
         << endl;
    throw "Config Error";
  }

  wc.saveparameter(root);
  cfg.writeFile(filename.c_str());
  save_seed("end_randi.dat");
}

#ifndef WIN32
#include <execinfo.h>

void handler(int sig) {
  void *array[10];
  size_t size;

  // get void*'s for all entries on the stack
  size = backtrace(array, 10);

  // print out all the frames to stderr
  cerr << "Error: signal: " << sig << endl;
  backtrace_symbols_fd(array, size, STDERR_FILENO);
  exit(1);
}
#endif

int main(int argc, char *argv[]) {
#ifndef WIN32
  signal(SIGSEGV, handler);
#endif
  if (argc != 2) {
    cerr << "usage " << argv[0] << " <configfile>" << endl;
    return 1;
  }

#ifdef PARALLEL
  int max_threads = omp_get_max_threads();
  cout << max_threads << " ";
#else
  int max_threads = 1;
#endif
  randi = new MTRand *[max_threads];
  for (int i = 0; i < max_threads; ++i) {
    randi[i] = new MTRand();
  }
#ifdef DEBUG
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);
#endif
  fflush(stdout);
  From_cfgFile(argv[1]);

  /*
      randi.save(randi_init);
      ofstream out(ranfile,ios::out|ios::binary);
      out.write(reinterpret_cast<const char*>(randi_init), sizeof randi_init);
      out.close();
  */
}
