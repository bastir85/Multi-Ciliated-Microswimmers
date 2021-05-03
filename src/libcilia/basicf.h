#ifndef BASICF
#define BASICF

#include "MersenneTwister.h"
#include "abbrev.cpp"
#include <fstream>
#include <iostream>
#include <math.h>
#include <omp.h>
#include <string>

using namespace std;
extern MTRand **randi;
// MTRand randi;
//-------------------------------------------------------
//----------------RND----------------------------
//-------------------------------------------------------
inline double get_rand_close() {
#ifdef PARALLEL
  int tid = omp_get_thread_num();
#else
  int tid = 0;
#endif
  return randi[tid]->rand();
}

inline double rnd() {
#ifdef PARALLEL
  int tid = omp_get_thread_num();
#else
  int tid = 0;
#endif
  return randi[tid]->randDblExc();
}

//-------------------------------------------------------
//----------------RNDint----------------------------
//-------------------------------------------------------
inline int ranint() {
// return (rand()%6);
#ifdef PARALLEL
  int tid = omp_get_thread_num();
#else
  int tid = 0;
#endif
  return randi[tid]->randInt(5);
}

//-------------------------------------------------------
//----------------NOrmal Rand------------------------
//-------------------------------------------------------
inline double normalrand() {
  double x1, x2;
#ifdef PARALLEL
  int tid = omp_get_thread_num();
#else
  int tid = 0;
#endif
  x1 = randi[tid]->randDblExc();
  x2 = randi[tid]->randDblExc();
  return sqrt(-2.0 * log(x1)) * cos(2.0 * PI * x2);
}

inline double angle2PI(double angle) {
  static const double twoPi = 2.0 * PI;
  return angle - twoPi * floor(angle / twoPi);
}

inline double normalrand(const double &mean, const double &variance) {
#ifdef PARALLEL
  int tid = omp_get_thread_num();
#else
  int tid = 0;
#endif
  return randi[tid]->randNorm(mean, variance);
}

inline double asin_trunk(double a) {

  if (a <= -1.0) {
    return -PI / 2.0;
  }
  if (a >= 1.0) {
    return PI / 2.0;
  }
  return asinf(a);
}
#endif
