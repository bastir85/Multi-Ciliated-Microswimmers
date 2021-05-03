#ifndef BASICF
#define BASICF

#include "MersenneTwister.h"
#include "abbrev.cpp"
#include <assert.h>
#include <fstream>
#include <iostream>
#include <math.h>
#include <omp.h>
#include <string>

using namespace std;
MTRand **randi;
#ifdef DEBUG
long unsigned int saat = 0;
// MTRand randi(saat);
#else
// MTRand randi[];
#endif
//-------------------------------------------------------
//----------------Squared distance------------------------
//-------------------------------------------------------
// function calculates squared distance
// between particles i and j
double sqdis(int i, int j) {
  assert(0);
  return 0;
  /*  return ( (xparticle[i]-xparticle[j])*(xparticle[i]-xparticle[j])
             +(yparticle[i]-yparticle[j])*(yparticle[i]-yparticle[j])
             +(zparticle[i]-zparticle[j])*(zparticle[i]-zparticle[j]));  */
}
//-------------------------------------------------------
//----------------min Squared distance------------------------
//-------------------------------------------------------
// function calculates squared distance
// between particles i and j, periodic in x and y
double xysqdis(int i, int j) {
  assert(0);
  /*  double
    dx=(xparticle[i]-xparticle[j])-sizex*round((xparticle[i]-xparticle[j])/sizex);
    double
    dy=(yparticle[i]-yparticle[j])-sizey*round((yparticle[i]-yparticle[j])/sizey);
    double dz=(zparticle[i]-zparticle[j]);
    return dx*dx+dy*dy+dz*dz;*/
  return 0;
}

//-------------------------------------------------------
//----------------RND----------------------------
//-------------------------------------------------------
inline double rnd() {
  int tid = omp_get_thread_num();
  return randi[tid]->randDblExc();
}

//-------------------------------------------------------
//----------------RNDint----------------------------
//-------------------------------------------------------
inline int ranint() {
  // return (rand()%6);
  int tid = omp_get_thread_num();
  return randi[tid]->randInt(5);
}

//-------------------------------------------------------
//----------------NOrmal Rand------------------------
//-------------------------------------------------------
inline double normalrand() {
  double x1, x2;
  int tid = omp_get_thread_num();
  x1 = randi[tid]->randDblExc();
  x2 = randi[tid]->randDblExc();
  return sqrt(-2.0 * log(x1)) * cos(2.0 * PI * x2);
}

#endif
