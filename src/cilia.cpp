//#define DEBUG
#include "cilia.h"
#include "MatrixIO.h"
#include "basicf.h"
#include "global.h"
#include <assert.h>
#include <fstream>
using namespace Numeric_lib;
using namespace std;

Cilia::Cilia(int rodlength, double beatst, double bending, double min_curv,
             double max_curv, double m_speed, double springconstant, double wx,
             double wy, double wz, double thetaX, double thetaY, double thetaZ,
             MPCD &mpcd, GlobalForces *gf, string name, int startidx,
             int num_anker, bool write_data)
    : MDObject(name, 3, rodlength, beatst, mpcd, gf, startidx, write_data),
      newlength(rodlength), oldlength(rodlength), start(num_anker),
      bending(bending),

      min_bend(min_curv), max_bend(max_curv), m_speed(m_speed),
      springconstant(springconstant), num_anker(num_anker)

{
  midpnt = start;
#ifdef CILIA_PHASE_POWER
  cout << "T dependant phase with T=" << CILIA_PHASE_POWER << endl;
#endif
#ifdef CILIA_PHASE_POSITION
  cout << "Postition dependend phase" << endl;
#endif
  cout << "num anker" << num_anker << endl;
  const double bl = bondlength;
  oldlength = bl;
  if (not open_datfile()) {
    create(wx, wy, wz, thetaX, thetaY, thetaZ);
  }

  newlength = bl;
}

Cilia::~Cilia() {
  cout << "save cilia" << endl;
  MDObject::save();
  fflush(stdout);
};

void Cilia::save(fstream &out) {
  MDObject::save(out);
  out << "#oldlength" << endl;
  out << oldlength;
  out << "#midpnt" << endl;
  out << midpnt << endl;
  out << "#tstroke" << endl;
  out << Tstroke;
}

void Cilia::load(istream &in) {
  string garb;
  MDObject::load(in);
  in >> garb;
  in >> oldlength;
  in >> garb;
  in >> midpnt;
  in >> garb;
  in >> Tstroke;
}

double Cilia::localbend(int a, int b) {
  double l = 0;

  for (int i = a; i < b - 1; i++) {
    l += 2.0 * sqrt(positions[0][i].sqdistance(positions[0][i + 1]));
    l -= sqrt(positions[1][i].sqdistance(positions[1][i + 1]));
    l -= sqrt(positions[2][i].sqdistance(positions[2][i + 1]));
  }
  return l;
}

double Cilia::length(int i) {
  double l = bondlength;
  double tl = sqrt(positions[0][i].sqdistance(positions[0][i + 1]));
  // TODO: improve and use both
  if (i < midpnt - 2.0) {
    // recovery stroke (lower part)
    // l=bondlength+beatstrength*bondlength*pow((-i+rodlength-1.0+num_anker)/(rodlength-num_anker),2.5);
    // l=bondlength+beatstrength*pow(1.0 -
    // (i-num_anker+3.0)/(rodlength-num_anker+2.0),2.5);
    l = bondlength +
        beatstrength *
            pow(1.0 - (i - num_anker - 1.0) / (rodlength - num_anker - 1.0),
                2.5) *
            0.5 *
            (1.0 -
             1.0 / (midpnt - i - 1.0)); // add to smooth movement even more
  } else {
    // power stroke (all parts, since MP==0)
    // was 2
    l = bondlength -
        2.0 * beatstrength /
            (pow((i - midpnt), 2.0) +
             1.0); // CHECK?? has been one but does not make sense was 4
    // l=bondlength-2.0*beatstrength/pow((i-midpnt)+1,1); //CHECK?? has been one
    // but does not make sense was 4
  }
  // stall force
  // if ((fabs(oldlength[i]-l)>0.0001)&&(fabs(l-tl)>bending)){
  if ((fabs(l - tl) > bending)) {
    tl = tl + bending - 2.0 * bending * (l < tl);
    if (((l - oldlength[i]) * (tl - oldlength[i])) > 0) {
      l = tl;
    } else {
      // Keep bondlength when going backward or in case of no activity.
      // stall when no bl change to allow proper fluctuations
      l = oldlength[i];
    }
  } else {
    // cout << "no stall force" << endl;
  }
  // assert(fabs(l-bondlength)< 1e-5);
  return l;
}

//-------------------------------------------------------
//--------------Forces     -------------
//-------------------------------------------------------
//
void Cilia::addcross() {
  // double diagonal_active=2.0*bondlength*sin(PI/3.0);
  // contacting rods to each other
  //#pragma omp for schedule (guided) nowait
  for (int i = 0; i < rodlength; ++i) {
    spring(i, i + rodlength, bondlength, springconstant);
    spring(i, i + 2 * rodlength, bondlength, springconstant);
    spring(i + rodlength, i, bondlength, springconstant);
    spring(i + rodlength, i + 2 * rodlength, bondlength, springconstant);
    spring(i + 2 * rodlength, i, bondlength, springconstant);
    spring(i + 2 * rodlength, i + rodlength, bondlength, springconstant);
    // add 2*2 new ones to rod Nr.3!
    /*spring(i+2*rodlength,i+3*rodlength,bondlength,crossstrength*springconstant);
    spring(i+3*rodlength,i+rodlength,bondlength,crossstrength*springconstant);
    spring(i+rodlength,i+3*rodlength,bondlength,crossstrength*springconstant);
    spring(i+3*rodlength,i+2*rodlength,bondlength,crossstrength*springconstant);*/
    // 0->3 diagonal
    // spring(i,i+3*rodlength,diagonal_active,crossstrength*springconstant);
    // spring(i+3*rodlength,i,diagonal_active,crossstrength*springconstant);
  }

  // foreward diagonals
  // now with adapted crosslenghts
  // double lengthi = crosslength;
  //#pragma omp for schedule (guided) nowait
  for (int i = 0; i < rodlength - 1; ++i) {
    double lengthi = sqrt(bondlength * bondlength + bondlength * newlength[i]);
    spring(i, i + rodlength + 1, lengthi, crossstrength * springconstant);
    spring(i, i + 2 * rodlength + 1, lengthi, crossstrength * springconstant);
    spring(i + rodlength, i + 1, lengthi, crossstrength * springconstant);
    spring(i + 2 * rodlength, i + 1, lengthi, crossstrength * springconstant);

    // Do we need and want this ones cross 1->2 lets use them for stability!
    spring(i + rodlength, i + 2 * rodlength + 1, crosslength,
           crossstrength * springconstant);
    spring(i + 2 * rodlength, i + rodlength + 1, crosslength,
           crossstrength * springconstant);
    // rod Nr.3
    /*    lengthi=sqrt(bondlength*bondlength+bondlength*newlN); //TODO: check
       symmetry?

        spring(i+3*rodlength,i+rodlength+1,lengthi,crossstrength*springconstant);
        spring(i+3*rodlength,i+2*rodlength+1,lengthi,crossstrength*springconstant);

        spring(i+rodlength,i+3*rodlength+1,lengthi,crossstrength*springconstant);
        spring(i+2*rodlength,i+3*rodlength+1,lengthi,crossstrength*springconstant);*/
  }

  // backward diagonals
  //#pragma omp for schedule (guided) nowait
  for (int i = 1; i < rodlength; ++i) {
    double lengthi =
        sqrt(bondlength * bondlength + bondlength * newlength[i - 1]);
    spring(i, i + rodlength - 1, lengthi, crossstrength * springconstant);
    spring(i, i + 2 * rodlength - 1, lengthi, crossstrength * springconstant);

    spring(i + rodlength, i - 1, lengthi, crossstrength * springconstant);
    spring(i + 2 * rodlength, i - 1, lengthi, crossstrength * springconstant);

    // Do we need and want this ones cross 1->2 lets use them for stability!
    spring(i + rodlength, i + 2 * rodlength - 1, crosslength,
           crossstrength * springconstant);
    spring(i + 2 * rodlength, i + rodlength - 1, crosslength,
           crossstrength * springconstant);

    // rod Nr.3
    /*    lengthi=sqrt(bondlength*bondlength+bondlength*newlN);
        spring(i+3*rodlength,i+rodlength-1,lengthi,crossstrength*springconstant);
        spring(i+3*rodlength,i+2*rodlength-1,lengthi,crossstrength*springconstant);

        spring(i+rodlength,i+3*rodlength-1,lengthi,crossstrength*springconstant);
        spring(i+2*rodlength,i+3*rodlength-1,lengthi,crossstrength*springconstant);*/
  }

} // ende von addcross

// newlength gibt spring eq länge, d.h. aktive force depends on midpoint
// => springs @first cilium and (during recovery) 2nd first and 7th
void Cilia::recalcforce(void) {
  if (fmod(wc.globalt + wc.ststep / 10., 0.01) <
      wc.ststep) { // sonst abhaengigkeit des Schlagmusters isb periode von
                   // ststep
    // newlengthcalc
    //        #pragma omp for schedule (guided) nowait
    for (int i = num_anker - 1; i < rodlength - 1; ++i) {
      newlength[i] = length(i);
      oldlength[i] = newlength[i];
    }
  }

  // Reset forces
  forces = 0.0;
  // active forces

  //    #pragma omp barrier

  //    #pragma omp for schedule (guided) nowait
  for (int i = 0; i < rodlength - 1; i++) {
    spring(i, i + 1, newlength[i], springconstant);
    spring(i + 1, i, newlength[i], springconstant);
    // no 4th rod
    /*spring(3*rodlength+i,3*rodlength+i+1,newlN,springconstant);
    spring(3*rodlength+i+1,3*rodlength+i,newlN,springconstant);*/
    spring(rodlength + i, rodlength + i + 1, bondlength, springconstant);
    spring(rodlength + i + 1, rodlength + i, bondlength, springconstant);
    spring(2 * rodlength + i, 2 * rodlength + i + 1, bondlength,
           springconstant);
    spring(2 * rodlength + i + 1, 2 * rodlength + i, bondlength,
           springconstant);
  }

  addcross();
}

void Cilia::switch_logic() {
  // if ((localbend(0,rodlength)<min_bend)&& pow(midpnt-start,2) < 0.0001){
  // midpnt+=1;
  if ((midpnt == start) && wc.globalt - Tstroke > 100.0) {
    midpnt += 1;
    cout << "switch to recover mode at " << wc.globalt << endl;
    cout << "m/n: " << midpnt << "/" << localbend(0, rodlength) << endl;
  }
  if ((localbend(0, rodlength) > max_bend) && (midpnt > start)) {
    cout << "switch to normal mode at " << wc.globalt << endl;
    cout << "m/n: " << midpnt << "/" << localbend(0, rodlength) << endl;
    Tstroke = wc.globalt;
    midpnt = start;
  }
  if (midpnt > start) { //&&(fmod(wc.globalt+0.1*wc.tstep,m_speed)<wc.tstep)){
    midpnt += wc.tstep / m_speed;
  }
}

void Cilia::debug() {
  // PR(localbend(0,rodlength));
  cout << "C=" << localbend(0, rodlength) << endl;
  cout << "midpnt:" << midpnt << endl;
}

istream &operator>>(istream &in, Cilia &cilia) {
  /*
   *  usage is >> cilia;
   */
  cilia.load(in);

  return in;
}

void Cilia::create(double wx, double wy, double wz, double thetaX,
                   double thetaY, double thetaZ) {
  Matrix<double, 1> vec(3);
  vec[0] = wx;
  vec[1] = wy;
  vec[2] = wz;
  positions = 0.0;
  velocities = 0.0;

  double l0 = bondlength * sin(PI / 3.0) * 2. / 3.;

  // inserting 3 rods
  for (int i = 0; i < rodlength; ++i) {
    positions[0][i][2] = bondlength * i;
    positions[0][i][0] -= l0;

    positions[1][i][2] = bondlength * i;
    positions[1][i][1] = bondlength / 2.0;
    positions[1][i][0] = bondlength * sin(PI / 3.0) - l0;

    positions[2][i][2] = bondlength * i;
    positions[2][i][1] = -bondlength / 2.0;
    positions[2][i][0] = bondlength * sin(PI / 3.0) - l0;
  }

  for (int j = 0; j < num_rods; ++j) {
    for (int i = 0; i < rodlength; ++i) {
      double xi = positions[j][i][0];
      double yi = positions[j][i][1];
      double zi = positions[j][i][2];
      positions[j][i][0] = zi * sin(thetaY) +
                           cos(thetaY) * (xi * cos(thetaZ) - yi * sin(thetaZ));
      positions[j][i][1] =
          cos(thetaX) * (yi * cos(thetaZ) + xi * sin(thetaZ)) -
          sin(thetaX) * (zi * cos(thetaY) + sin(thetaY) * (-(xi * cos(thetaZ)) +
                                                           yi * sin(thetaZ)));
      positions[j][i][2] =
          sin(thetaX) * (yi * cos(thetaZ) + xi * sin(thetaZ)) +
          cos(thetaX) * (zi * cos(thetaY) + sin(thetaY) * (-(xi * cos(thetaZ)) +
                                                           yi * sin(thetaZ)));
    }
  }

  for (int i = 0; i < num_rods; ++i) {
    for (int j = 0; j < rodlength; ++j) {
      positions[i][j][0] += wx;
      positions[i][j][1] += wy;
      positions[i][j][2] += wz;
    }
  }
}

#ifdef CILIA_PHASE_POSITION
#include <complex>

complex<double> Sn[] = {{1.0, 0.0},
                        {-0.015715716596, 0.0527879904963},
                        {-0.320133792036, -0.374825972894},
                        {0.244078112964, -0.145080891406},
                        {0.000801816544107, 0.228011094098},
                        {-0.290283207694, -0.135744155591},
                        {0.0600450128892, -0.20206341677},
                        {0.0232793738293, 0.238827200559},
                        {-0.228832824602, 0.116494753182},
                        {0.0811502695887, -0.130468260393},
                        {0.237124773952, 0.1390778565},
                        {-0.078692522666, 0.0947805797355},
                        {-0.0105884721343, -0.207543342063},
                        {0.162623349652, -0.0380313930705},
                        {-0.0800284157403, 0.128087721476},
                        {-0.148001105874, -0.0969006924129},
                        {0.122234390426, -0.0617788976376},
                        {0.0241496355132, 0.158723446476},
                        {-0.131802852521, -0.0487363459862},
                        {0.104078291599, -0.130444889418},
                        {0.0629746068851, 0.104763230756},
                        {-0.138288148555, 0.020506187496},
                        {-0.0093425029036, -0.134997387238},
                        {0.10283616895, 0.0658602470066},
                        {-0.0899719706342, 0.108149292689},
                        {-0.0581434785702, -0.0990394730391},
                        {0.134825995379, 0.0099807039945},
                        {-0.0119772041612, 0.106262099841},
                        {-0.0738890450486, -0.0886708510424},
                        {0.0942835502419, -0.0926134339427},
                        {0.017554557549, 0.101699839981},
                        {-0.145925273238, -0.0134375282806},
                        {0.0227122640155, -0.0973257818848},
                        {0.0641685881163, 0.126967055477},
                        {-0.0727021305514, 0.0626567094682},
                        {0.0215488934567, -0.117433932411},
                        {0.110407926355, 0.0222485731158},
                        {-0.0485391951886, 0.0475553689104},
                        {-0.0419395008576, -0.106633349876},
                        {0.0613870287221, -0.0364282568586},
                        {-0.0383274133929, 0.0939010824999},
                        {-0.0869814785589, -0.0415797529204},
                        {0.0763747542522, -0.0368808551222},
                        {0.00403604960112, 0.0874290304071},
                        {-0.0898125253082, 0.0275724518805},
                        {0.0556038318401, -0.0585787132544},
                        {0.0899281541357, 0.0362779008726},
                        {-0.0744716286443, 0.00145116175692},
                        {-0.045004813419, -0.0773885878563},
                        {0.0752265131062, 0.042975563151},
                        {-0.0221918062612, 0.0898910399115}};
double Cilia::getPhase() {
  complex<double> sum = 0.0;
  double psi0, psi1, phi;
  const std::complex<double> j(0, 1);
  double avapos[3][2] = {0, 0, 0, 0, 0, 0}; // 0,24,49 xz
  for (int i = 0; i < num_rods; ++i) {
    avapos[0][0] += positions[i][0][0];
    avapos[0][1] += positions[i][0][2];
    avapos[1][0] += positions[i][24][0];
    avapos[1][1] += positions[i][24][2];
    avapos[2][0] += positions[i][49][0];
    avapos[2][1] += positions[i][49][2];
  }

  avapos[0][0] /= num_rods;
  avapos[0][1] /= num_rods;
  avapos[1][0] /= num_rods;
  avapos[1][1] /= num_rods;
  avapos[2][0] /= num_rods;
  avapos[2][1] /= num_rods;
  psi0 =
      angle2PI(atan2(avapos[1][1] - avapos[0][1], avapos[1][0] - avapos[0][0]));
  //- 1.4;
  /*    cout << avapos[1][1] - avapos[0][1] << endl;
      cout << avapos[1][0] - avapos[0][0] << endl;
      cout << psi0 << " "<< psi0 -1.4 << endl;*/

  psi1 =
      angle2PI(atan2(avapos[2][1] - avapos[1][1], avapos[2][0] - avapos[1][0]));
  //- 1.8;
  /*cout << avapos[2][1] - avapos[1][1] << endl;
  cout << avapos[2][0] - avapos[1][0] << endl;
  cout << psi1 << " "<< psi1 -1.8 <<endl;
  cout << endl;*/
  phi = atan2(psi1 - 1.8, psi0 - 1.4);
  for (int i = 1; i < 51; ++i) {
    sum += Sn[i] / j / complex<double>(i, 0) *
           (exp(j * complex<double>(i, 0) * phi) - 1.0);
  }
  return angle2PI(phi + 2.0 * real(sum));
  // Map to 0,2PI intervall since we want the current rel phase.
}
#endif

#ifdef CILIA_PHASE_POWER
double Cilia::getPhase() {
  const double T = CILIA_PHASE_POWER;
  double phase = (wc.globalt - Tstroke) / T;
  if (phase > 1.0)
    return 2 * PI;
  return 2 * PI * phase;
}
#endif

Matrix<double, 2> Cilia::getBoundingBox() {
  Matrix<double, 2> minmax(2, 3);
  minmax[0] = 999999; // min
  minmax[1] = 0.0;    // max
  Matrix<double, 2> ava_rod = get_avarage_rod();
  for (int i = 0; i < rodlength; ++i) {
    Row<double, 1> ava = ava_rod[i];
    if (ava[0] > minmax[1][0]) {
      minmax[1][0] = ava[0];
    }
    if (ava[1] > minmax[1][1]) {
      minmax[1][1] = ava[1];
    }
    if (ava[2] > minmax[1][2]) {
      minmax[1][2] = ava[2];
    }

    if (ava[0] < minmax[0][0]) {
      minmax[0][0] = ava[0];
    }
    if (ava[1] < minmax[0][1]) {
      minmax[0][1] = ava[1];
    }
    if (ava[2] < minmax[0][2]) {
      minmax[0][2] = ava[2];
    }
  }

  return minmax;
}
