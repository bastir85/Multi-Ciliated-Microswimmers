#include "globalforces.h"
#include "MatrixIO.h"
#include "basicf.h"
#include "global.h"

GlobalForces::GlobalForces():
    num_objects(0),
    max_integration_error(0.0)
{
    #ifdef REPELLING
    cout << "activate vdw force with eps=" << epsilon << endl;
    #endif
}

void GlobalForces::set_md_objects(MDObject** md_objects, int num_objects){

    this->md_objects = md_objects;
    this->num_objects = num_objects;

}

#ifdef REPELLING
void inline GlobalForces::vdw_force(int cil0, int cil1, int beat0, int beat1) {
    //cout << "add force " << beat0 << " " << beat1 << endl;
    //for(int k=beat0-N/2; k < beat0+N/2 and k <= beat1;++k) {
    for(int k=beat0-N/2; k < beat0+N/2;++k) {
        Matrix<double,1> x0= cilias[cil0]->positions[0][k]+ cilias[cil0]->positions[1][k]+ cilias[cil0]->positions[2][k];
        for(int l=beat1+N/2; l >= beat1-N/2 and l >= k ;--l) {
            Matrix<double,1> delta = (x0 - (cilias[cil1]->positions[0][l]+
                                         cilias[cil1]->positions[1][l]+ cilias[cil1]->positions[2][l])
                                  )/3.0;

            double dis = sqrt(dot_product(delta,delta));

            /*double force = 24.*(2*pow(sigma,12)/pow(dis,13)-pow(sigma,6)/pow(dis,7))
                         /dis/3;*/

            //double force = 3*epsilon*(pow(sigma,12)/pow(dis,14));
            // just use harmonic
            double force =-2.0*epsilon*(1.0-cutoff/dis) *(dis < cutoff);
            //cout << "F " << k << "<->" << l << ":"<< force << endl;
            for(int rod=0;rod<cilias[cil0]->num_rods; ++rod) {
                //TODO: make parallel Watch out, race conditions possible!!!
                #pragma omp critical
                cilias[cil0]->forces[rod][k] += delta*force;
            }

            for(int rod=0;rod<cilias[cil1]->num_rods; ++rod) {
                //TODO: make parallel Watch out, race conditions possible!!!
                #pragma omp critical
                cilias[cil1]->forces[rod][l] -= delta*force;
            }
        }
    }
}

#endif

void GlobalForces::recalc_forces() {
    #ifdef REPELLING
    //boxing for short repulsive potential between MD particles
    #pragma omp for schedule (guided) nowait
    for(int i=0; i< num_objects; ++i) {
       Matrix<double,2> bb1 = cilias[i]->getBoundingBox();
       //cout << bb1 << endl;
        for(int j=num_objects-1; j>i; --j) {
            //cout << "cil:" << i << "->" << j << endl;
            Matrix<double,2> bb2 = cilias[j]->getBoundingBox();
            //cout << bb2 <<endl ;
            if(    (bb2[0][0] - bb1[1][0] < cutoff and bb1[0][0] - bb2[1][0] < cutoff)  //x1_min > x2_max and // x2_min > x1_max
               and (bb2[0][1] - bb1[1][1] < cutoff and bb1[0][1] - bb2[1][1] < cutoff)  // " for y
               and (bb2[0][2] - bb1[1][2] < cutoff and bb1[0][2] - bb2[1][2] < cutoff)  // " for z
               ){  // " for z
                //cout << "BB overlap " << endl;
                //x1_min  - x2_max >- Thresh. x2_min - x1_max > -Thresh
                //
                for(int k=1+N/2; k < cilias[i]->rodlength;k+=N) {
                    for(int l=cilias[j]->rodlength-1-N/2; l >= k; l-=N) {
                    //cout << "sp dis \t " << k << "<->" << l << endl;
                    //cout << sqrt(cilias[i]->positions[0][k].sqdistance(cilias[j]->positions[0][l]))<< endl;
                        if(sqrt(cilias[i]->positions[0][k].sqdistance(cilias[j]->positions[0][l])) < N*cutoff) {
                            //calulate rep force between cilia i,j
                            vdw_force(i,j,k,l);
                        }
                    }

                }
            }
        }

    }
    //cout << "end recaclc F" << endl <<endl;
    //fflush(stdout);
    #endif
}

double GlobalForces::get_max_sqforce() {
    double max_force = 0.0 ;
    for(int i=0; i< num_objects; ++i) {
        MDObject* obj = md_objects[i];
        for(int j=0; j< obj->num_rods; ++j) {
            for(int k=0; k< obj->rodlength; ++k) {
                double tmp = dot_product(obj->forces[j][k],obj->forces[j][k]);
                if (max_force < tmp) {
                    max_force = tmp;
                }
            }
        }
    }
    return max_force;
}



/*
 * This is only working serial, for parallel exceution think! Use a list of errors and reduce afterwards!
 */
void GlobalForces::set_integration_error(Vector3d& rel_error){
    for(int i=0; i<3;++i){
        double rel_err_abs = abs(rel_error[i]);
        if(max_integration_error < rel_err_abs){
            max_integration_error = rel_err_abs;
   //         cout << "max err= " << max_integration_error << endl;
        }
    }
}
#define TOL 1e-4
void GlobalForces::adapt_timestep(){
    double new_tstep = 2.0*0.9*TOL/max_integration_error;
    if (new_tstep > wc.tstep) new_tstep = wc.tstep;
    if (new_tstep/wc.ststep > 2.) new_tstep = wc.ststep*2.;

    wc.ststep = new_tstep;

    /*else if (wc.ststep < 1e-7) {
    //    cout << "ERROR" << endl;
        wc.ststep = 1e-7;
    }*/

//    cout << scientific << max_integration_error << endl;
//    cout << "Adapted timestep is " << scientific << wc.ststep << fixed<<endl;
    max_integration_error = 0.0;

}

