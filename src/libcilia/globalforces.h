#ifndef GLOBAL_FORCES_H
#define GLOBAL_FORCES_H
#include "Matrix.h"
#include "Vector.h"
#include <string.h>
#include "mdobject.h"
using namespace Numeric_lib;
using namespace std;

class GlobalForces {
    public:
        GlobalForces();
        void set_md_objects(MDObject** cilias, int num);
        void recalc_forces();
        Vector3d flowField(const Matrix<double,1>& pos);
        double get_max_sqforce();
        void set_integration_error(Vector3d& rel_error);
        void adapt_timestep();
    private:
        MDObject ** md_objects;
        int num_objects;
        #ifdef REPELLING
        constexpr static double cutoff=1.5;
        constexpr static int N = 4.0; // step size for spheres
        constexpr static double sigma=cutoff /2.5;
        constexpr static double epsilon=REPELLING;

        void inline vdw_force(int cil0, int cil1, int beat0, int beat1);
        #endif

        double max_integration_error;

};
#endif
