//
// Created by rode on 12/10/15.
//

#ifndef SIM_UTILS_H
#define SIM_UTILS_H

class SpringCfg {
public:
  SpringCfg(const Row<double, 1> &pos0, const Row<double, 1> &pos1,
            const Row<double, 1> &f0, const Row<double, 1> &f1, double l0)
      : pos0(pos0), pos1(pos1), force0(f0), force1(f1),
        l0(l0){
            //        cout << "create spring:" << pos0[0] << " <-> "<< pos1[0]
            //        << "("<<l0<<")"<<endl;
        };
  SpringCfg(const SpringCfg &&v)
      : pos0(move(v.pos0)), pos1(move(v.pos1)), force0(move(v.force0)),
        force1(move(v.force1)), l0(v.l0){
                                    // cout << "move spring:" << pos0[0] << "
                                    // <-> "<< pos1[0] << "("<<l0<<")"<<endl;
                                };
  void status() {
    double *tmp = pos0.data();
    cout << "pos0.data:" << tmp << "(" << &tmp[0] << ") - " << tmp[0] << endl;
  }

  Matrix<double, 1> pos0;
  Matrix<double, 1> pos1;
  Matrix<double, 1> force0;
  Matrix<double, 1> force1;
  double l0;
};
#endif // SIM_UTILS_H
