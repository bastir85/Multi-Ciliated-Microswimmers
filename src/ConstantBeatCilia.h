//
// Created by basti on 31.08.2016.
//

#ifndef SIM_CONSTANTBEATCILIA_H
#define SIM_CONSTANTBEATCILIA_H

#include "cilia.h"

class ConstantBeatCilia : public Cilia {
  friend class Clami;

public:
  ConstantBeatCilia(const Setting &cfg, MPCD &mpcd, GlobalForces *gf)
      : Cilia(cfg, mpcd, gf),
        power_stroke_duration(cfg["power_stroke_duration"]),
        recovery_stroke_duration(cfg["recovery_stroke_duration"]) {
    if (!loaded) {
      Tstroke = cfg["start_delay"];
    }
  }
  void switch_logic();

private:
  const double power_stroke_duration;
  const double recovery_stroke_duration;
};

#endif // SIM_CONSTANTBEATCILIA_H
