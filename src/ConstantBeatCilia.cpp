//
// Created by basti on 31.08.2016.
//

#include "ConstantBeatCilia.h"
#include "global.h"

void ConstantBeatCilia::switch_logic() {
  if (wc.globalt - Tstroke >= power_stroke_duration && midpnt == start) {
    midpnt += 1;
    cout << "switch to recover mode at " << wc.globalt << endl;
  }
  if (wc.globalt - Tstroke - power_stroke_duration >=
          recovery_stroke_duration &&
      midpnt > start) { // is it always thesame?
    cout << "switch to normal mode at " << wc.globalt << endl;
    Tstroke = wc.globalt;
    midpnt = start;
  }
  if (midpnt > start) { //&&(fmod(wc.globalt+0.1*wc.tstep,m_speed)<wc.tstep)){
    midpnt += wc.tstep / m_speed;
  }
}
