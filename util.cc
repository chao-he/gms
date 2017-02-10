#include <math.h>
#include <fstream>
#include "util.h"

const int EXPTABLE_SIZE = 90000;
static const double EPS = 1e-4;

static double ExpTable[EXPTABLE_SIZE];

static bool is_exptable_init = false;

inline void InitExpTable() {
  if (is_exptable_init) return;
  for(int i = 0; i < EXPTABLE_SIZE; ++ i) {
    ExpTable[i] = exp(-EPS*i);
  }
  is_exptable_init = true;
}

double Exp(double e) {
  InitExpTable();
  return e > 0 ? 1 : ExpTable[std::min<int>(EXPTABLE_SIZE - 1, -e/EPS)];
}
