#include "rpy.h"

namespace dashmm {

std::unique_ptr<RPYTable> rpy_table_; 

void update_rpy_table(double a, double kb, double T, double eta) {
  if (rpy_table_ == nullptr) {
    rpy_table_ = 
      std::unique_ptr<RPYTable>{new RPYTable{a, kb, T, eta}}; 
  } 
}

} // namespace dashmm
