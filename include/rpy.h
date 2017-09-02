#ifndef __RPY_EXPANSION_H__
#define __RPY_EXPANSION_H__

#include <cmath>
#include <complex>
#include <map>
#include <memory>
#include <vector>

#include "dashmm/index.h"
#include "builtins/laplace.h"
#include "dashmm/point.h"
#include "dashmm/types.h"
#include "dashmm/viewset.h"

namespace dashmm {

class RPYTable {
public:
  RPYTable(double a, double kb, double T, double eta) 
    : a_{a}, kb_{kb}, T_{T}, eta_{eta} {
      c0_ = kb * T / 6 / M_PI / a / eta; 
      c1_ = kb * T / 8 / M_PI / eta; 
      c2_ = kb * T  * a * a / 12 / M_PI;
    }

  double a() const {return a_;}  
  double c0() const {return c0_;}
  double c1() const {return c1_;}
  double c2() const {return c2_;}
  bool update(double a, double kb, double T, double eta) const {
    return (a != a_) || (kb * T / eta != kb_ * T_ / eta_); 
  }

private: 
  double a_;     // Radius of the bead 
  double kb_;    // Boltzmann constant 
  double T_;     // Absolute temperate
  double eta_;   // Solvent viscosity
  double c0_; 
  double c1_; 
  double c2_;
};

extern std::unique_ptr<RPYTable> rpy_table_; 

void update_rpy_table(double a, double kb = 1.0, 
                      double T = 1.0, double eta = 1.0); 

struct Bead {
  Bead() { }
  int index;                // Index of the bead
  dashmm::Point position;   // Position of the bead
  double q[3] = {0.0};      // "charges" carried by the bead
  double value[3] = {0.0}; 
}; 

template <typename Source, typename Target> 
class RPY {
public: 
  using source_t = Source; 
  using target_t = Target; 
  using expansion_t = RPY<Source, Target>; 

  RPY(ExpansionRole role, double scale = 1.0, Point center = Point{}) 
    : views_{ViewSet{role, center, scale}} {
    // The RPY class applies the Laplace kernel to create four views
    // Views 0, 1, 2 are generated using q[0], q[1], and q[2] as the "charge"
    // Views 3 is generated using -c1 * (q[0] * x + q[1] * y + q[2] * z) 

    // View size for each spherical harmonic expansion
    int p = builtin_laplace_table_->p(); 
    int nsh = (p + 1) * (p + 2) / 2; 

    // View size for each exponential expansion
    int nexp = builtin_laplace_table_->nexp(); 

    if (role == kSourcePrimary || role == kTargetPrimary) {
      size_t bytes = sizeof(dcomplex_t) * nsh; 
      for (int i = 0; i < 4; ++i) {
        char *data = new char[bytes](); 
        views_.add_view(i, bytes, data); 
      }
    } else if (role == kSourceIntermediate) {
      size_t bytes = sizeof(dcomplex_t) * nexp; 
      for (int i = 0; i < 24; ++i) {
        char *data = new char[bytes](); 
        views_.add_view(i, bytes, data);
      }
    } else if (role == kTargetIntermediate) {
      size_t bytes = sizeof(dcomplex_t) * nexp; 
      for (int i = 0; i < 112; ++i) {
        char *data = new char[bytes](); 
        views_.add_view(i, bytes, data);
      }
    }
  }

  RPY(const ViewSet &views) : views_{views} { }

  ~RPY() {
    int count = views_.count();
    if (count) {
      for (int i = 0; i < count; ++i) {
        delete [] views_.view_data(i);
      }
    }
  }

  void release() {views_.clear();}

  bool valid(const ViewSet &view) const {
    // \p view is assumed to be a subset of \p views_ (no range checking
    // performed). The function returns true if and only if each entry in the
    // required subset is associated with some data.
    bool is_valid = true;
    int count = view.count();
    for (int i = 0; i < count; ++i) {
      int idx = view.view_index(i);
      if (views_.view_data(idx) == nullptr) {
        is_valid = false;
        break;
      }
    }
    return is_valid;
  }
  
  int view_count() const { return views_.count(); }
  
  ViewSet get_all_views() const {return views_;}
  
  ExpansionRole role() const {return views_.role();}

  Point center() const {return views_.center();}

  size_t view_size(int view) const {
    return views_.view_bytes(view) / sizeof(dcomplex_t);
  }
  
  dcomplex_t view_term(int view, size_t i) const {
    dcomplex_t *data = reinterpret_cast<dcomplex_t *>(views_.view_data(view));
    return data[i];
  }  

  std::unique_ptr<expansion_t> S_to_M(const Source *first, 
                                      const Source *last) const {
    Point center = views_.center(); 
    double scale = views_.scale(); 
    double c1 = rpy_table_->c1(); 
    expansion_t *ret{new expansion_t{kSourcePrimary, scale, center}}; 
    dcomplex_t *M1 = reinterpret_cast<dcomplex_t *>(ret->views_.view_data(0)); 
    dcomplex_t *M2 = reinterpret_cast<dcomplex_t *>(ret->views_.view_data(1)); 
    dcomplex_t *M3 = reinterpret_cast<dcomplex_t *>(ret->views_.view_data(2)); 
    dcomplex_t *M4 = reinterpret_cast<dcomplex_t *>(ret->views_.view_data(3)); 

    for (auto i = first; i != last; ++i) {
      Point dist = point_sub(i->position, center); 
      lap_s_to_m(dist, i->q[0], scale, M1); 
      lap_s_to_m(dist, i->q[1], scale, M2); 
      lap_s_to_m(dist, i->q[2], scale, M3); 
      double q4 = -c1 * (i->q[0] * i->position.x() + 
                         i->q[1] * i->position.y() + 
                         i->q[2] * i->position.z()); 
      lap_s_to_m(dist, q4, scale, M4); 
    }

    return std::unique_ptr<expansion_t>{ret}; 
  }


  std::unique_ptr<expansion_t> S_to_L(const Source *first, 
                                      const Source *last) const {
    Point center = views_.center(); 
    double scale = views_.scale(); 
    double c1 = rpy_table_->c1();
    expansion_t *ret{new expansion_t{kTargetPrimary}}; 
    dcomplex_t *L1 = reinterpret_cast<dcomplex_t *>(ret->views_.view_data(0)); 
    dcomplex_t *L2 = reinterpret_cast<dcomplex_t *>(ret->views_.view_data(1)); 
    dcomplex_t *L3 = reinterpret_cast<dcomplex_t *>(ret->views_.view_data(2)); 
    dcomplex_t *L4 = reinterpret_cast<dcomplex_t *>(ret->views_.view_data(3)); 

    for (auto i = first; i != last; ++i) {
      Point dist = point_sub(i->position, center); 
      lap_s_to_l(dist, i->q[0], scale, L1); 
      lap_s_to_l(dist, i->q[1], scale, L2); 
      lap_s_to_l(dist, i->q[2], scale, L3); 
      double q4 = -c1 * (i->q[0] * i->position.x() + 
                         i->q[1] * i->position.y() + 
                         i->q[2] * i->position.z()); 
      lap_s_to_l(dist, q4, scale, L4); 
    }

    return std::unique_ptr<expansion_t>{ret};
  }

  
  std::unique_ptr<expansion_t> M_to_M(int from_child) const {
    dcomplex_t *M1 = reinterpret_cast<dcomplex_t *>(views_.view_data(0));
    dcomplex_t *M2 = reinterpret_cast<dcomplex_t *>(views_.view_data(0));
    dcomplex_t *M3 = reinterpret_cast<dcomplex_t *>(views_.view_data(0));
    dcomplex_t *M4 = reinterpret_cast<dcomplex_t *>(views_.view_data(0));

    expansion_t *ret{new expansion_t{kSourcePrimary}};
    dcomplex_t *W1 = reinterpret_cast<dcomplex_t *>(ret->views_.view_data(0));
    dcomplex_t *W2 = reinterpret_cast<dcomplex_t *>(ret->views_.view_data(0));
    dcomplex_t *W3 = reinterpret_cast<dcomplex_t *>(ret->views_.view_data(0));
    dcomplex_t *W4 = reinterpret_cast<dcomplex_t *>(ret->views_.view_data(0));

    lap_m_to_m(from_child, M1, W1);
    lap_m_to_m(from_child, M2, W2);
    lap_m_to_m(from_child, M3, W3);
    lap_m_to_m(from_child, M4, W4);

    return std::unique_ptr<expansion_t>{ret};
  }
 
  std::unique_ptr<expansion_t> M_to_L(Index s_index, Index t_index) const {
    return std::unique_ptr<expansion_t>{nullptr};
  }

  std::unique_ptr<expansion_t> L_to_L(int to_child) const {
    dcomplex_t *L1 = reinterpret_cast<dcomplex_t *>(views_.view_data(0));
    dcomplex_t *L2 = reinterpret_cast<dcomplex_t *>(views_.view_data(1));
    dcomplex_t *L3 = reinterpret_cast<dcomplex_t *>(views_.view_data(2));
    dcomplex_t *L4 = reinterpret_cast<dcomplex_t *>(views_.view_data(3));

    expansion_t *ret{new expansion_t{kTargetPrimary}};
    dcomplex_t *W1 = reinterpret_cast<dcomplex_t *>(ret->views_.view_data(0));
    dcomplex_t *W2 = reinterpret_cast<dcomplex_t *>(ret->views_.view_data(1));
    dcomplex_t *W3 = reinterpret_cast<dcomplex_t *>(ret->views_.view_data(2));
    dcomplex_t *W4 = reinterpret_cast<dcomplex_t *>(ret->views_.view_data(3));

    lap_l_to_l(to_child, L1, W1);
    lap_l_to_l(to_child, L2, W2);
    lap_l_to_l(to_child, L3, W3);
    lap_l_to_l(to_child, L4, W4);

    return std::unique_ptr<expansion_t>{ret};
  }  

  void M_to_T(Target *first, Target *last) const {
    double scale = views_.scale();
    dcomplex_t *M1 = reinterpret_cast<dcomplex_t *>(views_.view_data(0));
    dcomplex_t *M2 = reinterpret_cast<dcomplex_t *>(views_.view_data(0));
    dcomplex_t *M3 = reinterpret_cast<dcomplex_t *>(views_.view_data(0));
    dcomplex_t *M4 = reinterpret_cast<dcomplex_t *>(views_.view_data(0));

    for (auto i = first; i != last; ++i) {
      Point dist = point_sub(i->position, views_.center()); 
      auto result = rpy_m_to_t(dist, scale, M1, M2, M3, M4); 

      i->value[0] += result[0]; 
      i->value[1] += result[1];
      i->value[2] += result[2];
    }
  }

  void L_to_T(Target *first, Target *last) const {
    double scale = views_.scale();
    dcomplex_t *L1 = reinterpret_cast<dcomplex_t *>(views_.view_data(0));
    dcomplex_t *L2 = reinterpret_cast<dcomplex_t *>(views_.view_data(1));
    dcomplex_t *L3 = reinterpret_cast<dcomplex_t *>(views_.view_data(2));
    dcomplex_t *L4 = reinterpret_cast<dcomplex_t *>(views_.view_data(3));
    
    for (auto i = first; i != last; ++i) {
      Point dist = point_sub(i->position, views_.center()); 
      auto result = rpy_l_to_t(dist, scale, L1, L2, L3, L4); 
      
      i->value[0] += result[0]; 
      i->value[1] += result[1]; 
      i->value[2] += result[2]; 
    }
  }

  void S_to_T(const Source *s_first, const Source *s_last, 
              Target *t_first, Target *t_last) const {
    double a = rpy_table_->a(); 
    double c0 = rpy_table_->c0(); 
    double c1 = rpy_table_->c1();
    double c2 = rpy_table_->c2(); 

    for (auto i = t_first; i != t_last; ++i) {
      double v0 = 0.0, v1 = 0.0, v2 = 0.0; 
      for (auto j = s_first; j != s_last; ++j) {
        double x = i->position.x() - j->position.x(); 
        double y = i->position.y() - j->position.y(); 
        double z = i->position.z() - j->position.z(); 
        double q0 = j->q[0]; 
        double q1 = j->q[1];
        double q2 = j->q[2]; 
        double r = sqrt(x * x + y * y + z * z); 

        if (r >= 2 * a) {
          double t0 = (q0 * x + q1 * y + q2 * z) / pow(r, 3); 
          double t1 = c1 / r; 
          double t2 = c2 / pow(r, 3); 
          double t3 = 3 * c2 / pow(r, 2); 

          v0 += t1 * q0 + x * t0 + t2 * q0 - t3 * t0 * x;
          v1 += t1 * q1 + y * t0 + t2 * q1 - t3 * t0 * y;
          v2 += t1 * q2 + z * t0 + t2 * q2 - t3 * t0 * z;
        } else {
          double c3 = 9.0 / 32.0 / a; 
          double c4 = 3.0 / 32.0 / a;
          double A1 = c0 * (1 - c3 * r); 
          double A2 = c4 * (q0 * x + q1 * y + q2 * z) / r;
          
          v0 += A1 * q0 + A2 * x; 
          v1 += A1 * q1 + A2 * y; 
          v2 += A1 * q2 + A2 * z;
        }
      }
      i->value[0] += v0; 
      i->value[1] += v1; 
      i->value[2] += v2; 
    }      
  }

  std::unique_ptr<expansion_t> M_to_I() const {
    dcomplex_t *M1 = reinterpret_cast<dcomplex_t *>(views_.view_data(0));
    dcomplex_t *M2 = reinterpret_cast<dcomplex_t *>(views_.view_data(1));
    dcomplex_t *M3 = reinterpret_cast<dcomplex_t *>(views_.view_data(2));
    dcomplex_t *M4 = reinterpret_cast<dcomplex_t *>(views_.view_data(3));

    expansion_t *ret{new expansion_t{kSourceIntermediate}};

    lap_m_to_i(M1, ret->views_, 0); 
    lap_m_to_i(M2, ret->views_, 6); 
    lap_m_to_i(M3, ret->views_, 12); 
    lap_m_to_i(M4, ret->views_, 18); 

    return std::unique_ptr<expansion_t>(ret);
  }

  std::unique_ptr<expansion_t> I_to_I(Index s_index, Index t_index) const {
    ViewSet views{kTargetIntermediate};

    lap_i_to_i(s_index, t_index, views_, 0, 0, views); 
    lap_i_to_i(s_index, t_index, views_, 6, 28, views); 
    lap_i_to_i(s_index, t_index, views_, 12, 56, views); 
    lap_i_to_i(s_index, t_index, views_, 18, 84, views); 

    expansion_t *ret = new expansion_t{views};
    return std::unique_ptr<expansion_t>{ret};
  }
  
  std::unique_ptr<expansion_t> I_to_L(Index t_index) const {
    // t_index is the index of the child
    expansion_t *ret{new expansion_t{kTargetPrimary}};
    dcomplex_t *L1 = reinterpret_cast<dcomplex_t *>(ret->views_.view_data(0));
    dcomplex_t *L2 = reinterpret_cast<dcomplex_t *>(ret->views_.view_data(1));
    dcomplex_t *L3 = reinterpret_cast<dcomplex_t *>(ret->views_.view_data(2));
    dcomplex_t *L4 = reinterpret_cast<dcomplex_t *>(ret->views_.view_data(3));

    double scale = views_.scale() * 2;
    lap_i_to_l(views_, 0, t_index, scale, L1);
    lap_i_to_l(views_, 28, t_index, scale, L2); 
    lap_i_to_l(views_, 56, t_index, scale, L3); 
    lap_i_to_l(views_, 84, t_index, scale, L4); 

    return std::unique_ptr<expansion_t>(ret);
  }  

  void add_expansion(const expansion_t *temp1) {
    // This operation assumes that the views included in \p temp1 is a subset of
    // \p views_. No range checking performed.
    int count = temp1->views_.count();
    for (int i = 0; i < count; ++i) {
      int idx = temp1->views_.view_index(i);
      int size = temp1->views_.view_bytes(i) / sizeof(dcomplex_t);
      dcomplex_t *lhs = reinterpret_cast<dcomplex_t *>(views_.view_data(idx));
      dcomplex_t *rhs =
        reinterpret_cast<dcomplex_t *>(temp1->views_.view_data(i));
      
      for (int j = 0; j < size; ++j) {
        lhs[j] += rhs[j];
      }
    }
  }
  
  static void update_table(int n_digits, double domain_size, 
                           const std::vector<double> &kernel_params) {
    update_laplace_table(n_digits, domain_size); 
    update_rpy_table(kernel_params[0], kernel_params[1], 
                     kernel_params[2], kernel_params[3]); 
  }

  static void delete_table() { } 

  static double compute_scale(Index index) {
    return builtin_laplace_table_->scale(index.level()); 
  }

  static int weight_estimate(Operation op,
                             Index s = Index{}, Index t = Index{}) {
    int weight = 0;
    if (op == Operation::MtoI) {
      weight = 6;
    } else if (op == Operation::ItoI) {
      int dx = s.x() - 2 * t.x();
      int dy = s.y() - 2 * t.y();
      int dz = s.z() - 2 * t.z();
      for (int i = 0; i < 3; ++i) {
        int tag = merge_and_shift_table[dx + 2][dy + 2][dz + 2][i];
        if (tag == -1) {
          break;
        }
        weight++;
      }
    } else {
      weight = 1;
    }
    return weight;
  }  

private: 
  ViewSet views_; 
};    
       

} // namespace dashmm

#endif  // __RPY_EXPANSION_H__
