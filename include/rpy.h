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
  RPYTable(double a, double kb = 1.0, double T = 1.0, double eta = 1.0) 
    : a_{a}, kb_{kb}, T_{T}, eta_{eta} {
      c0_ = kb * T / 6 / M_PI / a / eta; 
      c1_ = kb * T / 8 / M_PI / eta; 
      c2_ = kb * T  * a * a / 12 / M_PI;
    }
  
  double c0() const {return c0_;}
  double c1() const {return c1_;}
  double c2() const {return c2_;}

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

void update_rpy_table(double a, double kb, double T, double eta); 

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
    double scale = views_.center(); 
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
    update_rpy_table(kernel_params[0], kernel_parms[1], 
                     kernel_params[2], kernel_parms[3]); 
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

  std::vector<double> rpy_m_to_t(Point t, double scale, 
                                 const dcomplex_t *M1, 
                                 const dcomplex_t *M2, 
                                 const dcomplex_t *M3, 
                                 const dcomplex_t *M4) const {
    std::vector<double> retval(3, 0.0); 

    double c1 = rpy_table_->c1(); 
    double c2 = rpy_table_->c2(); 
    int p = builtin_laplace_table_->p();
    const double *sqf = builtin_laplace_table_->sqf();

    std::vector<double> legendre((p + 2) * (p + 3) / 2); 
    std::vector<double> powers_r(p + 2); 
    std::vector<dcomplex_t> powers_ephi(p + 2); 
    
    powers_r[0] = 1.0; 
    powers_ephi[0] = dcomplex_t{1.0, 0.0}; 

    double x = t.x(); 
    double y = t.y(); 
    double z = t.z(); 
    double proj = sqrt(x * x + y * y); 
    double r = t.norm(); 
    double scale2 = scale * scale; 

    // Compute cosine of the polar angle theta
    double ctheta = (r <= 1e-14 ? 1.0 : z / r); 

    // Compute exp(i * phi) for the azimuthal angle phi
    dcomplex_t ephi = (proj / r <= 1e-14 ? dcomplex_t{1.0, 0.0} : 
                       dcomplex_t{x / proj, y / proj}); 

    // Compute powers of r 
    r *= scale; 
    for (int j = 1; j <= p + 1; ++j) 
      powers_r[j] = powers_r[j - 1] * r; 

    // Compute powers of exp(i * phi)
    for (int j = 1; j <= p + 1; ++j) 
      powers_ephi[j] = powers_ephi[j - 1] * ephi; 

    // Compute value of the Legendre polynomial
    legendre_Plm(p + 1, ctheta, legendre.data()); 

    // Process L1
    {
      auto s1 = eval_L(p, L1, sqf, powers_r, legendre, powers_ephi); 
      retval[0] += c1 * s1; 

      auto s2 = gradient0_L(p, L1, sqf, powers_r, legendre, powers_ephi);
      retval[2] -= c1 * x * s2 * scale; 

      auto s3 = gradientp_L(p, L1, sqf, powers_r, legendre, powers_ephi);
      retval[0] -= c1 * x * real(s3) * scale; 
      retval[1] -= c1 * x * imag(s3) * scale; 
      
      auto s4 = gradientp0_L(p, L1, sqf, powers_r, legendre, powers_ephi); 
      retval[2] += c2 * real(s4) * scale2; 

      auto s5 = gradientpp_L(p, L1, sqf, powers_r, legendre, powers_ephi); 
      auto s6 = gradientmp_L(p, L1, sqf, powers_r, legendre, powers_ephi); 

      retval[0] += c2 * real(s5 + s6) / 2 * scale2; 
      retval[1] += c2 * imag(s5) / 2 * scale2; 
    }

    // Process L2 
    {
      auto s1 = eval_L(p, L2, sqf, powers_r, legendre, powers_ephi); 
      retval[1] += c1 * s1; 

      auto s2 = gradient0_L(p, L2, sqf, powers_r, legendre, powers_ephi);
      retval[2] -= c1 * y * s2 * scale; 

      auto s3 = gradientp_L(p, L2, sqf, powers_r, legendre, powers_ephi); 
      retval[0] -= c1 * y * real(s3) * scale; 
      retval[1] -= c1 * y * imag(s3) * scale;

      auto s4 = gradientp0_L(p, L2, sqf, powers_r, legendre, powers_ephi); 
      retval[2] += c2 * imag(s4) * scale; 

      auto s5 = gradientpp_L(p, L2, sqf, powers_r, legendre, powers_ephi); 
      auto s6 = gradientmp_L(p, L2, sqf, powers_r, legendre, powers_ephi); 

      retval[0] += c2 * imag(s5) / 2 * scale2; 
      retval[1] += c2 * real(s6 - s5) / 2 * scale2; 
    }

    // Process L3 
    {
      auto s1 = eval_L(p, L3, sqf, powers_r, legendre, powers_ephi); 
      retval[2] += c1 * s1; 

      auto s2 = gradient0_L(p, L3, sqf, powers_r, legendre, powers_ephi); 
      retval[2] -= c1 * z * s2 * scale; 

      auto s3 = gradientp_L(p, L3, sqf, powers_r, legendre, powers_ephi); 
      retval[0] -= c1 * z * real(s3) * scale; 
      retval[1] -= c1 * z * imag(s3) * scale; 

      auto s4 = gradientp0_L(p, L3, sqf, powers_r, legendre, powers_ephi);
      retval[0] += c2 * real(s4) * scale2;
      retval[1] += c2 * imag(s4) * scale2; 

      auto s5 = gradient00_L(p, L3, sqf, powers_r, legendre, powers_ephi); 
      retval[2] += c2 * s5 * scale2;
    }

    // Process L4 
    {
      auto s1 = gradient0_L(p, L4, sqf, powers_r, legendre, powers_ephi); 
      retval[2] += s1 * scale; 

      auto s2 = gradientp_L(p, L4, sqf, powers_r, legendre, powers_ephi); 
      retval[0] += real(s2) * scale; 
      retval[1] += imag(s2) * scale;
    }
        
    return retval; 
  }



  std::vector<double> rpy_l_to_t(Point t, double scale, 
                                 const dcomplex_t *L1, 
                                 const dcomplex_t *L2, 
                                 const dcomplex_t *L3, 
                                 const dcomplex_t *L4) const {
    std::vector<double> retval(3, 0.0); 

    double c1 = rpy_table_->c1(); 
    double c2 = rpy_table_->c2(); 
    int p = builtin_laplace_table_->p();
    const double *sqf = builtin_laplace_table_->sqf();

    std::vector<double> legendre((p + 2) * (p + 3) / 2); 
    std::vector<double> powers_r(p + 2); 
    std::vector<dcomplex_t> powers_ephi(p + 2); 
    
    powers_r[0] = 1.0; 
    powers_ephi[0] = dcomplex_t{1.0, 0.0}; 

    double x = t.x(); 
    double y = t.y(); 
    double z = t.z(); 
    double proj = sqrt(x * x + y * y); 
    double r = t.norm(); 
    double scale2 = scale * scale; 

    // Compute cosine of the polar angle theta
    double ctheta = (r <= 1e-14 ? 1.0 : z / r); 

    // Compute exp(i * phi) for the azimuthal angle phi
    dcomplex_t ephi = (proj / r <= 1e-14 ? dcomplex_t{1.0, 0.0} : 
                       dcomplex_t{x / proj, y / proj}); 

    // Compute powers of r 
    r *= scale; 
    for (int j = 1; j <= p + 1; ++j) 
      powers_r[j] = powers_r[j - 1] * r; 

    // Compute powers of exp(i * phi)
    for (int j = 1; j <= p + 1; ++j) 
      powers_ephi[j] = powers_ephi[j - 1] * ephi; 

    // Compute value of the Legendre polynomial
    legendre_Plm(p + 1, ctheta, legendre.data()); 

    // Process L1
    {
      auto s1 = eval_L(p, L1, sqf, powers_r, legendre, powers_ephi); 
      retval[0] += c1 * s1; 

      auto s2 = gradient0_L(p, L1, sqf, powers_r, legendre, powers_ephi);
      retval[2] -= c1 * x * s2 * scale; 

      auto s3 = gradientp_L(p, L1, sqf, powers_r, legendre, powers_ephi);
      retval[0] -= c1 * x * real(s3) * scale; 
      retval[1] -= c1 * x * imag(s3) * scale; 
      
      auto s4 = gradientp0_L(p, L1, sqf, powers_r, legendre, powers_ephi); 
      retval[2] += c2 * real(s4) * scale2; 

      auto s5 = gradientpp_L(p, L1, sqf, powers_r, legendre, powers_ephi); 
      auto s6 = gradientmp_L(p, L1, sqf, powers_r, legendre, powers_ephi); 

      retval[0] += c2 * real(s5 + s6) / 2 * scale2; 
      retval[1] += c2 * imag(s5) / 2 * scale2; 
    }

    // Process L2 
    {
      auto s1 = eval_L(p, L2, sqf, powers_r, legendre, powers_ephi); 
      retval[1] += c1 * s1; 

      auto s2 = gradient0_L(p, L2, sqf, powers_r, legendre, powers_ephi);
      retval[2] -= c1 * y * s2 * scale; 

      auto s3 = gradientp_L(p, L2, sqf, powers_r, legendre, powers_ephi); 
      retval[0] -= c1 * y * real(s3) * scale; 
      retval[1] -= c1 * y * imag(s3) * scale;

      auto s4 = gradientp0_L(p, L2, sqf, powers_r, legendre, powers_ephi); 
      retval[2] += c2 * imag(s4) * scale; 

      auto s5 = gradientpp_L(p, L2, sqf, powers_r, legendre, powers_ephi); 
      auto s6 = gradientmp_L(p, L2, sqf, powers_r, legendre, powers_ephi); 

      retval[0] += c2 * imag(s5) / 2 * scale2; 
      retval[1] += c2 * real(s6 - s5) / 2 * scale2; 
    }

    // Process L3 
    {
      auto s1 = eval_L(p, L3, sqf, powers_r, legendre, powers_ephi); 
      retval[2] += c1 * s1; 

      auto s2 = gradient0_L(p, L3, sqf, powers_r, legendre, powers_ephi); 
      retval[2] -= c1 * z * s2 * scale; 

      auto s3 = gradientp_L(p, L3, sqf, powers_r, legendre, powers_ephi); 
      retval[0] -= c1 * z * real(s3) * scale; 
      retval[1] -= c1 * z * imag(s3) * scale; 

      auto s4 = gradientp0_L(p, L3, sqf, powers_r, legendre, powers_ephi);
      retval[0] += c2 * real(s4) * scale2;
      retval[1] += c2 * imag(s4) * scale2; 

      auto s5 = gradient00_L(p, L3, sqf, powers_r, legendre, powers_ephi); 
      retval[2] += c2 * s5 * scale2;
    }

    // Process L4 
    {
      auto s1 = gradient0_L(p, L4, sqf, powers_r, legendre, powers_ephi); 
      retval[2] += s1 * scale; 

      auto s2 = gradientp_L(p, L4, sqf, powers_r, legendre, powers_ephi); 
      retval[0] += real(s2) * scale; 
      retval[1] += imag(s2) * scale;
    }
        
    return retval; 
  }
  
  // Compute a multipole expansion
  double eval_M(int p, const dcomplex_t *M, const double *sqf, 
                const std::vector<double> &powers_r, 
                const std::vector<double> &legendre, 
                const std::vector<dcomplex_t> &powers_ephi) const {
    dcomplex_t s1{0.0, 0.0}; 

    for (int n = 0; n <= p; ++n) 
      s1 += M[midx(n, 0)] * powers_r[n] * legendre[midx(n, 0)]; 

    for (int n = 1; n <= p; ++n) {
      for (int m = 1; m <= n; ++m) {
        s1 += 2.0 * real(M[midx(n, m)]) * powers_ephi[m] * 
          powers_r[n] * legendre[midx(n, m)] * sqf[n - m] / sqf[n + m];
      }
    }

    return real(s1); 
  }

  // Compute d/dz of a multipole expansion
  double gradient0_M(int p, const dcomplex_t *M, const double *sqf, 
                     const std::vector<double> &powers_r, 
                     const std::vector<double> &legendre, 
                     const std::vector<dcomplex_t> &powers_ephi) const {
    dcomplex_t s1{0.0, 0.0}, s2{0.0, 0.0}; 

    for (int n = 0; n <= p; ++n) 
      s1 += M[midx(n, 0)] * (n + 1) * legendre[midx(n + 1, 0)] * 
        powers_r[n + 1]; 

    for (int n = 1; n <= p; ++n) {
      for (int m = 1; m <= n; ++m) {
        s2 += M[midx(n, m)] * powers_r[n + 1] * legendre[midx(n + 1, m)] * 
          sqf[n - m + 1] / sqf[n - m] * sqf[n - m + 1] / sqf[n + m] * 
          powers_ephi[m]; 
      }
    }
      
    return -real(s1 + 2 * s2); 
  }

  // Compute d/dx + i * d/dy of a multipole expansion
  dcomplex_t gradientp_M(int p, const dcomplex_t *M, const double *sqf, 
                         const std::vector<double> &powers_r, 
                         const std::vector<double> &legendre, 
                         const std::vector<dcomplex_t> &powers_ephi) const {
    dcomplex_t s1{0.0, 0.0}, s2{0.0, 0.0}; 

    for (int n = 0; n <= p; ++n) {
      for (int m = 0; m <= n; ++m) {
        s1 += M[midx(n, m)] * legendre[midx(n + 1, m + 1)] * powers_r[n + 1] *
          sqf[n - m] / sqf[n + m] * powers_ephi[m + 1]; 
      }

      for (int m = 1; m <= n; ++m) {
        s2 += M[midx(n, m)] * legendre[midx(n + 1, m - 1)] * powers_r[n + 1] *
          sqf[n - m + 2] / sqf[n - m] * sqf[n - m + 2] / sqf[n + m] * 
          powers_ephi[m - 1];
      }
    }

    return s1 - conj(s2);     
  }

  // Compute (d/dx + i * d/dy) * d/dz of a multipole expansion
  dcomplex_t gradientp0_M(int p, const dcomplex_t *M, const double *sqf,
                          const std::vector<double> &powers_r, 
                          const std::vector<double> &legendre, 
                          const std::vector<dcomplex_t> &powers_ephi) const {
    dcomplex_t s1{0.0, 0.0}, s2{0.0, 0.0}; 

    for (int n = 0; n <= p; ++n) {
      for (int m = 0; m <= n; ++m) {
        s1 += M[midx(n, m)] * powers_r[n + 2] * legendre[midx(n + 2, m + 1)] *
          sqf[n - m + 1] / sqf[n - m] * sqf[n - m + 1] / sqf[n + m] * 
          powers_ephi[m + 1]; 
      }

      for (int m = 1; m <= n; ++m) {
        s2 += M[midx(n, m)] * powers_r[n + 2] * legendre[midx(n + 2, m - 1)] * 
          sqf[n - m + 3] / sqf[n - m] * sqf[n - m + 3] / sqf[n + m] * 
          powers_ephi[m - 1];
      }
    }

    return -s1 + conj(s2); 
  }

  // Compute (d/dx + i * d/dy)^2 of a multipole expansion
  dcomplex_t gradientpp_M(int p, const dcomplex_t *M, const double *sqf, 
                          const std::vector<double> &powers_r, 
                          const std::vector<double> &legendre, 
                          const std::vector<dcomplex_t> &powers_ephi) const {
    dcomplex_t s1{0.0, 0.0}, s2{0.0, 0.0}; 

    for (int n = 0; n <= p; ++n) {
      for (int m = 0; m <= n; ++m) {
        s1 += M[midx(n, m)] * powers_r[n + 2] * legendre[midx(n + 2, m + 2)] * 
          sqf[n - m] / sqf[n + m] * powers_ephi[m + 2]; 
      }
    }

    for (int n = 2; n <= p; ++n) {
      for (int m = 2; m <= n; ++m) {
        s2 += M[midx(n ,m)] * powers_r[n + 2] * legendre[midx(n + 2, m - 2)] * 
          sqf[n - m + 4] / sqf[n - m] * sqf[n - m + 4] / sqf[n + m] *
          powers_ephi[m - 2]; 
      }
    }

    for (int n = 1; n <=p; ++n) {
      // m = 1
      


  }

  // Compute (d/dx - i * d/dy) * (d/dx + i * d/dy) of a local expansion
  double gradientmp_M(int p, const dcomplex_t *M, const double *sqf, 
                      const std::vector<double> &powers_r, 
                      const std::vector<double> &legendre, 
                      const std::vector<dcomplex_t> &powers_ephi) const {
    dcomplex_t s1{0.0, 0.0}, s2{0.0, 0.0}; 
    
    for (int n = 0; n <= p; ++n) {
      s1 += M[midx(n, 0)] * powers_r[n + 2] * legendre[midx(n + 2, 0)] * 
        (n + 1) * (n + 2); 

      for (int m = 1; m <= n; ++m) {
        s2 += M[midx(n, m)] * powers_r[n + 2] * legendre[midx(n + 2, m)] * 
          sqf[n - m + 2] / sqf[n - m] * sqf[n - m + 2] / sqf[n + m] * 
          powers_ephi[m]; 
      }
    }

    return -real(s1 + 2 * s2); 
  }

  // Compute d^2/dz^2 of a multipole expansion
  double gradient00_M(int p, const dcomplex_t *M, const double *sqf, 
                      const std::vector<double> &powers_r, 
                      const std::vector<double> &legendre, 
                      const std::vector<dcomplex_t> &powers_ephi) const {
    dcomplex_t s1{0.0, 0.0}, s2{0.0, 0.0}; 

    for (int n = 0; n <= p; ++n) {
      s1 += M[midx(n, 0)] * powers_r[n + 2] * legendre[midx(n + 2, 0)] * 
        (n + 1) * (n + 2);

      for (int m = 1; m <= n; ++m) {
        s2 += M[midx(n, m)] * powers_r[n + 2] * legendre[midx(n + 2, m)] * 
          sqf[n - m + 2] / sqf[n - m] * sqf[n - m + 2] / sqf[n + m] * 
          powers_ephi[m];
      }
    }

    return real(s1 + 2 * s2); 
  }

  // Compute a local expansion
  double eval_L(int p, const dcomplex_t *L, const double *sqf, 
                const std::vector<double> &powers_r, 
                const std::vector<double> &legendre, 
                const std::vector<dcomplex_t> &powers_ephi) const {
    dcomplex_t s1{0.0, 0.0}; 

    for (int n = 0; n <= p; ++n) {
      s1 += L[midx(n, 0)] * powers_r[n] * legendre[midx(n, 0)]; 
    }

    for (int n = 1; n <= p; ++n) {
      for (int m = 1; m <= n; ++m) {
        s1 += 2.0 * real(L[midx(n, m)]) * powers_ephi[m] * 
          powers_r[n] * legendre[midx(n, m)] * sqf[n - m] / sqf[n + m];
      }
    }

    return real(s1);
  }


  // Compute d/dz of a local expansion
  double gradient0_L(int p, const dcomplex_t *L, const double *sqf, 
                     const std::vector<double> &powers_r, 
                     const std::vector<double> &legendre, 
                     const std::vector<dcomplex_t> &powers_ephi) const {
    dcomplex_t s1{0.0, 0.0}, s2{0.0, 0.0}; 

    for (int n = 1; n <= p; ++n) {
      s1 += L[midx(n, 0)] * n * powers_r[n - 1] * legendre[midx(n - 1, 0)]; 

      for (int m = 1; m <= n - 1; ++m) {
        s2 += L[midx(n, m)] * powers_r[n - 1] * legendre[midx(n - 1, m)] * 
          sqf[n + m] / sqf[n + m - 1] * sqf[n - m] / sqf[n + m - 1] * 
          powers_ephi[m]; 

        // Simplified from 
        // sqrt((n + m) * (n - m)) * sqf[n - 1 - m] / sqf[n - 1 + m]
      }
    }

    return real(s1 + 2 * s2);
  }

  // Compute d/dx + i * d/dy of a local expansion
  dcomplex_t gradientp_L(int p, const dcomplex_t *L, const double *sqf,
                         const std::vector<double> &powers_r, 
                         const std::vector<double> &legendre,
                         const std::vector<dcomplex_t> &powers_ephi) const {
    dcomplex_t s1{0.0, 0.0}, s2{0.0, 0.0}; 

    for (int n = 2; n <= p; ++n) {
      for (int m = 0; m <= n - 2; ++m) {
        s1 += L[midx(n, m)] * powers_r[n - 1] * legendre[midx(n - 1, m + 1)] *
          sqf[n - m] / sqf[n + m] * powers_ephi[m + 1]; 

        // Simplified from 
        // sqrt((n - m) * (n - m - 1)) * 
        // sqf[n - 1 - (m + 1)] / sqf[n - 1 + m + 1]
      }
    }

    for (int n = 1; n <= p ++n) {
      for (int m = 1; m <= n; ++m) {
        s2 += L[midx(n, m)] * powers_r[n - 1] * legendre[midx(n - 1, m - 1)] *
          sqf[n + m] / sqf[n + m - 2] * sqf[n - m] / sqf[n + m - 2] * 
          powers_ephi[m - 1]; 

        // Simplified from 
        // sqrt((n + m) * (n + m - 1)) * sqf[n - 1 - (m - 1)] / 
        // sqf[n - 1 + m - 1]
      }
    }

    return s1 - conj(s2); 
  }


  // Compute (d/dx + i * d/dy) * d/dz of a local expansion
  dcomplex_t gradientp0_L(int p, const dcomplex_t *L, const double *sqf,
                          const std::vector<double> &powers_r, 
                          const std::vector<double> &legendre, 
                          const std::vector<dcomplex_t> &powers_ephi) const {
    dcomplex_t s1{0.0, 0.0}, s2{0.0, 0.0}; 

    for (int n = 3; n <= p; ++n) {
      for (int m = 0; m <= n - 3; ++m) {
        s1 += L[midx(n, m)] * powers_r[n - 2] * legendre[midx(n - 2, m + 1)] * 
          sqf[n - m] / sqf[n + m - 1] * sqf[n + m] / sqf[n + m - 1] * 
          powers_ephi[m + 1]; 

        // Simplified from 
        // sqrt((n + m) * (n - m) * (n - m - 1) * (n - m - 2)) * 
        // sqf[n - 2 - (m + 1)] / sqf[n - 2 + (m + 1)] 
      }
    }

    for (int n = 2; n <= p; ++n) {
      for (int m = 1; m <= n; ++m) {
        s2 += L[midx(n, m)] * powers_r[n - 2] * legendre[midx(n - 2, m - 1)] * 
          sqf[n - m] / sqf[n + m - 3] * sqf[n + m] / sqf[n + m - 3] * 
          powers_ephi[m + 1]; 

        // Simplified from
        // sqrt((n + m) * (n - m) * (n + m - 1) * (n + m - 2)) * 
        // sqf[n - 2 - (m - 1)] / sqf[n - 2 + (m - 1)]
      }
    }

    return s1 - conj(s2); 
  }

  // Compute (d/dx + i * d/dy)^2 of a local expansion
  dcomplex_t gradientpp_L(int p, const dcomplex_t *L, const double *sqf, 
                          const std::vector<double> &powers_r, 
                          const std::vector<double> &legendre, 
                          const std::vector<dcomplex_t> &powers_ephi) const {
    dcomplex_t s1{0.0, 0.0}, s2{0.0, 0.0}; 

    for (int n = 4; n <= p; ++n) {
      for (int m = 0; m <= n - 4; ++m) {
        s1 += L[midx(n, m)] * powers_r[n - 2] * legendre[midx(n - 2, m + 1)] 
          * sqf[n - m] / sqf[n + m] * powers_ephi[m + 2]; 

        // Simplified from 
        // sqrt((n - m)) * (n - m - 1) * (n - m - 2) * (n - m - 3)) * 
        // sqf[n - 2 - (m + 2)] / sqf[n - 2 + (m + 2)]
      }
    }

    for (int n = 2; n <= p; ++n) {
      for (int m = 2; m <= n; ++m) {
        s2 += L[midx(n, m)] * powers_r[n - 2] * legendre[midx(n - 2, m - 2)] * 
          * sqf[n + m] / sqf[n + m - 4] * sqf[n - m] / sqf[n + m - 4] 
          * powers_ephi[m - 2]; 

        // Simplified from 
        // sqrt((n + m) * (n + m - 1) * (n + m - 2) * (n + m - 3)) * 
        // sqf[n - 2 - (m - 2)] / sqf[n - 2 + (m - 2)]
      }
    }

    return s1 + conj(s2); 
  }

  // Compute (d/dx - i * d/dy) * (d/dx + i * d/dy) of a local expansion
  double gradientmp_L(int p, const dcomplex_t *L, const double *sqf, 
                      const std::vector<double> &powers_r, 
                      const std::vector<double> &legendre, 
                      const std::vector<dcomplex_t> &powers_ephi) const {
    dcomplex_t s1{0.0, 0.0}, s2{0.0, 0.0}; 
    
    for (int n = 2; n <= p; ++n) {
      s1 -= L[midx(n, 0)] * n * (n - 1) * powers_r[n - 2] * 
        legendre[midx(n - 2, 0)]; 

      for (int m = 1; m <= n; ++m) {
        s2 -= L[midx(n, m)] * powers_r[n - 2] * legendre[midx(n - 2, m)] * 
          sqf[n + m] / sqf[n + m - 2] * sqf[n - m] / sqf[n + m - 2] *
          powers_ephi[m]; 

        // Simplified from 
        // sqrt((n - m) * (n - m - 1) * (n + m) * (n + m - 1)) * 
        // sqf[n - 2 - m] / sqf[n - 2 + m]
      }
    }

    return real(s3) + 2 * real(s4);
  }

  // Compute d^2/dz^2 of a local expansion
  double gradient00_L(int p, const dcomplex_t *L, const double *sqf, 
                      const std::vector<double> &powers_r, 
                      const std::vector<double> &legendre, 
                      const std::vector<dcomplex_t> &powers_ephi) const {
    dcomplex_t s1{0.0, 0.0}, s2{0.0, 0.0}; 

    for (int n = 2; n <= p; ++n) {
      s1 += L3[midx(n, 0)] * n * (n - 1) * legendre[midx(n - 2, 0)] * 
        powers_r[n - 2]; 

      for (int m = 1; m <= n - 2; ++m) {
        s2 += L[midx(n, m)] * powers_r[n - 2] * legendre[midx(n - 2, m)] * 
          sqf[n + m] / sqf[n + m - 2] * sqf[n - m] / sqf[n + m - 2]; 

        // Simplified from 
        // sqrt((n + m) * (n - m) * (n - 1 + m) * (n - 1 - m)) * 
        // sqf[n - 2 - m] / sqf[n - 2 + m]
      }
    }

    return real(s1) + 2 * real(s2); 
  }

};    
       

} // namespace dashmm

#endif  // __RPY_EXPANSION_H__
