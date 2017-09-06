#include "rpy.h"

namespace dashmm {

std::unique_ptr<RPYTable> rpy_table_; 

void update_rpy_table(double a, double kb, double T, double eta) {
  if (rpy_table_ == nullptr) {
    // Create the table if it does not exist. 
    rpy_table_ = 
      std::unique_ptr<RPYTable>{new RPYTable{a, kb, T, eta}}; 
  } else if (rpy_table_->update(a, kb, T, eta)) {
    // Replace the existing table if it requires update
    rpy_table_.reset(new RPYTable{a, kb, T, eta});
  }
}

// Compute a multipole expansion
double eval_M(int p, const dcomplex_t *M, const double *sqf, 
              const double *powers_r, const double *legendre,
              const dcomplex_t *powers_ephi) {
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
                   const double *powers_r, const double *legendre,
                   const dcomplex_t *powers_ephi) {
  dcomplex_t s1{0.0, 0.0}, s2{0.0, 0.0}; 
  
  for (int n = 0; n <= p; ++n) 
    s1 += M[midx(n, 0)] * ((double) (n + 1)) * legendre[midx(n + 1, 0)] * 
      powers_r[n + 1]; 
  
  for (int n = 1; n <= p; ++n) {
    for (int m = 1; m <= n; ++m) {
      s2 += M[midx(n, m)] * powers_r[n + 1] * legendre[midx(n + 1, m)] * 
        sqf[n - m + 1] / sqf[n - m] * sqf[n - m + 1] / sqf[n + m] * 
        powers_ephi[m]; 
    }
  }
  
  return -real(s1 + 2.0 * s2); 
}

// Compute d/dx + i * d/dy of a multipole expansion
dcomplex_t gradientp_M(int p, const dcomplex_t *M, const double *sqf, 
                       const double *powers_r, const double *legendre,
                       const dcomplex_t *powers_ephi) {
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
                        const double *powers_r, const double *legendre,
                        const dcomplex_t *powers_ephi) {
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
                        const double *powers_r, const double *legendre,
                        const dcomplex_t *powers_ephi) {
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
    s2 -= M[midx(n, 1)] * powers_r[n + 2] * legendre[midx(n + 2, 1)] * 
      sqf[n + 3] / sqf[n - 1] * sqf[n + 3] / sqf[n + 1] * 
      conj(powers_ephi[1]); 
  }
  
  return s1 + conj(s2); 
}

// Compute (d/dx - i * d/dy) * (d/dx + i * d/dy) of a local expansion
double gradientmp_M(int p, const dcomplex_t *M, const double *sqf, 
                    const double *powers_r, const double *legendre,
                    const dcomplex_t *powers_ephi) {
  dcomplex_t s1{0.0, 0.0}, s2{0.0, 0.0}; 
  
  for (int n = 0; n <= p; ++n) {
    s1 += M[midx(n, 0)] * powers_r[n + 2] * legendre[midx(n + 2, 0)] * 
      ((double) (n + 1) * (n + 2)); 
    
    for (int m = 1; m <= n; ++m) {
      s2 += M[midx(n, m)] * powers_r[n + 2] * legendre[midx(n + 2, m)] * 
        sqf[n - m + 2] / sqf[n - m] * sqf[n - m + 2] / sqf[n + m] * 
        powers_ephi[m]; 
    }
  }
  
  return -real(s1 + 2.0 * s2); 
}

// Compute d^2/dz^2 of a multipole expansion
double gradient00_M(int p, const dcomplex_t *M, const double *sqf, 
                    const double *powers_r, const double *legendre,
                    const dcomplex_t *powers_ephi) {
  dcomplex_t s1{0.0, 0.0}, s2{0.0, 0.0}; 
  
  for (int n = 0; n <= p; ++n) {
    s1 += M[midx(n, 0)] * powers_r[n + 2] * legendre[midx(n + 2, 0)] * 
      ((double) (n + 1) * (n + 2));
    
    for (int m = 1; m <= n; ++m) {
      s2 += M[midx(n, m)] * powers_r[n + 2] * legendre[midx(n + 2, m)] * 
        sqf[n - m + 2] / sqf[n - m] * sqf[n - m + 2] / sqf[n + m] * 
        powers_ephi[m];
    }
  }
  
  return real(s1 + 2.0 * s2); 
}

std::vector<double> rpy_m_to_t(Point t, double scale, 
                               const dcomplex_t *M1, const dcomplex_t *M2, 
                               const dcomplex_t *M3, const dcomplex_t *M4) {
  std::vector<double> retval(3, 0.0); 
  
  double c1 = rpy_table_->c1(); 
  double c2 = rpy_table_->c2(); 
  int p = builtin_laplace_table_->p();
  const double *sqf = builtin_laplace_table_->sqf();
  
  double *legendre = new double[(p + 3) * (p + 4) / 2]; 
  double *powers_r = new double[p + 3]; 
  dcomplex_t *powers_ephi = new dcomplex_t[p + 3]; 
  
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
  
  // Compute powers of 1 / r
  powers_r[0] = 1.0 / r; 
  r *= scale; 
  for (int j = 1; j <= p + 2; ++j) 
    powers_r[j] = powers_r[j - 1] / r; 
  
  // Compute powers of exp(i * phi)
  powers_ephi[0] = dcomplex_t{1.0, 0.0}; 
  for (int j = 1; j <= p + 2; ++j) 
    powers_ephi[j] = powers_ephi[j - 1] * ephi; 
  
  // Compute value of the Legendre polynomial
  legendre_Plm(p + 2, ctheta, legendre); 
  
  // Process M1
  {
    auto s1 = eval_M(p, M1, sqf, powers_r, legendre, powers_ephi); 
    retval[0] += c1 * s1; 
    
    auto s2 = gradient0_M(p, M1, sqf, powers_r, legendre, powers_ephi);
    retval[2] -= c1 * x * s2 / scale; 
    
    auto s3 = gradientp_M(p, M1, sqf, powers_r, legendre, powers_ephi);
    retval[0] -= c1 * x * real(s3) / scale; 
    retval[1] -= c1 * x * imag(s3) / scale; 
    
    auto s4 = gradientp0_M(p, M1, sqf, powers_r, legendre, powers_ephi); 
    retval[2] += c2 * real(s4) / scale2; 
    
    auto s5 = gradientpp_M(p, M1, sqf, powers_r, legendre, powers_ephi); 
    auto s6 = gradientmp_M(p, M1, sqf, powers_r, legendre, powers_ephi); 
    
    retval[0] += c2 * real(s5 + s6) / 2 / scale2; 
    retval[1] += c2 * imag(s5) / 2 / scale2; 
  }

  // Process M2
  {
    auto s1 = eval_M(p, M2, sqf, powers_r, legendre, powers_ephi); 
    retval[1] += c1 * s1; 
    
    auto s2 = gradient0_M(p, M2, sqf, powers_r, legendre, powers_ephi);
    retval[2] -= c1 * y * s2 / scale; 
    
    auto s3 = gradientp_M(p, M2, sqf, powers_r, legendre, powers_ephi); 
    retval[0] -= c1 * y * real(s3) / scale; 
    retval[1] -= c1 * y * imag(s3) / scale;
    
    auto s4 = gradientp0_M(p, M2, sqf, powers_r, legendre, powers_ephi); 
    retval[2] += c2 * imag(s4) / scale2; 
    
    auto s5 = gradientpp_M(p, M2, sqf, powers_r, legendre, powers_ephi); 
    auto s6 = gradientmp_M(p, M2, sqf, powers_r, legendre, powers_ephi); 
    
    retval[0] += c2 * imag(s5) / 2 / scale2; 
    retval[1] += c2 * real(s6 - s5) / 2 / scale2; 
  }
  
  // Process M3 
  {
    auto s1 = eval_M(p, M3, sqf, powers_r, legendre, powers_ephi); 
    retval[2] += c1 * s1; 
    
    auto s2 = gradient0_M(p, M3, sqf, powers_r, legendre, powers_ephi); 
    retval[2] -= c1 * z * s2 / scale; 
    
    auto s3 = gradientp_M(p, M3, sqf, powers_r, legendre, powers_ephi); 
    retval[0] -= c1 * z * real(s3) / scale; 
    retval[1] -= c1 * z * imag(s3) / scale; 
    
    auto s4 = gradientp0_M(p, M3, sqf, powers_r, legendre, powers_ephi);
    retval[0] += c2 * real(s4) / scale2;
    retval[1] += c2 * imag(s4) / scale2; 
    
    auto s5 = gradient00_M(p, M3, sqf, powers_r, legendre, powers_ephi); 
    retval[2] += c2 * s5 / scale2;
  }

  // Process L4 
  {
    auto s1 = gradient0_M(p, M4, sqf, powers_r, legendre, powers_ephi); 
    retval[2] += s1 / scale; 
    
    auto s2 = gradientp_M(p, M4, sqf, powers_r, legendre, powers_ephi); 
    retval[0] += real(s2) / scale; 
    retval[1] += imag(s2) / scale;
  }
  
  delete [] powers_r; 
  delete [] legendre;
  delete [] powers_ephi; 

  return retval; 
}


// Compute a local expansion
double eval_L(int p, const dcomplex_t *L, const double *sqf, 
              const double *powers_r, const double *legendre,
              const dcomplex_t *powers_ephi) {
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
                   const double *powers_r, const double *legendre,
                   const dcomplex_t *powers_ephi) {
  dcomplex_t s1{0.0, 0.0}, s2{0.0, 0.0}; 
  
  for (int n = 1; n <= p; ++n) {
    s1 += L[midx(n, 0)] * ((double) n) * powers_r[n - 1] * 
      legendre[midx(n - 1, 0)]; 
    
    for (int m = 1; m <= n - 1; ++m) {
      s2 += L[midx(n, m)] * powers_r[n - 1] * legendre[midx(n - 1, m)] * 
        sqf[n + m] / sqf[n + m - 1] * sqf[n - m] / sqf[n + m - 1] * 
        powers_ephi[m]; 
    }
  }
  
  return real(s1 + 2.0 * s2);
}

// Compute d/dx + i * d/dy of a local expansion
dcomplex_t gradientp_L(int p, const dcomplex_t *L, const double *sqf,
                       const double *powers_r, const double *legendre,
                       const dcomplex_t *powers_ephi) {
  dcomplex_t s1{0.0, 0.0}, s2{0.0, 0.0}; 
  
  for (int n = 2; n <= p; ++n) {
    for (int m = 0; m <= n - 2; ++m) {
      s1 += L[midx(n, m)] * powers_r[n - 1] * legendre[midx(n - 1, m + 1)] *
        sqf[n - m] / sqf[n + m] * powers_ephi[m + 1]; 
    }
  }
  
  for (int n = 1; n <= p; ++n) {
    for (int m = 1; m <= n; ++m) {
      s2 += L[midx(n, m)] * powers_r[n - 1] * legendre[midx(n - 1, m - 1)] *
        sqf[n + m] / sqf[n + m - 2] * sqf[n - m] / sqf[n + m - 2] * 
        powers_ephi[m - 1]; 
    }
  }
  
  return s1 - conj(s2); 
}


// Compute (d/dx + i * d/dy) * d/dz of a local expansion
dcomplex_t gradientp0_L(int p, const dcomplex_t *L, const double *sqf,
                        const double *powers_r, const double *legendre,
                        const dcomplex_t *powers_ephi) {
  dcomplex_t s1{0.0, 0.0}, s2{0.0, 0.0}; 
  
  for (int n = 3; n <= p; ++n) {
    for (int m = 0; m <= n - 3; ++m) {
      s1 += L[midx(n, m)] * powers_r[n - 2] * legendre[midx(n - 2, m + 1)] * 
        sqf[n - m] / sqf[n + m - 1] * sqf[n + m] / sqf[n + m - 1] * 
        powers_ephi[m + 1]; 
    }
  }
  
  for (int n = 2; n <= p; ++n) {
    for (int m = 1; m <= n; ++m) {
      s2 += L[midx(n, m)] * powers_r[n - 2] * legendre[midx(n - 2, m - 1)] * 
        sqf[n - m] / sqf[n + m - 3] * sqf[n + m] / sqf[n + m - 3] * 
        powers_ephi[m + 1]; 
    }
  }
  
  return s1 - conj(s2); 
}

// Compute (d/dx + i * d/dy)^2 of a local expansion
dcomplex_t gradientpp_L(int p, const dcomplex_t *L, const double *sqf, 
                        const double *powers_r, const double *legendre,
                        const dcomplex_t *powers_ephi) {
  dcomplex_t s1{0.0, 0.0}, s2{0.0, 0.0}; 
  
  for (int n = 4; n <= p; ++n) {
    for (int m = 0; m <= n - 4; ++m) {
      s1 += L[midx(n, m)] * powers_r[n - 2] * legendre[midx(n - 2, m + 1)] *
        sqf[n - m] / sqf[n + m] * powers_ephi[m + 2]; 
    }
  }
  
  for (int n = 2; n <= p; ++n) {
    for (int m = 2; m <= n; ++m) {
      s2 += L[midx(n, m)] * powers_r[n - 2] * legendre[midx(n - 2, m - 2)] * 
        sqf[n + m] / sqf[n + m - 4] * sqf[n - m] / sqf[n + m - 4] *
        powers_ephi[m - 2]; 
    }
  }
  
  return s1 + conj(s2); 
}

// Compute (d/dx - i * d/dy) * (d/dx + i * d/dy) of a local expansion
double gradientmp_L(int p, const dcomplex_t *L, const double *sqf, 
                    const double *powers_r, const double *legendre, 
                    const dcomplex_t *powers_ephi) {
  dcomplex_t s1{0.0, 0.0}, s2{0.0, 0.0}; 
  
  for (int n = 2; n <= p; ++n) {
    s1 -= L[midx(n, 0)] * ((double) n * (n - 1)) * powers_r[n - 2] * 
      legendre[midx(n - 2, 0)]; 
    
    for (int m = 1; m <= n; ++m) {
      s2 -= L[midx(n, m)] * powers_r[n - 2] * legendre[midx(n - 2, m)] * 
        sqf[n + m] / sqf[n + m - 2] * sqf[n - m] / sqf[n + m - 2] *
        powers_ephi[m]; 
    }
  }
  
  return real(s1 + 2.0 * s2);
}

// Compute d^2/dz^2 of a local expansion
double gradient00_L(int p, const dcomplex_t *L, const double *sqf, 
                    const double *powers_r, const double *legendre, 
                    const dcomplex_t *powers_ephi) {
  dcomplex_t s1{0.0, 0.0}, s2{0.0, 0.0}; 
  
  for (int n = 2; n <= p; ++n) {
    s1 += L[midx(n, 0)] * ((double) n * (n - 1)) * legendre[midx(n - 2, 0)] * 
        powers_r[n - 2]; 
    
    for (int m = 1; m <= n - 2; ++m) {
      s2 += L[midx(n, m)] * powers_r[n - 2] * legendre[midx(n - 2, m)] * 
        sqf[n + m] / sqf[n + m - 2] * sqf[n - m] / sqf[n + m - 2]; 
    }
  }
  
  return real(s1 + 2.0 * s2); 
}

std::vector<double> rpy_l_to_t(Point t, double scale, 
                               const dcomplex_t *L1, const dcomplex_t *L2,
                               const dcomplex_t *L3, const dcomplex_t *L4) {
  std::vector<double> retval(3, 0.0); 
  
  double c1 = rpy_table_->c1(); 
  double c2 = rpy_table_->c2(); 
  int p = builtin_laplace_table_->p();
  const double *sqf = builtin_laplace_table_->sqf();
  
  double *legendre = new double[(p + 2) * (p + 3) / 2]; 
  double *powers_r = new double[p + 2]; 
  dcomplex_t *powers_ephi = new dcomplex_t[p + 2]; 
  
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
  legendre_Plm(p + 1, ctheta, legendre); 
  
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
    retval[2] += c2 * imag(s4) * scale2; 
    
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

  delete [] powers_r; 
  delete [] legendre; 
  delete [] powers_ephi;
  
  return retval; 
}
  
} // namespace dashmm
