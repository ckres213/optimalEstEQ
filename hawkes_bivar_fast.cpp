// hawkes_bivar_fast.cpp
// Fast kernels for a focused bivariate compact-support Hawkes simulation.
//
// Model:
//   lambda_i(t) = mu_i + sum_{j=1}^2 alpha_{ij} X_j(t; beta),
//   X_j(t; beta) = sum_{s in N_j, s<t, t-s<=A} k_beta(t-s),
//   k_beta(u) = beta exp(-beta u) / (1 - exp(-beta A)), 0 <= u <= A.
//
// Parameter order:
//   (mu1, mu2, alpha11, alpha12, alpha21, alpha22, beta)
// where row i receives excitation from source column j.
//
// Implements cluster simulation, exact conditional log-likelihood/score,
// exact D_theta non-score martingale moment, derivative-enriched
// overidentified bounded-inverse GMM moments, and time-average population
// matrices for the Godambe comparison.

// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>
using namespace Rcpp;

struct Event {
  double t;
  int type; // 1 or 2
};

static inline double den_texp(double beta, double A) {
  return 1.0 - std::exp(-beta * A);
}

static inline double k_texp(double u, double beta, double A) {
  if (!(u >= 0.0 && u <= A)) return 0.0;
  const double den = den_texp(beta, A);
  return beta * std::exp(-beta * u) / den;
}

static inline double dk_dbeta(double u, double beta, double A) {
  if (!(u >= 0.0 && u <= A)) return 0.0;
  const double den = den_texp(beta, A);
  const double eA = std::exp(-beta * A);
  const double kval = beta * std::exp(-beta * u) / den;
  const double q = 1.0 / beta - u - A * eA / den;
  return kval * q;
}

static inline double K_texp(double lo, double hi, double beta, double A) {
  lo = std::max(0.0, std::min(A, lo));
  hi = std::max(0.0, std::min(A, hi));
  if (!(hi > lo)) return 0.0;
  const double den = den_texp(beta, A);
  return (std::exp(-beta * lo) - std::exp(-beta * hi)) / den;
}

static inline double dK_dbeta(double lo, double hi, double beta, double A) {
  lo = std::max(0.0, std::min(A, lo));
  hi = std::max(0.0, std::min(A, hi));
  if (!(hi > lo)) return 0.0;
  const double eA = std::exp(-beta * A);
  const double den = 1.0 - eA;
  const double num = std::exp(-beta * lo) - std::exp(-beta * hi);
  const double dnum = -lo * std::exp(-beta * lo) + hi * std::exp(-beta * hi);
  const double dden = A * eA;
  return (dnum * den - num * dden) / (den * den);
}

static inline double spectral_radius_2x2(double a11, double a12, double a21, double a22) {
  const double tr = a11 + a22;
  const double disc = (a11 - a22) * (a11 - a22) + 4.0 * a12 * a21;
  return 0.5 * (tr + std::sqrt(std::max(0.0, disc)));
}

static inline bool parse_theta(const NumericVector& theta,
                               double& mu1, double& mu2,
                               double& a11, double& a12,
                               double& a21, double& a22,
                               double& beta,
                               double rho_max) {
  if (theta.size() != 7) return false;
  mu1 = theta[0]; mu2 = theta[1];
  a11 = theta[2]; a12 = theta[3]; a21 = theta[4]; a22 = theta[5];
  beta = theta[6];
  if (!R_finite(mu1) || !R_finite(mu2) || !R_finite(a11) || !R_finite(a12) ||
      !R_finite(a21) || !R_finite(a22) || !R_finite(beta)) return false;
  if (mu1 <= 0.0 || mu2 <= 0.0 || a11 <= 0.0 || a12 <= 0.0 ||
      a21 <= 0.0 || a22 <= 0.0 || beta <= 0.0) return false;
  if (!(spectral_radius_2x2(a11, a12, a21, a22) < rho_max)) return false;
  return true;
}

static inline void fill_events(const NumericVector& time, const IntegerVector& type,
                               std::vector<Event>& events,
                               std::vector<double>& t1,
                               std::vector<double>& t2) {
  events.clear(); t1.clear(); t2.clear();
  const int n = std::min(time.size(), type.size());
  events.reserve(n); t1.reserve(n / 2 + 1); t2.reserve(n / 2 + 1);
  for (int i = 0; i < n; ++i) {
    if (!R_finite(time[i])) continue;
    const int tp = type[i];
    if (tp != 1 && tp != 2) continue;
    events.push_back(Event{time[i], tp});
    if (tp == 1) t1.push_back(time[i]); else t2.push_back(time[i]);
  }
  std::sort(events.begin(), events.end(), [](const Event& a, const Event& b) {
    if (a.t == b.t) return a.type < b.type;
    return a.t < b.t;
  });
  std::sort(t1.begin(), t1.end());
  std::sort(t2.begin(), t2.end());
}

static inline void x_dx_at(double t, const std::vector<double>& times,
                           double beta, double A,
                           double& X, double& dX) {
  X = 0.0; dX = 0.0;
  const double lo = t - A;
  auto it0 = std::lower_bound(times.begin(), times.end(), lo);
  auto it1 = std::lower_bound(times.begin(), times.end(), t); // strictly before t
  for (auto it = it0; it != it1; ++it) {
    const double u = t - *it;
    X += k_texp(u, beta, A);
    dX += dk_dbeta(u, beta, A);
  }
}

static inline double sum_K_source(const std::vector<double>& times,
                                  double beta, double A, double T) {
  double ans = 0.0;
  for (double s : times) {
    if (!(s < T && s + A > 0.0)) continue;
    const double lo = std::max(0.0, -s);
    const double hi = std::min(A, T - s);
    ans += K_texp(lo, hi, beta, A);
  }
  return ans;
}

static inline double sum_dK_source(const std::vector<double>& times,
                                   double beta, double A, double T) {
  double ans = 0.0;
  for (double s : times) {
    if (!(s < T && s + A > 0.0)) continue;
    const double lo = std::max(0.0, -s);
    const double hi = std::min(A, T - s);
    ans += dK_dbeta(lo, hi, beta, A);
  }
  return ans;
}

static inline double prod_k_k_one(double s, double r, double beta1, double beta2,
                                  double A, double T) {
  const double t0 = std::max(0.0, std::max(s, r));
  const double t1 = std::min(T, std::min(s + A, r + A));
  if (!(t1 > t0)) return 0.0;

  const double B = beta1 + beta2;
  const double L = t1 - t0;
  const double den1 = den_texp(beta1, A);
  const double den2 = den_texp(beta2, A);

  // Stable form:
  // exp(beta1*s + beta2*r) * exp(-B*t)
  // = exp(-beta1*(t-s) - beta2*(t-r)).
  // Since t >= max(s,r), the exponent is non-positive and cannot overflow.
  const double lag1 = t0 - s;
  const double lag2 = t0 - r;
  const double c = beta1 * beta2 *
                   std::exp(-beta1 * lag1 - beta2 * lag2) /
                   (den1 * den2);

  return c * (-std::expm1(-B * L)) / B;
}

static inline double prod_dk_k_one(double s, double r, double beta1, double beta2,
                                   double A, double T) {
  const double t0 = std::max(0.0, std::max(s, r));
  const double t1 = std::min(T, std::min(s + A, r + A));
  if (!(t1 > t0)) return 0.0;

  const double B = beta1 + beta2;
  const double L = t1 - t0;
  const double den1 = den_texp(beta1, A);
  const double den2 = den_texp(beta2, A);

  const double lag1 = t0 - s;
  const double lag2 = t0 - r;

  const double c = beta1 * beta2 *
                   std::exp(-beta1 * lag1 - beta2 * lag2) /
                   (den1 * den2);

  // dk/dbeta = k(u) * [1/beta - u - A exp(-beta A)/(1-exp(-beta A))],
  // with u = t - s.  Write t = t0 + w and integrate over w in [0,L].
  const double eA = std::exp(-beta1 * A);
  const double C = 1.0 / beta1 - A * eA / den1;
  const double q0 = C - lag1;

  const double expBL = std::exp(-B * L);
  const double int_exp = (-std::expm1(-B * L)) / B;
  const double int_w_exp = (1.0 - expBL * (1.0 + B * L)) / (B * B);

  return c * (q0 * int_exp - int_w_exp);
}

static inline double product_sources(const std::vector<double>& t1s,
                                     const std::vector<double>& t2s,
                                     double beta1, double beta2,
                                     double A, double T,
                                     bool diff_first) {
  double ans = 0.0;
  for (double s : t1s) {
    if (!(s < T && s + A > 0.0)) continue;
    auto it0 = std::lower_bound(t2s.begin(), t2s.end(), std::max(0.0, s) - A);
    auto it1 = std::lower_bound(t2s.begin(), t2s.end(), std::min(T, s + A));
    for (auto it = it0; it != it1; ++it) {
      const double r = *it;
      if (!(r < T && r + A > 0.0)) continue;
      ans += diff_first ? prod_dk_k_one(s, r, beta1, beta2, A, T)
                        : prod_k_k_one(s, r, beta1, beta2, A, T);
    }
  }
  return ans;
}

static inline void deriv_rows(double X1, double X2, double dX1, double dX2,
                              double a11, double a12, double a21, double a22,
                              double* D1, double* D2) {
  for (int k = 0; k < 7; ++k) { D1[k] = 0.0; D2[k] = 0.0; }
  D1[0] = 1.0;
  D1[2] = X1;
  D1[3] = X2;
  D1[6] = a11 * dX1 + a12 * dX2;
  D2[1] = 1.0;
  D2[4] = X1;
  D2[5] = X2;
  D2[6] = a21 * dX1 + a22 * dX2;
}


static inline double aug_phi(double lam, double stabilizer) {
  return stabilizer / (stabilizer + lam);
}

static inline void aug_factors(double lam, double stabilizer, double* fac) {
  const double ph = aug_phi(lam, stabilizer);
  fac[0] = 1.0;
  fac[1] = ph;
  fac[2] = ph * ph;
}

static inline void ls_comp_exact(const std::vector<double>& t1,
                                 const std::vector<double>& t2,
                                 double mu1, double mu2,
                                 double a11, double a12,
                                 double a21, double a22,
                                 double beta, double A, double T,
                                 NumericVector& comp) {
  const double S1 = sum_K_source(t1, beta, A, T);
  const double S2 = sum_K_source(t2, beta, A, T);
  const double D1int = sum_dK_source(t1, beta, A, T);
  const double D2int = sum_dK_source(t2, beta, A, T);
  const double P11 = product_sources(t1, t1, beta, beta, A, T, false);
  const double P12 = product_sources(t1, t2, beta, beta, A, T, false);
  const double P22 = product_sources(t2, t2, beta, beta, A, T, false);
  const double Q11 = product_sources(t1, t1, beta, beta, A, T, true);
  const double Q12 = product_sources(t1, t2, beta, beta, A, T, true);
  const double Q21 = product_sources(t2, t1, beta, beta, A, T, true);
  const double Q22 = product_sources(t2, t2, beta, beta, A, T, true);

  comp[0] = mu1 * T + a11 * S1 + a12 * S2;
  comp[1] = mu2 * T + a21 * S1 + a22 * S2;
  comp[2] = mu1 * S1 + a11 * P11 + a12 * P12;
  comp[3] = mu1 * S2 + a11 * P12 + a12 * P22;
  comp[4] = mu2 * S1 + a21 * P11 + a22 * P12;
  comp[5] = mu2 * S2 + a21 * P12 + a22 * P22;
  comp[6] =
    (a11 * mu1 + a21 * mu2) * D1int +
    (a12 * mu1 + a22 * mu2) * D2int +
    (a11 * a11 + a21 * a21) * Q11 +
    (a11 * a12 + a21 * a22) * Q12 +
    (a12 * a11 + a22 * a21) * Q21 +
    (a12 * a12 + a22 * a22) * Q22;
}

// [[Rcpp::export]]
DataFrame simulate_hawkes_bivar_cluster_cpp(double T, double A,
                                            NumericVector mu,
                                            NumericMatrix alpha,
                                            double beta,
                                            double burnin) {
  Rcpp::RNGScope scope;
  if (mu.size() != 2 || alpha.nrow() != 2 || alpha.ncol() != 2) {
    stop("mu must have length 2 and alpha must be 2x2");
  }
  const double tmin = -burnin - A;
  std::vector<Event> events;
  const int n_imm1 = R::rpois(mu[0] * (T - tmin));
  const int n_imm2 = R::rpois(mu[1] * (T - tmin));
  events.reserve(static_cast<size_t>(std::max(100, 2 * (n_imm1 + n_imm2))));
  for (int n = 0; n < n_imm1; ++n) events.push_back(Event{tmin + R::runif(0.0, 1.0) * (T - tmin), 1});
  for (int n = 0; n < n_imm2; ++n) events.push_back(Event{tmin + R::runif(0.0, 1.0) * (T - tmin), 2});

  size_t q = 0;
  while (q < events.size()) {
    const Event parent = events[q];
    if (parent.t < T) {
      const double L = std::min(A, T - parent.t);
      if (L > 0.0) {
        const double trunc = 1.0 - std::exp(-beta * L);
        for (int child_type = 1; child_type <= 2; ++child_type) {
          const double a = alpha(child_type - 1, parent.type - 1);
          if (a <= 0.0) continue;
          const double mean_child = a * trunc / den_texp(beta, A);
          const int nc = R::rpois(mean_child);
          for (int c = 0; c < nc; ++c) {
            const double u = R::runif(0.0, 1.0);
            const double lag = -std::log(1.0 - u * trunc) / beta;
            events.push_back(Event{parent.t + lag, child_type});
          }
        }
      }
    }
    ++q;
  }

  std::vector<Event> keep;
  keep.reserve(events.size());
  for (const auto& e : events) if (e.t >= -A && e.t <= T) keep.push_back(e);
  std::sort(keep.begin(), keep.end(), [](const Event& a, const Event& b) {
    if (a.t == b.t) return a.type < b.type;
    return a.t < b.t;
  });
  NumericVector out_time(keep.size());
  IntegerVector out_type(keep.size());
  for (size_t i = 0; i < keep.size(); ++i) {
    out_time[i] = keep[i].t;
    out_type[i] = keep[i].type;
  }
  return DataFrame::create(_["time"] = out_time, _["type"] = out_type);
}

// [[Rcpp::export]]
List hawkes_bivar_loglik_score_cpp(NumericVector theta,
                                   NumericVector time,
                                   IntegerVector type,
                                   double T, double A,
                                   double rho_max = 0.98) {
  double mu1, mu2, a11, a12, a21, a22, beta;
  NumericVector zero(7);
  if (!parse_theta(theta, mu1, mu2, a11, a12, a21, a22, beta, rho_max)) {
    return List::create(_["valid"] = false, _["loglik"] = R_NegInf, _["score"] = zero);
  }
  std::vector<Event> events;
  std::vector<double> t1, t2;
  fill_events(time, type, events, t1, t2);

  double ll = 0.0;
  NumericVector score(7);
  for (const auto& e : events) {
    if (!(e.t >= 0.0 && e.t <= T)) continue;
    double X1, dX1, X2, dX2;
    x_dx_at(e.t, t1, beta, A, X1, dX1);
    x_dx_at(e.t, t2, beta, A, X2, dX2);
    const double lam1 = mu1 + a11 * X1 + a12 * X2;
    const double lam2 = mu2 + a21 * X1 + a22 * X2;
    if (!(lam1 > 0.0) || !(lam2 > 0.0) || !R_finite(lam1) || !R_finite(lam2)) {
      return List::create(_["valid"] = false, _["loglik"] = R_NegInf, _["score"] = zero);
    }
    double D1[7], D2[7];
    deriv_rows(X1, X2, dX1, dX2, a11, a12, a21, a22, D1, D2);
    if (e.type == 1) {
      ll += std::log(lam1);
      for (int k = 0; k < 7; ++k) score[k] += D1[k] / lam1;
    } else {
      ll += std::log(lam2);
      for (int k = 0; k < 7; ++k) score[k] += D2[k] / lam2;
    }
  }
  const double S1 = sum_K_source(t1, beta, A, T);
  const double S2 = sum_K_source(t2, beta, A, T);
  const double D1int = sum_dK_source(t1, beta, A, T);
  const double D2int = sum_dK_source(t2, beta, A, T);
  ll -= mu1 * T + mu2 * T + (a11 + a21) * S1 + (a12 + a22) * S2;
  score[0] -= T;
  score[1] -= T;
  score[2] -= S1;
  score[3] -= S2;
  score[4] -= S1;
  score[5] -= S2;
  score[6] -= (a11 + a21) * D1int + (a12 + a22) * D2int;
  return List::create(_["valid"] = true, _["loglik"] = ll, _["score"] = score);
}

// [[Rcpp::export]]
NumericVector hawkes_bivar_ls_moment_exact_cpp(NumericVector theta,
                                               NumericVector time,
                                               IntegerVector type,
                                               double T, double A,
                                               double rho_max = 0.98) {
  double mu1, mu2, a11, a12, a21, a22, beta;
  NumericVector bad(7); for (int k = 0; k < 7; ++k) bad[k] = 1e50;
  if (!parse_theta(theta, mu1, mu2, a11, a12, a21, a22, beta, rho_max)) return bad;

  std::vector<Event> events;
  std::vector<double> t1, t2;
  fill_events(time, type, events, t1, t2);

  NumericVector event(7), comp(7);
  for (const auto& e : events) {
    if (!(e.t >= 0.0 && e.t <= T)) continue;
    double X1, dX1, X2, dX2;
    x_dx_at(e.t, t1, beta, A, X1, dX1);
    x_dx_at(e.t, t2, beta, A, X2, dX2);
    double Dr1[7], Dr2[7];
    deriv_rows(X1, X2, dX1, dX2, a11, a12, a21, a22, Dr1, Dr2);
    if (e.type == 1) for (int k = 0; k < 7; ++k) event[k] += Dr1[k];
    else for (int k = 0; k < 7; ++k) event[k] += Dr2[k];
  }

  const double S1 = sum_K_source(t1, beta, A, T);
  const double S2 = sum_K_source(t2, beta, A, T);
  const double D1int = sum_dK_source(t1, beta, A, T);
  const double D2int = sum_dK_source(t2, beta, A, T);
  const double P11 = product_sources(t1, t1, beta, beta, A, T, false);
  const double P12 = product_sources(t1, t2, beta, beta, A, T, false);
  const double P22 = product_sources(t2, t2, beta, beta, A, T, false);
  const double Q11 = product_sources(t1, t1, beta, beta, A, T, true); // integral dX1 * X1
  const double Q12 = product_sources(t1, t2, beta, beta, A, T, true); // integral dX1 * X2
  const double Q21 = product_sources(t2, t1, beta, beta, A, T, true); // integral dX2 * X1
  const double Q22 = product_sources(t2, t2, beta, beta, A, T, true); // integral dX2 * X2

  comp[0] = mu1 * T + a11 * S1 + a12 * S2;
  comp[1] = mu2 * T + a21 * S1 + a22 * S2;
  comp[2] = mu1 * S1 + a11 * P11 + a12 * P12;
  comp[3] = mu1 * S2 + a11 * P12 + a12 * P22;
  comp[4] = mu2 * S1 + a21 * P11 + a22 * P12;
  comp[5] = mu2 * S2 + a21 * P12 + a22 * P22;
  comp[6] =
    (a11 * mu1 + a21 * mu2) * D1int +
    (a12 * mu1 + a22 * mu2) * D2int +
    (a11 * a11 + a21 * a21) * Q11 +
    (a11 * a12 + a21 * a22) * Q12 +
    (a12 * a11 + a22 * a21) * Q21 +
    (a12 * a12 + a22 * a22) * Q22;

  NumericVector out(7);
  for (int k = 0; k < 7; ++k) out[k] = (event[k] - comp[k]) / T;
  return out;
}



static inline int sanitize_aug_degree(int degree) {
  if (degree < 1) return 1;
  if (degree > 2) return 2;
  return degree;
}

static inline void build_aug_breaks(const std::vector<Event>& events,
                                    double A, double T,
                                    std::vector<double>& breaks) {
  breaks.clear();
  breaks.reserve(2 * events.size() + 2);
  breaks.push_back(0.0);
  breaks.push_back(T);
  for (const auto& e : events) {
    if (R_finite(e.t) && e.t > 0.0 && e.t < T) breaks.push_back(e.t);
    const double u = e.t + A;
    if (R_finite(u) && u > 0.0 && u < T) breaks.push_back(u);
  }
  std::sort(breaks.begin(), breaks.end());
  std::vector<double> uniq;
  uniq.reserve(breaks.size());
  const double eps = 1e-11 * std::max(1.0, T);
  for (double x : breaks) {
    if (!R_finite(x)) continue;
    x = std::max(0.0, std::min(T, x));
    if (uniq.empty() || std::fabs(x - uniq.back()) > eps) uniq.push_back(x);
  }
  if (uniq.empty() || uniq.front() > 0.0) uniq.insert(uniq.begin(), 0.0);
  if (uniq.back() < T) uniq.push_back(T);
  breaks.swap(uniq);
}

static inline bool eval_aug_state(double t,
                                  const std::vector<double>& t1,
                                  const std::vector<double>& t2,
                                  double mu1, double mu2,
                                  double a11, double a12,
                                  double a21, double a22,
                                  double beta, double A,
                                  double s1, double s2,
                                  double& lam1, double& lam2,
                                  double* D1, double* D2,
                                  double* f1, double* f2) {
  double X1, dX1, X2, dX2;
  x_dx_at(t, t1, beta, A, X1, dX1);
  x_dx_at(t, t2, beta, A, X2, dX2);
  lam1 = mu1 + a11 * X1 + a12 * X2;
  lam2 = mu2 + a21 * X1 + a22 * X2;
  if (!(lam1 > 0.0) || !(lam2 > 0.0) || !R_finite(lam1) || !R_finite(lam2)) return false;
  deriv_rows(X1, X2, dX1, dX2, a11, a12, a21, a22, D1, D2);
  aug_factors(lam1, s1, f1);
  aug_factors(lam2, s2, f2);
  return true;
}

static inline void integrate_aug_comp_event_adaptive(
    const std::vector<Event>& events,
    const std::vector<double>& t1,
    const std::vector<double>& t2,
    double mu1, double mu2,
    double a11, double a12,
    double a21, double a22,
    double beta, double A, double T,
    double s1, double s2,
    int degree,
    double hmax,
    NumericVector& comp,
    int& n_eval,
    int& n_intervals) {

  static const double nodes[5] = {
    -0.90617984593866399280,
    -0.53846931010568309104,
     0.0,
     0.53846931010568309104,
     0.90617984593866399280
  };
  static const double weights[5] = {
    0.23692688505618908751,
    0.47862867049936646804,
    0.56888888888888888889,
    0.47862867049936646804,
    0.23692688505618908751
  };

  std::vector<double> breaks;
  build_aug_breaks(events, A, T, breaks);
  if (!(hmax > 0.0) || !R_finite(hmax)) hmax = T;
  hmax = std::max(hmax, 1e-8);

  n_eval = 0;
  n_intervals = 0;
  double D1[7], D2[7], f1[3], f2[3];
  for (size_t ii = 0; ii + 1 < breaks.size(); ++ii) {
    const double left0 = breaks[ii];
    const double right0 = breaks[ii + 1];
    if (!(right0 > left0)) continue;
    const int nsub = std::max(1, static_cast<int>(std::ceil((right0 - left0) / hmax)));
    const double step = (right0 - left0) / static_cast<double>(nsub);
    for (int ss = 0; ss < nsub; ++ss) {
      const double left = left0 + ss * step;
      const double right = (ss == nsub - 1) ? right0 : (left + step);
      if (!(right > left)) continue;
      const double mid = 0.5 * (left + right);
      const double half = 0.5 * (right - left);
      ++n_intervals;
      for (int g = 0; g < 5; ++g) {
        const double tt = mid + half * nodes[g];
        double lam1, lam2;
        if (!eval_aug_state(tt, t1, t2, mu1, mu2, a11, a12, a21, a22,
                            beta, A, s1, s2, lam1, lam2, D1, D2, f1, f2)) {
          continue;
        }
        const double ww = half * weights[g];
        ++n_eval;
        for (int r = 1; r <= degree; ++r) {
          for (int k = 0; k < 7; ++k) {
            comp[7 * r + k] += ww * (lam1 * f1[r] * D1[k] + lam2 * f2[r] * D2[k]);
          }
        }
      }
    }
  }
}

static inline void integrate_aug_omega_event_adaptive(
    const std::vector<Event>& events,
    const std::vector<double>& t1,
    const std::vector<double>& t2,
    double mu1, double mu2,
    double a11, double a12,
    double a21, double a22,
    double beta, double A, double T,
    double s1, double s2,
    int degree,
    double hmax,
    NumericMatrix& Omega_aug,
    int& n_eval,
    int& n_intervals) {

  static const double nodes[5] = {
    -0.90617984593866399280,
    -0.53846931010568309104,
     0.0,
     0.53846931010568309104,
     0.90617984593866399280
  };
  static const double weights[5] = {
    0.23692688505618908751,
    0.47862867049936646804,
    0.56888888888888888889,
    0.47862867049936646804,
    0.23692688505618908751
  };

  std::vector<double> breaks;
  build_aug_breaks(events, A, T, breaks);
  if (!(hmax > 0.0) || !R_finite(hmax)) hmax = T;
  hmax = std::max(hmax, 1e-8);

  n_eval = 0;
  n_intervals = 0;
  double D1[7], D2[7], f1[3], f2[3];
  const int q_aug = 7 * (degree + 1);
  for (size_t ii = 0; ii + 1 < breaks.size(); ++ii) {
    const double left0 = breaks[ii];
    const double right0 = breaks[ii + 1];
    if (!(right0 > left0)) continue;
    const int nsub = std::max(1, static_cast<int>(std::ceil((right0 - left0) / hmax)));
    const double step = (right0 - left0) / static_cast<double>(nsub);
    for (int ss = 0; ss < nsub; ++ss) {
      const double left = left0 + ss * step;
      const double right = (ss == nsub - 1) ? right0 : (left + step);
      if (!(right > left)) continue;
      const double mid = 0.5 * (left + right);
      const double half = 0.5 * (right - left);
      ++n_intervals;
      for (int g = 0; g < 5; ++g) {
        const double tt = mid + half * nodes[g];
        double lam1, lam2;
        if (!eval_aug_state(tt, t1, t2, mu1, mu2, a11, a12, a21, a22,
                            beta, A, s1, s2, lam1, lam2, D1, D2, f1, f2)) {
          continue;
        }
        const double ww = half * weights[g] / T;
        ++n_eval;
        for (int r = 0; r <= degree; ++r) {
          for (int a = 0; a < 7; ++a) {
            const int row = 7 * r + a;
            for (int rp = 0; rp <= degree; ++rp) {
              for (int b = 0; b < 7; ++b) {
                const int col = 7 * rp + b;
                if (row < q_aug && col < q_aug) {
                  Omega_aug(row, col) += ww *
                    (lam1 * f1[r] * f1[rp] * D1[a] * D1[b] +
                     lam2 * f2[r] * f2[rp] * D2[a] * D2[b]);
                }
              }
            }
          }
        }
      }
    }
  }
}

// [[Rcpp::export]]
NumericVector hawkes_bivar_aug_moment_adaptive_cpp(NumericVector theta,
                                                    NumericVector time,
                                                    IntegerVector type,
                                                    double T, double A,
                                                    double rho_max = 0.98,
                                                    double s1 = 0.4,
                                                    double s2 = 0.4,
                                                    int aug_degree = 1,
                                                    double hmax = 0.5) {
  aug_degree = sanitize_aug_degree(aug_degree);
  const int q_aug = 7 * (aug_degree + 1);
  double mu1, mu2, a11, a12, a21, a22, beta;
  NumericVector bad(q_aug); for (int k = 0; k < q_aug; ++k) bad[k] = 1e50;
  if (!parse_theta(theta, mu1, mu2, a11, a12, a21, a22, beta, rho_max)) return bad;
  if (!(T > 0.0) || !(s1 > 0.0) || !(s2 > 0.0)) return bad;

  std::vector<Event> events;
  std::vector<double> t1, t2;
  fill_events(time, type, events, t1, t2);

  NumericVector event(q_aug), comp(q_aug);

  double D1[7], D2[7], f1[3], f2[3];
  for (const auto& e : events) {
    if (!(e.t >= 0.0 && e.t <= T)) continue;
    double lam1, lam2;
    if (!eval_aug_state(e.t, t1, t2, mu1, mu2, a11, a12, a21, a22,
                        beta, A, s1, s2, lam1, lam2, D1, D2, f1, f2)) return bad;
    if (e.type == 1) {
      for (int r = 0; r <= aug_degree; ++r) for (int k = 0; k < 7; ++k) event[7 * r + k] += f1[r] * D1[k];
    } else {
      for (int r = 0; r <= aug_degree; ++r) for (int k = 0; k < 7; ++k) event[7 * r + k] += f2[r] * D2[k];
    }
  }

  NumericVector comp0(7);
  ls_comp_exact(t1, t2, mu1, mu2, a11, a12, a21, a22, beta, A, T, comp0);
  for (int k = 0; k < 7; ++k) comp[k] = comp0[k];

  int n_eval = 0, n_intervals = 0;
  integrate_aug_comp_event_adaptive(events, t1, t2, mu1, mu2, a11, a12, a21, a22,
                                    beta, A, T, s1, s2, aug_degree, hmax,
                                    comp, n_eval, n_intervals);
  if (n_eval <= 0) return bad;

  NumericVector out(q_aug);
  for (int k = 0; k < q_aug; ++k) out[k] = (event[k] - comp[k]) / T;
  return out;
}

// [[Rcpp::export]]
List hawkes_bivar_aug_omega_adaptive_cpp(NumericVector theta,
                                          NumericVector time,
                                          IntegerVector type,
                                          double T, double A,
                                          double rho_max = 0.98,
                                          double s1 = 0.4,
                                          double s2 = 0.4,
                                          int aug_degree = 1,
                                          double hmax = 0.5) {
  aug_degree = sanitize_aug_degree(aug_degree);
  const int q_aug = 7 * (aug_degree + 1);
  double mu1, mu2, a11, a12, a21, a22, beta;
  NumericMatrix Omega_aug(q_aug, q_aug);
  if (!parse_theta(theta, mu1, mu2, a11, a12, a21, a22, beta, rho_max) ||
      !(s1 > 0.0) || !(s2 > 0.0) || !(T > 0.0)) {
    return List::create(_["valid"] = false, _["Omega_aug"] = Omega_aug,
                        _["n_eval"] = 0, _["n_intervals"] = 0,
                        _["q_aug"] = q_aug, _["aug_degree"] = aug_degree);
  }
  std::vector<Event> events;
  std::vector<double> t1, t2;
  fill_events(time, type, events, t1, t2);
  int n_eval = 0, n_intervals = 0;
  integrate_aug_omega_event_adaptive(events, t1, t2, mu1, mu2, a11, a12, a21, a22,
                                     beta, A, T, s1, s2, aug_degree, hmax,
                                     Omega_aug, n_eval, n_intervals);
  return List::create(_["valid"] = n_eval > 0, _["Omega_aug"] = Omega_aug,
                      _["n_eval"] = n_eval, _["n_intervals"] = n_intervals,
                      _["q_aug"] = q_aug, _["aug_degree"] = aug_degree,
                      _["hmax"] = hmax);
}


// [[Rcpp::export]]
List hawkes_bivar_time_average_matrices_cpp(NumericVector theta,
                                            NumericVector time,
                                            IntegerVector type,
                                            NumericVector t_eval,
                                            double A,
                                            double rho_max = 0.98,
                                            double s1 = 0.4,
                                            double s2 = 0.4,
                                            int aug_degree = 1) {
  aug_degree = sanitize_aug_degree(aug_degree);
  const int q_aug = 7 * (aug_degree + 1);
  double mu1, mu2, a11, a12, a21, a22, beta;
  NumericMatrix I(7, 7), A_ls(7, 7), Omega_ls(7, 7), A_aug(q_aug, 7), Omega_aug(q_aug, q_aug);
  if (!parse_theta(theta, mu1, mu2, a11, a12, a21, a22, beta, rho_max) || !(s1 > 0.0) || !(s2 > 0.0)) {
    return List::create(_["valid"] = false, _["I"] = I, _["A_ls"] = A_ls,
                        _["Omega_ls"] = Omega_ls, _["A_aug"] = A_aug,
                        _["Omega_aug"] = Omega_aug, _["n_eval"] = 0,
                        _["q_aug"] = q_aug, _["aug_degree"] = aug_degree);
  }
  std::vector<Event> events;
  std::vector<double> t1, t2;
  fill_events(time, type, events, t1, t2);
  int n_used = 0;
  double D1[7], D2[7], f1[3], f2[3];
  for (int m = 0; m < t_eval.size(); ++m) {
    const double tt = t_eval[m];
    if (!R_finite(tt)) continue;
    double lam1, lam2;
    if (!eval_aug_state(tt, t1, t2, mu1, mu2, a11, a12, a21, a22,
                        beta, A, s1, s2, lam1, lam2, D1, D2, f1, f2)) continue;
    ++n_used;
    for (int a = 0; a < 7; ++a) {
      for (int b = 0; b < 7; ++b) {
        const double o1 = D1[a] * D1[b];
        const double o2 = D2[a] * D2[b];
        I(a, b) += o1 / lam1 + o2 / lam2;
        A_ls(a, b) += o1 + o2;
        Omega_ls(a, b) += lam1 * o1 + lam2 * o2;
      }
    }
    for (int r = 0; r <= aug_degree; ++r) {
      for (int a = 0; a < 7; ++a) {
        const int row = 7 * r + a;
        for (int b = 0; b < 7; ++b) {
          A_aug(row, b) += f1[r] * D1[a] * D1[b] + f2[r] * D2[a] * D2[b];
        }
        for (int rp = 0; rp <= aug_degree; ++rp) {
          for (int b = 0; b < 7; ++b) {
            const int col = 7 * rp + b;
            Omega_aug(row, col) +=
              lam1 * f1[r] * f1[rp] * D1[a] * D1[b] +
              lam2 * f2[r] * f2[rp] * D2[a] * D2[b];
          }
        }
      }
    }
  }
  if (n_used > 0) {
    const double inv_n = 1.0 / static_cast<double>(n_used);
    for (int a = 0; a < 7; ++a) for (int b = 0; b < 7; ++b) {
      I(a, b) *= inv_n;
      A_ls(a, b) *= inv_n;
      Omega_ls(a, b) *= inv_n;
    }
    for (int a = 0; a < q_aug; ++a) {
      for (int b = 0; b < 7; ++b) A_aug(a, b) *= inv_n;
      for (int b = 0; b < q_aug; ++b) Omega_aug(a, b) *= inv_n;
    }
  }
  return List::create(_["valid"] = true, _["I"] = I, _["A_ls"] = A_ls,
                      _["Omega_ls"] = Omega_ls, _["A_aug"] = A_aug,
                      _["Omega_aug"] = Omega_aug, _["n_eval"] = n_used,
                      _["q_aug"] = q_aug, _["aug_degree"] = aug_degree);
}
