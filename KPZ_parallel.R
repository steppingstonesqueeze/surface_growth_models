# KPZ parallelised for large lattices --
# KPZ_surface_growth.R
# Fast KPZ integrator (1D & 2D) with periodic BCs, explicit Euler time-step.
# Parallelized finite-difference stencils via RcppParallel; noise generated in R each step.
suppressPackageStartupMessages({
  library(ggplot2)
  library(reshape2)
  library(Rcpp)
  library(RcppParallel)
})

# ---- Roughness ----
roughness <- function(h) {
  mu <- mean(h)
  sqrt(mean((h - mu)^2))
}

# ---- Two-point correlation (2D) ----
two_point_corr_2d <- function(h) {
  L <- nrow(h)
  maxr <- floor(L/2)
  cr <- numeric(maxr)
  for (r in 1:maxr) {
    hx <- h - h[, c((r+1):L, 1:r)]
    hy <- h - h[c((r+1):L, 1:r), ]
    cr[r] <- mean(hx^2 + hy^2) / 2
  }
  data.frame(r = 1:maxr, C_r = cr)
}

# ---- Parallel KPZ stencil (1D & 2D), Euler step ----
# h <- h + dt*(nu*lap + (lambda/2)*|grad|^2) + sqrt(2*D*dt) * eta
# Here, "eta" (noise) is passed in from R to keep RNG thread-safe.

Rcpp::cppFunction(code = '
#include <Rcpp.h>
#include <RcppParallel.h>
using namespace Rcpp;
using namespace RcppParallel;

inline int modp(int i, int L){ int m = i % L; return (m < 0) ? m + L : m; }

// ---------------- 1D KPZ ----------------
struct KPZ1DStep : public Worker {
  const RVector<double> h_in;
  const RVector<double> noise;
  RVector<double> h_out;
  const int L;
  const double nu, lambda, dt, noise_scale;

  KPZ1DStep(NumericVector h_in_, NumericVector noise_, NumericVector h_out_,
            int L_, double nu_, double lambda_, double dt_, double ns_)
    : h_in(h_in_), noise(noise_), h_out(h_out_),
      L(L_), nu(nu_), lambda(lambda_), dt(dt_), noise_scale(ns_) {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; ++i) {
      int il = (i == 0 ? L-1 : (int)i-1);
      int ir = (i == (size_t)L-1 ? 0 : (int)i+1);
      double hi = h_in[i];
      double hl = h_in[il];
      double hr = h_in[ir];

      double lap = hl + hr - 2.0*hi;
      double grad = (hr - hl) * 0.5;             // central difference
      double det = nu * lap + 0.5 * lambda * grad * grad;

      h_out[i] = hi + dt * det + noise_scale * noise[i];
    }
  }
};

// [[Rcpp::export]]
NumericVector kpz1d_step_parallel(NumericVector h, NumericVector noise,
                                  double nu, double lambda, double dt, double noise_scale) {
  int L = h.size();
  NumericVector out(L);
  KPZ1DStep worker(h, noise, out, L, nu, lambda, dt, noise_scale);
  parallelFor(0, (size_t)L, worker);
  return out;
}

// ---------------- 2D KPZ ----------------
struct KPZ2DStep : public Worker {
  const RMatrix<double> Hin;
  const RMatrix<double> Noise;
  RMatrix<double> Hout;
  const int L;
  const double nu, lambda, dt, noise_scale;

  KPZ2DStep(NumericMatrix Hin_, NumericMatrix Noise_, NumericMatrix Hout_,
            int L_, double nu_, double lambda_, double dt_, double ns_)
    : Hin(Hin_), Noise(Noise_), Hout(Hout_),
      L(L_), nu(nu_), lambda(lambda_), dt(dt_), noise_scale(ns_) {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t k = begin; k < end; ++k) {
      int i = (int)(k / L);
      int j = (int)(k % L);

      int iu = (i == 0 ? L-1 : i-1);
      int id = (i == L-1 ? 0 : i+1);
      int jl = (j == 0 ? L-1 : j-1);
      int jr = (j == L-1 ? 0 : j+1);

      double h  = Hin(i,j);
      double hu = Hin(iu,j);
      double hd = Hin(id,j);
      double hl = Hin(i,jl);
      double hr = Hin(i,jr);

      double lap = hu + hd + hl + hr - 4.0*h;
      double gx  = (hr - hl) * 0.5;
      double gy  = (hd - hu) * 0.5;
      double det = nu * lap + 0.5 * lambda * (gx*gx + gy*gy);

      Hout(i,j) = h + dt * det + noise_scale * Noise(i,j);
    }
  }
};

// [[Rcpp::export]]
NumericMatrix kpz2d_step_parallel(NumericMatrix H, NumericMatrix Noise,
                                  double nu, double lambda, double dt, double noise_scale) {
  int L = H.nrow();
  NumericMatrix out(L, L);
  size_t N = (size_t)L * (size_t)L;
  KPZ2DStep worker(H, Noise, out, L, nu, lambda, dt, noise_scale);
  parallelFor(0, N, worker);
  return out;
}
', plugins = c("cpp11"))

# ---- Runner ----
# Note: noise_scale = sqrt(2*D*dt)
run_KPZ <- function(L = 256, T_max = 2000, dim = 2,
                    nu = 1.0, lambda = 1.0, D = 1.0, dt = 0.01,
                    seed = 123) {
  set.seed(seed)
  noise_scale <- sqrt(2 * D * dt)
  
  if (dim == 1) {
    h <- numeric(L)
    Wt <- numeric(T_max)
    for (t in 1:T_max) {
      eta <- rnorm(L)
      h <- kpz1d_step_parallel(h, eta, nu, lambda, dt, noise_scale)
      Wt[t] <- roughness(h)
    }
    df <- data.frame(x = 1:L, h = h)
    p1 <- ggplot(df, aes(x, h)) +
      geom_line(linewidth = 0.5) +
      labs(title = "KPZ 1D Surface (final)", x = "x", y = "h") +
      theme_minimal()
    
    p2 <- ggplot(data.frame(t = 1:T_max, W = Wt), aes(t, W)) +
      geom_line(linewidth = 0.5) +
      labs(title = "KPZ Roughness W(t)", x = "t", y = "W") +
      theme_minimal()
    
    print(p1); print(p2)
    invisible(list(h = h, Wt = Wt))
  } else {
    H <- matrix(0.0, nrow = L, ncol = L)
    Wt <- numeric(T_max)
    for (t in 1:T_max) {
      Eta <- matrix(rnorm(L*L), nrow = L, ncol = L)
      H <- kpz2d_step_parallel(H, Eta, nu, lambda, dt, noise_scale)
      Wt[t] <- roughness(as.numeric(H))
    }
    df <- melt(H)
    names(df) <- c("x","y","h")
    p1 <- ggplot(df, aes(x, y, fill = h)) +
      geom_raster(interpolate = FALSE) +
      scale_fill_gradient(name = "h", limits = range(df$h)) +
      coord_fixed() +
      labs(title = "KPZ 2D Surface (final)", x = "x", y = "y") +
      theme_minimal()
    
    p2 <- ggplot(data.frame(t = 1:T_max, W = Wt), aes(t, W)) +
      geom_line(linewidth = 0.5) +
      labs(title = "KPZ Roughness W(t)", x = "t", y = "W") +
      theme_minimal()
    
    corr_df <- two_point_corr_2d(H)
    p3 <- ggplot(corr_df, aes(r, C_r)) +
      geom_line(linewidth = 0.6) +
      labs(title = "KPZ 2D Two-point Roughness Correlation C(r)",
           x = "r (lattice units)", y = "C(r)") +
      theme_minimal()
    
    print(p1); print(p2); print(p3)
    invisible(list(H = H, Wt = Wt, corr = corr_df))
  }
}
