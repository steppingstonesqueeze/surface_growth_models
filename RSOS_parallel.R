# Parallelised RSOS for running large 2D lattices -- 

# RSOS_surface_growth.R
# Fast RSOS growth with parallel red–black sweeps (growth-only).
# Optional slower symmetric (+/-1) RSOS fallback is provided (sequential).
suppressPackageStartupMessages({
  library(ggplot2)
  library(reshape2)
  library(Rcpp)
  library(RcppParallel)
})

# ---- Roughness (width) ----
roughness <- function(h) {
  mu <- mean(h)
  sqrt(mean((h - mu)^2))
}

# ---- Two-point correlation (2D): axis-averaged height-difference structure function ----
# C(r) = average over x,y of (h(x,y) - h(x+r,y))^2 and (h(x,y) - h(x,y+r))^2 (periodic), averaged.
two_point_corr_2d <- function(h) {
  L <- nrow(h)
  maxr <- floor(L/2)
  cr <- numeric(maxr)
  for (r in 1:maxr) {
    # shift in x
    hx <- h - h[, c((r+1):L, 1:r)]
    # shift in y
    hy <- h - h[c((r+1):L, 1:r), ]
    cr[r] <- mean(hx^2 + hy^2) / 2
  }
  data.frame(r = 1:maxr, C_r = cr)
}

# ---- Parallel RSOS (growth-only) via RcppParallel ----
# Red–black checkerboard: each sweep does black then red sublattice.
# Constraint: |h_i - h_neighbor| <= 1 after deposition attempt; if violates, reject.
# We implement separate 1D and 2D workers.

Rcpp::cppFunction(code = '
#include <Rcpp.h>
#include <RcppParallel.h>
using namespace Rcpp;
using namespace RcppParallel;

// Periodic helpers
inline int modp(int i, int L){ int m = i % L; return (m < 0) ? m + L : m; }

// ---------------- 1D RSOS growth-only (red-black) ----------------
struct RSOS1DPass : public Worker {
  const RVector<int> h_in;
  RVector<int> h_out;
  const int L;
  const int parity; // 0 for even, 1 for odd

  RSOS1DPass(IntegerVector h_in_, IntegerVector h_out_, int L_, int parity_)
    : h_in(h_in_), h_out(h_out_), L(L_), parity(parity_) {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t idx = begin; idx < end; ++idx) {
      int i = (int)idx;
      // Only update indices with i % 2 == parity
      if ((i & 1) != parity) {
        h_out[i] = h_in[i];
        continue;
      }
      int il = modp(i-1, L), ir = modp(i+1, L);
      int hi = h_in[i];
      int hl = h_in[il];
      int hr = h_in[ir];

      // propose growth-only: hi_new = hi + 1
      int hi_new = hi + 1;
      // RSOS constraint: |hi_new - neighbor| <= 1
      if (std::abs(hi_new - hl) <= 1 && std::abs(hi_new - hr) <= 1) {
        h_out[i] = hi_new;
      } else {
        h_out[i] = hi;
      }
    }
  }
};

// ---------------- 2D RSOS growth-only (red-black) ----------------
struct RSOS2DPass : public Worker {
  const RMatrix<int> H_in;
  RMatrix<int> H_out;
  const int L;
  const int parity; // (i+j)%2 parity

  RSOS2DPass(IntegerMatrix Hin, IntegerMatrix Hout, int L_, int parity_)
    : H_in(Hin), H_out(Hout), L(L_), parity(parity_) {}

  void operator()(std::size_t begin, std::size_t end) {
    // Flattened index => (i,j)
    for (std::size_t k = begin; k < end; ++k) {
      int i = (int)(k / L);
      int j = (int)(k % L);
      if (((i + j) & 1) != parity) {
        H_out(i,j) = H_in(i,j);
        continue;
      }
      int iu = (i == 0 ? L-1 : i-1);
      int id = (i == L-1 ? 0 : i+1);
      int jl = (j == 0 ? L-1 : j-1);
      int jr = (j == L-1 ? 0 : j+1);

      int h  = H_in(i,j);
      int hu = H_in(iu,j);
      int hd = H_in(id,j);
      int hl = H_in(i,jl);
      int hr = H_in(i,jr);

      int hnew = h + 1; // growth-only
      if (std::abs(hnew - hu) <= 1 &&
          std::abs(hnew - hd) <= 1 &&
          std::abs(hnew - hl) <= 1 &&
          std::abs(hnew - hr) <= 1) {
        H_out(i,j) = hnew;
      } else {
        H_out(i,j) = h;
      }
    }
  }
};

// [[Rcpp::export]]
IntegerVector rsos1d_growth_parallel(IntegerVector h, int sweeps) {
  int L = h.size();
  IntegerVector a = clone(h), b = clone(h);
  for (int s = 0; s < sweeps; ++s) {
    // even sublattice
    RSOS1DPass pass_even(a, b, L, 0);
    parallelFor(0, (size_t)L, pass_even);
    // odd sublattice (swap roles)
    RSOS1DPass pass_odd(b, a, L, 1);
    parallelFor(0, (size_t)L, pass_odd);
  }
  return a;
}

// [[Rcpp::export]]
IntegerMatrix rsos2d_growth_parallel(IntegerMatrix H, int sweeps) {
  int L = H.nrow();
  IntegerMatrix A = clone(H), B = clone(H);
  size_t N = (size_t)L * (size_t)L;

  for (int s = 0; s < sweeps; ++s) {
    RSOS2DPass pass_black(A, B, L, 0);
    parallelFor(0, N, pass_black);
    RSOS2DPass pass_red(B, A, L, 1);
    parallelFor(0, N, pass_red);
  }
  return A;
}
', plugins = c("cpp11"))

# ---- SLOW symmetric RSOS (±1) fallback (sequential random sequential updates) ----
# You probably won’t need this; it’s here for completeness.
rsos_update_seq_1d <- function(h) {
  L <- length(h)
  i <- sample.int(L, 1)
  il <- if (i == 1) L else i - 1
  ir <- if (i == L) 1 else i + 1
  step <- sample(c(-1L, 1L), 1)
  h_new <- h[i] + step
  if (abs(h_new - h[il]) <= 1 && abs(h_new - h[ir]) <= 1) h[i] <- h_new
  h
}
rsos_update_seq_2d <- function(H) {
  L <- nrow(H)
  i <- sample.int(L, 1); j <- sample.int(L, 1)
  iu <- if (i==1) L else i-1; id <- if (i==L) 1 else i+1
  jl <- if (j==1) L else j-1; jr <- if (j==L) 1 else j+1
  step <- sample(c(-1L, 1L), 1)
  h_new <- H[i,j] + step
  nbrs <- c(H[iu,j], H[id,j], H[i,jl], H[i,jr])
  if (all(abs(h_new - nbrs) <= 1)) H[i,j] <- h_new
  H
}

# ---- Runner ----
# mode = "growth" (fast parallel deposition-only) or "symmetric" (slow ±1 sequential)
run_RSOS <- function(L = 256, T_max = 2000, dim = 2, mode = "growth",
                     sweeps_per_frame = 1, seed = 123) {
  set.seed(seed)
  if (dim == 1) {
    h <- integer(L)
    Wt <- numeric(T_max)
    for (t in 1:T_max) {
      if (mode == "growth") {
        h <- rsos1d_growth_parallel(h, sweeps = sweeps_per_frame)
      } else {
        # ~L random sequential attempts per "sweep"
        for (k in 1:(sweeps_per_frame * L)) h <- rsos_update_seq_1d(h)
      }
      Wt[t] <- roughness(h)
    }
    df <- data.frame(x = 1:L, h = as.numeric(h))
    p1 <- ggplot(df, aes(x, h)) +
      geom_line(linewidth = 0.5) +
      labs(title = "RSOS 1D Surface (final)", x = "x", y = "h") +
      theme_minimal()
    
    p2 <- ggplot(data.frame(t = 1:T_max, W = Wt), aes(t, W)) +
      geom_line(linewidth = 0.5) +
      labs(title = "RSOS Roughness W(t)", x = "t", y = "W") +
      theme_minimal()
    
    print(p1); print(p2)
    invisible(list(h = h, Wt = Wt))
  } else {
    H <- matrix(0L, nrow = L, ncol = L)
    Wt <- numeric(T_max)
    for (t in 1:T_max) {
      if (mode == "growth") {
        H <- rsos2d_growth_parallel(H, sweeps = sweeps_per_frame)
      } else {
        # ~L^2 random sequential attempts per "sweep"
        for (k in 1:(sweeps_per_frame * L * L)) H <- rsos_update_seq_2d(H)
      }
      Wt[t] <- roughness(as.numeric(H))
    }
    df <- melt(H)
    names(df) <- c("x","y","h")
    p1 <- ggplot(df, aes(x, y, fill = h)) +
      geom_raster(interpolate = FALSE) +
      scale_fill_gradient(name = "h", limits = range(df$h)) +
      coord_fixed() +
      labs(title = "RSOS 2D Surface (final)", x = "x", y = "y") +
      theme_minimal()
    
    p2 <- ggplot(data.frame(t = 1:T_max, W = Wt), aes(t, W)) +
      geom_line(linewidth = 0.5) +
      labs(title = "RSOS Roughness W(t)", x = "t", y = "W") +
      theme_minimal()
    
    corr_df <- two_point_corr_2d(H)
    p3 <- ggplot(corr_df, aes(r, C_r)) +
      geom_line(linewidth = 0.6) +
      labs(title = "RSOS 2D Two-point Roughness Correlation C(r)",
           x = "r (lattice units)", y = "C(r)") +
      theme_minimal()
    
    print(p1); print(p2); print(p3)
    invisible(list(H = H, Wt = Wt, corr = corr_df))
  }
}
