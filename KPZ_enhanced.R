# KPZ 1D and 2D

# KPZ model in 1 and 2D

# KPZ_surface_growth.R
library(ggplot2)
library(reshape2)

# Roughness
roughness <- function(h) {
  mean_val <- mean(h)
  sqrt(mean((h - mean_val)^2))
}

# Two-point correlation for 2D - simplified via averaging and presenting on 1 axis
two_point_corr <- function(h) {
  L <- nrow(h)
  corr <- numeric(L %/% 2)
  for (r in 1:(L %/% 2)) {
    diffs <- (h - cbind(h[, (r+1):L], h[, 1:r]))
    corr[r] <- mean(diffs^2)
  }
  data.frame(r = 1:(L %/% 2), C_r = corr)
}

# Two point radially averaged correlation function

# ---- True isotropic radial 2D two-point correlation ----
# Returns height–difference correlation G(r) = < [h(x) - h(x+r)]^2 >
# Radially averaged over all directions with same |r|
two_point_corr_2d_radial <- function(h) {
  L <- nrow(h)
  # Precompute coordinate grid
  coords <- expand.grid(x = 1:L, y = 1:L)
  center_vals <- as.numeric(h[cbind(coords$x, coords$y)])
  
  # All displacements
  dx <- outer(1:L, 1:L, "-")
  dx <- ifelse(dx >  L/2, dx - L,
               ifelse(dx < -L/2, dx + L, dx))
  dy <- t(dx)  # same shape
  
  # Flatten all pairs into displacement vectors
  disps <- expand.grid(dx = 0:(L-1), dy = 0:(L-1))
  disps$dx <- ifelse(disps$dx > L/2, disps$dx - L, disps$dx)
  disps$dy <- ifelse(disps$dy > L/2, disps$dy - L, disps$dy)
  
  # Remove duplicate directions by symmetry (keep unique displacement vectors)
  disps <- unique(disps)
  
  # Bin by radial distance
  disps$r <- sqrt(disps$dx^2 + disps$dy^2)
  disps$rbin <- round(disps$r)  # integer binning
  
  max_r <- floor(sqrt(2) * L / 2)
  G_r <- numeric(max_r)
  counts <- numeric(max_r)
  
  for (rb in 1:max_r) {
    # Find displacement vectors in this radial bin
    vecs <- subset(disps, rbin == rb)
    if (nrow(vecs) == 0) next
    
    diffsq <- c()
    for (k in 1:nrow(vecs)) {
      dx <- vecs$dx[k]; dy <- vecs$dy[k]
      shifted <- h[((1:L - 1 + dx) %% L) + 1, ((1:L - 1 + dy) %% L) + 1]
      diffsq <- c(diffsq, (as.numeric(h) - as.numeric(shifted))^2)
    }
    G_r[rb] <- mean(diffsq)
    counts[rb] <- length(diffsq)
  }
  
  data.frame(r = 1:max_r, C_r = G_r, count = counts)
}

# ===================== ADDITIONS: spectra + exponents =====================

# 1D structure factor S(k)
structure_factor1d <- function(h1d, detrend = TRUE) {
  N <- length(h1d)
  x <- if (detrend) as.numeric(h1d) - mean(h1d) else as.numeric(h1d)
  Hk <- fft(x)
  S  <- (Mod(Hk)^2) / N
  kraw <- 0:(N-1)
  kabs <- pmin(kraw, N - kraw)
  kmax <- floor(N/2)
  out_k <- 0:kmax
  out_S <- numeric(length(out_k))
  for (i in seq_along(out_k)) {
    mask <- (kabs == out_k[i])
    out_S[i] <- mean(S[mask])
  }
  data.frame(k = out_k, S_k = out_S)
}

# 2D structure factor S(k), radially binned
structure_factor2d <- function(h2d, detrend = TRUE) {
  nr <- nrow(h2d); nc <- ncol(h2d)
  X <- if (detrend) h2d - mean(h2d) else h2d
  Hk <- fft(X)
  S  <- (Mod(Hk)^2) / (nr * nc)
  kx <- 0:(nc-1); kx <- pmin(kx, nc - kx)
  ky <- 0:(nr-1); ky <- pmin(ky, nr - ky)
  K  <- outer(ky, kx, function(a,b) sqrt(a^2 + b^2))
  kmax <- floor(min(nr, nc)/2)
  bins <- 0:kmax
  Sk   <- numeric(kmax + 1)
  for (i in seq_along(bins)) {
    mask <- K >= (bins[i] - 0.5) & K < (bins[i] + 0.5)
    if (!any(mask)) { Sk[i] <- NA_real_; next }
    Sk[i] <- mean(S[mask])
  }
  data.frame(k = bins, S_k = Sk)
}

# Fit alpha from S(k) ~ k^{-(2α + d)}
fit_alpha_from_Sk <- function(Sk_df, d, k_min = NULL, k_max = NULL) {
  df <- Sk_df
  df <- df[is.finite(df$k) & is.finite(df$S_k) & df$k > 0 & df$S_k > 0, ]
  if (nrow(df) < 5) return(list(alpha = NA_real_, slope = NA_real_, r2 = NA_real_, window = c(NA, NA)))
  if (is.null(k_min)) k_min <- as.numeric(quantile(df$k, 0.2))
  if (is.null(k_max)) k_max <- as.numeric(quantile(df$k, 0.6))
  dfw <- df[df$k >= k_min & df$k <= k_max, ]
  if (nrow(dfw) < 3) return(list(alpha = NA_real_, slope = NA_real_, r2 = NA_real_, window = c(k_min, k_max)))
  fit <- lm(log(S_k) ~ log(k), data = dfw)
  slope <- as.numeric(coef(fit)[2])
  alpha <- -(slope + d) / 2
  r2 <- summary(fit)$r.squared
  list(alpha = alpha, slope = slope, r2 = r2, window = c(k_min, k_max))
}

# Growth exponent beta from roughness vs time (early-time window)
estimate_beta <- function(times, widths, t_min_frac = 0.1, t_max_frac = 0.5) {
  n <- length(times)
  if (n < 5) return(list(beta = NA_real_, r2 = NA_real_, idx = integer(0)))
  i1 <- max(1L, floor(n * t_min_frac))
  i2 <- min(n,  ceiling(n * t_max_frac))
  idx <- i1:i2
  xs <- log(pmax(1e-9, times[idx]))
  ys <- log(pmax(1e-12, widths[idx]))
  fit <- lm(ys ~ xs)
  list(beta = as.numeric(coef(fit)[2]),
       r2 = summary(fit)$r.squared,
       idx = idx)
}

estimate_z <- function(alpha, beta) {
  if (!is.finite(alpha) || !is.finite(beta) || beta == 0) return(NA_real_)
  alpha / beta
}

# ===================== END ADDITIONS =====================

# Discrete KPZ update
kpz_update <- function(h, dim, nu = 1, lambda = 1, D = 1, dt = 0.01) {
  if (dim == 1) {
    L <- length(h)
    lap <- c(h[L], h[1:(L-1)]) + c(h[2:L], h[1]) - 2*h
    grad <- (c(h[2:L], h[1]) - c(h[L], h[1:(L-1)])) / 2
    h <- h + dt * (nu * lap + (lambda/2) * grad^2) +
      sqrt(2 * D * dt) * rnorm(L)
  } else {
    L <- nrow(h)
    lap <- (rbind(h[L,], h[1:(L-1),]) +
              rbind(h[2:L,], h[1,]) +
              cbind(h[,L], h[,(1:(L-1))]) +
              cbind(h[,2:L], h[,1]) - 4*h)
    gradx <- (cbind(h[,2:L], h[,1]) - cbind(h[,L], h[,(1:(L-1))])) / 2
    grady <- (rbind(h[2:L,], h[1,]) - rbind(h[L,], h[1:(L-1),])) / 2
    h <- h + dt * (nu * lap + (lambda/2) * (gradx^2 + grady^2)) +
      sqrt(2 * D * dt) * matrix(rnorm(L*L), L, L)
  }
  h
}

# Main KPZ runner
run_KPZ <- function(L = 50, T_max = 5000, dim = 1, nu = 1, lambda = 1, D = 1, dt = 0.01) {
  if (dim == 1) {
    h <- rep(0, L)
    W <- numeric(T_max)
    for (t in 1:T_max) {
      h <- kpz_update(h, 1, nu, lambda, D, dt)
      W[t] <- roughness(h)
    }
    df <- data.frame(x = 1:L, h = h)
    p1 <- ggplot(df, aes(x, h)) + geom_line() +
      ggtitle("KPZ 1D Surface") + theme_minimal()
    p2 <- ggplot(data.frame(t = 1:T_max, W = W),
                 aes(t, W)) + geom_line() +
      ggtitle("Roughness W(t)") + theme_minimal()
    print(p1)
    Sys.sleep(15)
    print(p2)
    
    # -------- Exponent summary (1D) --------
    Sk <- structure_factor1d(h)
    fit_a <- fit_alpha_from_Sk(Sk, d = 1)
    times <- 1:T_max
    fit_b <- estimate_beta(times, W)
    z_est <- estimate_z(fit_a$alpha, fit_b$beta)
    cat("\n[KPZ 1D] alpha(spec) =", signif(fit_a$alpha,4),
        " (slope=", signif(fit_a$slope,4), ", R2=", signif(fit_a$r2,4), ", k-window=", paste(signif(fit_a$window,4), collapse=".."), ")\n",
        "[KPZ 1D] beta =", signif(fit_b$beta,4), " (R2=", signif(fit_b$r2,4), ")\n",
        "[KPZ 1D] z = alpha/beta =", signif(z_est,4), "\n", sep="")
    
  } else {
    h <- matrix(0, L, L)
    W <- numeric(T_max)
    for (t in 1:T_max) {
      h <- kpz_update(h, 2, nu, lambda, D, dt)
      W[t] <- roughness(as.vector(h))
    }
    df <- melt(h)
    p1 <- ggplot(df, aes(Var1, Var2, fill = value)) +
      geom_tile() + scale_fill_viridis_c() +
      ggtitle("KPZ 2D Surface") + theme_minimal()
    p2 <- ggplot(data.frame(t = 1:T_max, W = W),
                 aes(t, W)) + geom_line() +
      ggtitle("Roughness W(t)") + theme_minimal()
    
    corr_df <- two_point_corr(h)
    corr_df_radial_avg <- two_point_corr_2d_radial(h)
    
    p3 <- ggplot(corr_df, aes(r, C_r)) +
      geom_line() + ggtitle("Two-point Correlation C(r)") + theme_minimal()
    p4 <- ggplot(corr_df_radial_avg, aes(r, C_r)) +
      geom_line(linewidth = 0.6) +
      labs(title = "KPZ 2D Isotropic Two-point Correlation",
           x = "r (lattice units)", y = "C(r)") +
      theme_minimal()
    
    print(p1)
    Sys.sleep(15)
    print(p2)
    Sys.sleep(15)
    print(p3)
    Sys.sleep(15)
    print(p4)
    
    # -------- Exponent summary (2D) --------
    Sk <- structure_factor2d(h)
    fit_a <- fit_alpha_from_Sk(Sk, d = 2)
    times <- 1:T_max
    fit_b <- estimate_beta(times, W)
    z_est <- estimate_z(fit_a$alpha, fit_b$beta)
    cat("\n[KPZ 2D] alpha(spec) =", signif(fit_a$alpha,4),
        " (slope=", signif(fit_a$slope,4), ", R2=", signif(fit_a$r2,4), ", k-window=", paste(signif(fit_a$window,4), collapse=".."), ")\n",
        "[KPZ 2D] beta =", signif(fit_b$beta,4), " (R2=", signif(fit_b$r2,4), ")\n",
        "[KPZ 2D] z = alpha/beta =", signif(z_est,4), "\n", sep="")
  }
}

# 1D run example
run_KPZ()

# 2D run example
run_KPZ(dim = 2)
