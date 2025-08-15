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
  }
}

# 1D run example
run_KPZ()

# 2D run example
run_KPZ(dim = 2)
