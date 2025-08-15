# Restricted Solid on Solid (RSOS) model in 1 and 2D

# RSOS_surface_growth.R
library(ggplot2)
library(reshape2)

# Compute roughness
roughness <- function(h) {
  mean_val <- mean(h)
  sqrt(mean((h - mean_val)^2))
}

# Two-point correlation function for 2D - simplified averaged 1D axis
two_point_corr <- function(h) {
  L <- nrow(h)
  corr <- numeric(L %/% 2)
  for (r in 1:(L %/% 2)) {
    diffs <- (h - cbind(h[, (r+1):L], h[, 1:r]))
    corr[r] <- mean(diffs^2)
  }
  data.frame(r = 1:(L %/% 2), C_r = corr)
}

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


# RSOS update rule
rsos_update <- function(h, dim) {
  if (dim == 1) {
    i <- sample(1:length(h), 1)
    left <- h[ifelse(i == 1, length(h), i-1)]
    right <- h[ifelse(i == length(h), 1, i+1)]
    if (abs(h[i] - left) < 2 && abs(h[i] - right) < 2) {
      h[i] <- h[i] + sample(c(-1, 1), 1)
    }
  } else {
    L <- nrow(h)
    i <- sample(1:L, 1)
    j <- sample(1:L, 1)
    neigh <- c(h[ifelse(i == 1, L, i-1), j],
               h[ifelse(i == L, 1, i+1), j],
               h[i, ifelse(j == 1, L, j-1)],
               h[i, ifelse(j == L, 1, j+1)])
    if (all(abs(h[i, j] - neigh) < 2)) {
      h[i, j] <- h[i, j] + sample(c(-1, 1), 1)
    }
  }
  h
}

# Main RSOS runner
run_RSOS <- function(L = 50, T_max = 5000, dim = 1) {
  if (dim == 1) {
    h <- rep(0, L)
    W <- numeric(T_max)
    for (t in 1:T_max) {
      h <- rsos_update(h, 1)
      W[t] <- roughness(h)
    }
    df <- data.frame(x = 1:L, h = h)
    p1 <- ggplot(df, aes(x, h)) + geom_line() +
      ggtitle("RSOS 1D Surface") + theme_minimal()
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
      h <- rsos_update(h, 2)
      W[t] <- roughness(as.vector(h))
    }
    df <- melt(h)
    
    p1 <- ggplot(df, aes(Var1, Var2, fill = value)) +
      geom_tile() + scale_fill_viridis_c() +
      ggtitle("RSOS 2D Surface") + theme_minimal()
    p2 <- ggplot(data.frame(t = 1:T_max, W = W),
                 aes(t, W)) + geom_line() +
      ggtitle("Roughness W(t)") + theme_minimal()
    
    corr_df <- two_point_corr(h)
    p3 <- ggplot(corr_df, aes(r, C_r)) +
      geom_line() + ggtitle("Two-point Correlation C(r)") + theme_minimal()
    
    #2D radially averaged correlation function
    corr_df_radial_avg <- two_point_corr_2d_radial(h)
    
    p4 <- ggplot(corr_df_radial_avg, aes(r, C_r)) +
      geom_line(linewidth = 0.6) +
      labs(title = "RSOS 2D Isotropic Two-point Correlation",
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

# run it
run_RSOS() #default params

# 2D RSOS
run_RSOS(dim = 2) # 2D 
