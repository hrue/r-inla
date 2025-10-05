rprec.rw2 <- function(n, cyclic = FALSE, noise = 1e-8) {
  # Define diagonal elements
  main_diag <- rep(6, n)
  off1_diag <- rep(-4, n - 1)
  off2_diag <- rep(1, n - 2)
  # Boundary corrections (non-cyclic)
  main_diag[1:2] <- c(1, 5)
  main_diag[(n - 1):n] <- c(5, 1)
  off1_diag[c(1, n - 1)] <- -2
  # Build sparse matrix
  Q <- INLA:::inla.as.sparse(
    bandSparse(n, k = c(0, 1, 2, -1, -2),
               diagonals = 
                 list(main_diag, off1_diag, off2_diag, off1_diag, off2_diag)
               )
  )
  # Cyclic adjustment
  if (cyclic) {
    # main diagonals
    Q[1, 1]   <- Q[n, n]     <- 6
    Q[2, 2]   <- Q[n-1, n-1] <- 6
    # 1-step neighbours
    Q[1, 2]   <- Q[2, 1]     <- -4
    Q[n-1, n] <- Q[n, n-1]   <- -4
    # 2-step neighbours
    Q[1, 3]   <- Q[3, 1]     <- 1
    Q[n-2, n] <- Q[n, n-2]   <- 1
    # cyclic part
    Q[1, n]   <- Q[n, 1]     <- -4
    Q[1, n-1] <- Q[n-1, 1]   <- 1
    Q[2, n]   <- Q[n, 2]     <- 1
  }
  # Add noise to diagonal
  diag(Q) <- diag(Q) + noise
  return(Q)
}
