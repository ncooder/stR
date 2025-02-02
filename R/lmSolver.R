# Solves system || y - Xb || -> min against b.
# Environment env is used to pass some useful data
# (such as last solution and Cholesky decomposition)
# between different calls of lmSolver.
lmSolver <- function(X,
                     y,
                     type = "Matrix",
                     method = "cholesky",
                     iterControl = list(maxiter = 20, tol = 1e-6),
                     trace = TRUE) {
  if (length(y) > nrow(X)) {
    stop("y is too long in lmSolver...")
  }

  y <- c(y, rep(0, nrow(X) - length(y)))

  y <- as.numeric(y)

  if (type == "Matrix") {
    if (method == "cholesky") {
      result <- .solve.dgC.chol(t(X), y)
      return(result$coef)
    }
    if (method == "qr") {
      b <- solve(qr(X), y)
      return(b)
    }
    stop("Unknown 'method' for type='Matrix' in lmSolver...")
  }

  if (type == "iterative") {
    if (method == "cg") {
      f <- function(z) crossprod(X, as.numeric(X %*% z))
      invm <- 1 / colSums(X^2)
      invf <- function(z) invm * z

      Xty <- crossprod(X, y)
      b0 <- rep(0, length(Xty))

      result <- olscg(FUN = f,
                      y = Xty,
                      b = b0,
                      invFUN = invf,
                      iterControl = iterControl,
                      trace = trace)
      return(result$b)
    }

    if (method == "cg-chol") {
      cholResult <- .solve.dgC.chol(t(X), y)
      eL <- Matrix::expand(cholResult$L)
      L <- eL$L
      Lt <- t(eL$L)
      P <- eL$P
      Pt <- t(eL$P)

      f <- function(z) crossprod(X, as.numeric(X %*% z))
      invf <- function(z) solve(P, solve(Lt, solve(L, solve(Pt, z))))

      Xty <- crossprod(X, y)
      b0 <- rep(0, length(Xty))

      result <- olscg(FUN = f,
                      y = Xty,
                      b = b0,
                      invFUN = invf,
                      iterControl = iterControl,
                      trace = trace)

      return(as.numeric(result$b))
    }

    if (method == "lsmr-chol") {
      cholResult <- .solve.dgC.chol(t(X), y)
      eL <- Matrix::expand(cholResult$L)
      L <- eL$L
      Lt <- t(L)
      P <- eL$P
      Pt <- t(P)

      invL <- function(x) as.numeric(solve(L, x, system = "L"))
      invLt <- function(x) as.numeric(solve(Lt, x, system = "U"))

      invP <- function(x) as.numeric(t(P) %*% x)
      invPt <- function(x) as.numeric(t(Pt) %*% x)

      A <- function(x, k) {
        if (k == 1) {
          return(X %*% invP(invLt(x)))
        } else {
          tmp <- crossprod(X, x)
          tmp <- invPt(tmp)
          return(invL(tmp))
        }
      }

      x0 <- 0
      bVec <- y

      result <- lsmr(A = A,
                     b = bVec,
                     atol = iterControl$tol,
                     btol = iterControl$tol,
                     itnlim = iterControl$maxiter)

      tmp <- invLt(result$x)
      tmp <- invP(tmp)
      return(as.numeric(tmp + x0))
    }

    if (method == "lsmr") {
      D <- sqrt(colSums(X^2))
      invD <- 1 / D

      A <- function(x, k) {
        if (k == 1) {
          return(X %*% (invD * x))
        } else {
          return(invD * crossprod(X, x))
        }
      }

      x0 <- rep(0, ncol(X))
      bVec <- y

      result <- lsmr(A = A,
                     b = bVec,
                     atol = iterControl$tol,
                     btol = iterControl$tol,
                     itnlim = iterControl$maxiter)

      if (trace) {
        cat("\nIterations:", result$itn, "  ")
      }
      return(invD * result$x + x0)
    }

    stop("Unknown 'method' for type='iterative' in lmSolver...")
  }

  stop("Unknown type in lmSolver...")
}

