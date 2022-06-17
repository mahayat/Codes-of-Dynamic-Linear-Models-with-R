dlmFilterDF <- function (y, mod, simplify = FALSE, DF) 
{
    ## storage.mode(y) <- "double"
    mod1 <- mod
    yAttr <- attributes(y)
    ytsp <- tsp(y)
    y <- as.matrix(y)
    timeNames <- dimnames(y)[[1]]
    stateNames <- names(mod$m0)
    m <- rbind(mod$m0, matrix(0, nr = nrow(y), nc = length(mod$m0)))
    a <- matrix(0, nr = nrow(y), nc = length(mod$m0))
    f <- matrix(0, nr = nrow(y), nc = ncol(y))
    U.C <- vector(1 + nrow(y), mode = "list")
    D.C <- matrix(0, 1 + nrow(y), length(mod$m0))
    U.R <- vector(nrow(y), mode = "list")
    D.R <- matrix(0, nrow(y), length(mod$m0))
    U.W <- vector(nrow(y), mode = "list")
    D.W <- matrix(0, nrow(y), length(mod$m0))
    Wliste <- vector(nrow(y), mode = "list")
    P <- vector(nrow(y), mode = "list")
    tmp <- La.svd(mod$V, nu = 0)
    Uv <- t(tmp$vt)
    Dv <- sqrt(tmp$d)
    Dv.inv <- 1/Dv
    Dv.inv[abs(Dv.inv) == Inf] <- 0
    sqrtVinv <- Dv.inv * t(Uv)
    sqrtV <- Dv * Uv
    tmp <- La.svd(mod$C0, nu = 0)
    U.C[[1]] <- t(tmp$vt)
    D.C[1, ] <- sqrt(tmp$d)
    for (i in seq(length = nrow(y))) {
      tF.Vinv <- t(mod$FF) %*% crossprod(sqrtVinv)
      a[i, ] <- mod$GG %*% m[i, ]
      P[[i]] <- mod$GG %*% crossprod(D.C[i,] * t(U.C[[i]])) %*% t(mod$GG)
      Wliste[[i]] <- P[[i]]* ((1-DF)/DF)
      svdW <- La.svd( Wliste[[i]] , nu = 0)
      sqrtW <- sqrt(svdW$d) * svdW$vt
      U.W[[i]] <- t(svdW$vt)
      D.W[i, ] <- sqrt(svdW$d)
      tmp <- La.svd(rbind(D.C[i, ] * t(mod$GG %*% U.C[[i]]), 
            sqrtW), nu = 0)
      U.R[[i]] <- t(tmp$vt)
      D.R[i, ] <- tmp$d
      f[i, ] <- mod$FF %*% a[i, ]
      D.Rinv <- 1/D.R[i, ]
      D.Rinv[abs(D.Rinv) == Inf] <- 0
      tmp <- La.svd(rbind(sqrtVinv %*% mod$FF %*% U.R[[i]], 
      diag(x = D.Rinv, nrow = length(D.Rinv))), nu = 0)
      U.C[[i + 1]] <- U.R[[i]] %*% t(tmp$vt)
      foo <- 1/tmp$d
      foo[abs(foo) == Inf] <- 0
      D.C[i + 1, ] <- foo
      m[i + 1, ] <- a[i, ] + crossprod(D.C[i + 1, ] * t(U.C[[i
          + 1]])) %*% tF.Vinv %*% as.matrix(y[i, ] - f[i,])
    }        
     m <- drop(m)
     a <- drop(a)
     f <- drop(f)
     attributes(f) <- yAttr
     ans <- list(m = m, U.C = U.C, D.C = D.C, a = a, U.R = U.R, 
         D.R = D.R, f = f, U.W=U.W, D.W=D.W)
     ans$m <- drop(ans$m)
     ans$a <- drop(ans$a)
     ans$f <- drop(ans$f)
     attributes(ans$f) <- yAttr
     if (!is.null(ytsp)) {
        tsp(ans$a) <- ytsp
        tsp(ans$m) <- c(ytsp[1] - 1/ytsp[3], ytsp[2:3])
        class(ans$a) <- class(ans$m) <- if (length(mod$m0) > 
            1) 
            c("mts", "ts")
        else "ts"
    }
    if (!(is.null(timeNames) && is.null(stateNames))) {
        dimnames(ans$a) <- list(timeNames, stateNames)
        dimnames(ans$m) <- list(if (is.null(timeNames)) NULL else c("", 
            timeNames), stateNames)
    }
    if (simplify) 
        return(c(mod = list(mod1), ans))
    else {
        attributes(y) <- yAttr
        ans <- c(y = list(y), mod = list(mod1), ans)
        class(ans) <- "dlmFiltered"
        return(ans)
    }
}
