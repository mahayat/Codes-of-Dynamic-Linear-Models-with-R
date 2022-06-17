### Gibbs sampler for a Bayesian Analysis of the Output gap

gdpGibbs <- function(y, a.theta, b.theta, shape.theta, rate.theta,
                     dV = 1e-7, m0 = rep(0,4),
                     C0 = diag(x=c(rep(1e7,2), rep(1,2))),
                     n.sample = 1, thin = 0, save.states = FALSE)
{
    mod <- dlmModPoly(2, dV = dV, dW = rep(1,2)) +
        dlmModARMA(ar = rep(0,2), sigma2 = 1)
    mod$m0 <- m0
    mod$C0 <- C0
    p <- 4 # dim of state space
    r <- 3 # number of unknown variances
    nobs <- NROW(y)
    if ( is.numeric(thin) && (thin <- as.integer(thin)) >= 0 )
    {
        every <- thin + 1
        mcmc <- n.sample * every
    }
    else
        stop("\"thin\" must be a nonnegative integer")
    ## check hyperpriors for precision(s) of 'theta'
    msg4 <- paste("Either \"a.theta\" and \"b.theta\" or \"shape.theta\"",
                  "and \"rate.theta\" must be specified")
    msg5 <- "Unexpected length of \"shape.theta\" and/or \"rate.theta\""
    msg6 <- "Unexpected length of \"a.theta\" and/or \"b.theta\""
    if (is.null(a.theta))
        if (is.null(shape.theta)) stop(msg4)
        else
            if (is.null(rate.theta)) stop(msg4)
            else
            {
                ## check length of shape.theta and rate.theta
                if (!all(c(length(shape.theta), length(rate.theta)) %in% c(1,r)))
                    warning(msg5)
            }
    else
        if (is.null(b.theta)) stop(msg4)
        else
        {
            if (!all(c(length(a.theta), length(b.theta)) %in% c(1,r)))
                warning(msg6)
            shape.theta <- a.theta^2 / b.theta
            rate.theta <- a.theta / b.theta
        }
    shape.theta <- shape.theta + 0.5 * nobs
    theta <- matrix(0, nobs + 1, p)
    if ( save.states )
        gibbsTheta <- array(0, dim = c(nobs + 1, p, n.sample))
    gibbsPhi <- matrix(0, nrow = n.sample, ncol = 2)
    gibbsVars <- matrix(0, nrow = n.sample, ncol = r)
    AR2support <- function(u)
    {
        ## stationarity region for AR(2) parameters
        (sum(u) < 1) && (diff(u) < 1) && (abs(u[2]) < 1)
    }
    ARfullCond <- function(u)
    {
        ## log full conditional density for AR(2) parameters
        mod$GG[3:4,3] <- u
        -dlmLL(y, mod) + sum(dnorm(u, sd = c(2,1) * 0.33, log=TRUE))
    }
    it.save <- 0
    for (it in 1:mcmc)
    {
        ## generate AR parameters
        mod$GG[3:4,3] <- arms(mod$GG[3:4,3],
                              ARfullCond, AR2support, 1)
        ## generate states - FFBS
        modFilt <- dlmFilter(y, mod, simplify=TRUE)
        theta[] <- dlmBSample(modFilt)
        ## generate W
        theta.center <- theta[-1,-4,drop=FALSE] -
            (theta[-(nobs + 1),,drop=FALSE] %*% t(mod$GG))[,-4]
        SStheta <- drop(sapply( 1 : 3, function(i)
                               crossprod(theta.center[,i])))
        diag(mod$W)[1:3] <-
            1 / rgamma(3, shape = shape.theta,
                       rate = rate.theta + 0.5 * SStheta)
        ## save current iteration, if appropriate
        if ( !(it %% every) )
        {
            it.save <- it.save + 1
            if ( save.states )
                gibbsTheta[,,it.save] <- theta
            gibbsPhi[it.save,] <- mod$GG[3:4,3]
            gibbsVars[it.save,] <- diag(mod$W)[1:3]
        }
    }
    if ( save.states )
        return(list(phi = gibbsPhi, vars = gibbsVars, theta = gibbsTheta))
    else
        return(list(phi = gibbsPhi, vars = gibbsVars))
}
