## quantile function for weighted particle clouds
weighted.quantile <- function(x, w, probs)
{
    ## Make sure 'w' is a probability vector
    if ((s <- sum(w)) != 1)
        w <- w / s
    ## Sort 'x' values
    ord <- order(x)
    x <- x[ord]
    w <- w[ord]
    ## Evaluate cdf
    W <- cumsum(w)
    ## Invert cdf
    tmp <- outer(probs, W, "<=")
    n <- length(x)
    quantInd <- apply(tmp, 1, function(x) (1 : n)[x][1])
    ## Return
    ret <- x[quantInd]
    ret[is.na(ret)] <- x[n]
    return(ret)
}
