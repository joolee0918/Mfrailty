# ========================================================= 
# Gauss Legendre Integration Integrate from x1 to x2.  Reference : Numerical Recipes in Fortran
# =========================================================

gauleg.f <- function(n, x1, x2) {
    EPS <- 3e-14
    
    m <- (n + 1)/2
    xm <- 0.5 * (x2 + x1)
    xl <- 0.5 * (x2 - x1)
    
    x <- rep(0, n)
    w <- rep(0, n)
    for (i in 1:m) {
        z <- cos(pi * (i - 0.25)/(n + 0.5))
        
        tol <- 9999
        while (tol > EPS) {
            p1 <- 1
            p2 <- 0
            for (j in 1:n) {
                p3 <- p2
                p2 <- p1
                p1 <- ((2 * j - 1) * z * p2 - (j - 1) * p3)/j
            }
            
            pp <- n * (z * p1 - p2)/(z * z - 1)
            z1 <- z
            z <- z1 - p1/pp
            
            tol <- abs(z - z1)
            if (tol <= EPS) {
                break
            }
        }
        
        x[i] = xm - xl * z
        x[n + 1 - i] = xm + xl * z
        w[i] = (2 * xl)/((1 - z * z) * pp * pp)
        w[n + 1 - i] = w[i]
    }
    
    return(data.frame(location = x, weight = w))
}
