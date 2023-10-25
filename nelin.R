library(ggplot2)
library(Ryacas)
check_device <- function()
{
    while (!is.null(dev.list())) Sys.sleep(1)
}

phi <- function(t, ei_vals, ei_vecs) {
    phi_1 <- exp(Re(ei_vals[1])*t)*(cos(Im(ei_vals[1])*t)*Re(ei_vecs[,1]) - sin(Im(ei_vals[1])*t)*Im(ei_vecs[,1]))
    phi_2 <- exp(Re(ei_vals[2])*t)*(cos(Im(ei_vals[2])*t)*Im(ei_vecs[,2]) + sin(Im(ei_vals[2])*t)*Re(ei_vecs[,2]))
    a <- cbind(phi_1, phi_2)
    return(a)
}

print("A:")
A_tilde <- matrix(c(2*sqrt(2), 2*sqrt(2), -sqrt(2), 0), 2)
delta <- 0.04
A_tilde

print("eigen normalized:")
ei <- eigen(A_tilde)
ei_vals <- ei$values
ei_vecs <- ei$vectors
print("ei_vals:")
ei_vals
print("ei_vecs:")
ei_vecs

print("TEST FOR Im(ei_vecs[1, 1]):")
Im(ei_vecs[1, 1])
print("TEST FOR Im(ei_vecs[1, 2]):")
Im(ei_vecs[1, 2])

print("phi(0): ")
phi(0, ei_vals, ei_vecs)
print("phi(delta): ")
phi(delta, ei_vals, ei_vecs)

print("PHI(delta)PHI^{-1}(0):")
A <- phi(delta, ei_vals, ei_vecs)%*%solve(phi(0, ei_vals, ei_vecs))
print(A)

print("B:")
B <- (A - diag(2))%*%solve(A_tilde)
B
print("V:")
V <- matrix(c(1, 1, -1, 1, 1, -1, -1, -1), 2)
V
print("U:")
U <- B%*%V
U
print("H:")
H <- t(solve(B))%*%diag(2)%*%solve(B)
H

# print("TEST (x1, x2)%*%H%*%t((x1, x2)):")
# x1 <- Sym("x_1")
# x2 <- Sym("x_2")
# ellipse_function <- matrix(c("x1", "x2"), 1)%*%H

show <- yac_str(paste0("Plot2D({x^2 + 2*y + 3})"))
