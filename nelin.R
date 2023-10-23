library(ggplot2)

phi <- function(t) {
    phi_1 <- exp(Re(ei_vals[1])*t)*(cos(Im(ei_vals[1])*t)*Re(ei_vecs[,1]) - sin(Im(ei_vals[1])*t)*Im(ei_vecs[,1]))
    phi_2 <- exp(Re(ei_vals[2])*t)*(cos(Im(ei_vals[2])*t)*Re(ei_vecs[,2]) - sin(Im(ei_vals[2])*t)*Im(ei_vecs[,2]))
    a <- cbind(phi_1, phi_2)
    return(a)
}

print("A:")
A_tilde <- matrix(c(2*sqrt(2), 2*sqrt(2), -sqrt(2), 0), 2)
A_tilde

print("eigen normalized:")
ei <- eigen(A_tilde)
ei_vals <- ei$values
ei_vecs <- ei$vectors
ei_vals
ei_vecs

print("CHECK :")
A_check <- ei_vecs%*%diag(ei_vals)%*%solve(ei_vecs)
print(A_check)

print(Im(ei_vals[1]))
ei_vecs[1]

print("vecs")
print(ei_vecs)

print("phi(0): ")
phi(0)
print("phi(0.04): ")
phi(0.04)

A <- phi(0.04)%*%solve(phi(0))
print(A)
