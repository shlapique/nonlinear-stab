library(ggplot2)

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

pplot <- function(M, Color) {
    tmp <- M[1:nrow(M), chull(t(M))]
    data <- as.data.frame(t(cbind(tmp, tmp[,1])))
    colnames(data) <- c("x", "y")
    return(geom_path(data, mapping=aes(x=x, y=y, z=NULL), color=Color))
}

xsum <- function(h, vec) {
    return((solve(h)%*%vec) / c(sqrt(vec%*%(solve(h)%*%vec))))
}


print("A:")
A_tilde <- matrix(c(2*sqrt(2), 2*sqrt(2), -sqrt(2), 0), 2)
delta <- 0.04
A_tilde

print("eigen (normalized)...")
ei <- eigen(A_tilde)
ei_vals <- ei$values
ei_vecs <- ei$vectors
print("ei_vals:")
ei_vals
print("ei_vecs:")
ei_vecs

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

x<-seq(-0.1, 0.1, length=1000)
y<-seq(-0.1, 0.1, length=1000)
z<-outer(x, y, function(x,y) H[1, 1]*x^2 + H[2, 2]*y^2 + H[1, 2]*2*x*y -1)

df <- data.frame(expand.grid(x = x, y = y), z = c(z)) 
window <- ggplot(df, aes(x = x, y = y, z = z)) +
  geom_contour(aes(z = z), breaks = 0)

# outer 4
window <- window + pplot(U, "green")

ei <- eigen(H)
ei_vals <- ei$values
ei_vecs <- ei$vectors

print("p1:")
p1 <- Re(ei_vecs[,1])
p1
print("p2:")
p2 <- Re(ei_vecs[,2])
p2

print("xp1:")
xp1 <- xsum(H, p1)
xp1
print("xp2:")
xp2 <- xsum(H, p2)
xp2

X <- matrix(c(xp1, xp2, -xp1, -xp2), 2)
X
# inner 4
window <- window + pplot(X, "purple")

P <- matrix(c(X, xsum(H, X[,2] + X[,1]), xsum(H, X[,3] + X[,2]), xsum(H, X[,4] + X[,3]),
              xsum(H, X[,4] + X[,1])), 2)
P
# inner 8
window <- window + pplot(P, "red")

print("p:")
p <- matrix(c(p1, p2, -p1, -p2, p1), 2)
p

print("new U")
U <- cbind(X, X[,1])
U
U_out <- matrix(c(0, 0), 2)
U_out

for(i in 1:(ncol(U)-1))
{
    print(paste("Iteration(i) =", i))
    tmp <- solve(t(P[1:2, i:(i+1)]))%*%c(t(P[,i])%*%U[,i], t(P[,i+1])%*%U[,i+1])
    print(tmp)
    U_out <- cbind(U_out, tmp)
    print(U_out)
}
U_out
window <- window + pplot(U_out, "yellow")

print("U: new")
U <- U[, -ncol(U)]
U <- cbind(U, xsum(H, P[, 1] + P[, 2]))
U

print("TEST U_5")
U_5 <- U[1:nrow(U), chull(t(U))]
U_5 <- cbind(U_5, U_5[, 1])
U_5

print("p: new")
p <- p[, -ncol(p)]
p <- cbind(p, p[, 1] + p[, 2])
p
print("p_5")
p_5 <- p[1:nrow(p), chull(t(p))]
p_5 <- cbind(p_5, p_5[, 1])
p_5

U_out_5 <- matrix(c(0, 0), 2)
U_out_5

for(i in 1:(ncol(U)))
{
    print(paste("Iteration(i) =", i))
    tmp <- solve(t(p_5[1:2, i:(i+1)]))%*%c(t(p_5[,i])%*%U_5[,i], t(p_5[,i+1])%*%U_5[,i+1])
    print(tmp)
    U_out_5 <- cbind(U_out_5, tmp)
    print(U_out_5)
}
U_out_5
# outer 5
window <- window + pplot(U_out_5, "pink")    

# inner 5
window <- window + pplot(U_5[, -ncol(U_5)], "black")

X11()
window
check_device()
