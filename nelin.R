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
    convex <- chull(t(M))
    tmp <- M[1:2, convex]
    data <- as.data.frame(t(cbind(tmp, tmp[,1])))
    colnames(data) <- c("x", "f")
    # data <- head(data, -1)
    plot(data, type="b", col=Color)
}

xsum <- function(h, vec) {
    return((solve(h)%*%vec) / c(sqrt(vec%*%(solve(h)%*%vec))))
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

x<-seq(-0.1, 0.1, length=1000)
y<-seq(-0.1, 0.1, length=1000)
z<-outer(x, y, function(x,y) 558.84381*x^2 + 624.81490*y^2 -1)
# z<-outer(x, y, function(x,y) 370.94*x^2 + 352.55*y^2 -1 )

# X11()
# contour(x, y, z, levels=0)
# check_device()

# X11()
# pplot(U, "green")
# check_device()

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
xp1 <- (solve(H)%*%p1) / c(sqrt(p1%*%(solve(H)%*%p1)))
xp1
print("xp2:")
xp2 <- (solve(H)%*%p2) / c(sqrt(p2%*%(solve(H)%*%p2)))
xp2

X <- matrix(c(xp1, xp2, -xp1, -xp2), 2)
X

# X11()
# pplot(X, "green")
# check_device()

a <- xsum(H, X[,2] + X[,1])
a
P <- matrix(c(X, xsum(H, X[,2] + X[,1]), xsum(H, X[,3] + X[,2]), xsum(H, X[,4] + X[,3]),
              xsum(H, X[,4] + X[,1])), 2)
P

# X11()
# pplot(P, "green")
# check_device()

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
X11()
pplot(U_out, "green")
check_device()
