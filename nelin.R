library(ggplot2)

NullSpace <- function (A) {
  m <- dim(A)[1]; n <- dim(A)[2]
  ## QR factorization and rank detection
  QR <- base::qr.default(A)
  r <- QR$rank
  ## cases 2 to 4
  if ((r < min(m, n)) || (m < n)) {
    R <- QR$qr[1:r, , drop = FALSE]
    P <- QR$pivot
    F <- R[, (r + 1):n, drop = FALSE]
    I <- base::diag(1, n - r)
    B <- -1.0 * base::backsolve(R, F, r)
    Y <- base::rbind(B, I)
    X <- Y[base::order(P), , drop = FALSE]
    return(X)
    }
  ## case 1
  return(base::matrix(0, n, 1))
}

eigen_vectors <- function(A, eigen_values){
    ans <- NULL
    for(i in eigen_values)
    {
        print(A - diag(i, nrow(A)))
        cbind(ans, NullSpace(A - diag(i, nrow(A))))
    }
    return(ans)
}

phi <- function(t){
    phi_1 <- exp(Re(ei_vals[1])*t)*(cos(Im(ei_vals[1])*t)*Re(ei_vecs[,1]) - sin(Im(ei_vals[1])*t)*Im(ei_vecs[,1]))
    phi_2 <- exp(Re(ei_vals[2])*t)*(cos(Im(ei_vals[2])*t)*Re(ei_vecs[,2]) - sin(Im(ei_vals[2])*t)*Im(ei_vecs[,2]))
    a <- cbind(phi_1, phi_2)
    return(a)
}

A_tilde <- cbind(c(2*sqrt(2), 2*sqrt(2)), c(-sqrt(2), 0))
A_tilde

print(eigen(A_tilde))
ei <- eigen(A_tilde)
ei_vals <- ei$values
ei_vecs <- eigen_vectors(A_tilde, ei_vals)
print("EIVAEEKAJ:FKJA:KFKJ:A")
print(ei_vecs)

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
