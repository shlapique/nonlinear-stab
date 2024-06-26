library(ggplot2)
library(lpSolve)
library(Ryacas)
library(stringr)
library(nloptr)

check_device <- function()
{
    while (!is.null(dev.list())) Sys.sleep(1)
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }
 if (numPlots==1) {
    print(plots[[1]])
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

phi <- function(t, ei_vals, ei_vecs) {
    phi_1 <- exp(Re(ei_vals[1])*t)*(cos(Im(ei_vals[1])*t)*Re(ei_vecs[,1]) - sin(Im(ei_vals[1])*t)*Im(ei_vecs[,1]))
    phi_2 <- exp(Re(ei_vals[2])*t)*(cos(Im(ei_vals[2])*t)*Im(ei_vecs[,2]) + sin(Im(ei_vals[2])*t)*Re(ei_vecs[,2]))
    a <- cbind(phi_1, phi_2)
    return(a)
}

pplot <- function(M, Color) {
    tmp <- M[, chull(t(M))]
    data <- as.data.frame(t(cbind(tmp, tmp[,1])))
    colnames(data) <- c("x", "y")
    return(geom_path(data, mapping=aes(x=x, y=y, z=NULL), color=Color))
}

xsum <- function(h, vec) {
    return((solve(h)%*%vec) / c(sqrt(vec%*%(solve(h)%*%vec))))
}

minkowski_sum <- function(mat1, mat2)
{
    temp <- matrix(0, nrow(mat1))
	for(i in 1:ncol(mat1)){
        for(j in 1:ncol(mat2)){
            temp <- cbind(temp, mat1[, i] + mat2[, j])
        }
    }
    temp <- temp[, chull(t(temp))]
    return(temp)
}

# bulids polyhedral approximations m in {4, 5, 6, 7, 8}
poly_approximations <- function(P, pp, Uin_prev, Xu)
{
    a <- 5 # bc we have 4p fig before
    b <- 8
    plots <- list()
    for(j in 1:(b-a+1))
    {
        P[, ncol(P)] <- pp[, j]
        Uin <- cbind(Uin_prev[, -ncol(Uin_prev)], Xu[, j+1])
        order <- chull(t(Uin))
        Uin <- Uin[, order]
        Uin <- cbind(Uin, Uin[, 1])
        P <- P[, order]
        P <- cbind(P, P[, 1])
        Uout <- matrix(c(0, 0), 2)
        for(i in 1:(ncol(Uin)-1))
        {
            tmp <- solve(t(P[1:2, i:(i+1)]))%*%c(t(P[,i])%*%Uin[,i], t(P[,i+1])%*%Uin[,i+1])
            # print(tmp)
            Uout <- cbind(Uout, tmp)
            # print(Uout)
        }
        plots[[j]] <- list(Uin, Uout)
        Uin_prev <- Uin
    }
    return(plots)
}

get_x_trajectory <- function(X, Nmin, A, U, control, trajectory)
{
    ans <- list()
    lambda <- ysym(c("l1", "l2", "l3", "l4"))
    for(k in 1:(Nmin-1)) 
    {
        print("PRINT FROM FUNC FOR TESTING trajectory:")
        print(paste("K=", k))
        print(trajectory[, k])
        mu <- (c("mu1"))
        for(i in 2:ncol(X[[Nmin-k]]))
        {
            mu[i] <- paste0("mu", i)
        }
        mu <- ysym(mu)
        print(mu)
        temp <- ysym(X[[Nmin-k]])%*%mu - ysym(U)%*%lambda - A%*%trajectory[, k]
        print(temp)
        constr1 <- as.numeric(str_extract_all(temp[1], "[+-]?[0-9]*\\.[0-9]*e?[+-]?[0-9]*")[[1]])
        constr1_rhs <- -constr1[length(constr1)]
        constr1 <- constr1[-length(constr1)]
        constr2 <- as.numeric(str_extract_all(temp[2], "[+-]?[0-9]*\\.[0-9]*e?[+-]?[0-9]*")[[1]])
        constr2_rhs <- -constr2[length(constr2)]
        constr2 <- constr2[-length(constr2)]

        constr3 <- c(replicate(ncol(X[[Nmin-k]]), 0))
        constr3 <- append(constr3, c(replicate(4, 1)))

        objective_coefs <- replicate(ncol(X[[Nmin-k]]), 1)
        objective_coefs <- append(objective_coefs, c(replicate(4, 0)))
        print(objective_coefs)
        constraints_matrix <- matrix(c(constr1, constr2, constr3), nrow=3, byrow=TRUE)
        print("CONSTRAINTS MATRIX:")
        print(constraints_matrix)

        constraints_rhs <- c(constr1_rhs, constr2_rhs, 1)
        print("CONSTR RHS:")
        print(constraints_rhs)
        constraints_dir <- c("=", "=", "=")
        res <- lp("min", objective_coefs, constraints_matrix, constraints_dir, constraints_rhs)
        print(paste("K=", k))
        print(res$solution)
        # FIXME
        sol <- res$solution
        sol <- c(sol[!sol == 0], 0)
        # FIXME
        print(sol)
        control <- cbind(control, U%*%sol)
        print(control)
        trajectory <- cbind(trajectory, A%*%trajectory[, k] + control[, k + 1])
        print(control)
        print(trajectory)
    }
    return(trajectory)
}

Norm <- function(vec, q) {
    s <- sum(abs(vec)^q)
    return(s^(1/q))
}

# solve nonlinear system of equations for maximization problem
optim_el <- function(x1_max, x2_max, P_rot, alpha, rr) {
    s1 <- 100
    s2 <- 100
    mesh_x_1 <- seq(0, x1_max, length=s1)
    mesh_x_2 <- seq(0, x2_max, length=s2)

    # mesh_bulok <- matrix(list(), ncol=s1, nrow=s2)
    mesh_bulok <- list()
    # m <- matrix(list(), ncol=s1, nrow=s2)
    for(i in 1:length(mesh_x_1)) {
        for(j in 1:length(mesh_x_2)) {
            # check constraints for every node
            tmp <- TRUE
            for(nn in 1:ncol(P_rot)) {
                tmp <- Norm(diag(c(mesh_x_1[i], mesh_x_2[j]))%*%P_rot[, nn], rr) <= alpha[1]
                if(tmp == FALSE) {
                    break
                }
            }
            if(tmp == TRUE) {
                mesh_bulok <- cbind(mesh_bulok, c(mesh_x_1[i], mesh_x_2[j],
                                                   mesh_x_1[i] * mesh_x_2[j]))
            }
        }
    }
    mmm <- 0
    a <- c(0, 0) # for a1 a2
    for(l in 1:ncol(mesh_bulok)) {

        if(mesh_bulok[[3, l]][1] >= mmm) {
            mmm <- mesh_bulok[[3, l]][1]
            a <- c(mesh_bulok[[1, l]][1], mesh_bulok[[2, l]][1])
        }
    }
    return(a)
}

# mera hyper-ellips 
# a_ast -- list()
Mera <- function(a_ast, r_i) {
    n <- 2
    return(a_ast[[1]] * a_ast[[2]] * ((2 * gamma(1/r_i + 1))^n) / (gamma(n/r_i + 1)))
}

print("A_tilde:")
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

print("A = PHI(delta)PHI^{-1}(0):")
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
U <- U[, chull(t(U))]
# U <- matrix(c(-1, 1, 3, 3, 1, -1, -3, -3), 2)

newU <- cbind(U[, ncol(U)], U[, 2:ncol(U)-1])
P <- matrix(c(0, 1, -1, 0), 2)%*%(U - newU)
alpha <- matrix(replicate(ncol(P), 0))
for(i in 1:ncol(P)) {
    alpha[i] <- P[, i]%*%U[, i]
}

x1_max <- max(U[1, ])
x2_max <- max(U[2, ])
N1 <- 100
N2 <- 100
I_11 <- 0
I_22 <- 0
I_12 <- 0
for(i in 1:N1) {
    for(j in 1:N2) {
        x <- matrix(c(i * 2 * x1_max / N1 - x1_max, j * 2 * x2_max / N2 -
                      x2_max), 2)
        if(all(t(P)%*%x - alpha < 0) == TRUE) {
            I_11 <- I_11 + x[2]^2 * 2 * x1_max * 2 * x2_max / (N1 * N2)
            I_22 <- I_22 + x[1]^2 * 2 * x1_max * 2 * x2_max / (N1 * N2)
            I_12 <- I_12 - x[1] * x[2] * 2 * x1_max * 2 * x2_max / (N1 * N2)
        }
    }
}

J <- matrix(c(I_11, I_12, I_12, I_22), 2)
lambda <- eigen(J)$values
S <- eigen(J)$vectors

U_rot <- solve(S)%*%U

window <- ggplot() + pplot(U, "green") + pplot(U_rot, "purple")

r_arr <- c(8/7, 6/5, 4/3, 2, 4, 6, 8)
P_rot <- solve(S)%*%P

constraints_matrix <- matrix()
# ||diag(a)*p^k||^q <= alpha_k 

aas <- list()
for(ri in 1:length(r_arr)) {
    aas <- cbind(aas, optim_el(x1_max, x2_max, P_rot, alpha, r_arr[ri]))
}

# compute the best `mera` for hyper-ellipse
best_mera <- 0
best_mera_index <- 0
for(i in 1:ncol(aas)) {
    if(Mera(aas[, i], r_arr[i]) > best_mera) {
        best_mera_index <- i
        best_mera <- Mera(aas[, i], r_arr[i])
    }
}

# find ellipse equation
y1 <- ysym("y1")
y2 <- ysym("y2")
Y <- ysym(S) %*% c(y1, y2)

x <- seq(-15, 15, length=500)
y <- seq(-15, 15, length=500)

z <- outer(x, y, function(x, y) abs(x / aas[, best_mera_index][[1]])^r_arr[best_mera_index] + abs(y / aas[, best_mera_index][[2]])^r_arr[best_mera_index] -1)

df <- data.frame(expand.grid(x = x, y = y), z = c(z)) 
el <- ggplot(df, aes(x = x, y = y, z = z)) + geom_contour(aes(z = z), breaks = 0, colour="red")

window2 <- el + pplot(U_rot, "purple")

