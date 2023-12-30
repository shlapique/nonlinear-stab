library(ggplot2)
library(lpSolve)
library(Ryacas)
library(stringr)

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
U

# 0-controll sets
X <- list(-solve(A)%*%U)
X[[1]]
window1 <- ggplot() + pplot(X[[1]], "purple")
for(i in 2:6){
    X[[i]] <- minkowski_sum(solve(A)%*%X[[i-1]], X[[1]])
    window1 <- window1 + pplot(X[[i]], "red")
}
print("X:")
X

# N_min
Nmin <- 6
x0 <- c(0.1, 0.24)
point <- data.frame(x=x0[1], y=x0[2])
window1 <- window1 + geom_point(point, mapping=aes(x=x, y=y))

# LP Nmin 
objective_coefs <- replicate(ncol(X[[Nmin]]), 1)
print(objective_coefs)
constraints_matrix <- matrix(c(X[[Nmin]][1, ],
                               X[[Nmin]][2, ]), nrow=2, byrow=TRUE)
constraints_rhs <- c(x0[1], x0[2])
constraints_dir <- c("=", "=")
problem_type <- "min"
res <- lp("min", objective_coefs, constraints_matrix, constraints_dir, constraints_rhs)
print("ZLP (N_min):")
res
# LP Nmin-1
objective_coefs <- replicate(ncol(X[[Nmin-1]]), 1)
print(objective_coefs)
constraints_matrix <- matrix(c(X[[Nmin-1]][1, ],
                               X[[Nmin-1]][2, ]), nrow=2, byrow=TRUE)
constraints_rhs <- c(x0[1], x0[2])
constraints_dir <- c("=", "=")
problem_type <- "min"
res2 <- lp("min", objective_coefs, constraints_matrix, constraints_dir, constraints_rhs)
print("ZLP (N_min - 1):")
res2

# LP for trajectory
print("TETS FOR TRAJECTORY")

mu <- (c("mu1"))
for(i in 2:ncol(X[[Nmin-1]]))
{
    mu[i] <- paste0("mu",i)
}
mu <- ysym(mu)
mu
lambda <- ysym(c("l1", "l2", "l3", "l4"))
lambda  
temp <- ysym(X[[Nmin-1]])%*%mu - ysym(U)%*%lambda - A%*%x0
temp

constr1 <- as.numeric(str_extract_all(temp[1], "[+-]?[0-9]*\\.[0-9]*e?[+-]?[0-9]*")[[1]])
constr1_rhs <- -constr1[length(constr1)]
constr1 <- constr1[-length(constr1)]
constr2 <- as.numeric(str_extract_all(temp[2], "[+-]?[0-9]*\\.[0-9]*e?[+-]?[0-9]*")[[1]])
constr2_rhs <- -constr2[length(constr2)]
constr2 <- constr2[-length(constr2)]

constr3 <- c(replicate(ncol(X[[Nmin-1]]), 0))
constr3 <- append(constr3, c(replicate(4, 1)))

objective_coefs <- replicate(ncol(X[[Nmin-1]]), 1)
objective_coefs <- append(objective_coefs, c(replicate(4, 0)))
print(objective_coefs)
constraints_matrix <- matrix(c(constr1, constr2, constr3), nrow=3, byrow=TRUE)
print("CONSTRAINTS MATRIX:")
constraints_matrix

constraints_rhs <- c(constr1_rhs, constr2_rhs, 1)
print("CONSTR RHS:")
constraints_rhs
constraints_dir <- c("=", "=", "=")
res3 <- lp("min", objective_coefs, constraints_matrix, constraints_dir, constraints_rhs)
res3
res3$solution
str(res3)

###########
print("TESTING FUNCTIONS..>")
if(res3$objval < 1)
{
    control <- matrix(c(0, 0), nrow=2, byrow=TRUE)
    print(control)
    trajectory <- matrix(x0, nrow=2, byrow=TRUE)
    print(trajectory)
    print("TEST FOR K")
    out <- get_x_trajectory(X, Nmin, A, U, control, trajectory)
    print(out)
}
out <- data.frame(x=out[1, ], y=out[2, ])
window1 <- window1 + geom_point(out, mapping=aes(x=x, y=y))

print("H:")
H <- t(solve(B))%*%diag(2)%*%solve(B)
H

# ellipse
x<-seq(-0.1, 0.1, length=100)
y<-seq(-0.1, 0.1, length=100)
z<-outer(x, y, function(x,y) H[1, 1]*x^2 + H[2, 2]*y^2 + H[1, 2]*2*x*y -1)
df <- data.frame(expand.grid(x = x, y = y), z = c(z)) 
ellipse <- ggplot(df, aes(x = x, y = y, z = z)) +
  geom_contour(aes(z = z), breaks = 0, colour="green")

# window <- window + pplot(U, "green")

ei <- eigen(H)
ei_vals <- ei$values
ei_vecs <- ei$vectors

print("p1:")
p1 <- Re(ei_vecs[,1])
p1
print("p2:")
p2 <- Re(ei_vecs[,2])
p2
print("p3:")
p3 <- -p1
p3
print("p4:")
p4 <- -p2
p4

print("P:")
P <- matrix(c(p1, p2, p3, p4, p1), 2)
P

print("Xu1:")
Xu1 <- xsum(H, p1)
Xu1
print("Xu2:")
Xu2 <- xsum(H, p2)
Xu2
print("Xu3:")
Xu3 <- xsum(H, p3)
Xu3
print("Xu4:")
Xu4 <- xsum(H, p4)
Xu4

print("pp:")
pp <- matrix(c(p1 + p2, p2 + p3, p3 + p4, p4 + p1), 2)
pp

Xu <- matrix(c(0, 0), 2)
for(i in 1:4)
    Xu <- cbind(Xu, xsum(H, pp[, i]))
print("Xu:")
Xu

# 4 ###
Uin4 <- matrix(c(Xu1, Xu2, Xu3, Xu4, Xu1), 2)
Uout4 <- matrix(c(0, 0), 2)
for(i in 1:(ncol(Uin4)-1))
{
    tmp <- solve(t(P[1:2, i:(i+1)]))%*%c(t(P[,i])%*%Uin4[,i], t(P[,i+1])%*%Uin4[,i+1])
    # print(tmp)
    Uout4 <- cbind(Uout4, tmp)
    # print(Uout4)
}
poly4 <- ellipse + pplot(Uin4, "red") + pplot(Uout4, "blue")
#######

l <- poly_approximations(P, pp, Uin4, Xu)
poly <- list(poly4)
for(i in 1:length(l))
{
    p <- ellipse + pplot(l[[i]][[1]], "red") + pplot(l[[i]][[2]], "blue")
    poly[[i+1]] <- p
}


X11()
window1
check_device()

# X11()
# multiplot(plotlist=poly, cols=2)
# check_device()
