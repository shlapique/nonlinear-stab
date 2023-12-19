library(ggplot2)
library(expm)

check_device <- function()
{
    while (!is.null(dev.list())) Sys.sleep(1)
}

pplot <- function(M, Color) {
    tmp <- M[1:nrow(M), chull(t(M))]
    data <- as.data.frame(t(cbind(tmp, tmp[,1])))
    colnames(data) <- c("x", "y")
    return(geom_path(data, mapping=aes(x=x, y=y, z=NULL), color=Color))
}

minkowski_sum <- function(mat1, mat2)
{
    temp <- matrix(0, nrow(mat1))
	for(i in 1:ncol(mat1)){
        for(j in 1:ncol(mat2)){
            temp <- cbind(temp, mat1[, i] + mat2[, j])
        }
    }
    temp <- temp[1:nrow(temp), chull(t(temp))]
    return(temp)
}

find_max <- function(mat){
    max_col <- mat[, 1]
    val <- sqrt(mat[1, 1]^2 + mat[2, 1]^2)
    for(i in 2:ncol(mat)){
        tt <- sqrt(mat[1, i]^2 + mat[2, i]^2)
        if(val <= tt){
            max_col <- mat[, i]
            val <- tt
            print(paste("val:", tt))
            print(paste("i: ", i))
        }
    }
    return(max_col)
}

circleFun <- function(center = c(0,0), diameter = 1, npoints = 100){
    r = diameter / 2
    tt <- seq(0,2*pi,length.out = npoints)
    xx <- center[1] + r * cos(tt)
    yy <- center[2] + r * sin(tt)
    return(data.frame(x = xx, y = yy))
}

A <- matrix(c(1.1163353, 0.1196577, -0.05982886, 0.99667761), 2)
print("A from 2nd lab:")
A
U <- matrix(c(0.04113075, 0.04230539, -0.04348004, 0.03760683, 0.04348004, -0.03760683, -0.04113075, -0.04230539), 2)
print("U:")
U

Ad <- A
print("Ad:")
Ad


SV <-eigen(Ad)
SV

lambda <- SV$values
lambda <- sqrt(Re(lambda[1]^2 + Im(lambda[1])^2))
print("labmda:")
lambda

# print("S:")
S <- matrix(c(Re(SV$vectors[, 1]), Im(SV$vectors[, 1])), 2)
S

new_U <- Re(solve(S)%*%U)
print("new_U:")
new_U

R_max <- find_max(new_U)
print("RADIUS:")
R_max <- sqrt(R_max[1]^2 + R_max[2]^2)

# draw a circle
dat <- circleFun(c(0, 0), diameter=2*R_max, npoints=100)
window <- ggplot(dat, aes(x, y)) + geom_path()
###
window <- window + pplot(new_U, "green")

R <- R_max / (lambda - 1)
R

print("H:")
H <- t(solve(S))%*%solve(S)*(1/R^2)
H

# draw an ellipse
x<-seq(-R, R, length=1000)
y<-seq(-R, R, length=1000)
z<-outer(x, y, function(x,y) H[1, 1]*x^2 + H[2, 2]*y^2 + H[1, 2]*2*x*y -1)

df <- data.frame(expand.grid(x = x, y = y), z = c(z)) 
window_el <- ggplot(df, aes(x = x, y = y, z = z)) +
  geom_contour(aes(z = z), breaks = 0, colour="red")

X <- list(-solve(A)%*%U)
X[[1]]
window_el <- window_el + pplot(X[[1]], "purple")
for(i in 2:(100)){
    X[[i]] <- minkowski_sum(solve(A)%*%X[[i-1]], X[[1]])
    window_el <- window_el + pplot(X[[i]], "yellow")
}

X11()
window
check_device()
X11()
window_el
check_device()
