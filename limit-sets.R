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

A <- matrix(c(1.1163353, 0.1196577, -0.05982886, 0.99667761), 2)
# A <- matrix(c(1.105391960, -0.0221048920, -0.0221048920, 1.105391960), 2)
print("A from 2nd lab:")
A
U <- matrix(c(0.04113075, 0.04230539, -0.04348004, 0.03760683, 0.04348004, -0.03760683, -0.04113075, -0.04230539), 2)
print("U:")
U

Ad <- A
print("Ad:")
# Ad <- (A - diag(2))%^%5
Ad


SV <-eigen(Ad)
SV

lambda <- SV$values
lambda

print("S:")
S <- SV$vectors
S

tmpU <- Re(solve(S)%*%U)
print("tmpU:")
tmpU

# u1_max <- max(tmpU[1, ])
# u1_max
# u2_max <- max(tmpU[2, ])
# u2_max

