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

find_norm <- function(mat, m){
    # tmp <- matrix(0, nrow(mat))
    tmp <- list()
    for(i in 1:ncol(mat[m])){
        # tmp <- cbind(tmp, abs(mat[1, i]) + abs(mat[2, i]))
        # append(tmp, abs(mat[1, i]) + abs(mat[2, i]))
        tmp[[i]] <- abs(mat[1, i]) + abs(mat[2, i])
    }
    return(max(tmp))
}

A <- matrix(c(1.1163353, -0.05982886, 0.1196577, 0.99667761), 2)
U <- matrix(c(0.04113075, 0.04230539, -0.04348004, 0.03760683, 0.04348004, -0.03760683, -0.04113075, -0.04230539), 2)
U

Ad <- A
Ad

# Ad
# eigen(Ad)

Temp <- solve(Ad)%^%1
Temp
alpha_1 <- max((c(abs(Temp[1, 1]) + abs(Temp[2, 1]), abs(Temp[1, 2]) + abs(Temp[2, 2]))))
alpha_1
m <- 1
while(alpha_1 > 1) {
    m <- m + 1
    Temp <- solve(Ad)%^%m
    alpha_1 <- max((c(abs(Temp[1, 1]) + abs(Temp[2, 1]), abs(Temp[1, 2]) + abs(Temp[2, 2]))))
}
m
alpha_1

X <- list(-solve(A)%*%U)
X[[1]]
window <- ggplot() + pplot(X[[1]], "green")
for(i in 2:(m*20)){
    print(paste("i: ", i))
    print("solve(A):")
    print(A)
    X[[i]] <- minkowski_sum(solve(A)%*%X[[i-1]], X[[1]])
    print(X[[i]])
    window <- window + pplot(X[[i]], "red")
}

print(paste("class:", class(X)))
print("X:")
print(X)

ma <- find_norm(X)
mainf <- max(X)
ma
mainf

X11()
window
check_device()
