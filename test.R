library(ggplot2)
library(expm)
library(animation)

check_device <- function() {
  while (!is.null(dev.list())) Sys.sleep(1)
}

pplot <- function(M, Color) {
  tmp <- M[1:nrow(M), chull(t(M))]
  data <- as.data.frame(t(cbind(tmp, tmp[, 1])))
  colnames(data) <- c("x", "y")
  return(geom_path(data, mapping = aes(x = x, y = y, z = NULL), color = Color))
}

minkowski_sum <- function(mat1, mat2) {
  temp <- matrix(0, nrow(mat1))
  for (i in 1:ncol(mat1)) {
    for (j in 1:ncol(mat2)) {
      temp <- cbind(temp, mat1[, i] + mat2[, j])
    }
  }
  temp <- temp[1:nrow(temp), chull(t(temp))]
  return(temp)
}

find_norm <- function(mat, m) {
  tmp <- list()
  tt <- mat[[m]]
  for (i in 1:ncol(mat[[m]])) {
    tmp[[i]] <- abs(tt[1, i]) + abs(tt[2, i])
  }
  return(max(unlist(tmp)))
}

Rn <- function(alpha, m, N) {
  return(alpha^N * m / (1 - alpha))
}

# X,Y,Z -- lists of matrices
# N -- number of plots in gif
animate <- function(X, Y, Z, N) {
  p <- ggplot() +
    pplot(X[[1]], "green")
  saveGIF(
    {
      for (i in 1:N) {
        p <- p + pplot(X[[i]], "red") + pplot(Y[[i]], "coral") + pplot(Z[[i]], "darkblue")
        print(p)
        Sys.sleep(1) # Pause for 1 sec between each plot
      }
    },
    movie.name = "animation.gif",
    interval = 0.5,
    ani.width = 800,
    ani.height = 600
  )
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

Temp <- solve(Ad) %^% 1
print("Ad^-1:")
Temp

alpha_1 <- max((c(abs(Temp[1, 1]) + abs(Temp[2, 1]), abs(Temp[1, 2]) + abs(Temp[2, 2]))))
alpha_1
m <- 1
while (alpha_1 > 0.8) {
  m <- m + 1
  Temp <- solve(Ad) %^% m
  alpha_1 <- max((c(abs(Temp[1, 1]) + abs(Temp[2, 1]), abs(Temp[1, 2]) + abs(Temp[2, 2]))))
}
m
alpha_1
print(paste("m for alpha_1:", m, " | alpha_1:", alpha_1))

# compute alpha_inf
Temp22 <- solve(Ad) %^% 1
Temp22
alpha_inf <- max((c(abs(Temp22[1, 1]) + abs(Temp22[1, 2]), abs(Temp22[2, 1]) + abs(Temp22[2, 2]))))
alpha_inf
mtmp <- 1
while (alpha_inf > 0.8) {
  mtmp <- mtmp + 1
  Temp22 <- solve(Ad) %^% mtmp
  alpha_inf <- max((c(abs(Temp22[1, 1]) + abs(Temp22[1, 2]), abs(Temp22[2, 1]) + abs(Temp22[2, 2]))))
}
mtmp
alpha_inf
print(paste("m for alpha_inf:", mtmp, " | alpha_inf:", alpha_inf))

X <- list(-solve(A) %*% U)
X[[1]]
window <- ggplot() +
  pplot(X[[1]], "green")
for (i in 2:(m * 20)) {
  X[[i]] <- minkowski_sum(solve(A) %*% X[[i - 1]], X[[1]])
  window <- window + pplot(X[[i]], "red")
}

ma <- find_norm(X, m)
mainf <- max(unlist(X))
ma
mainf

N <- 10
r1 <- Rn(alpha_1, ma, N)
rinf <- Rn(alpha_inf, mainf, N)
print("R1 and Rinf")
r1
rinf

Ball <- matrix(c(1, 0, 0, 1, -1, 0, 0, -1), 2)
Ball_inf <- matrix(c(1, 1, -1, 1, -1, -1, 1, -1), 2)
print("BALL AND BALL_INF:")
Ball
Ball_inf

Num <- 20
xest <- list()
for (i in 1:Num) {
  xest[[i]] <- minkowski_sum(X[[i * m]], Rn(alpha_1, ma, i) * Ball)
  window <- window + pplot(xest[[i]], "coral")
}

xest_inf <- list()
for (i in 1:Num) {
  xest_inf[[i]] <- minkowski_sum(X[[i * m]], Rn(alpha_inf, mainf, i) * Ball_inf)
  window <- window + pplot(xest_inf[[i]], "darkblue")
}

# animation
# generates a gif
animate(X, xest, xest_inf, Num)

X11()
window
check_device()
