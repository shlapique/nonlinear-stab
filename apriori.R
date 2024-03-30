library(ggplot2)
library(Ryacas)

check_device <- function()
{
    while (!is.null(dev.list())) Sys.sleep(1)
}

my_f <- function(x, y){
    return(abs(2*x + 3*y)^4 + abs(x + 2*y)^4 -1)
}
circleFun <- function(center = c(0,0), r = 1, npoints = 100){
    tt <- seq(0,2*pi,length.out = npoints)
    xx <- center[1] + r * cos(tt)
    yy <- center[2] + r * sin(tt)
    return(data.frame(x = xx, y = yy))
}

find_points <- function(x, y) {
    for(i in x) {
        for(j in y) {
            if(all.equal(my_f(i, j), 0, 1e-6)){
                print("YESAH")
            } else {
               print("NOOOO")
            }
        }
    }
}

print("A_1:")
A_1 <- matrix(c(-1.3, -1.8, 0.3, 0.2), 2)
A_1

A_2 <- matrix(c(-0.1991, 1.2842, -0.1976, -1.1869), 2)
A_2

print("eigen1...")
ei1 <- eigen(A_1)
abs(ei1$values)

print("eigen2...")
ei2 <- eigen(A_2)
abs(ei2$values)

print("TEST VECTORS")
ei2$vectors

x<-seq(-5, 5, length=100)
y<-seq(-5, 5, length=100)
z<-outer(x, y, )
df <- data.frame(expand.grid(x = x, y = y), z = c(z)) 

ellipse <- ggplot(df, aes(x = x, y = y, z = z)) +
  geom_contour(aes(z = z), breaks = 0, colour="green")


print("s2")
S2 <- matrix(c(Re(ei2$vectors[, 1]), Im(ei2$vectors[, 1])), ncol=2)
S2

# find ellipse equation
y1 <- ysym("y1")
y2 <- ysym("y2")
Y <- ysym(S2) %*% c(y1, y2)

c1 <- as.numeric(yac_str(paste0("Coef(", Y[1],", y1, 1)")))
c1[2] <- as.numeric(yac_str(paste0("Coef(", Y[1],", y2, 1)")))
c2 <- as.numeric(yac_str(paste0("Coef(", Y[2],", y1, 1)")))


x <- seq(-15, 15, length=500)*c1[1] + seq(-15, 15, length=500)*c1[2]
y <- seq(-15, 15, length=500)*c2

z <- outer(x, y, function(x,y) abs(2*x + 3*y)^4 + abs(x + 2*y)^4 -1)
df <- data.frame(expand.grid(x = x, y = y), z = c(z)) 
el <- ggplot(df, aes(x = x, y = y, z = z)) + geom_contour(aes(z = z), breaks = 0, colour="red")

# to extract data from implicit ggplot methods
pg_data_x <- ggplot_build(el)$data[[1]]["x"]
pg_data_y <- ggplot_build(el)$data[[1]]["y"]
el_data <- cbind(pg_data_x, pg_data_y)

maxim <- max(sqrt(abs(el_data[, 1])^2 + abs(el_data[, 2])^2))
minim <- min(sqrt(abs(el_data[, 1])^2 + abs(el_data[, 2])^2))

circ_max <- geom_path(circleFun(c(0, 0), r=maxim, npoints=500), mapping=aes(x, y, z=NULL), color="purple")
circ_min <- geom_path(circleFun(c(0, 0), r=minim, npoints=500), mapping=aes(x, y, z=NULL), color="green")

el + circ_max + circ_min
