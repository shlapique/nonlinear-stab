library(ggplot2)
library(Ryacas)
library(stringr)

check_device <- function() {
  while (!is.null(dev.list())) Sys.sleep(1)
}

A <- matrix(c(1.1163353, 0.1196577, -0.05982886, 0.99667761), 2)
print("A from 1st lab:")
A
H <- matrix(c(558.84381, -18.29773, -18.29773, 624.81490), 2)
print("H:")
H
print("x0:")
x0 <- ysym(c(0.1, 0.2))
x0

Nmin <- 9
h <- list(t(A) %*% H %*% A)
for (i in 2:Nmin) {
  h[[i]] <- t(A) %*% h[[i - 1]] %*% A
}
print("harr:")
h

psi_yac <- ysym(c("psi0", "psi1"))
psi_yac

fun <- ysym("0")
for (i in 1:Nmin) {
  fun <- fun + (solve(ysym(h[[i]])) %*% psi_yac) / sqrt(psi_yac %*% solve(ysym(h[[i]])) %*% psi_yac)
}
fun
# text <- fun[1]
# print("THEN THE TEXT ITSELF")
# text <- text$yacas_cmd
# print(class(text))
# text


# Matching pattern
pattern <- "[+-]?[0-9]*\\.[0-9]*e?[+-]?[0-9]*"

# Extracting matches
# matches <- grep(pattern, text, value = TRUE, perl = TRUE)
# matches <- sub("[+-]?[0-9]*\\.[0-9]*e?[+-]?[0-9]*", text)
# matches <- str_extract_all(text, "[+-]?[0-9]*\\.[0-9]*e?[+-]?[0-9]*")
# matches <- str_extract_all(text, "\\d+")

# Output
print("ANSWER:")
# print(matches)

# test <- y_eval(ysym("Simplify(-x0[1]/a)"), x0=fun)
# test
sol <- y_eval(ysym("OldSolve({-x0[1]/a == y0[1], -x0[2]/a == y0[2], psi0^2 + psi1^2 == 1}, {a, psi0, psi1})"), y0 = fun, psi0 = psi_yac[1], psi1 = psi_yac[2], x0 = x0)
# sol <- y_eval(ysym("OldSolve({-x0[1]/a == y0[1], -x0[2]/a == y0[2], psi0^2 + psi1^2 == 1}, {a, psi0, psi1})"), y0=fun, psi0=psi_yac[1], psi1=psi_yac[2], x0=x0)
# sol <- y_eval(ysym("OldSolve({(x*y)/2-3*a == 0, x^2/4-(3*a)/2 == 0, 45-(3*x+(3*y)/2) == 0}, {x, y, a})"))
sol
