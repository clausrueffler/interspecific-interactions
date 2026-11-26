# parameters
r <- 0.5
e <- 0.05
a <- 0.05
d <- 0.05
k <- 150

# isocline plot
R <- 0:k
Riso <- function(R, r, k, a, h, g, e){
  (((1 + a * h * R) * (k - R) * r) / (a * k - a * g * (k - R) * r))
}

h <- c(0.0001, 0.1, 0.2, 0.3)
g <- c(0.0001, 0.1, 0.2, 0.3, 0.4, 0.5)


par(mfrow=c(6, 4), mar = c(1, 1, 1.5, 0.2))
for (i in seq_along(g)) {
  for (j in seq_along(h)) {
    curve(Riso(x, r, k, a, h[j], g[i], e), from = 0, to = k, type = "l", lty =2, col = "red", ylim = c(0, 15),
          xlab = "", ylab = "", main = paste0("h=", h[j], ", g=", g[i]), xaxt = "n", yaxt = "n")
    abline(a = -1/(a*g[i]), b = (e - d * h[j]) / (d * g[i]), lty = 2, col = "blue") 
  }
}

# equilibrium values
Req <- matrix(NA, nrow = length(g), ncol = length(h), dimnames = list(paste0("g=", g), paste0("h=", h)))
Ceq <- Req

for(i in seq_along(g)) {
  for(j in seq_along(h)) {
    f <- function(R){
      Riso(R, r, k, a, h[j], g[i], e) -
        (-1/(a * g[i]) +
           ((e - d * h[j]) / (d * g[i])) * R)
      }
    xgrid <- seq(0, k, length.out = 5000)
    ygrid <- sapply(xgrid, f)
    idx <- which(
      ygrid[-1] * ygrid[-length(ygrid)] < 0
    )
    if(length(idx) > 0) {
      Req[i, j] <- uniroot( # use uniroot to solve for the intersection point
        f,
        lower = xgrid[idx[1]],
        upper = xgrid[idx[1] + 1]
      )$root
      Ceq[i, j] <- Riso(
        Req[i, j],
        r, k, a,
        h[j], g[i], e
      )
    }
  }
}

round(Req, 3) # matrix of resource equilibrium values
round(Ceq, 3) # matrix of consumer equilibrium values


# adding population dynamics to phase plane
R0 <- 20
C0 <- 2

MaxTime <- 1000
Time <- seq(0, MaxTime, by = 0.1)

par(mfrow=c(6, 4), mar = c(1, 1, 1.5, 0.2))
for (i in 1:6) {
  for (j in 1:4) {
    
    params <- c(r = r, k = k, a = a, h = h[j], g = g[i], e = e, d = d)
    interference.out <- ode(c(R0, C0), Time, predprey.interference, params)
    plot(interference.out[, 2], interference.out[, 3], type = "l",  lty = 1, col = 1,  ylim = c(0, 20), ylab = "", xlab = "",
         main = paste0("h=", h[j], ", g=", g[i]), xaxt = "n", yaxt = "n")
    
    curve(Riso(x, r, k, a, h[j], g[i], e), from = 0, to = k, type = "l", lty =2, col = "red", add = TRUE)
    abline(a = -1/(a*g[i]), b = (e - d * h[j]) / (d * g[i]), lty = 2, col = "blue") 
  }
}