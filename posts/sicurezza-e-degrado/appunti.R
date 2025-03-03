library(BiasedUrn)

data |> 
  select(popolazione_residente_in_eta_0_14_anni,
         popolazione_residente_in_eta_15_64_anni,
         popolazione_residente_over64,
         cittadini_residenti,
         stranieri_residenti) |> 
  slice_sample(n = 1)

z1 <- 7894
z2 <- 34515
z3 <- 8016
M <- 40622
N <- 9803

z <- c(z1, z2, z3)

sim <- matrix(nrow = 10000, ncol = 6)

colnames(sim) <- str_c(rep(c( "p", "q"), c(3, 3)), rep(1:3, 2))

LP <- double(10000)

p <- sim[1, 1:3] <- z / (sum(M + N))
q <- sim[1, 4:6] <- z / (sum(M + N))

x <- round(p * M)
y <- z - x

LP[1] <- ddirichlet(p, rep(1, 3), log = TRUE) +
  ddirichlet(q, rep(1, 3), log = TRUE) + 
  dmultinom(x, M, p, log = TRUE) +
  dmultinom(y, N, q, log = TRUE)

accept <- double(10000)


for(i in 2:10000) {
  
  x_prop <- y_prop <- double(3)
  
  prop <- rMFNCHypergeo(1, c(M, N), z[1], c(p[1] * p[3], q[1] * q[3]))
  
  x_prop[1] <- prop[1]
  y_prop[1] <- prop[2]
  
  prop <- rMFNCHypergeo(1, c(M - x_prop[1], N + x_prop[1] - z[1]), z[2], c(p[2] * p[3], q[2] * q[3]))
  
  x_prop[2] <- prop[1]
  y_prop[2] <- prop[2]
  
  x_prop[3] <- M - sum(x_prop)
  y_prop[3] <- N - sum(y_prop)
  
  log_obj <- dmultinom(x_prop, M, p, log = TRUE) + dmultinom(y_prop, N, q, log = TRUE) - 
    dmultinom(x, M, p, log = TRUE) - dmultinom(y, N, q, log = TRUE)
  
  log_proposal <- log(dMFNCHypergeo(c(x[1], y[1]), c(M, N), z[1], c(p[1] * p[3], q[1] * q[3]))) +
    log(dMFNCHypergeo(c(x[2], y[2]), c(M - x[1], N + x[1] - z[1]), z[2], c(p[2] * p[3], q[2] * q[3]))) -
    log(dMFNCHypergeo(c(x_prop[1], y_prop[1]), c(M, N), z[1], c(p[1] * p[3], q[1] * q[3]))) -
    log(dMFNCHypergeo(c(x_prop[2], y_prop[2]), c(M - x_prop[1], N + x_prop[1] - z[1]), z[2], c(p[2] * p[3], q[2] * q[3])))
  
  
  u <- runif(1)
  
  if(log(u) <= log_obj + log_proposal) {
    
    x <- x_prop
    y <- y_prop
    
    accept[i] <- 1
    
  }
  
  
  p <- sim[i, 1:3] <- as.vector(rdirichlet(1, x + 1))
  q <- sim[i, 4:6] <- as.vector(rdirichlet(1, y + 1))
  
  LP[i] <- ddirichlet(p, rep(1, 3), log = TRUE) +
    ddirichlet(q, rep(1, 3), log = TRUE) + 
    dmultinom(x, M, p, log = TRUE) +
    dmultinom(y, N, q, log = TRUE)
  
  if(i %% 100 == 0) {
    
    print(paste(i %/% 100, "%"))
  }
  
}

LML <- LML(theta = sim[, -c(3, 6)], LL = LP)

list(
  posterior = sim,
  LML = LML
)


sim[, -c(3, 6)] |> 
  acf()

sim[, 1] |> 
  plot(type = 'l')

mean(accept)
