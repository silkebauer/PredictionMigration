
# the model for investigating how predictability along a migratory route
# affects migratory behaviour
#
# codes for the definition of the model

# functions {{{
## model {{{
f.norm <- function(x, m, v) {#{{{
  exp(-(x-m)^2/(2*v))
}#}}}


read.param <- function(file) {#{{{
  p <- scan(file=file, what=character(), sep="\n", quiet=TRUE)
  pp <- strsplit(p, "\\s+")
  np <- unlist(lapply(pp, function(z) z[1]))
  pp <- lapply(pp, function(z) as.numeric(z[-1]))
  names(pp) <- np
  # rescaling to grid time
  pp$beta <- pp$beta/pp$De
  pp$R.wide <- pp$R.wide/pp$De
  pp$tau <- pp$tau/pp$De
  pp$sigma <- pp$sigma/pp$De
  pp$mort.max <- pp$mort.max - pp$mort.min
  pp
}#}}}


init.param <- function(P=NULL, ...) {#{{{
  # create a structure which contains the input values for the program
  # 'P <- init.param()' creates P with default values
  # 'P <- init.param(P, De=0.2)' change the existing value of De in P
  n.valt <- names(valt <- list(...));
  # list of real timed variables
  v.real <- c("beta", "R.wide", "tau", "sigma")
  if(is.null(P)) { # P is not given, create it with default values
    K = 10  # max site index
    P <- list(
        De = 0.1,   # \Delta, grid step
        N = 1000, # max time in time grid
        K = K,
        Q = 100,  # max grid value for the advancement of spring
        tau = seq(1, 80, length.out=(K+1)), # the average time spring starts at site k
        sigma = rep(10, (K+1)), # the K+1 number sigma values for the sites
        rho = c(NA, rep(0.3, K)),  # the K number rho values for the sites k 
        alpha = 0.01, # slope of the sigmoid function at the inflection point
        beta = 0, # position of the inflection point
        mort.max = 0.05, # maximal mortality
        mort.min = 0.005, # minimal mortality
        R.wide = 10, # the spread of the breeding season at site K
        R.future = 0, # the future reproductive value
        clutch.size = 1, # the value of current brood
        last=NA  # just for programing convinience
        )
    # rescaling to grid time
    P$beta <- P$beta/P$De
    P$R.wide <- P$R.wide/P$De
    P$tau <- P$tau/P$De
    P$sigma <- P$sigma/P$De
    P$mort.max <- P$mort.max - P$mort.min
  }
  # be carefull when setting individual parameters, NO checking here
  mort.max.set <- FALSE
  for(i in n.valt) {
    if(i %in% names(P)) {
      P[[i]] <- valt[[i]];
      if(i %in% v.real) {P[[i]] <- P[[i]]/P$De}
      if(i == "mort.max") mort.max.set <- TRUE
    } else {
      stop(paste("Setting unused variable: ",i));
    }
  }
  if(mort.max.set) P$mort.max <- P$mort.max - P$mort.min
  if(length(P$tau) != P$K+1 || length(P$sigma) != P$K+1 ||
      length(P$rho) != P$K+1) {
    cat("length of tau:", length(P$tau), "\n")
    cat("length of sigma:", length(P$sigma), "\n")
    cat("length of rho:", length(P$rho), "\n")
    cat("these lengthes should be", P$K+1, "\n")
    stop("Not enough values for 'tau', 'sigma' or 'rho' are given!")
  }
  P
}#}}}


fill.Ak <- function(i, j, k, v, P=P) {#{{{
  # calculate a_{i,j}(k)
  m <- i * P$De * P$rho[k] * (P$sigma[k] / P$sigma[k-1])
  x <- j * P$De
  #cat("m:", m, "x:", x, "v:", v, "\n")
  f.norm(x, m, v)
}#}}}


init.A <- function(P) {#{{{
  # initialise the transition matricies which are actually an array
  ii <- jj <- seq(-P$Q, P$Q, 1)
  kk <- seq(0, P$K, 1)
  A <- array(NA, dim=c(length(ii), length(jj), length(kk)),
      dimnames=list(i=ii, j=jj, k=kk))

  v <- P$sigma^2 * (1 - P$rho^2)
  for(k in 2:length(kk)) {
    z <- outer(ii, jj, fill.Ak, k=k, v=v[k], P=P) 
    z <- apply(z, 1, function(y) y/sum(y))
    A[,,k] <- t(z)
  }
  invisible(A)
}#}}}


surv.breeding <- function(y, tt, P) {
	# this function calculates the probability of survival from tt to T on
	# the breeding ground
	1 - M(y + tt, P$K+1, P)
}

R <- function(y, tt, P) {#{{{
  # the terminal reward as a function of 'y' and 't'

  # original version: no effect of mortality on site K
  f.norm(y + tt, P$tau[P$K+1], P$R.wide) 

  # mortality on site K affects the rewards
  #(P$R.future + f.norm(y + tt, P$tau[P$K+1], P$R.wide)) * (1 - M(y + tt,
                                                               #P$K+1,
                                                               #P))
}#}}}

M <- function(ty, k, P=P) {#{{{
  # mortality function, it is a sigmoid
  # it starts from mort.max and decreases to mort.min, the inflection
  # point is at beta and the slope is alpha
  # 'ty' is t+y
  x <- (ty - P$tau[k])
  P$mort.min[k] + P$mort.max[k]/(1 + exp(-P$alpha[k]*(P$beta[k] - x)))
}#}}}

init.RV <- function(P) {#{{{
  # create and initilise RV, the array of reproductive values from
  # reproduction
  ii <- seq(-P$Q, P$Q, 1)
  nn <- seq(0, P$N, 1)
  RV <- array(0, dim=c(length(ii), length(nn)),
      dimnames=list(i=ii, n=nn))
  RV <- outer(ii, nn, R, P=P)
  RV <- t(apply(RV, 1, function(z) z/sum(z)))
	RV <- P$clutch.size * RV/max(RV)
  S <- array(0, dim=c(length(ii), length(nn)),
      dimnames=list(i=ii, n=nn))
  S <- outer(ii, nn, surv.breeding, P=P)
	s <- t(apply(S, 1, function(z) {
						 zz <- cumprod(rev(z)) * P$R.future
						 rev(zz)
			}))
	RV <- RV + s
  invisible(RV)
}#}}}

init.V <- function(P) {#{{{
  # create and initilise V, the array of reproductive values
  ii <- seq(-P$Q, P$Q, 1)
  kk <- seq(0, P$K, 1)
  nn <- seq(0, P$N, 1)
  V <- array(0, dim=c(length(ii), length(kk), length(nn)),
      dimnames=list(i=ii, k=kk, n=nn))
  #V[,P$K+1,] <- outer(ii, nn, R, P=P)
  V[,P$K+1,] <- init.RV(P)
  invisible(V)
}#}}}

init.Action <- function(P, init.val=NA) {#{{{
  # create and initilise Action the array of optimal behaviour
  ii <- seq(-P$Q, P$Q, 1)
  kk <- seq(0, P$K, 1)
  nn <- seq(0, P$N, 1)
  Action <- array(init.val, dim=c(length(ii), length(kk), length(nn)),
      dimnames=list(i=ii, k=kk, n=nn))
  invisible(Action)
}#}}}


backward <- function(P=P) {#{{{
  # the backward computations
  Res <- list( #{{{
        A = NA, # store the value of init.A
        V = NA, # store the value of init.V
        #RV = NA, # store the value of init.RV
        Action = NA, # store the value of init.Action
        Diff = NA, # store the value of init.Action
        H.L = NA, # store the value of init.Action
        H.S = NA, # store the value of init.Action
        S = NA # store the value of init.Action
        ) #}}}
  # initialising arrays
  Res$A <- init.A(P)
  Res$V <- init.V(P)
  #Res$RV <- init.RV(P)
  Res$Action <- init.Action(P, 1)
  Res$Diff <- init.Action(P)
  Res$H.L <- init.Action(P)
  Res$H.S <- init.Action(P)
  Res$S <- init.Action(P)
  nn <- seq(0, P$N, 1)
  kk <- seq(0, P$K, 1)
  ii <- seq(-P$Q, P$Q, 1)
  
  for(n in rev(1:length(nn[-length(nn)]))) {
    for(k in rev(1:length(kk[-length(kk)]))) {
      for(i in 1:length(ii)) {
        # print(paste(i,k,n, sep=":"))
        # rescaling removed as the M function expects grid time
        # Question: how should we deal with the rescaling of the rate of
        # mortality? At the moment I have removed that rescaling too.
        # John's version
        #H.stay <- (1 - M((nn[n]+ii[i])*P$De, k, P)*P$De)*Res$V[i, k, n+1]
				H.stay <- (1 - M((nn[n]+ii[i]), k, P)) * Res$V[i, k, n+1]
				H.leave <- sum(Res$A[i,,k+1] * Res$V[,k+1,n+1])
        Res$H.S[i,k,n] <- H.stay
        Res$H.L[i,k,n] <- H.leave
        Res$V[i,k,n] <- max(H.stay, H.leave)
        # Adding 1e-6 to avoid ambiguity caused by imprecise machine
        # computation. This is not a very clean solution, we would need
        # errors-in-decision-making approach here.
        Res$Action[i,k,n] <- ifelse(H.stay+1e-6 > H.leave, 1, 0)
        Res$Diff[i,k,n] <- H.stay - H.leave
      }
    }
  }
  invisible(Res)
}#}}}

forward <- function(Res, P) {#{{{
  # the forward computations
  nn <- seq(0, P$N, 1)
  kk <- seq(0, P$K, 1)
  ii <- seq(-P$Q, P$Q, 1)
  for(n in 1:length(nn[-length(nn)])) {
    for(k in 1:length(kk)) {
    #for(k in 1:length(kk[-length(kk)])) {
      for(i in 1:length(ii)) {
        if(Res$S[i, k, n] > 0) {
          if(Res$Action[i, k, n] > 0.5) { # stay
            mort <- M(nn[n]+ii[i], k, P)
            Res$S[i, k, n+1] <- Res$S[i, k, n+1] + (1-mort) * Res$S[i, k, n]
          } else { # leave
            Res$S[, k+1, n+1] <- Res$S[, k+1, n+1] + Res$S[i, k, n] *
            Res$A[i, , k+1]
          }
        }
      }
    }
  }
  invisible(Res)
}#}}}
#}}}

## calculation on the result of the forward computation {{{
stay <- function(x, P=P, pop.size=1000) {#{{{
  l.i <- min(which(x > 1/pop.size))
  u.i <- max(which(x > 1/pop.size))
  Time <- P$De * (0:P$N)[l.i:u.i]
  pop.size <- 1-x[l.i:u.i]
  col.pal <- cut(pop.size, breaks=c(-0.01, seq(0.1,1,0.1)))
  list(Time=Time, pop.size=pop.size, col.pal=as.numeric(col.pal))
}#}}}

pl.stay <- function(S, P, xlim=c(0, P$N*P$De), ...) {#{{{
  # plot the staying time for each site in a composite plot
  hc <- heat.colors(10)
  z <- apply(S, c(2,3), sum)
  plot(1,1, type="n", ylim=c(0,nrow(z)-1), xlim=xlim, 
       xlab="time (real)", ylab="sites", ...)
  for(i in 1:nrow(z)) {
    s <- stay(z[i,], P)
    points(s$Time, rep(i-1, length(s$Time)), pch=16, col=hc[s$col.pal])
  }
}#}}}

m.arr <- function(x, P=P, pop.size=1000) {#{{{
  l.i <- min(which(x > 1/pop.size))
  u.i <- which.max(x)
  Time <- (0:P$N)[l.i:u.i]
  pop.size <- x[l.i:u.i]
  P$De*sum(Time*pop.size)/sum(pop.size)
}#}}}

m.arrival.time <- function(S) {#{{{
  # this function calculates the mean arrival time for each site
  # separatelly
  # S is returned by 'forward'
  # arrival time is the mean time during the period when the population
  # size in a given size increasing
  z <- apply(S, c(2,3), sum)
  apply(z, 1, m.arr, P=P, pop.size=1000)
}#}}}

pl.m.arr <- function(a, P=P) {#{{{
  # plot the mean arrival times for the sites
  # 'a' is returned by 'm.arrival.time'
  plot(1:P$K, a[-1], type="S")
  points(1:P$K, a[-1], pch=16, col=1)
}#}}}
#}}}


## graphics {{{
pl.mortality <- function(P) {#{{{
  # plot the mortality function for each site as a surface
  # 'i' is the advancement of spring
  n <- 0:P$N
  k <- 0:P$K
  op <- par(mar=c(2,4,0,1)+0.1, oma=c(2,0,2,0))
  layout(matrix(1:3, ncol=1))
  for(i in c(-P$Q, 0, P$Q)) {
    m <- sapply(k, function(kk) M(n+i, kk+1, P=P))
    #nn <- seq(1,1001, 10)
    #persp(m[nn,], theta=45, phi=45, xlab="n", ylab="k",
    #zlab="mortality", main=paste("i=", i))
    matplot(n*P$De, m, type="l", xlab=,
        ylab="mortality", lty=k+1, col=k+1)
    legend("topright", legend=k, lty=k+1, col=k+1, title="sites",
        bty="n")
    mtext(paste("i=", i), side=3, line=0, adj=0.5)
  }
  #axis(1, outer=TRUE)
  mtext("time (real time)", side=1, outer=TRUE, line=1, adj=0.55)
  layout(1)
  par(op)
}#}}}
pl.V <- function(P) {#{{{
  # plot the mortality function for each site as a surface
  # 'i' is the advancement of spring
  n <- 0:P$N
  k <- 0:P$K
  ii <- seq(-P$Q, P$Q, 1)
  layout(matrix(1:3, ncol=3))
  for(i in c(-P$Q, 0, P$Q)) {
    nn <- seq(1,1001, 10)
    j <- which(ii==i)
    persp(t(V[j,,nn]), theta=45, phi=50, xlab="n", ylab="k",
        zlab="reproductive value", main=paste("i=", i))
  }
  layout(1)
}#}}}
pl.Arr <- function(m, P) {#{{{
  # plot the values of an array created by 'init.Action'
  n <- 0:P$N
  k <- 0:P$K
  ii <- seq(-P$Q, P$Q, 1)
  layout(matrix(1:3, ncol=3))
  for(i in c(-P$Q, 0, P$Q)) {
    nn <- seq(1,1001, 10)
    j <- which(ii==i)
    persp(t(m[j,,nn]), theta=45, phi=50, xlab="n", ylab="k",
        zlab="value", main=paste("i=", i))
  }
  layout(1)
}#}}}
pl.Action <- function(P) {#{{{
  # plot the mortality function for each site as a surface
  # 'i' is the advancement of spring
  n <- 0:P$N
  k <- 0:P$K
  ii <- seq(-P$Q, P$Q, 1)
  layout(matrix(1:3, ncol=3))
  for(i in c(-P$Q, 0, P$Q)) {
    nn <- seq(1,1001, 10)
    j <- which(ii==i)
    persp(t(Action[j,,nn]), theta=45, phi=50, xlab="n", ylab="k",
        zlab="staying", main=paste("i=", i))
  }
  layout(1)
}#}}}
pl.Diff <- function(P) {#{{{
  # plot the mortality function for each site as a surface
  # 'i' is the advancement of spring
  n <- 0:P$N
  k <- 0:P$K
  ii <- seq(-P$Q, P$Q, 1)
  layout(matrix(1:3, ncol=3))
  for(i in c(-P$Q, 0, P$Q)) {
    nn <- seq(1,1001, 10)
    j <- which(ii==i)
    persp(t(Diff[j,,nn]), theta=45, phi=50, xlab="n", ylab="k",
        zlab="H.stay - H.leave", main=paste("i=", i))
  }
  layout(1)
}#}}}
pl.res <- function() {#{{{
layout(matrix(c(1,1,1,2,2,2,3,3,4,4,5,5), ncol=2))
matplot(P$De*seq(-P$Q, P$Q, 1), L[,-(P$K+1)], type="l", lty=0:P$K+1,
    col=0:P$K+1, xlab="advancement of spring (real time)", 
    ylab="leaving time (real time)", xlim=P$De*c(-P$Q, P$Q*1.05))
legend("topright", legend=0:(P$K-1), lty=1:(P$K), col=1:(P$K), bty="n",
    title="location", cex=0.55)

# Plot 2 - Terminal Reward on breeding site & onset of spring on al sites
plot(V[P$Q+1,P$K+1,] ~ c(1:(P$N+1)), 
    main="Terminal Reward on breeding site \n & onset of spring", 
    xlab="time (in grid time)", ylab="Reproductive value", type="l",
    lty="solid", col="grey", xlim=c(min(0,P$tau[1]), P$N),
    ylim=c(0,0.1))
abline(v=P$tau, lty=3, col="lightblue")
nn <- seq(-500, 1000, length=10000)
for(jj in 0:K+1) {
  hx <- dnorm(nn, mean=P$tau[jj], sd=P$sigma[jj])
  lines(nn, hx, col=jj)
}
# Plot 3 - Mortality functions
n <- 0:P$N
k <- 0:P$K
for(i in c(-P$Q, 0, P$Q)) {
  m <- sapply(k, function(kk) M(n+i, kk+1, P=P))
  matplot(n*P$De, m, type="l", xlab=, ylab="mortality", lty=k+1, col=k+1)
  legend("topright", legend=k, lty=k+1, col=k+1, title="sites", bty="n")
  mtext(paste("i=", i), side=3, line=0, adj=0.5)
  abline(v=L[(P$Q+1)+i,-ncol(L)], col=1:(P$K))
}
layout(1)
}#}}}
#}}}
#}}}

# vim: foldmethod=marker
