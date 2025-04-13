# Robert and Casella (2009), Introducing Monte Carlo Methods with R
X = c(91,504,557,609,693,727,803,857,929,970,1043,1089,1195,1384,1713)

# Parameters
xbar = mean(X); n = length(X)
theta0 = xbar; tau.sq = var(X); a = (n-1)/2; b = (n-1)*var(X)/2

N = 5000
theta = sig.sq = numeric(N)
theta[1] = rnorm(1, mean = theta0, sd = sqrt(tau.sq))
sig.sq[1] = 1/rgamma(1, a, b)

cond.mean = function(s2) s2*theta0 / (s2 + n*tau.sq) + n*tau.sq*xbar/(s2 + n*tau.sq)
cond.sd = function(s2) sqrt(s2*tau.sq / (s2 + n*tau.sq))
cond.a = n/2 + a
cond.b = function(t) 0.5 * sum((X-t)^2) + b

for (t in 1:(N-1)){
  theta[t+1] = rnorm(1, cond.mean(sig.sq[t]), cond.sd(sig.sq[t]))
  sig.sq[t+1] = 1/rgamma(1, cond.a, cond.b(theta[t+1]))
}

log(quantile(theta, c(0.05, 0.95)))
log(sqrt(quantile(sig.sq, c(0.05, 0.95))))

# Use last 1000 value
log(quantile(theta[4001:5000], c(0.05, 0.95)))
log(sqrt(quantile(sig.sq[4001:5000], c(0.05, 0.95))))

# Non-bayesian interval based on t-dist
log(t.test(X, conf.level=0.9)$conf.int)

# Non-Bayesian interval based on chi-square distribution
c(log(sqrt((n-1)*var(X)/qchisq(0.95, df=n-1))),
  log(sqrt((n-1)*var(X)/qchisq(0.05, df=n-1))))
