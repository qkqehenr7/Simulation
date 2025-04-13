# Robert and Casella (2009), Introducing Monte Carlo Methods With R
a = 2.7; b = 6.3;

# Sample size
N = 5000;

# trial
n.sim = 100

target.f = function(x) x^(a-1)*(1-x)^(b-1)
pval.1 = pval.2 = pval.3 = numeric(n.sim)

# Use Metropolis-Hastings Algorithm For Generating Beta Dist
for (i in 1:n.sim) {
  x = numeric(N)
  x[1] = runif(1) 
  for (n in 1:(N-1)) {
    y = runif(1)
    if (runif(1) < target.f(y) / target.f(x[n])) x[n+1] = y else x[n+1] = x[n]
  }
  
  # K-S Test
  pval.1[i] = ks.test(jitter(x[1:100]), "pbeta", a, b)$p.value
  pval.2[i] = ks.test(jitter(x[4901:5000]), "pbeta", a, b)$p.value
  pval.3[i] = ks.test(jitter(x[seq(4010, 5000, 10)]), "pbeta", a, b)$p.value
}

## Result
c(mean(pval.1), mean(pval.2), mean(pval.3))