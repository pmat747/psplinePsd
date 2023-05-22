library(psplinePsd)

set.seed(12345)

# Simulate AR(4) data
n = 2 ^ 7
ar.ex = c(0.9, -0.9, 0.9, -0.9)
data = arima.sim(n, model = list(ar = ar.ex))
data = data - mean(data)
data = data / stats::sd(data)


# Run MCMC and make plots
mcmc = gibbs_pspline(data, burnin=100, Ntotal=400, degree = 3,  eqSpacedKnots=TRUE)
plot(mcmc)
plot(mcmc, ylog = FALSE, main = "Estimate of PSD using the P-spline method")
