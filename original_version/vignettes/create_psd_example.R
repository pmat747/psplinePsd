library(psplinePsd)

set.seed(12345)

# Simulate AR(4) data
n = 2 ^ 7
ar.ex = c(0.9, -0.9, 0.9, -0.9)
data = arima.sim(n, model = list(ar = ar.ex))
data = data - mean(data)
data = data / stats::sd(data)


degree = 3
k = 5
FZ <- fast_ft(data);
pdgrm <- abs(FZ) ^ 2;
omega <- 2 * (1:(n / 2 + 1) - 1) / n;
v = get_initial_weights(pdgrm, k)
knots = knotLoc(data = data, k = k, degree = degree, eqSpaced = TRUE);
db.list <- dbspline(omega, knots, degree);
psd <- qpsd(omega, k, v, degree, db.list);

plot(psd)
