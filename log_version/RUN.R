setwd("G:/My Drive/ResearchExp/PsplinesExperiments/logpsplinePsd_AR")
library(Rcpp)
sourceCpp("internalFunctions.cpp")
sourceCpp("internal_beyondWhittle.cpp")
source("internal_gibbs_core.R")
source("internal_gibbs_s3.R")
source("gibbs_pspline.R")			
source("internal_gibbs_util.R")
source("gibbs_pspline_postProposal.R")	
source("internal_psdsPlot_s3.R")
source("gibbs_pspline_simple.R")		
source("otherFunctions.R")
#source("postprocess.R")
#source("spec_pspline.R")

#############
### AR(4) ###
#############

# Generate AR(4) data 
n = 256;
set.seed(1);
data = arima.sim(n, model = list(ar = c(0.9, -0.9, 0.9, -0.9)));
data = data - mean(data);

# Run MCMC (may take some time)
#mcmc  = gibbs_pspline(data, 30000, 5000, 10, eqSpacedKnots = FALSE)
mcmc = gibbs_pspline(data, 100000, 50000, 10, amh=TRUE)

omega  <- 2 * (1:(n / 2 + 1) - 1) / n;
lambda <- omega * pi;
freq = 2 * pi / n * (1:(n / 2 + 1) - 1)[-c(1, n / 2 + 1)]  # Remove first and last frequency
psd.trueAR4 = psd_arma(freq, ar = c(0.9, -0.9, 0.9, -0.9), ma =0,
                       sigma2 = 1)  # True PSD

burnin = 2500;
thin   = 1;
N      = length(mcmc$ll.trace); N;
index  = seq(burnin+1, N, by = thin)
ts.plot(mcmc$ll.trace[index])

ts.plot(mcmc$V[20,index])
ts.plot(mcmc$tau[index])
ts.plot(mcmc$phi[index])
ts.plot(mcmc$delta[index])
ts.plot(mcmc$V[25, index])
plot(mcmc$V[17,index],mcmc$V[21,index])

post = mcmc$fpsd.sample[, index];
m    = apply(post, 1, median);

# Transformed versions of these for uniform CI construction
log.fpsd.sample <- post
log.fpsd.s    <- apply(log.fpsd.sample, 1, stats::median);
log.fpsd.mad  <- apply(log.fpsd.sample, 1, stats::mad);
log.fpsd.help <- apply(log.fpsd.sample, 1, uniformmax);
log.Cvalue    <- stats::quantile(log.fpsd.help, 0.9);

# Compute Uniform CIs
psd.u95 <- log.fpsd.s + log.Cvalue * log.fpsd.mad;
psd.u05 <- log.fpsd.s - log.Cvalue * log.fpsd.mad;

#pdf("ar4.pdf")
plot(lambda[-c(1, length(lambda))], log(psd.trueAR4), type = "l",
     ylim = c(min(psd.u05[-c(1, length(lambda))]), 
              max(psd.u95[-c(1, length(lambda))])),
     xlab = "Frequency", ylab = "log psd",
     main = "AR(4) case - n=256", col = "gray");

lines(lambda[-c(1, length(lambda))], m[-c(1, length(m))]);
lines(lambda[-c(1, length(lambda))], psd.u05[-c(1, length(lambda))], lty = 2);
lines(lambda[-c(1, length(lambda))], psd.u95[-c(1, length(lambda))], lty = 2);
#dev.off()

#########
###   ###
#########

### 
data1 = scan("1second.txt");
R = gibbs_pspline(data1, 15000, 2000, amh= TRUE);
m = apply(R$fpsd.sample, 1, median);
n = length(data1);
omega  <- 2 * (1:(n / 2 + 1) - 1) / n;
lambda <- omega * pi;
plot(seq(0,pi,length = length(R$log.psd.median)), log(R$pdgrm), col = "gray", type = "l");
lines(seq(0,pi,length = length(R$log.psd.median)), R$log.psd.median);

plot(R);
plot(R$ll.trace, type = "l");


m = apply(R$fpsd.sample[,-c(1:4000)], 1, median);
plot(seq(0,pi,length = length(R$psd.median)), log(R$pdgrm), col = "gray", 
     type = "l");
lines(seq(0,pi,length = length(m)), m, type = "l");
abline(v=c(0.475,0.53,0.75,0.8, 1.45,1.6,2.19,2.35));



R1 = gibbs_pspline(data1, 75000, 40000, 10, 
                  nfreqbin=c(0.475,0.53,0.75,0.8, 1.45,1.6,2.19,2.35),
                  wfreqbin = c(10,10,5,10,5,10,5,10,3), amh=TRUE);
plot(R1$ll.trace, type = "l")
m1 = apply(R1$fpsd.sample, 1, median)
plot(seq(0,pi,length = length(R1$psd.median)), log(R1$pdgrm), col = "gray", type = "l")
lines(seq(0,pi,length = length(m1)), m1, type = "l")
abline(v=c(0.475,0.53,0.75,0.8, 1.45,1.6,2.19,2.35))
knots = knotLoc(data1, k = 40, degree = 3, eqSpaced = FALSE, 
                nfreqbin = c(0.475,0.53,0.75,0.8, 1.45,1.6,2.19,2.35), 
                wfreqbin = c(10,10,5,10,5,10,5,10,3));
knots*pi
points(knots*pi, rep(-110, length(knots)), col = "brown")

plot(seq(0,pi,length = length(R1$log.psd.median)), log(R1$pdgrm), col = "gray", type = "l");
lines(seq(0,pi,length = length(R1$log.psd.median)), R1$log.psd.median);


#############
### Plots ###
#############

# Comparision freq partition

pdf("1sec.pdf", width = 14)
plot(seq(0,pi,length = length(R1$psd.median)), log(R1$pdgrm), col = "gray", 
     type = "l", xlab = "Frequency", ylab = "log psd");
lines(seq(0,pi,length = length(R$psd.median)), m, col = "black");
lines(seq(0,pi,length = length(m1)), m1, type = "l", col = "red");
abline(v=c(0.475,0.53,0.75,0.8, 1.45,1.6,2.19,2.35), lty="dotted");
knots1 = knotLoc(data1, k = 40, degree = 3, eqSpaced = FALSE, 
                 nfreqbin = c(0.475,0.53,0.75,0.8, 1.45,1.6,2.19,2.35), 
                 wfreqbin = c(10,10,5,10,5,10,5,10,3));

knots1*pi;
points(knots1*pi, rep(-110, length(knots)), col = "red");

knots = knotLoc(data1, k = 40, degree = 3, eqSpaced = FALSE, 
                nfreqbin = 1);
points(knots*pi, rep(-108, length(knots)), col = "black");

legend("topright", legend = c("No Partition", "Partition"), 
       col = c("black","red"), lty =1)
dev.off()

###################################
### Posterior samples analysis ####
###################################

ts.plot(mcmc$tau, type="l")
plot(mcmc$V[1,],mcmc$V[2,])
ts.plot(mcmc$phi, type="l")
ts.plot(mcmc$delta, type="l")
ts.plot(mcmc$V[2,])
ess_ar4  = apply(mcmc$V,1,function(x)coda::effectiveSize(x));
apply(mcmc$V,1,function(x)summary(x));


ts.plot(R$ll.trace, type="l")
ts.plot(R1$ll.trace, type="l")
ts.plot(R1$tau, type="l")
ts.plot(R$tau, type="l")

coda::effectiveSize(R$tau)

plot(R$V[1,],R$V[2,])

ts.plot(R$phi, type="l")
ts.plot(R1$phi, type="l")

ts.plot(R$delta, type="l")
ts.plot(R1$delta, type="l")

dim(R$fpsd.sample)

ts.plot(R1$fpsd.sample[4,])

k   = length(knots) + 3 - 1
w = R$pdgrm / sum(R$pdgrm);
w = w[round(seq(1, length(w), length = k))];
w[which(w==0)] = 1e-50; # prevents errors when there are zeros
w = w/sum(w);
#w = w[-k]; # CHANGED
v = w; # log(w)

sum(w)

# weights
plot(R$V[25,],R$V[38,])
ts.plot(R$V[25,])
ts.plot(R$V[38,])

coda::effectiveSize(R$V[25,])
coda::effectiveSize(R$V[38,])

plot(R1$V[25,],R1$V[38,])
ts.plot(R1$V[25,])
ts.plot(R1$V[38,])

coda::effectiveSize(R1$V[25,])
coda::effectiveSize(R1$V[38,])

burnin = 7500;
thin   = 10;
N      = dim(R1$V)[2];
index  = seq(from = burnin + 1, to = N, by = thin);ts.plot(R1$V[25,])

ts.plot(R1$V[38,index])
coda::effectiveSize(R1$V[38,index])

###########
### ESS ###
###########

# ESS for log-psd estimates at specific frequencies
ess  = apply(R$fpsd.sample,1,function(x)coda::effectiveSize(x));
ess1 = apply(R1$fpsd.sample,1,function(x)coda::effectiveSize(x));

###################
### Correlation ###
###################

library(corrplot)
corrplot(cor(t(R1$V)))

burnin = 5;
thin   = 3;

###############
### Sunspot ###   
###############

setwd("G:/My Drive/ResearchExp/PsplinesExperiments/datasets")
sunspot = read.table("sunspot.txt")
sunspot = sunspot[,2];

Rsunspot  = gibbs_pspline(sunspot, 150000, 50000, 10, amh=TRUE);

### plots ###

burnin = 1000;
thin   = 2;
N      = length(Rsunspot$ll.trace); N;
index  = seq(burnin+1, N, by = thin)
ts.plot(Rsunspot$ll.trace[index])

ts.plot(Rsunspot$V[20,index])
ts.plot(Rsunspot$tau[index])
ts.plot(Rsunspot$phi[index])
ts.plot(Rsunspot$delta[index])
ts.plot(Rsunspot$V[25, index])
plot(Rsunspot$V[17,index],Rsunspot$V[21,index])

post = Rsunspot$fpsd.sample[, index];
m    = apply(post, 1, median);

# Transformed versions of these for uniform CI construction
log.fpsd.sample <- post
log.fpsd.s    <- apply(log.fpsd.sample, 1, stats::median);
log.fpsd.mad  <- apply(log.fpsd.sample, 1, stats::mad);
log.fpsd.help <- apply(log.fpsd.sample, 1, uniformmax);
log.Cvalue    <- stats::quantile(log.fpsd.help, 0.9);

# Compute Uniform CIs
psd.u95 <- log.fpsd.s + log.Cvalue * log.fpsd.mad;
psd.u05 <- log.fpsd.s - log.Cvalue * log.fpsd.mad;

n = length(sunspot);
omega  <- 2 * (1:(n / 2 + 1) - 1) / n;
lambda <- omega * pi;

aux = c(psd.u05[-c(1, length(lambda))], psd.u95[-c(1, length(lambda))],
        log(Rsunspot$pdgrm)[-1]); 

pdf("sunspot.pdf")
plot(seq(0,pi,length = length(Rsunspot$log.psd.median))[-1], 
     log(Rsunspot$pdgrm)[-1], type = "l",
     ylim = c(min(aux), max(aux)),
     xlab = "Frequency", ylab = "log psd",
     main = "1 second", col = "gray");

lines(lambda[-c(1, length(lambda))], m[-c(1, length(m))]);
lines(lambda[-c(1, length(lambda))], psd.u05[-c(1, length(lambda))], lty = 2);
lines(lambda[-c(1, length(lambda))], psd.u95[-c(1, length(lambda))], lty = 2);
dev.off()
#write.table(exp(m), "sunspot_psd.txt", row.names = FALSE, col.names = FALSE)

### COMPARED TO SPECTRM FUNCTION ###

p1 = spectrum(sunspot, method = "pgram", plot = "FALSE")
plot(p1$freq, log(2*p1$spec), type = "l", col = "gray", xlab = "Frequency", ylab = "log(psd)");

p2 = spectrum(sunspot, method = "ar", plot = "FALSE")
lines(p2$freq, log(2*p2$spec));

lines(seq(0, 0.5, length = length(Rsunspot$log.psd.median)), 
      Rsunspot$log.psd.median + log(2*2*pi), col = "red");

# The sampling frequency (fs) here is 1, otherwise the factor would be log(2*2*pi/fs)

###############
### Carinae ###   
###############

setwd("C:/Users/em14066/Desktop/PsplinesExperiments/datasets");
carinae = read.csv("Carinae.csv", header = FALSE);
index   = which(carinae == 0); # Finding missing data
aux     = mean(carinae[- index, 1]); # calculating mean
newdata = carinae; # replacing missing data with mean value
newdata[index, 1] = aux; 
newdata = unlist(sqrt(newdata));
newdata = newdata - mean(newdata);
carinae = unlist(newdata);

Rcarinae  = gibbs_pspline(carinae, 150000, 50000, 10, amh=TRUE);

burnin = 1000;
thin   = 2;
N      = length(Rcarinae$ll.trace); N;
index  = seq(burnin+1, N, by = thin)
ts.plot(Rcarinae$ll.trace[index])

ts.plot(Rcarinae$V[20,index])
ts.plot(Rcarinae$tau[index])
ts.plot(Rcarinae$phi[index])
ts.plot(Rcarinae$delta[index])
ts.plot(Rcarinae$V[25, index])
plot(Rcarinae$V[17,index],Rcarinae$V[21,index])

post = Rcarinae$fpsd.sample[, index];
m    = apply(post, 1, median);

# Transformed versions of these for uniform CI construction
log.fpsd.sample <- post
log.fpsd.s    <- apply(log.fpsd.sample, 1, stats::median);
log.fpsd.mad  <- apply(log.fpsd.sample, 1, stats::mad);
log.fpsd.help <- apply(log.fpsd.sample, 1, uniformmax);
log.Cvalue    <- stats::quantile(log.fpsd.help, 0.9);

# Compute Uniform CIs
psd.u95 <- log.fpsd.s + log.Cvalue * log.fpsd.mad;
psd.u05 <- log.fpsd.s - log.Cvalue * log.fpsd.mad;

n = length(carinae);
omega  <- 2 * (1:(n / 2 + 1) - 1) / n;
lambda <- omega * pi;

aux = c(psd.u05[-c(1, length(lambda))], psd.u95[-c(1, length(lambda))],
        log(Rcarinae$pdgrm)[-1]); 

pdf("carinae.pdf")
plot(seq(0,pi,length = length(Rcarinae$log.psd.median))[-1], 
     log(Rcarinae$pdgrm)[-1], type = "l",
     ylim = c(min(aux), max(aux)),
     xlab = "Frequency", ylab = "log psd",
     main = "1 second", col = "gray");

lines(lambda[-c(1, length(lambda))], m[-c(1, length(m))]);
lines(lambda[-c(1, length(lambda))], psd.u05[-c(1, length(lambda))], lty = 2);
lines(lambda[-c(1, length(lambda))], psd.u95[-c(1, length(lambda))], lty = 2); 
dev.off()

######################
### SGWB Power Law ###
######################

x = scan("sgwb_T100_fs25.txt")
R = gibbs_pspline(x, 25000, 5000, 10, amh=TRUE);

setwd("C:/Users/em14066/Desktop/PsplinesExperiments/PSDs")
source("functions.R")

burnin = 100;
thin   = 3;
N      = length(R$ll.trace); N;
index  = seq(burnin+1, N, by = thin); length(index);
ts.plot(R$ll.trace[index])

ts.plot(R$V[20,index]);
ts.plot(R$tau[index]);
ts.plot(R$phi[index]);
ts.plot(R$delta[index]);
ts.plot(R$V[25, index]);
plot(R$V[17,index], R$V[21,index]);

log.psd.est = apply(R$fpsd.sample[, index], 1, stats::median);

fs = 25;
freq = seq(0, fs/2, length = length(R$log.psd.median))
plot(freq, log.psd.est + log(2 * pi * 2 / fs), type = "l");

lines(freq, log(omega_GW(freq)), col = "red", lwd = 1.5);

##################
### LISA Noise ###
##################

y = scan("lisanoise_T100_fs25.txt")
R1 = gibbs_pspline(y, 100000, 10000, 10, amh=TRUE);

burnin = 1700;
thin   = 3;
N      = length(R1$ll.trace); N;
index  = seq(burnin+1, N, by = thin); length(index);
ts.plot(R1$ll.trace[index])

ts.plot(R1$V[20,index]);
ts.plot(R1$tau[index]);
ts.plot(R1$phi[index]);
ts.plot(R1$delta[index]);
ts.plot(R1$V[25, index]);
plot(R1$V[17,index], R1$V[21,index]);

log.psd.est = apply(R1$fpsd.sample[, index], 1, stats::median);

fs = 25;
freq = seq(0, fs/2, length = length(R1$log.psd.median))
plot(freq, log.psd.est + log(2 * pi * 2 / fs), type = "l");

lines(freq, log(Sn(freq)), col = "red", lwd = 1.5);

