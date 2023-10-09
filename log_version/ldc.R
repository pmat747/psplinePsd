# Implementation of Bayesian using Whittle likelihood and Pspline Prior  applied to SGWB A channel from noiseless MLDC
#library(beyondWhittle)
library(rjags)
library(rhdf5)
library(bspec)
library(beyondWhittle)
#library(psplinePsd)
library(ggplot2)
#SGWB 

setwd("C:/Users/em14066/Desktop/PsplinesExperiments/datasets")
Ha <- H5Fopen(name="LDC1-6_SGWB_v1_noiserand.hdf5")
ampl<-Ha$H5LISA$GWSources$`SGWB-0`$EnergyAmplitude
slope <- Ha$H5LISA$GWSources$`SGWB-0`$EnergySlope
FreqRef <- Ha$H5LISA$GWSources$`SGWB-0`$FrequencyRef
FreqShape <- Ha$H5LISA$GWSources$`SGWB-0`$FrequencyShape
FreqSky <- Ha$H5LISA$GWSources$`SGWB-0`$Sky
XYZ<- t(Ha$H5LISA$PreProcess$TDIdata)

file = H5Fopen(name="LDC1-6_SGWB_v1_noiserand.hdf5")
dataset = list.datasets(Ha)
for(i in 1:length(dataset)){
  print(i)
  print(dataset[i])
}

h5closeAll()

setwd("C:/Users/em14066/Desktop/PsplinesExperiments/logpsplinePsd_AR")

Ha$H5LISA$Observatory$DataSources$pm3s$PowerSpectralDensity

#########################################################################

set.seed(1)
start<-100000
n <- 2^13
N <- n/2
X <- XYZ[(start+1):(start+n),2]
Y <- XYZ[(start+1):(start+n),3]
Z <- XYZ[(start+1):(start+n),4]
deltat<- XYZ[2,1]-XYZ[1,1] # sampled every 5 seconds
fs <- 1/deltat            #sampling frequency,

A<-(Z-X)/sqrt(2)
E<-(X-2*Y+Z)/sqrt(6)
T<-(X+Y+Z)/sqrt(3)
#par(mfrow=c(3,1))
#plot.ts(A)
#plot.ts(E)
#plot.ts(T)



hw = hannwindow(n)  # Requires bspec package, first Hann window, then standardize
A <- hw*A
E <- hw*T
T <- hw*T
#plot.ts(A)
#plot.ts(E)
#plot.ts(T)


#Plotting of Periodogram and true PSD
#####################################

pdgrm<-rep(0,N)
pdgrm <- 2*deltat/N*((abs(fft(A)))^2 )[2:(N-1)]
omega<-seq(0,fs/2,length=N)[-c(1,N)] # in Hz

# ggplots
d1<-data.frame(omega,pdgrm)

#response function for Channel A, E
H0 <- 2.175*10^(-18)# Hubble constant
L <- 2.5*10^6 # arm length in km
c <- 299792 # speed of light in km/s
fstar <- c/(2*pi*L)
fq <- omega/fstar
RAA <- 4*(3/10 + 169/1680*fq^2 + 85/6048*fq^4 - 17873/15667200*fq^6 + 19121/2476656000*fq^8) # Guillaume (50)
RA <- RAA*16/9*2/pi*fq^4

Omega_GW <- ampl*(omega/FreqRef)^slope  # true power law for Omega_GW(f)
d2<-data.frame(omega,Omega_GW)

factor <- 2*pi^2/(3*H0^2)*omega^3/RA
PL <- Omega_GW/factor   # true power law transformed to PSD
Omega_pdgrm <- factor*pdgrm  #periodogram transformed to GW energy density Omega_GW


#plot pdgrm in Hz on log scale overlaid by true power law
#p1 <- ggplot(d1,aes(x=omega, y=pdgrm)) + scale_x_continuous(trans = "log10") + scale_y_continuous(trans = "log10") + geom_line(color="gray") + geom_line(d2,mapping=aes(x=omega, y=PL, color="red")) + ggtitle("Periodogram of A Channel")
#p1

#plot energy pdgrm in Hz on log scale overlaid by true power law of Omega_GW
#p2 <- ggplot(d1,aes(x=omega, y=Omega_pdgrm)) + scale_x_continuous(trans = "log10") + scale_y_continuous(trans = "log10") + geom_line(color="gray") + geom_line(d2,mapping=aes(x=omega, y=Omega_GW, color="red"))+ ggtitle("Energy PSD of A Channel")
#p2


#nonparametric spectral density, data standardized

# standardize the time series to avoid overflow
AS <- (A-mean(A))/sd(A)
pdgrmAS<-rep(0,N)
pdgrmAS <- 2*deltat/N*((abs(fft(AS)))^2 )[2:(N-1)]
pdgrmA <- sd(A)^2*pdgrmAS

pilotmcmc = gibbs_pspline(data=AS, 100000, 50000, 10, amh= TRUE)
mcmcP1 = mcmcP2 =  pilotmcmc;

post = mcmcP1$fpsd.sample
# Transformed versions of these for uniform CI construction
log.fpsd.sample <- post
log.fpsd.s      <- apply(log.fpsd.sample, 1, stats::median);
log.fpsd.mad    <- apply(log.fpsd.sample, 1, stats::mad);
log.fpsd.help   <- apply(log.fpsd.sample, 1, uniformmax);
log.Cvalue      <- stats::quantile(log.fpsd.help, 0.9);

# Compute Uniform CIs
psd.u95 <- log.fpsd.s + log.Cvalue * log.fpsd.mad;
psd.u05 <- log.fpsd.s - log.Cvalue * log.fpsd.mad;

psd.u95 <- psd.u95[2:(N-1)];
psd.u05 <- psd.u05[2:(N-1)];


#########################################
mcmcpsdmedian <- log.fpsd.s[2:(N-1)];
n = length(AS);
omega  <- 2 * (1:(n / 2 + 1) - 1) / n;
length(omega)
length(pilotmcmc$pdgrm)
length(log.fpsd.s)
length(log.fpsd.mad)
length(pdgrmA)
length(psd.u05)
###########################################

dP2 <- data.frame(omega,pilotmcmc$pdgrm, mcmcpsdmedian, psd.u05,psd.u95);

#pdf("sgwb.pdf", width = 14)
plotP2 <- ggplot(dP2,aes(x=omega, y=pdgrmA)) + scale_x_continuous(trans = "log10") + scale_y_continuous(trans = "log10") + geom_line(color="gray") + 
  geom_line(dP2,mapping=aes(x=omega, y=sd(A)^2*2*pi*2*deltat*exp(mcmcpsdmedian), color="Whittle", linetype="median"))  +
  geom_line(dP2,mapping=aes(x=omega, y=sd(A)^2*2*pi*2*deltat*exp(psd.u05), color="Whittle", linetype="uniform")) + 
  geom_line(dP2,mapping=aes(x=omega, y=sd(A)^2*2*pi*2*deltat*exp(psd.u95), color="Whittle", linetype="uniform")) +
  ggtitle("Periodogram of A Channel, Whittle Likelihood")
plotP2
#dev.off()

stats::spectrum(AS, n.freq = 100, method = "ar", 
                plot = FALSE)$method;

############
### Plot ###
############

burnin = 1000;
thin   = 2;
N      = length(mcmcP1$ll.trace); N;
index  = seq(burnin+1, N, by = thin);
length(index)
ts.plot(mcmcP1$ll.trace[index])

ts.plot(mcmcP1$V[20,index])
ts.plot(mcmcP1$tau[index])
ts.plot(mcmcP1$phi[index])
ts.plot(mcmcP1$delta[index])
ts.plot(mcmcP1$V[25, index])
plot(mcmcP1$V[17,index],mcmcP1$V[21,index])

post = mcmcP1$fpsd.sample[, index];
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

n = length(AS);
omega  <- 2 * (1:(n / 2 + 1) - 1) / n;
lambda <- omega * pi;

aux = c(psd.u05[-c(1, length(lambda))], psd.u95[-c(1, length(lambda))],
        log(mcmcP1$pdgrm)[-1]); 

#pdf("ldc.pdf")
plot(seq(0,pi,length = length(mcmcP1$log.psd.median))[-1], 
     log(mcmcP1$pdgrm)[-1], type = "l",
     ylim = c(min(aux), max(aux)),
     xlab = "Frequency", ylab = "log psd",
     main = "LDC", col = "gray");

lines(lambda[-c(1, length(lambda))], m[-c(1, length(m))]);
lines(lambda[-c(1, length(lambda))], psd.u05[-c(1, length(lambda))], lty = 2);
lines(lambda[-c(1, length(lambda))], psd.u95[-c(1, length(lambda))], lty = 2); 
#dev.off()
