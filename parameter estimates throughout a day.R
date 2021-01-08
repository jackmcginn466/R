#Peak frequency estimate throughout the day
#Parameter estimates of real world data

library(R.matlab) # for loading mat files
library(resample) # for std of columns
#######################################################

buoyName = "Scilly1"
year = 2017
month = 6
day = 27
folderName = sprintf("%s}%04d",buoyName,year)
fileName = sprintf("%s/%s}%04d-%02d-%02d}raw.mat",folderName,buoyName,year,month,day)
waveData = readMat(fileName)$timeSeries

##Function gives spectral density of JONSWAP for a vector input of angular frequencies and vector of parameters. As it spectral density it has been halved for positive and negative frequencies
jonswap <- function(omega,theta){
  sigmaV <- 0.07 + 0.02 * ifelse(omega > theta[[2]],1,0)
  g <- ifelse(omega>0,(theta[[1]]*(omega^(-theta[[4]]))*exp(-(theta[[4]]/4)*((omega/theta[[2]])^(-4)))*(theta[[3]]^(exp(-(1/2*sigmaV^2)*((omega/theta[[2]])-1)^2))))/2,0)
  return(g)
}

##Function that converts a continuous spectral density function to discrete approximation for certain values of omega, taking in to account aliasing
DiscreteSpectralDensityCalculation <- function(w,samplingInterval,summationMax,theta){
  discreteSpecDense <- c()
  for(i in 1:length(w))
  {
    angularfrequencies <- seq(w[[i]] - (2*pi*summationMax/samplingInterval),w[[i]] + (2*pi*summationMax/samplingInterval),2*pi/samplingInterval)
    spectraldensityfunction <- jonswap(abs(angularfrequencies),theta)
    discreteSpecDense[[i]] <- sum(unlist(spectraldensityfunction))
  }
  return(discreteSpecDense)
}

#######################################################
### data exploration
dim(waveData)
Delta = 1/1.28 # this is 1.28 Hz data
N = 2304
Imatrix <- matrix(rep(0,48*N/2 +1),nrow = N/2 +1, ncol = 48)


## periodogram
omega = seq(0,pi/Delta,by = 2*pi/Delta/N)
for(i in 1:48)
  {
  I = abs(fft(waveData[,i]))^2*(Delta/2/pi/N)
  Imatrix[,i] = I[1:(N/2+1)]
}

par(mfrow = c(2,1))
plot(omega,I,type = 'l')
abline(v=0.35,col = 'red')
abline(v=3.5,col = 'red')

debiasedwhittleLiklihood <- function(thetaList,periodogramInput = I, AngularFreq = omega, TimeInt = Delta, TotalNumberSample = N, lowerFreq = 0.35, upperFreq = 3.5) {
  omegashrunk <- c()
  NotFilled <- TRUE
  startpos <-0
  endpos <- Inf
  for(i in 1:length(AngularFreq)){
    if(AngularFreq[[i]] > lowerFreq && AngularFreq[[i]]<upperFreq){
      omegashrunk <- c(omegashrunk,AngularFreq[[i]])
      NotFilled <- FALSE
    }else{
      if(i >startpos && NotFilled){
        startpos <- i
      }else if(i < endpos && !NotFilled){
        endpos <- i
      }
    }
  }
  AngularFreq <- omegashrunk
  DiscreteSpectralDensityWhittle <- DiscreteSpectralDensityCalculation(AngularFreq,TimeInt, 5, thetaList)
  autoCov <- Re(fft(unlist(DiscreteSpectralDensityWhittle),inverse = TRUE))*(2*pi/(TimeInt*length(AngularFreq)))
  samplePoints <- length(autoCov)
  discretespecdensefbar <- c() 
  numbersequence <- seq(1,1-((samplePoints-1)/samplePoints),-1/samplePoints)
  discretespecdensefbar <- numbersequence*autoCov
  fbar <- (1/(2*pi))*(Re(fft(unlist(discretespecdensefbar)))*(2*TimeInt) - (TimeInt*autoCov[[1]]))
  whittleL <- -(log(fbar)+(periodogramInput[(startpos + 1):(endpos -1)]/fbar))
  return(-sum(unlist(whittleL))) ##CHANGE BACK TO NEGATIVE FOR OPTIM
}

deBiasedbestomegap <- c()
Imatrix[is.na(Imatrix)] <- 0 #fix how this is done

for(j in 24:48){
deBiasedbestomegap <- c(deBiasedbestomegap,((optim(c(0.5,omega[which.max(Imatrix[,j])],3,5), debiasedwhittleLiklihood,periodogramInput = Imatrix[,j] , lower=c(0.001,0.001,1,1), upper=c(0.99,0.99,10,10), method="L-BFGS-B"))[[1]])[[2]])
print(j)
}

times <- seq(0.5,24,0.5)
deBiasedbestomegap[24] <- NA

plot(x = NULL, y =NULL, main = "Estimate for omegaP throughout a day", xlab = "Time of day", ylab = "OmegaP (Hz)", xlim = c(0,24), ylim=c(0,1))
lines(times, deBiasedbestomegap, xlim=range(times), ylim=range(deBiasedbestomegap), pch=16, col = "red")

deBiasedbestomegap



