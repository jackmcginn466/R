#Removing swell
#23rd june a good day and 15th
#june 8th has good 2 peak split
#Parameter estimates of real world data
#install.packages("bayestestR")
#install.packages("ggpmisc")
#install.packages("pracma")
library(pracma)
library(ggpmisc)
library(bayestestR)
library(R.matlab) # for loading mat files
library(resample) # for std of columns
#######################################################
### example of loading a record
## make sure the working directory is set to the correct location (if it is different you may need to add folder names etc)
# note you don't have to use sprintf, but this will help you to do this in a loop if you want

##CURRENT DATES AND TIME RETURN REALLY GOOD FREQUENCY SELECTION GRAPHS

buoyName = "Wavehub"
year = 2017
month = 6
day = 14
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
N = 2300
## periodogram
DiscreteFrequenciesFirstHalf <- seq(0,((2*pi)*((N/2)-1))/(Delta*N),(2*pi)/(Delta*N))
DiscreteFrequenciesSecondHalf <- seq(pi/Delta, (2*pi)/(Delta*N), -(2*pi)/(Delta*N))
omega <- c(DiscreteFrequenciesFirstHalf, DiscreteFrequenciesSecondHalf)
#omega = seq(0,pi/Delta,by = 2*pi/Delta/N)
y = (waveData[,5])[1:N]
omega = seq(0,pi/Delta,by = 2*pi/Delta/N)
I = abs(fft(y))^2*(Delta/2/pi/N)
I = I[1:(N/2+1)]

lowerfreq = 0.35
upperfreq = 3.5
length(y)
##Amount periodogram is split in to to be averaged over
amountSplit <- 25
##Splitting the wave in to differnt sections, calculating the periodogram for each and averaging
ySplit <-matrix(y,nrow = amountSplit, ncol = N/amountSplit, byrow = TRUE)
periodogramSplit <-matrix(rep(0,N),nrow = amountSplit, ncol = (N/2+1)/amountSplit)
for(i in 1:nrow(ySplit)){
  Isplit <- (1/(2*pi*ncol(ySplit)))*Mod((fft(ySplit[i,])))^2 ##ADDED 
  length(Isplit[1:((N/2+1)/amountSplit)])
  periodogramSplit[i,] <- Isplit[1:((N/2+1)/amountSplit)]
}
averagedPeriodogram <- colMeans(periodogramSplit)
averagedPeriodogram <- averagedPeriodogram*Delta

omegaSplit <- seq(0,pi/Delta,by =(2*pi)/(Delta*((N-1)/amountSplit)))

peakpositions <- findpeaks(averagedPeriodogram, minpeakheight = 0.4*max(averagedPeriodogram))
findpeaks(averagedPeriodogram, minpeakheight = 0.5*max(averagedPeriodogram))

if(nrow(peakpositions)>1){
  lowerfreq = omegaSplit[peakpositions[1,4]]
  if(nrow(peakpositions)>2){
    upperfreq = omegaSplit[peakpositions[2,4]]
  }
}

# 0.008365604 0.746540503 3.172130018 3.084209930

plot(omega,I,type = 'l')
abline(v=lowerfreq,col = 'red')
abline(v=upperfreq,col = 'red')

times <- seq(1,2300,1)

plot(times,y,'l')

omegashrunk <- c()
NotFilled <- TRUE
startPos <-0
endPos <- Inf
for(i in 1:length(omega)){
  if(omega[[i]] > lowerfreq && omega[[i]]<upperfreq){
    omegashrunk <- c(omegashrunk,omega[[i]])
    NotFilled <- FALSE
  }else{
     if(i >startPos && NotFilled){
      startPos <- i
    }else if(i < endPos && !NotFilled){
      endPos <- i
    }
  }
}

debiasedwhittleLiklihood <- function(thetaList,periodogramInput = I, AngularFreq = omegashrunk, TimeInt = Delta, TotalNumberSample = N,endpos = endPos,startpos = startPos) {
  DiscreteSpectralDensityWhittle <- DiscreteSpectralDensityCalculation(AngularFreq,TimeInt, 5, thetaList)
  autoCov <- Re(fft(unlist(DiscreteSpectralDensityWhittle),inverse = TRUE))*(2*pi/(TimeInt*length(AngularFreq)))
  samplePoints <- length(autoCov)
  discretespecdensefbar <- c() 
  numbersequence <- seq(1,1-((samplePoints-1)/samplePoints),-1/samplePoints)
  discretespecdensefbar <- numbersequence*autoCov
  fbar <- (1/(2*pi))*(Re(fft(unlist(discretespecdensefbar)))*(2*TimeInt) - (TimeInt*autoCov[[1]]))
 whittleL <- -(log(fbar)+(periodogramInput[(startpos + 1):(endpos -1)] /fbar))
  return(-sum(unlist(whittleL))) ##CHANGE BACK TO NEGATIVE FOR OPTIM
}

SwellwhittleLiklihood <- function(thetaList,periodogramInput = I, AngularFreq = omega, TimeInt = Delta, TotalNumberSample = N) {
  DiscreteSpectralDensityWhittle <- DiscreteSpectralDensityCalculation(AngularFreq,TimeInt, 5, thetaList)
  autoCov <- Re(fft(unlist(DiscreteSpectralDensityWhittle),inverse = TRUE))*(2*pi/(TimeInt*length(AngularFreq)))
  samplePoints <- length(autoCov)
  discretespecdensefbar <- c() 
  numbersequence <- seq(1,1-((samplePoints-1)/samplePoints),-1/samplePoints)
  discretespecdensefbar <- numbersequence*autoCov
  fbar <- (1/(2*pi))*(Re(fft(unlist(discretespecdensefbar)))*(2*TimeInt) - (TimeInt*autoCov[[1]]))
  whittleL <- -(log(fbar)+(periodogramInput /fbar))
  return(-sum(unlist(whittleL))) ##CHANGE BACK TO NEGATIVE FOR OPTIM
}

initial_alpha = area_under_curve(x = omegashrunk, y = I[(startPos+1):(endPos-1)], method = "trapezoid" )/area_under_curve(x = omegashrunk, y =jonswap(omegashrunk,c(1,omega[which.max(I)],3,5)),method="trapezoid")

deBiasedbestTheta <- optim(c(0.7,omega[which.max(I)],3,4), debiasedwhittleLiklihood , lower=c(0.00001,0.001,1,1), upper=c(0.99,0.99,10,10), method="L-BFGS-B")
#swellTheta <- optim(c(0.7,omega[which.max(I)],3,4), SwellwhittleLiklihood , lower=c(0.00001,0.001,1,1), upper=c(0.99,0.99,10,10), method="L-BFGS-B")

estimatedjonswap <- jonswap(omega,deBiasedbestTheta[[1]])
#swelljonswap <- jonswap(omega,swellTheta[[1]])

plot(x = NULL, y = NULL, main = "Periodogram and estimated spectral density", xlab = "Omega", ylab = "Estimated spectral density", xlim = c(0,3.5), ylim=c(-0,max(I)))
lines(omega, I, xlim=range(omega), ylim=range(I), pch=16, col = "black")
#lines(omega, swelljonswap, xlim=range(omega), ylim=range(swelljonswap), pch=16, col = "purple")
lines(omega, estimatedjonswap, xlim=range(omega), ylim=range(estimatedjonswap), pch=16, col = "red")
lines(omegaSplit, averagedPeriodogram, xlim=range(omegaSplit), ylim=range(averagedPeriodogram), pch=16, col = "blue")
abline(v=lowerfreq,col = 'green')
abline(v=upperfreq,col = 'green')
legend(2.2, max(I), legend=c( "Swell predicted JONSWAP", "De-biased Whittle Predicted JONSWAP","Bartletts periodogram"),col=c("purple", "red","blue"), lty=1:2, cex=0.8)

deBiasedbestTheta[[1]]


