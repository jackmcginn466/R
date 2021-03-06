
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

##Key Variables
TimeInterval <- 1
TotalNumberSamplePoints <- 10000
##Parameters for JONSWAP
theta <- c(0.7,0.7,3.3,4)

##Generate fourier frequencies that can be present for given time interval and total number of sample points
DiscreteFrequenciesFirstHalf <- seq(0,((2*pi)*((TotalNumberSamplePoints/2)-1))/(TimeInterval*TotalNumberSamplePoints),(2*pi)/(TimeInterval*TotalNumberSamplePoints))
DiscreteFrequenciesSecondHalf <- seq(pi/TimeInterval, (2*pi)/(TimeInterval*TotalNumberSamplePoints), -(2*pi)/(TimeInterval*TotalNumberSamplePoints))
AngularFrequencies <- c(DiscreteFrequenciesFirstHalf, DiscreteFrequenciesSecondHalf)


DiscreteSpectralDensity <- DiscreteSpectralDensityCalculation(AngularFrequencies,TimeInterval,5,theta)


autocovarianceFunction <-Re(fft(unlist(DiscreteSpectralDensity),inverse = TRUE))*(2*pi/(TimeInterval*TotalNumberSamplePoints))

sigma <- toeplitz(autocovarianceFunction)

matrixL <- chol.default(sigma)

x <- rnorm(TotalNumberSamplePoints,0,1)

y <- t(matrixL) %*% x   



times <- seq(0,(TimeInterval*TotalNumberSamplePoints)-TimeInterval,TimeInterval)

##Plot wave
plot(x = NULL,y = NULL,main="Simulated Guassian wave JONSWAP spectral density", xlab="Time(s)", ylab="Height(m)", xlim=c(0,TimeInterval*TotalNumberSamplePoints +1), ylim = c(-4,4))
points(times, y,pch = "." ,col = "red")
lines(times, y, xlim=range(times), ylim=range(y), pch=16)

##Calculate noisy periodogram and JONSWAP and compare them by plotting them together
angfreqforjonswapplot <-seq(0,pi/TimeInterval,(2*pi)/(100*(TimeInterval*TotalNumberSamplePoints)))

periodogram <- (1/(2*pi*TotalNumberSamplePoints))*Mod((fft(y)))^2 ##ADDED IN HALF FOR NEGATIVE, CHECK THIS, I'VE DONE DIVIDE BY 4 BECAUSE I THOUGHT THE HALF WOULD BE SQUARED WITH THE FOURIER TRANSFORM AND IT FIT THE JONSWAP CURVE BETTER WHEN I DID THIS
periodogram <- periodogram *TimeInterval ##Dont think need TO MULTIPLY IT SO COMMENT OUT ##CHECKED DO NEED

jonswapPoints <-jonswap(angfreqforjonswapplot,theta)

plot(x = NULL, y =NULL, main = "Spectral density of simulated wave with JONSWAP overlaid", xlab = "Angular Frequency (radians/s)", ylab = "Spectral Density", xlim = c(0,pi/TimeInterval), ylim = c(0,4))
lines(angfreqforjonswapplot, jonswapPoints, xlim=range(angfreqforjonswapplot), ylim=range(jonswapPoints), pch=16, col = "red")
lines(AngularFrequencies, periodogram, xlim=range(AngularFrequencies), ylim=range(periodogram), pch=16, col ="blue")

##Whittle Liklihood
##de biased whittle liklihood - get autocovariance function for each theta from spectral density equation

##Negative so can get max in optim
debiasedwhittleLiklihood <- function(thetaList,periodogramInput = periodogram, AngularFreq = AngularFrequencies, TimeInt = TimeInterval, TotalNumberSample = TotalNumberSamplePoints ){
  whittle <- c()
  DiscreteSpectralDensityWhittle <- DiscreteSpectralDensityCalculation(AngularFreq,TimeInt, 5, thetaList)
  autoCov <- Re(fft(unlist(DiscreteSpectralDensityWhittle),inverse = TRUE))*(2*pi/(TimeInt*TotalNumberSample))
  samplePoints <- length(autoCov)
  discretespecdensefbar <- c() 
  numbersequence <- seq(1,1-((samplePoints-1)/samplePoints),-1/samplePoints)
  discretespecdensefbar <- numbersequence*autoCov
  fbar <- (1/(2*pi))*(Re(fft(unlist(discretespecdensefbar)))*(2*TimeInt) - (TimeInt*autoCov[[1]]))
  whittleL <- -(log(fbar)+(periodogramInput/fbar))
  return(-sum(unlist(whittleL))) ##CHANGE BACK TO NEGATIVE FOR OPTIM
}

#Function to calculate Whittle liklihood for given parameters, ##returns minus whittle so it finds max value in optim not min
#THIS IS THE ALIASED WHITTLE LIKLIHOOD
whittleLiklihood <- function(thetaList,periodogramInput = periodogram, AngularFreq = AngularFrequencies, TimeInt = TimeInterval )
{
  
  discreteSpecDenseW <- DiscreteSpectralDensityCalculation(AngularFreq,TimeInt, 5, thetaList)
  whittleL <- -(log(unlist(discreteSpecDenseW))+(periodogramInput/unlist(discreteSpecDenseW)))
  return(-sum(whittleL))
}

##Look in to using optim function to find the best values for whittle parameters - only use periodogram not bartletts
bestTheta <- optim(c(0.5,0.5,5,5), whittleLiklihood, lower=c(0.1,0.1,1.01,1.01), upper=c(0.99,0.99,5,5), method="L-BFGS-B")
deBiasedbestTheta <- optim(c(0.5,0.5,5,5), debiasedwhittleLiklihood , lower=c(0.1,0.1,1.01,1.01), upper=c(0.99,0.99,10,10), method="L-BFGS-B")

jonswapPoints3 <- jonswap(angfreqforjonswapplot,bestTheta[[1]])
jonswapPoints4 <- jonswap(angfreqforjonswapplot,deBiasedbestTheta[[1]])

plot(x = NULL, y =NULL, main = "Whittle predicted parameters of simulated wave and spectral density it was simulated from", xlab = "Angular Frequency (radians/s)", ylab = "Spectral Density", xlim = c(0,pi/TimeInterval), ylim = c(0,4))
lines(angfreqforjonswapplot, jonswapPoints, xlim=range(angfreqforjonswapplot), ylim=range(jonswapPoints), pch=16, col = "red")
lines(angfreqforjonswapplot, jonswapPoints3, xlim=range(angfreqforjonswapplot), ylim=range(jonswapPoints3), pch=16, col = "blue")
lines(angfreqforjonswapplot, jonswapPoints4, xlim=range(angfreqforjonswapplot), ylim=range(jonswapPoints4), pch=16, col = "green")
legend(2.2, 4, legend=c( "JONSWAP", "Whittle Predicted JONSWAP", "De-biased Whittle"),col=c("red", "blue", "green"), lty=1:2, cex=0.8)


amountSplit <- 100

##Splitting the wave in to differnt sections, calculating the periodogram for each and averaging
ySplit <-matrix(y,nrow = amountSplit, ncol = TotalNumberSamplePoints/amountSplit, byrow = TRUE)
periodogramSplit <-matrix(rep(0,TotalNumberSamplePoints),nrow = amountSplit, ncol = TotalNumberSamplePoints/amountSplit)
for(i in 1:nrow(ySplit)){
  periodogramSplit[i,] <- (1/(2*pi*ncol(ySplit)))*Mod((fft(ySplit[i,])))^2 ##ADDED IN HALF FOR NEGATIVE, CHECK THIS, I'VE DONE DIVIDE BY 4 BECAUSE I THOUGHT THE HALF WOULD BE SQUARED WITH THE FOURIER TRANSFORM AND IT FIT THE JONSWAP CURVE BETTER WHEN I DID THIS
}
averagedPeriodogram <- colMeans(periodogramSplit)
averagedPeriodogram <- averagedPeriodogram*TimeInterval

##Calculating new frequencies the periodogram is for now its been split to represent smaller ranges in time of the wave
DiscreteFrequenciesFirstHalfSplit <- seq(0,((2*pi)*(((TotalNumberSamplePoints/amountSplit)/2)-1))/(TimeInterval*(TotalNumberSamplePoints/amountSplit)),(2*pi)/(TimeInterval*(TotalNumberSamplePoints/amountSplit)))
DiscreteFrequenciesSecondHalfSplit <- seq(pi/TimeInterval, (2*pi)/(TimeInterval*(TotalNumberSamplePoints/amountSplit)), -(2*pi)/(TimeInterval*(TotalNumberSamplePoints/amountSplit)))
AngularFrequenciesSplit <- c(DiscreteFrequenciesFirstHalfSplit, DiscreteFrequenciesSecondHalfSplit)


##Plotting the averaged periodogram and JONSWAP on the same plot
plot(x = NULL, y =NULL, main = "Bartletts periodogram of simulated wave and spectral density it was simulated from", xlab = "Angular Frequency (radians/s)", ylab = "Spectral Density", xlim = c(0,pi/TimeInterval), ylim = c(0,2))
lines(angfreqforjonswapplot, jonswapPoints, xlim=range(angfreqforjonswapplot), ylim=range(jonswapPoints), pch=16, col = "red")
lines(AngularFrequenciesSplit, averagedPeriodogram, xlim=range(AngularFrequenciesSplit), ylim=range(averagedPeriodogram), pch=16, col ="black")
lines(angfreqforjonswapplot, jonswapPoints3, xlim=range(angfreqforjonswapplot), ylim=range(jonswapPoints3), pch=16, col = "blue")
lines(angfreqforjonswapplot, jonswapPoints4, xlim=range(angfreqforjonswapplot), ylim=range(jonswapPoints4), pch=16, col = "green")
legend(2.2, 2, legend=c( "JONSWAP", "Whittle Predicted JONSWAP", "De-biased Whittle","Bartletss Periodogram"),col=c("red", "blue", "green", "black"), lty=1:2, cex=0.8)





bestTheta[[1]]
deBiasedbestTheta[[1]]
