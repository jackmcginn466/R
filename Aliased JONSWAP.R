##DOING BARTLETTS METHOD FOR SMOOTHING THE PERIODOGRAM AND COMPARING IT WITH JONSWAP
##Simulated wave from a Gaussian process, autocovariance matrix for this wave is applied from the JONSWAP equation with given parameters


##Function gives spectral density of JONSWAP for a vector input of angular frequencies and vector of parameters. As it spectral density it has been halved for positive and negative frequencies
jonswap <- function(w,theta){
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
TimeInterval <- 0.5
TotalNumberSamplePoints <- 100
##Parameters for JONSWAP
theta <- c(0.7,0.7,3.3,4)

##Generate fourier frequencies that can be present for given time interval and total number of sample points
DiscreteFrequenciesFirstHalf <- seq(0,((2*pi)*((TotalNumberSamplePoints/2)-1))/(TimeInterval*TotalNumberSamplePoints),(2*pi)/(TimeInterval*TotalNumberSamplePoints))
DiscreteFrequenciesSecondHalf <- seq(pi/TimeInterval, (2*pi)/(TimeInterval*TotalNumberSamplePoints), -(2*pi)/(TimeInterval*TotalNumberSamplePoints))
AngularFrequencies <- c(DiscreteFrequenciesFirstHalf, DiscreteFrequenciesSecondHalf)
##CHECK THIS IS DONE CORRECTLY, I have done 0 to ((2*pi)*((M/2)-1))/(delta*M) for the first half and then the frequencies go back down again from pi/delta to 2*pi/M*delta (goes to being one interval off 0 again)

##TO DO THE OVERSAMPLING
angfreqforautocov <- seq(pi/TimeInterval,20*pi/TimeInterval, 2*pi/TimeInterval*TotalNumberSamplePoints)
angfreqforautocov <- c(AngularFrequencies, angfreqforautocov)

DiscreteSpectralDensity <- DiscreteSpectralDensityCalculation(angfreqforautocov,TimeInterval,5,theta)
DiscreteSpectralDensity <- DiscreteSpectralDensity[1:TotalNumberSamplePoints]

##TOOK THE REAL PART OF FFT, ALL COMPLEX PARTS WERE 0 JUST WOULDN'T WORK UNLESS I DID
autocovarianceFunction <-Re(fft(unlist(DiscreteSpectralDensity),inverse = TRUE))*(2*pi/(TimeInterval*TotalNumberSamplePoints))
autocovarianceFunction <- autocovarianceFunction[1:TotalNumberSamplePoints]

sigma <- toeplitz(autocovarianceFunction)

matrixL <- chol(sigma)##THIS ONLY WORKS FOR SPECIFIC TIME INTERVAL AND SAMPLE POINTS COMBINATIONS - CHECK WHY

x <- rnorm(TotalNumberSamplePoints,0,1)

y <- t(matrixL) %*% x   ## times transpose of matrixL as t(matrixL) * matrixL gave sigma not the other way around

##CHECK RANGE ON TIM, MINUSING ONE TIME INTERVAL AT THE END
times <- seq(0,(TimeInterval*TotalNumberSamplePoints)-TimeInterval,TimeInterval)

##Plot wave
plot(x = NULL,y = NULL,main="Simulated Guassian wave JONSWAP spectral density", xlab="Time(s)", ylab="Height(m)", xlim=c(0,TimeInterval*TotalNumberSamplePoints +1), ylim = c(-4,4))
points(times, y,pch = "." ,col = "red")
lines(times, y, xlim=range(times), ylim=range(y), pch=16)
length(y)
length(times)
##Generate JONSWAP for plot
angfreqforjonswapplot <-seq(0,pi/TimeInterval,(2*pi)/(100*(TimeInterval*TotalNumberSamplePoints)))
jonswapPointsaliased <-DiscreteSpectralDensityCalculation(angfreqforjonswapplot,TimeInterval,6,theta)
jonswapPoints <- jonswap(angfreqforjonswapplot,theta)

jonswap(AngularFrequencies,c(0.7,0.7,3.3,4.5))

##Amount periodogram is split in to to be averaged over
amountSplit <- 10

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
plot(x = NULL, y =NULL, main = "Bartletts periodogram of simulated wave and spectral density it was simulated from", xlab = "Angular Frequency (radians/s)", ylab = "Spectral Density", xlim = c(0,pi/TimeInterval), ylim = c(0,4))
lines(angfreqforjonswapplot, jonswapPointsaliased, xlim=range(angfreqforjonswapplot), ylim=range(jonswapPointsaliased), pch=16, col = "red")
lines(AngularFrequenciesSplit, averagedPeriodogram, xlim=range(AngularFrequenciesSplit), ylim=range(averagedPeriodogram), pch=16, col ="blue")
lines(angfreqforjonswapplot, jonswapPoints, xlim=range(angfreqforjonswapplot), ylim=range(jonswapPoints), pch=16, col = "green")
legend(2.2, 4, legend=c( "Bartletts Periodogram", "Alliased JONSWAP", "JONSWAP"),col=c("blue", "red","green"), lty=1:2, cex=0.8)

