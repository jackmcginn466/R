#box plots and histograms for whittle method

##Function gives spectral density of JONSWAP for a vector input of angular frequencies and vector of parameters. As it spectral density it has been halved for positive and negative frequencies
jonswap <- function(omega,theta){
  sigmaV <- 0.07 + 0.02 * ifelse(omega > theta[[2]],1,0)
  g <- ifelse(omega>0,(theta[[1]]*(omega^(-theta[[4]]))*exp(-(theta[[4]]/4)*((omega/theta[[2]])^(-4)))*(theta[[3]]^(exp(-(1/2*sigmaV^2)*((omega/theta[[2]])-1)^2))))/2,0)
  return(g)
}==

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
TotalNumberSamplePoints <- 1000
##Parameters for JONSWAP
theta <- c(0.7,0.7,3.3,4)

##Generate fourier frequencies that can be present for given time interval and total number of sample points
DiscreteFrequenciesFirstHalf <- seq(0,((2*pi)*((TotalNumberSamplePoints/2)-1))/(TimeInterval*TotalNumberSamplePoints),(2*pi)/(TimeInterval*TotalNumberSamplePoints))
DiscreteFrequenciesSecondHalf <- seq(pi/TimeInterval, (2*pi)/(TimeInterval*TotalNumberSamplePoints), -(2*pi)/(TimeInterval*TotalNumberSamplePoints))
AngularFrequencies <- c(DiscreteFrequenciesFirstHalf, DiscreteFrequenciesSecondHalf)
##CHECK THIS IS DONE CORRECTLY, I have done 0 to ((2*pi)*((M/2)-1))/(delta*M) for the first half and then the frequencies go back down again from pi/delta to 2*pi/M*delta (goes to being one interval off 0 again)

##EACH COLUMN OF X AND Y

DiscreteSpectralDensity <- DiscreteSpectralDensityCalculation(AngularFrequencies,TimeInterval,5,theta)

autocovarianceFunction <-Re(fft(unlist(DiscreteSpectralDensity),inverse = TRUE))*(2*pi/(TimeInterval*TotalNumberSamplePoints))

sigma <- toeplitz(autocovarianceFunction)

matrixL <- chol.default(sigma)##THIS ONLY WORKS FOR SPECIFIC TIME INTERVAL AND SAMPLE POINTS COMBINATIONS - CHECK WHY

x <- matrix(rnorm(TotalNumberSamplePoints*TotalNumberSamplePoints,0,1), nrow = TotalNumberSamplePoints, ncol = TotalNumberSamplePoints)

y <- matrix(rep(0,TotalNumberSamplePoints*TotalNumberSamplePoints), nrow = TotalNumberSamplePoints, ncol = TotalNumberSamplePoints)

for(l in 1:TotalNumberSamplePoints){
  y[,l] <- t(matrixL) %*% x[,l]
}

periodogram <- matrix(rep(0,TotalNumberSamplePoints*TotalNumberSamplePoints), nrow = TotalNumberSamplePoints, ncol = TotalNumberSamplePoints)
for(k in 1:TotalNumberSamplePoints){
  periodogram[,k] <-  TimeInterval*(1/(2*pi*TotalNumberSamplePoints))*Mod((fft(y[,k])))^2 ##ADDED IN HALF FOR NEGATIVE, CHECK THIS, I'VE DONE DIVIDE BY 4 BECAUSE I THOUGHT THE HALF WOULD BE SQUARED WITH THE FOURIER TRANSFORM AND IT FIT THE JONSWAP CURVE BETTER WHEN I DID THIS

}

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

#Function to calculate Whittle liklihood for given parameters, returns minus whittle so it finds max value in optim not min

#THIS IS THE ALIASED WHITTLE LIKLIHOOD
whittleLiklihood <- function(thetaList,periodogramInput = periodogram, AngularFreq = AngularFrequencies, TimeInt = TimeInterval )
{
  discreteSpecDenseW <- DiscreteSpectralDensityCalculation(AngularFreq,TimeInt, 5, thetaList)
  whittleL <- -(log(unlist(discreteSpecDenseW))+(periodogramInput/unlist(discreteSpecDenseW)))
  return(-sum(whittleL))
}

whittleestimatematrix <- matrix(rep(0,4*TotalNumberSamplePoints),nrow = TotalNumberSamplePoints, ncol = 4)
debiasedwhittleestimatematrix <- matrix(rep(0,4*TotalNumberSamplePoints),nrow = TotalNumberSamplePoints, ncol = 4)

for(i in 1:TotalNumberSamplePoints){
  whittleestimatematrix[i,] <- optim(c(0.5,0.5,5,5), whittleLiklihood, periodogramInput= periodogram[,i], lower=c(0.1,0.1,1.01,1.01), upper=c(0.99,0.99,5,5), method="L-BFGS-B")[[1]]
}

for(i in 1:TotalNumberSamplePoints){
  debiasedwhittleestimatematrix[i,] <- optim(c(0.5,0.5,5,5), debiasedwhittleLiklihood, periodogramInput= periodogram[,i],lower=c(0.1,0.1,1.01,1.01), upper=c(0.99,0.99,5,5), method="L-BFGS-B")[[1]]
}

boxplot(debiasedwhittleestimatematrix[,1], whittleestimatematrix[,1],main = "Whittle Parameter estimates of alpha",at = c(1,2),names = c("debiased whittle", "whittle"),las = 2,col = c("blue","red"),border = "brown",horizontal = TRUE,notch = FALSE)
boxplot(debiasedwhittleestimatematrix[,2], whittleestimatematrix[,2],main = "Whittle Parameter estimates of omega_p",at = c(1,2),names = c("debiased whittle", "whittle"),las = 2,col = c("blue","red"),border = "brown",horizontal = TRUE,notch = FALSE)
boxplot(debiasedwhittleestimatematrix[,3], whittleestimatematrix[,3],main = "Whittle Parameter estimates of gamma",at = c(1,2),names = c("debiased whittle", "whittle"),las = 2,col = c("blue","red"),border = "brown",horizontal = TRUE,notch = FALSE)
boxplot(debiasedwhittleestimatematrix[,4], whittleestimatematrix[,4],main = "Whittle Parameter estimates of r",at = c(1,2),names = c("debiased whittle", "whittle"),las = 2,col = c("blue","red"),border = "brown",horizontal = TRUE,notch = FALSE)








