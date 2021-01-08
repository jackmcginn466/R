#peak frequnecy significant wave height correlation

Delta = 1/1.28 # this is 1.28 Hz data
N = 2304
daysInMonth <- c(31,28,31,30,31,30,31,31,30,31,30,31)

Hs <- c()
peakfrequencies <- c()
weeklyAverageHs <- c()


## periodogram
omega = seq(0,pi/Delta,by = 2*pi/Delta/N)
buoyName = "Wavehub"
year = 2017
month = 1
day = 1
folderName = sprintf("%s}%04d",buoyName,year)
for(l in 1:12){
  day = 1
  for(i in 1:daysInMonth[[l]]){
    fileName = sprintf("%s/%s}%04d-%02d-%02d}raw.mat",folderName,buoyName,year,month,day)
    waveData = readMat(fileName)$timeSeries
    Hs <- c(Hs, mean(4*colStdevs(waveData, na.rm = TRUE),na.rm = TRUE))
    for(j in 1:48){
      I = abs(fft(waveData[,j]))^2*(Delta/2/pi/N)
      I = I[1:(N/2+1)]
      ImatrixDay[,j] <- I
    }
    Iday <- rowMeans(ImatrixDay, na.rm = TRUE)
    peakfrequencies <- c(peakfrequencies, omega[which.max(Iday)])
    day = day + 1
  }
  
  month = month + 1
}


plot(x = peakfrequencies, y =Hs, main = "Significant Wave Height against peak frequency", xlab = "Peak frequency (Hz)", ylab = "Significant wave height (m)", xlim = c(0,max(peakfrequencies)), ylim=c(0,max(Hs)))


cor(peakfrequencies,Hs)
