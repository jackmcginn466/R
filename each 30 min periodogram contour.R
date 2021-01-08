#Periodogram Throughout the year contour plot

library(R.matlab) # for loading mat files
library(resample) # for std of columns
#######################################################
### example of loading a record
## make sure the working directory is set to the correct location (if it is different you may need to add folder names etc)
# note you don't have to use sprintf, but this will help you to do this in a loop if you want


#dim(waveData)

Delta = 1/1.28 # this is 1.28 Hz data
N = 2304
daysInMonth <- c(31,28,31,30,31,30,31,31,30,31,30,31)

Imatrixyear <- matrix(rep(0,1153*sum(daysInMonth)*48),nrow = 1153,ncol =sum(daysInMonth)*48)

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
    for(j in 1:48){
      I = abs(fft(waveData[,j]))^2*(Delta/2/pi/N)
      I = I[1:(N/2+1)]
      Imatrixyear[,j+((i-1)*48)+(sum(daysInMonth[0:(l-1)])*48)] <- I
    }
    day = day + 1
  }
  month = month + 1
}



days = seq(0.5,365*24,0.5)



filled.contour( x = days , y = omega , z = t(Imatrixyear), zlim = c(0,0.02), nlevels = 100,col = cm.colors(99))
#zlim = range(t(Imatrixyear))
dim(Imatrixyear)
length(days)
length(omega)

