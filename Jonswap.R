

##theta = (alpha,omega p,gamma,r)
## for negative or positive frequncies

jonswap <- function(w,theta){
  g <-c()
  for(i in 1:length(w))
    {
    if(w[[i]]>theta[[2]])
      {
      sigma <-0.09
      }
    else
      {
      sigma<-0.07
      }
    ##ADDED IN DIVIDE BY 2++
    g[[i]]<-(theta[[1]]*w[[i]]^(-theta[[4]])*exp(-(theta[[4]]/4)*(w[[i]]/theta[[2]])^-4)*theta[[3]]^(exp(-(1/2*sigma^2)*((w[[i]]/theta[[2]])-1)^2)))/2
  }
  return(g)
}

x <- seq(-3,3,by=0.0001)
theta1 <- c(0.7,0.7,3.3,4) ##red
theta2<-c(0.6,0.7,3.3,4)  ##blue
theta3 <-c(0.8,0.7,3.3,4) ## green
y1 <- jonswap(x,theta1)
y2 <- jonswap(x,theta2)
y3 <- jonswap(x,theta3)

plot(x = NULL,y = NULL,main="Plotting JONSWAP for varied r", xlab="Frequency", ylab="Spectral density", xlim=c(-3,3), ylim = c(0,4))
points(x, y1,pch = "." ,col = "red")
points(x, y2,pch = "." ,col = "blue")
points(x, y3,pch = "." ,col = "green")
legend(3.3, 4, legend=c("r = 3", "r = 4", "r = 5"),col=c("blue", "red", "green"), lty=1:2, cex=0.8)

