BLpost <- function(R) {
leaveout  <- -(1:100)
print(paste('Leave out the first', tail(leaveout, n=1), "samples"))
Post <- list()

#Fluctuation1
fluct <- function(x) sds=mean(abs(x),na.rm=TRUE)
Fl <- sapply(as.data.frame(R$Innovation[leaveout,]),fluct)
Post$Fluctuation1 <- Fl

#Stability
stabi <- function(x) sds=sd((x),na.rm=TRUE)
St  <- sapply(as.data.frame(R$Prediction[leaveout,]),stabi)
Post$Stability <- St

#Fluctuation^2
#fluct2 <- function(x) sds=sd((x),na.rm=TRUE)
Fl <- sapply(as.data.frame((R$Innovation[leaveout,])),stabi)
Post$Fluctuation2 <- Fl

#Total time of outliers
Out <- R$Outliers[leaveout,]
Out[is.na(as.vector(R$Outliers[leaveout,]))] = 0
Post$OutlierSum <- sum(Out)

# Constancy
constancy   <- function(x) {
  x <- t(x)
  dx <- diff(x)
  sds = sapply(as.data.frame(t(dx)),sd,na.rm=TRUE)
  z=det(diag(sds))
  return(z)
}

rpred <-  as.matrix(R$Prediction[leaveout,])
C <- constancy(rpred)
Post$Constancy <- C

#Parameters
Post$parameters = R$parameters
return(Post)
}
