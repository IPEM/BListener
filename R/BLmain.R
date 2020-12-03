BLmain <- function(u, meter, tg = NULL, outt = 0.7, ssa = 100,
                   g = 0, V = 0.00001, W = 0.01) {

if (!is.numeric(u)) stop("u is not numeric.")
if (!is.numeric(meter)) stop("meter is not numeric.")

ml <- length(meter)
meter <- log2(meter)
if (is.null(tg)) {
  su <- u[which(u < 1200)]
  su <- su[which(su > 150)]
  x <- stats::kmeans(su, centers = ml, nstart = 25)$centers
  tg <- x[which.min(x)]}
v <- sqrt(V)
w <- sqrt(W)
Vt <- diag(ml) * v^2
Wt <- diag(ml) * w^2
str_var <- sprintf("r=%s, g=%s, o=%s",
                    toString(round(V/W, digits = 5)),
                    toString(round(g, digits = 5)),
                    toString(round(outt, digits = 5)))
IOIguess <- matrix( meter + log2(tg), nrow = ml, ncol = 1)
Meter0 <- matrix(meter, nrow = ml, ncol = 1)
Varguess <- diag(ml) * v
rsa <- ssa / 1000
uu <- round(u*rsa)
len <- tail(cumsum(uu),n = 1)
CtOt <- IOIguess
VarCtOt <- Varguess
Ot <- IOIguess
I <- diag(ml); A <- diag(ml); Phi <- diag(ml)

ObservedTrials  <- list()
ObservedTrialsCouples <- list()
Results <- list(
  Observed     = matrix(NA, nrow = len, ncol = ml),
  Prediction   = matrix(NA,nrow = len, ncol = ml),
  Variance     = matrix(NA,nrow = len, ncol = ml),
  Innovation   = matrix(NA,nrow = len, ncol = ml),
  AdaptiveGain = matrix(NA,nrow = len, ncol = ml),
  Narrative    = matrix(NA,nrow = length(u) , ncol = 1),
  Outliers     = matrix(NA,nrow = len, ncol = 1),
  Timeline     = c(1:len),
  String       = str_var)
ti <- 1

for (i in 1:(length(uu))) {
  Observe <- as.numeric(u[i])
  difference  <- abs(log2(Observe) - CtOt)
  mindifference <- min(difference)
  A  <- diag(ml) * 0
  if(mindifference < outt) {
    index <- which.min(difference)
    Results$Narrative[i] <- index
    Ot[index] <- log2(Observe)
    A <- diag(ml) * 0
    A[index, index] <- 1
    } else {
      Results$Outliers[ti] <- Observe
      print(paste('--------> Outlier detected, in ms:', Observe))
    }; if(Observe <= 0) print('--------> Observe is 0')

  for (Time in 1:uu[i]) {
    vt <- as.matrix(rnorm(ml,0,v), nrow = ml, ncol = 1)
    wt <- as.matrix(rnorm(ml,0,w), nrow = ml, ncol = 1)
    Ctm1Otm1 <- CtOt
    VarCtm1Otm1 <- VarCtOt
    x <- Ctm1Otm1 - Ctm1Otm1[length(ml)]
    d <- (Meter0 - x) - mean(Meter0 - x)
    exo <- g*d
    CtOtm1 <- Phi %*% Ctm1Otm1 + vt + exo
    VarCtOtm1 <- Phi %*%VarCtm1Otm1 %*% t(Phi) + Vt
      if(Time == 1) {
        OOt <- A %*% Ot + wt
        Et <- OOt - (A %*% CtOtm1)
        R <- (VarCtOtm1)
        Q <- (VarCtOtm1 + Wt)
        K <- (R %*% t(A)) / (A %*% Q %*% t(A))
        K[is.nan(K)] <- 0
        CtOt <- CtOtm1 + (K %*% Et)
        VarCtOt <- (I - (K %*% A)) %*% R
        Ctm1Otm1 <- CtOt
        VarCtm1Otm1 <- VarCtOt
        Results$Observed[ti,] <- t(A %*% Ot)
        Results$Innovation[ti,] <- t(A %*% Et)
        Results$AdaptiveGain[ti,] <- t(A %*% diag(K))
      } else {
        CtOt <- CtOtm1
        VarCtOt <- VarCtOtm1
      }
    Results$Prediction[ti,] <- t(CtOt)
    Results$Variance[ti,] <- t(diag(VarCtOt))
    ti <- ti + 1
  }
}

Results$Observed[as.vector(Results$Observed) == 0 ] <- NA
Results$Innovation[as.vector(Results$Innovation) == 0 ] <- NA
Results$AdaptiveGain[as.vector(Results$AdaptiveGain) == 0 ] <- NA
Results$Prediction[as.vector(Results$Prediction) == 0 ] <- NA
Results$Outliers[as.vector(Results$Outliers) == 0 ] <- NA
Results$Tempo  <- (rowMeans(sweep(Results$Prediction,2,t(Meter0) ), na.rm = TRUE))+1
Results$Tempo[Results$Tempo==0] <-NA
Results$parameters <- list(u, meter, tg, outt, g, V, W)
names(Results$parameters) <- c("u","meter", "tg", "outt", "g", "V", "W")
class(Results) <- "BL"
return(Results)
}
