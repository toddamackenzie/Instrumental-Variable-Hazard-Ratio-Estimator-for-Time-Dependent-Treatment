
library(survival)
library(MASS)
library(rootSolve)
library(Rcpp)

cppFunction('
double HRSums(NumericVector StartTime, NumericVector StopTime, IntegerVector Event,
NumericVector HR, NumericVector W) {
  // Assumes StopTime is sorted smallest to largest, and then Events =1 before 0
  int n = StartTime.size();
  double EE = 0, sumHR, sumWHR ;
  for (int i = 0; i < n; i++) {
    if (Event(i) == 1) {
      sumHR = 0;
      sumWHR = 0 ;
      for (int j = i; j < n; j++) {
        if (StartTime(j) < StopTime(i)) {
          sumHR += HR(j) ;
          sumWHR += W(j) * HR(j) ;
          }
      }
      EE += W(i) - sumWHR / sumHR ;
    }
  }
  return EE;
}')


cppFunction('
List ff(NumericMatrix mm) {
  List Out(2) ;
  Out["nrow"] = mm.nrow() ;
  Out["ncol"] = mm.ncol() ;
  return Out ;
}')

cppFunction('
List CoVarCalc(NumericVector StartTime, NumericVector StopTime, IntegerVector Event,
NumericVector HR, NumericVector X1, NumericVector W1, NumericVector X2, NumericVector W2) {
  // Assumes StopTime is sorted smallest to largest, and then Events =1 before 0
  int n = StartTime.size() ;
  double VarEE = 0, Deriv = 0, 
    sumHR, sumW1HR, sumW2HR, sumX2HR, sumW1X2HR ;
  List Out(2) ;
  for (int i = 0; i < n; i++) {
    if (Event(i) == 1) {
      sumHR = 0; sumW1HR = 0 ; sumW2HR = 0 ; sumX2HR = 0 ; sumW1X2HR = 0 ;
      for (int j = i; j < n; j++) {
        if (StartTime(j) < StopTime(i)) {
          sumHR += HR(j) ;
          sumW1HR += W1(j) * HR(j) ;
          sumW2HR += W2(j) * HR(j) ;
          sumX2HR += X2(j) * HR(j) ;
          sumW1X2HR += W1(j) * X2(j) * HR(j) ;
        }
      }
    VarEE += (W1(i) - sumW1HR/sumHR) * (W2(i) - sumW2HR/sumHR);
    Deriv -= sumW1X2HR / sumHR - sumW1HR * sumX2HR / pow(sumHR, 2) ;
    }
  }
  Out["varEE"] = VarEE ;
  Out["Deriv"] = Deriv  ;
  return Out ;
}')


IVHR.TDC <- function(Time, Status, X, IV, Z=NULL, Start.Time=NULL, Start.Est=NULL, Grid.Est=NULL) {
  # Start.Time is the start of the interval (i.e. counting process notation), by default it is set to zero
  # X is the treatment being recieved at beginning of interval, e.g. 1 for treatment, 0 for control
  # W is the instrument, e.g. 1 for randomized to treatment, 0 for randomized to control
  n <- length(Time)
  IV <- data.frame(IV)
  if (is.null(Z)) Z <- matrix(nrow=n, ncol=0)
  Z <- data.frame(Z)
  n.iv <- dim(IV)[2]
  n.z <- dim(Z)[2]
	if (is.null(Start.Time)) Start.Time <- rep(0, n)
	DF <- data.frame(Start.Time, Time, Status, X, IV=IV, Z=Z)
	ord <- order(Time)
	DF <- DF[ord, ]
	IV <- as.matrix(as.matrix(IV)[ord, ])
	Z <- as.matrix(as.matrix(Z)[ord, ])
	Est.Equat <- function(beta) {
	  XZ <- cbind(DF$X, Z)
	  DF$HR <- exp(as.matrix(XZ) %*% beta)
	  o <- c()
		for (i in 1:dim(IV)[2])
		  o <- c(o, HRSums(DF$Start.Time, DF$Time, DF$Status, DF$HR, IV[,i]))
	  if (dim(Z)[2] > 0) {
  	  for (i in 1:dim(Z)[2])
  	    o <- c(o, HRSums(DF$Start.Time, DF$Time, DF$Status, DF$HR, Z[,i]))
	  }
  	o
	}
	CovarEstMatrix <- function(beta) {
	  X.Z <- cbind(DF$X, Z)
	  IV.Z <- cbind(IV, Z)
	  Degree <- dim(IV)[2] + dim(Z)[2]
	  DF$HR <- exp(as.matrix(X.Z) %*% beta)
	  Cov <- Deriv <- matrix(nrow=Degree, ncol=Degree)
	  for (i in 1:Degree) {
	    for (j in 1:Degree) {
	      oCov <- CoVarCalc(DF$Start.Time, DF$Time, DF$Status, DF$HR, 
	                 X.Z[,i], IV.Z[,i], X.Z[,j], IV.Z[,j])
	      Cov[i,j] <- oCov$varEE
	      Deriv[i,j] <- oCov$Deriv
	    }
	  }
	  Inv.Deriv <- ginv(Deriv)
	  Inv.Deriv %*% Cov %*% t(Inv.Deriv)
	}
	if (is.null(Grid.Est)) {
  	if (is.null((Start.Est))) Start.Est <- rep(0,n.iv+n.z)
  	out.solution <- multiroot(Est.Equat, start=Start.Est)
  	beta.hat <- rep(NA, n.iv+n.z)
  	beta.Cov <- matrix(NA, nrow=n.iv+n.z, ncol=n.iv+n.z)
  	if(!is.na(out.solution$estim.precis)) {
  	  if(out.solution$estim.precis < 0.00001) {
  	    beta.hat <- out.solution$root
  	    beta.Cov <- CovarEstMatrix(beta.hat)
  	  }
  	}
  	o <- list(Est.log.HR=beta.hat, Cov.log.HR=beta.Cov) #, SE=SE, SE.2=SE.2)
	}
	if (!is.null(Grid.Est)) {
	  EE <- c()
	  for (est in Grid.Est) {
	    ee <- Est.Equat(est)
	    EE <- c(EE,ee)
	  }
	o <- cbind(Grid.Est, EE)  
	}
	o
}

