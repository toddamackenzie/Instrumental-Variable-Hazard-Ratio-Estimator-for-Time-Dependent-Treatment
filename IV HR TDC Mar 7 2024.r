
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


IVHR.TDC <- function(Time, Status, X, IV, Z=NULL, Start.Time=NULL, Start.Est=NULL) {
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
	#  to.min <- function(beta) {Est.Equat(beta)^2}
	  #  out.solution <- nlminb(start=0, to.min)
	  #  beta.hat <- ifelse(out.solution$convergence == 0, out.solution$par, NA)
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
	list(Est.log.HR=beta.hat, Cov.log.HR=beta.Cov) #, SE=SE, SE.2=SE.2)
}

######################################################################
# Simulations of Marginal TDCPH Cox model
# gaussian copula achieves association between Wait and 
# Time-to-event, T0, under no treatment (infinite wait)
######################################################################

rset <- function(x) {
  u <- runif(1)
  u <- ceiling(u * length(x))
  x[u]
}

SimMarginal <- function(Ns=c(500,2000), MedSurv0=10, rhoIVWaits=c(0.2,0.4,0.6),
    rhoT0Waits= c(-0.3, 0, +0.3), MedWaits = c(1,2,5), 
    HRs = c(1/3,0.50,2/3,1,1.5,2,3), CensPercs=c(0.50,0.75)) {
  # sum of squares of rho's must be 1 or less
  clock0 <- Sys.time()
  Out <- c()
  while(TRUE) {
    N <- rset(Ns)
    rhoIVWait <- rset(rhoIVWaits)
    rhoT0Wait <- rset(rhoT0Waits)
    MedWait <- rset(MedWaits)
    HR <- rset(HRs)
    CensPerc <- rset(CensPercs)
    # Variables
    IV <- rnorm(N)
    Z1 <- rnorm(N)
    Z2 <- rnorm(N)
    Z3 <- rnorm(N)
    UT0 <- runif(N)
    T0 <- -log(1-UT0) * MedSurv0 / exp(Z1 - Z3) / exp(1) # last part to standardize for log normal Zs
    UWait <- pnorm(rhoIVWait * IV + rhoT0Wait * qnorm(UT0) + sqrt(1 -rhoIVWait^2  -rhoT0Wait^2)*rnorm(N))
    Wait <- MedWait * pmax(0.0001, -log(1-UWait))
    #cor.test(T0, Wait, method="spearman")
    T1 <- (T0 - Wait) / HR # HR speeds up or slows down exp-distributed event (could adopt to Weibull)
    Cens <- runif(N)
    Time <- ifelse(Wait < T0, Wait + T1, T0)
    Scale <- quantile(Time/Cens, 1 - CensPerc)
    Cens <- Cens * Scale
    minTC <- pmin(Time, Cens)
    Event <- ifelse(Time < Cens, 1, 0)
    # IV-Wait association
    tWait <- pmin(Wait, Time, Cens)
    o <- coxph(Surv(tWait, tWait < pmin(Time,Cens)) ~ IV)
    chi.IV.Wait <- diff(o$loglik)
    # Final data step
    DF0 <- data.frame(IV, Treat=0, Start=0, Stop = pmin(Wait, minTC), Event = ifelse(minTC < Wait, Event, 0), Z1, Z2, Z3)
    DF1 <- data.frame(IV, Treat=1, Start=Wait, Stop = minTC, Event = Event, Z1, Z2, Z3)[Wait < minTC,]
    DF <- rbind(DF0, DF1)
    # Some statistics of relevance
    ne0 <- sum(DF0$Event)
    ne1 <- sum(DF1$Event)
    total.time.X0 <- sum(DF0$Stop - DF0$Start)
    total.time.X1 <- sum(DF1$Stop - DF1$Start)
    # Housekeeping
    DF$Start <- floor(DF$Start*10^4)/10^4
    DF$Stop <- ceiling(DF$Stop*10^4)/10^4
    DF <- DF[DF$Start < DF$Stop, ]
    # Estimation 
    o.Cox <- coxph(Surv(Start, Stop, Event, type="counting") ~ Treat + Z1 + Z2 + Z3, data=DF)
    Est.Cox.Treat <- o.Cox$coef[1]
    SE.Cox.Treat <- sqrt(o.Cox$var[1,1])
    DF <- DF[order(DF$Stop), ]
    # The estimator showed bias if chi < 10
    Est.IV.Treat <- SE.IV.Treat<- NA
    if(chi.IV.Wait >= 10) { 
      Est.IV <- IVHR.TDC(Time=DF$Stop, Status=DF$Event, X=DF$Treat, IV=DF$IV, Z=DF[, c("Z1", "Z2", "Z3")], Start.Time=DF$Start, Start.Est = o.Cox$coef)
      Est.IV.Treat <- Est.IV$Est.log.HR[1]
      SE.IV.Treat <- sqrt(Est.IV$Cov.log.HR[1,1])  
    }  
    NewRow <- c(Est.Cox.Treat, SE.Cox.Treat, Est.IV.Treat, SE.IV.Treat, ne0, ne1, chi.IV.Wait, N, rhoIVWait, rhoT0Wait, MedWait, HR, CensPerc)
    Out <- rbind(Out, NewRow)
    #
    clock1 <- Sys.time()
    if (as.numeric(clock1) - as.numeric(clock0) > 60) {
      filename <- paste("../Simulation Results/Out", Sys.Date(),".txt", sep="")
      write.table(Out, sep="\t", row.names=FALSE, file=filename)
      print(dim(Out)[1])
      print(Sys.time())
      clock0 <- clock1
    }
  }
}


SimMarginal()

Results <- read.delim(filename, row.names = NULL)
dimnames(Results)[[2]] <- c("Est.Cox", "Est.IV", "ne0", "ne1", "chi.IV.Wait", "N", "rhoIVWait", "rhoT0Wait", "MedWait", "HR", "CensPerc")
dimnames(Results)[[2]] <- c("Est.Cox", "SE.Cox.Treat", "Est.IV", "SE.IV.Treat", "ne0", "ne1", "chi.IV.Wait", "N", "rhoIVWait", "rhoT0Wait", "MedWait", "HR", "CensPerc")

table(Miss <- is.na(Results$Est.IV))
summary(glm(Miss ~ ne0 + ne1 + N + chi.IV.Wait + rhoT0Wait + MedWait + CensPerc, data=Results, family=binomial))$coef

summary(glm(Miss ~ ne0 + ne1 + chi.IV.Wait + rhoT0Wait + MedWait, data=Results, family=binomial))$coef

Figure <- function(Subset = NULL) {
  if (is.null(Subset)) Subset <- rep(TRUE, dim(Results)[1])
  Results <- Results[Subset,]
  par(mar=c(4,4,1,0.5))
  par(mfrow=c(3,3))
  
  HRs <- sort(unique(Results$HR))
  rhoIVWaits <- sort(unique(Results$rhoIVWait))
  rhoT0Waits <- sort(unique(Results$rhoT0Wait))
  for (rhoT0Wait in rhoT0Waits) {
    for (rhoIVWait in rhoIVWaits) {
      keep <- Results$rhoT0Wait == rhoT0Wait & Results$rhoIVWait == rhoIVWait
      #make a boxplot type
      dist.IV <- aggregate(Results$Est.IV[keep], list(Results$HR[keep]), quantile, c(0.25, 0.5,0.75), na.rm=TRUE)[,2]
      dist.Cox <- aggregate(Results$Est.Cox[keep], list(Results$HR[keep]), quantile, c(0.25, 0.5,0.75), na.rm=TRUE)[,2]
      plot(log(HRs), dist.IV[,"50%"], xlab="", ylab="Estimate", pch=16, col=3, ylim=c(-1.75,+1.75), axes=FALSE)
      axis(1, log(HRs), round(HRs,2))
      yy <- c(0.33, 0.50, 1, 2, 3)
      axis(2, log(yy), yy)
      text(0, -1.4, "True HR", adj=0)
      title(paste("Confounding ", rhoT0Wait, "IV strength", rhoIVWait))
      points(log(HRs), dist.IV[,"25%"], col=3)
      points(log(HRs), dist.IV[,"75%"], col=3)
      for (i in 1:length(HRs)) lines(rep(log(HRs[i]),2 ), dist.IV[i, c("25%","75%")],  col=3, lwd=2)
      points(log(HRs)+0.03, dist.Cox[,"50%"], col=2, pch=16)
      points(log(HRs)+0.03, dist.Cox[,"25%"], col=2)
      points(log(HRs)+0.03, dist.Cox[,"75%"], col=2)
      for (i in 1:length(HRs)) lines(rep(log(HRs[i])+0.03,2 ), dist.Cox[i, c("25%","75%")], col=2, lwd=2)
      abline(a=0,b=1)
    }  
  }
}  


png("sim125.png", width = 8, height = 5, units = 'in', res=300)
Figure(Results$N==500 & Results$CensPerc==0.75)
dev.off()

png("sim500.png", width = 8, height = 5, units = 'in', res=300)
Figure(Results$N==2000 & Results$CensPerc==0.50)
dev.off()


CI.IV.Coverage <- (Results$Est.IV - 1.96*Results$SE.IV.Treat < log(Results$HR)) & log(Results$HR) < (Results$Est.IV + 1.96*Results$SE.IV.Treat)

mean(CI.IV.Coverage, na.rm=TRUE)
summary(glm(CI.IV.Coverage ~ ne0 + ne1 + N + chi.IV.Wait + factor(rhoT0Wait) + MedWait + CensPerc, data=Results, family=binomial))$coef
summary(glm(CI.IV.Coverage ~ ne0 + ne1 + factor(rhoT0Wait) + MedWait + CensPerc, data=Results, family=binomial))$coef

tapply(CI.IV.Coverage, Results$rhoT0Wait, mean, na.rm=TRUE)
tapply(CI.IV.Coverage, Results$rhoIVWait, mean, na.rm=TRUE)
tapply(CI.IV.Coverage, Results$MedWait, mean, na.rm=TRUE)
tapply(CI.IV.Coverage, Results$CensPerc, mean, na.rm=TRUE)
tapply(CI.IV.Coverage, cut(Results$ne1, c(0, 50, 100, Inf)), mean, na.rm=TRUE)
tapply(CI.IV.Coverage, Results$N, mean, na.rm=TRUE)
tapply(CI.IV.Coverage, cut(Results$HR, c(-Inf, 0.99,1.01, Inf)), mean, na.rm=TRUE)

CI.Cox.Coverage <- (Results$Est.Cox - 1.96*Results$SE.Cox.Treat < log(Results$HR)) & log(Results$HR) < (Results$Est.Cox + 1.96*Results$SE.Cox.Treat)

mean(CI.Cox.Coverage, na.rm=TRUE)
summary(glm(CI.Cox.Coverage ~ ne0 + ne1 + N + chi.IV.Wait + factor(rhoT0Wait) + MedWait + CensPerc, data=Results, family=binomial))$coef
summary(glm(CI.Cox.Coverage ~ ne0 + ne1 + factor(rhoT0Wait) + MedWait + CensPerc, data=Results, family=binomial))$coef

tapply(CI.Cox.Coverage, Results$rhoT0Wait, mean, na.rm=TRUE)
tapply(CI.Cox.Coverage, Results$rhoIVWait, mean, na.rm=TRUE)
tapply(CI.Cox.Coverage, Results$MedWait, mean, na.rm=TRUE)
tapply(CI.Cox.Coverage, Results$CensPerc, mean, na.rm=TRUE)
tapply(CI.Cox.Coverage, cut(Results$ne1, c(0, 100, Inf)), mean, na.rm=TRUE)
tapply(CI.Cox.Coverage, cut(Results$HR, c(-Inf, 0.99,1.01, Inf)), mean, na.rm=TRUE)

