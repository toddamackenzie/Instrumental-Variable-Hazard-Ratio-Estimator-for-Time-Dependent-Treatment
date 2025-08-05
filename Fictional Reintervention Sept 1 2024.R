remove(list=ls())
library(survival)
library(cmprsk)
source("IV HR TDC Source Mar 29 2024.r")

# Settings to simulate dataset
N <- 20000
N.Ctr <- 150
rhoIVWait <- 0.5
rhoT0Wait <- 0.3
MedWait <- 10
MedSurv0 <- 10
HR <- 0.8
#CensPerc <- 0.5
set.seed(31032024)

# Variables
p.ctr <- (p <- rgamma(N.Ctr, shape=2)) /  sum(p)
N.Per.Ctr <- rmultinom(1, size=N, p=p.ctr)
Center <- rep(1:N.Ctr, times=N.Per.Ctr)
Ctr.Reint.Force <- rnorm(N.Ctr)
preIV <- Ctr.Reint.Force[Center]
Age <- 75 + 7 * rnorm(N)
Female <- runif(N) < 0.2
NotElective <- rnorm(N) < 0.1
UT0 <- runif(N)
T0 <- -log(1-UT0) * MedSurv0 / exp(.03*(Age-75) + 0.4*(NotElective-0.1)) # last part to standardize for log normal Zs
UWait <- pnorm(rhoIVWait * preIV + rhoT0Wait * qnorm(UT0) + sqrt(1 -rhoIVWait^2  -rhoT0Wait^2)*rnorm(N))
Wait <- MedWait * pmax(0.0001, -log(1-UWait))
#cor.test(T0, Wait, method="spearman")
T1 <- (T0 - Wait) / HR # HR speeds up or slows down exp-distributed event (could adopt to Weibull)
Cens <- 10 * runif(N)
Time <- ifelse(Wait < T0, Wait + T1, T0)
#Scale <- quantile(Time/Cens, 1 - CensPerc)
Cens <- Cens #* Scale
minTC <- pmin(Time, Cens)
Event <- ifelse(Time < Cens, 1, 0)
Reint <- Wait < pmin(Time,Cens)
tWait <- pmin(Wait, Time, Cens)

Ctr.Rate <- tapply(Reint, Center, sum) / tapply(tWait, Center, sum)
Ctr.Rate.SE <- sqrt(Ctr.Rate / tapply(tWait, Center, sum))
# IV at patient level
IV <- Ctr.Rate[Center]


# Start Stop Data Frame
DF0 <- data.frame(IV, Treat=0, Start=0, Stop = pmin(Wait, minTC), Event = ifelse(minTC < Wait, Event, 0), Age, Female, NotElective)
DF1 <- data.frame(IV, Treat=1, Start=Wait, Stop = minTC, Event = Event, Age, Female, NotElective)[Wait < minTC,]
DF <- rbind(DF0, DF1)



########################################
# Start of Graphic
########################################
png("Reintervention.png", width = 6, height = 6, units = 'in', res=300)

par(mfrow=c(2,2), mar=c(4,4,2.5,2))
oci <- cuminc(tWait, ifelse(Wait < minTC,1 , ifelse(Event==1,2,0)))
plot(oci$`1 1`$time, oci$`1 1`$est, type="l", axes=FALSE, xlab="Years of Followup", ylim=c(0,0.8), ylab="Cumulative %")
title("A. Mortality and\nReintervention", adj=0)
lines(survfit(Surv(minTC, Event) ~ 1), fun="F", )
axis(1)
axis(2, p <- c(0, 0.2,0.4,0.6,0.8), 100*p)
text(4, 0.20, "Reintervention", adj=0)
text(6, 0.55, "Mortality", adj=1)

# Calculate IV based on incidence of reintervnetion
Time.Cuts <- c(0:9)
Rate <- N.Ev <- Rate.SE <- matrix(nrow=length(Time.Cuts)-1, ncol=2)
for (i in 2:length(Time.Cuts)) {
  int.left <- pmax(DF$Start, Time.Cuts[i-1])
  int.right <- pmin(DF$Stop, Time.Cuts[i])
  int.Event <- ifelse(DF$Event == 1, int.left < DF$Stop & DF$Stop <= int.right, 0)
  fut <- pmax(0, int.right - int.left)
  N.Ev[i-1,] <- tapply(int.Event, DF$Treat, sum)
  Rate[i-1,] <- tapply(int.Event, DF$Treat, sum) / tapply(fut, DF$Treat, sum)
  Rate.SE[i-1,] <- sqrt(tapply(int.Event, DF$Treat, sum)) / tapply(fut, DF$Treat, sum)
}

int.midpoint <- 0.5*(Time.Cuts[-1] + Time.Cuts[-length(Time.Cuts)])
CI.Left <- Rate * exp(c(-2)/sqrt(N.Ev))
CI.Right <- Rate * exp(c(+2)/sqrt(N.Ev))
y.max <- quantile(CI.Right, 0.90)
y.min <- quantile(CI.Left, 0.10)
plot(int.midpoint, Rate[,1], axes=FALSE, type="l", log="y", xlab="Years of Followup", ylab="Mortality Rate (/yr)", pch=16, ylim=c(y.min,y.max), lwd=2)
title("B. Mortality Rate\n by Year and Treatment", adj=0)
lines(int.midpoint, Rate[,2], lwd=2, lty=3)
for (i in 1:2) 
  for(j in 2:length(Time.Cuts)) 
    lines(rep(int.midpoint[j-1],2), c(CI.Left[j-1,i], CI.Right[j-1,i]), lty=ifelse(i, 1, 3))
axis(1)
axis(2)
text(0.5, 0.085, "No Reint.", adj=0)
text(8.5, 0.15, "Reint.", adj=1)

# waterfall plot
ord <- order(Ctr.Rate)
Rts <- Ctr.Rate[ord]
plot(1:N.Ctr, Rts, axes=FALSE, xlab="Centers Ranked by Rate", ylab="Rate (#/pers-yr)", pch=16)
title("C. Reintervention Rate by Center", adj=0)
SEs <- Ctr.Rate.SE[ord]
for (i in 1:N.Ctr) lines(c(i,i), Rts[i]+c(-2,2)*SEs[i] )
axis(2)


# Center Mortality Rates vs Center Reintervention Rates
Mort.Rate <- tapply(Event, Center, sum) / tapply(minTC, Center, sum)
plot(Ctr.Rate, Mort.Rate, log="x", pch=16, cex=N.Per.Ctr/300, xlab="Reintervention Rate",
     ylab="Mortality Rate")
title("D. Mortality by Reint. Rate", adj=0)

dev.off()
########################################
# End of Graphic
########################################

ns <- table(Center)
summary(lm(Mort.Rate ~ Ctr.Rate, weight=ns))
summary(lm(log(Mort.Rate) ~ log(Ctr.Rate), weight=ns))

# IV-Wait association
(o <- coxph(Surv(tWait, tWait < pmin(Time,Cens)) ~ preIV))
(o <- coxph(Surv(tWait, tWait < pmin(Time,Cens)) ~ IV))
(chi.IV.Wait <- diff(o$loglik))
#summary(crr(tWait, ifelse(Wait < minTC,1 , ifelse(Event==1,2,0)), IV))
#(o <- coxph(Surv(tWait, tWait < pmin(Time,Cens)) ~ factor(Center)))
#diff(o$loglik)/N.Ctr

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
summary(o.Cox.0 <- coxph(Surv(Start, Stop, Event, type="counting") ~ Treat, data=DF))

summary(o.Cox.1 <- coxph(Surv(Start, Stop, Event, type="counting") ~ Treat + Age + Female + NotElective, data=DF))

DF <- DF[order(DF$Stop), ]

# Overall
(Est.IV <- IVHR.TDC(Time=DF$Stop, Status=DF$Event, X=DF$Treat, IV=DF$IV, Z=NULL, Start.Time=DF$Start, Start.Est = o.Cox.0$coef))
cbind(Est.IV$Est.log.HR, sqrt(diag(Est.IV$Cov.log.HR))) 

(t0<-Sys.time())
(Est.IV <- IVHR.TDC(Time=DF$Stop, Status=DF$Event, X=DF$Treat, IV=DF$IV, Z=DF[, c("Age", "Female", "NotElective")], Start.Time=DF$Start, Start.Est = o.Cox.1$coef))
cbind(Est.IV$Est.log.HR, sqrt(diag(Est.IV$Cov.log.HR))) 
Sys.time()-t0

# Less than 3 years
(Est.IV <- IVHR.TDC(Time=DF$Stop, Status=DF$Event & DF$Stop<3, X=DF$Treat, IV=DF$IV, Z=NULL, Start.Time=DF$Start, Start.Est = o.Cox.0$coef))
cbind(Est.IV$Est.log.HR, sqrt(diag(Est.IV$Cov.log.HR))) 

(Est.IV <- IVHR.TDC(Time=DF$Stop, Status=DF$Event & DF$Stop<3 X=DF$Treat, IV=DF$IV, Z=DF[, c("Age", "Female", "NotElective")], Start.Time=DF$Start, Start.Est = o.Cox.1$coef))
cbind(Est.IV$Est.log.HR, sqrt(diag(Est.IV$Cov.log.HR))) 

# 5 years +
(Est.IV <- IVHR.TDC(Time=DF$Stop, Status=DF$Event & DF$Stop>=3, X=DF$Treat, IV=DF$IV, Z=NULL, Start.Time=DF$Start, Start.Est = o.Cox.0$coef))
cbind(Est.IV$Est.log.HR, sqrt(diag(Est.IV$Cov.log.HR))) 

(Est.IV <- IVHR.TDC(Time=DF$Stop, Status=DF$Event & DF$Stop>=3, X=DF$Treat, IV=DF$IV, Z=DF[, c("Age", "Female", "NotElective")], Start.Time=DF$Start, Start.Est = o.Cox.1$coef))
cbind(Est.IV$Est.log.HR, sqrt(diag(Est.IV$Cov.log.HR))) 

# by 1st or >1 year of reinter
DF0 <- data.frame(IV, Treat1=0, Treat2=0, Start=0,      Stop = pmin(Wait, minTC), Event = ifelse(minTC < Wait, Event, 0), Age, Female, NotElective)
DF1 <- data.frame(IV, Treat1=1, Treat2=0, Start=Wait,   Stop = pmin(Wait+1, minTC), Event = Event, Age, Female, NotElective)[Wait < minTC,]
DF2 <- data.frame(IV, Treat1=0, Treat2=1, Start=Wait+1, Stop = minTC, Event = Event, Age, Female, NotElective)[Wait < minTC,]
DF <- rbind(DF0, DF1, DF2)

(Est.IV <- IVHR.TDC(Time=DF$Stop, Status=DF$Event, X=DF[,c("Treat1", "Treat2")], IV=cbind(DF$IV,DF$IV), Z=NULL, Start.Time=DF$Start, Start.Est = o.Cox.0$coef))
cbind(Est.IV$Est.log.HR, sqrt(diag(Est.IV$Cov.log.HR))) 

(Est.IV <- IVHR.TDC(Time=DF$Stop, Status=DF$Event, X=DF$Treat, IV=cbind(DF$IV, DF$IV), Z=DF[, c("Age", "Female", "NotElective")], Start.Time=DF$Start, Start.Est = o.Cox$coef))
cbind(Est.IV$Est.log.HR, sqrt(diag(Est.IV$Cov.log.HR))) 




####################################################################
# Basement




R <- 1000
N <- 300
beta <- 0.25
Z <- rep(NA, R)
for (r in 1:R) {
  X <- rnorm(N)
  Time0 <- rexp(N) / exp(beta * X)
  Cens <- runif(N)
  scale <- quantile(Time0/Cens, 0.46)
  Cens <- scale * Cens
  Time <- pmin(Time0, Cens)
  Status <- ifelse(Time0<Cens, 1, 0)
  ocs <- summary(coxph(Surv(Time, Status) ~ X))
  Z[r] <- ocs$coef[1,4]
}
mean(Z>= 1.96)



R <- 1000
N <- 300
beta <- log(1.8)
Z <- IZ <- rep(NA, R)
for (r in 1:R) {
  X <- rnorm(N)
  Fat <- runif(N) < 0.5
  Time0 <- rexp(N) / exp(beta * X)
  Cens <- runif(N)
  scale <- quantile(Time0/Cens, 0.25)
  Cens <- scale * Cens
  Time <- pmin(Time0, Cens)
  Status <- ifelse(Time0<Cens, 1, 0)
  ocs <- summary(coxph(Surv(Time, Status) ~ X))
  ocsi <- summary(coxph(Surv(Time, Status) ~ X*Fat))
  Z[r] <- ocs$coef[1,4]
  IZ[r] <- ocsi$coef[3,4]
}
mean(Z>= qnorm(1-0.05/700))
mean(IZ>= qnorm(1-0.05/700))

