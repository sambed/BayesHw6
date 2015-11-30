---
title: "bayeshw6"
author: "sambed-x"
date: "November 22, 2015"
output: html_document
---

##1a.

```{r}
par(mfrow=c(1,2))

# load("C:/Users/Sapana/Desktop/R/hw6.RData")
## y = faults
## x= fabric length

library(rstan)
mod3code='
  data{
    int<lower=0> n;
    int x[n];   //fabric length
    int y[n];    //faults
  }

  parameters{
    real beta0;
    real beta1;
  }

transformed parameters{
    real<lower=0> lambda[n];
    for(i in 1:n){
      lambda[i] <- exp(beta0 + beta1*(x[i]));
    }
}

  model{
    for(i in 1:n){
      y[i] ~ poisson(lambda[i]);
    }
    beta0 ~ normal(0,.99);
    beta1 ~ normal(0,.085);
  }

'

mod3data = with(fabric,list(n=nrow(fabric), x=fabric$len, y=fabric$faults))
mod3fit = stan(model_code=mod3code, data=mod3data, seed=889)

print(mod3fit,pars=c("beta0","beta1"),digits=3)
```

##1b.
```{r}
# Now let's extract the coefficient draws
par(mfrow=c(1,2))
b0.1 <- extract(mod3fit)$beta0
b1.1 <- extract(mod3fit)$beta1
plot(b0.1,b1.1, main = 'Beta1 vs. Beta0')
nsim <- length(b0.1)

ww <- seq(122,952,by=1)
n <-length(ww)

lambdas.1 <- matrix(numeric(n*length(ww)),nrow=n)


predlambda <- function(length,b0,b1){
  return(exp(b0+b1*(length)))
}

for(i in 1:n){
  lambdas.1[i,] <- apply(as.matrix(ww),1,predlambda,b0=b0.1[i],b1=b1.1[i])
}

# Obtain the conditional means and plot them with the data
cmean1 <- colMeans(lambdas.1)
plot(jitter(fabric$len,amount=.1),fabric$faults,pch=16,bty="l",xlab="Length",ylab="# of Faults", main ='Conditional Means/Prediciton Interval Plots')
lines(ww,cmean1,col='blue')

# Obtain 95% credible intervals for conditional means, and plot
low.cmean1 <- apply(lambdas.1,2,quantile,.025)
lines(ww,low.cmean1,lty=2,col='red')
high.cmean1 <- apply(lambdas.1,2,quantile,.975)
lines(ww,high.cmean1,lty=2, col= 'red')


### Now do the same for the posterior predictive intervals
predy <- function(length,b0,b1){
  return(rpois(1,exp(b0+b1*(length))))
}

# Construct a matrix to hold the predicted values
ys <- matrix(numeric(n*length(ww)),nrow=n)
for(i in 1:n){
  ys[i,] <- apply(as.matrix(ww),1,predy,b0=b0.1[i],b1=b1.1[i])
}

# Obtain 95% prediction intervals, and plot
low.y <- apply(ys,2,quantile,.025)
lines(ww,low.y,lty=1,col="green ")
high.y <- apply(ys,2,quantile,.975)
lines(ww,high.y,lty=1,col="green")


# 3 points lie outside of the PI i.e about 9%

###Checking Overdispersion


## I'll cut the
## data into 6 equal-length intervals based on length, and examine
## the means and variances in each group. If the variances far
## exceed the means, we have evidence of overdispersion relative
## to the Poisson model.

lengthgroups <- cut(fabric$len,breaks=6)
res <- cbind(by(fabric$faults,lengthgroups,mean),by(fabric$faults,lengthgroups,var))
colnames(res) <- c("Mean","Var")
res

```

In majority of the groups the variance is greater than mean. Therefore there is evidence of overdispersion.


##1c.

```{r}

### Bayes model with proper reference prior, with random intercept
mod5_code='
  data{
    int<lower=0> n;
    int length[n];
    int faults[n];
  }
  
  parameters{
    real b0;         // interceot
    real b1;         // slope 
    real tau[n];  //random intercept adjustment
    real<lower=0> sig; // std of tau, lower bounded by zero
  }
  
   transformed parameters{
    real<lower=0> lambda[n];
      for(i in 1:n){
        lambda[i] <- exp(b0 + b1*(length[i]) + tau[i]);
      }
  }
  model{
    for(i in 1:n){
      faults[i] ~ poisson(lambda[i]);
      tau[i] ~ normal(0,sig);
    }
    b0 ~ normal(0,.9);
    b1 ~ normal(0,.085);
    sig ~ gamma(0.1,0.1);
  }
'

mod5_data = with(fabric,list(n=nrow(fabric),length=fabric$len, faults=fabric$faults))

mod5_fit = stan(model_code=mod5_code, data=mod5_data,seed = 1000)

print(mod5_fit,pars=c("b0","b1",'sig'),digits=3)
```


##1d.

```{r}
# Now let's extract the coefficient draws
## Examine the fit
par(mfrow=c(1,2))
b0.2 <- extract(mod5_fit)$b0
b1.2 <- extract(mod5_fit)$b1
sig <- extract(mod5_fit)$sig
plot(b0.2,b1.2, main = 'Beta1 vs. Beta0(Int. Model)')

# Construct a matrix to hold the predicted means 
# based on each MCMC draw
n <- length(b0.2)
lambdas.2 <- matrix(numeric(n*length(ww)),nrow=n)

predlambda2 <- function(length,b0,b1,sig){
  return(exp(b0+b1*(length)))
}

for(i in 1:n){
  lambdas.2[i,] <- apply(as.matrix(ww),1,
                         predlambda2,b0=b0.2[i],b1=b1.2[i])                         
}

# Obtain the conditional means and plot them with the data
cmean2 <- colMeans(lambdas.2)
plot(jitter(fabric$len,amount=.1),fabric$faults,pch=16,bty="l",
     xlab="Length",ylab="Number of Faults",main='Cond Means/PI Plots for Intercept Model')
lines(ww,cmean2)

# Obtain 95% credible intervals for conditional means, and plot
low.cmean2 <- apply(lambdas.2,2,quantile,.025)
lines(ww,low.cmean2,lty=2,col='red')
high.cmean2 <- apply(lambdas.2,2,quantile,.975)
lines(ww,high.cmean2,lty=2,col='red')


### Now do the same for the posterior predictive intervals
predy2 <- function(width,b0,b1,sig){
  return(rpois(1,exp(b0+b1*(width-26)+rnorm(1,0,sig))))
}

# Construct a matrix to hold the predicted values
ys2 <- matrix(numeric(n*length(ww)),nrow=n)
for(i in 1:n){
  ys2[i,] <- apply(as.matrix(ww),1,predy2,b0=b0.2[i],b1=b1.2[i],sig=sig[i])
}


# Obtain 95% prediction intervals, and plot
low.y2 <- apply(ys2,2,quantile,.025)
lines(ww,low.y2,lty=3,col="green")
high.y2 <- apply(ys2,2,quantile,.975)
lines(ww,high.y2,lty=3,col="green")

```


Yes, adding random intercept improves the model. Our PI of the model with random intercept includes all most all data points except for one. Although, overfiting might be of a concern here.

##2a.

```{r}
##EDA
# load("C:/Users/Sapana/Desktop/R/crab.RData")
par(mfrow=c(1,1))
crab$pa <-ifelse(crab$satell>1,1,0)
head(crab)
with(crab,plot(satell~width, main = 'Satellite vs. Width'))
```

```{r}

invlogit <- function(x) { exp(x)/(1+exp(x))}
logit <- function(x) { log(x/(1-x))}

## "noninformative" priors

mod1code <- '
  data{
    int<lower=0> n;             // sample size
    int<lower=0,upper=1> y[n];  // presence or absence of satellite
    real x[n];                  // width
  }
  parameters{
    real beta0;
    real beta1;
  }
  model{
    for(i in 1:n){
      y[i] ~ bernoulli(inv_logit(beta0+beta1*x[i]));
    }
    beta0 ~ normal(0,10);
    beta1 ~ normal(0,10);
  }
'

mod1dat <- list(n=nrow(crab),y=crab$pa,x=crab$width)
mod1 <- stan(model_code=mod1code,data=mod1dat)
print(mod1,pars = c('beta0','beta1'),digits=3)
```

```{r}
par(mfrow=c(1,2))
plot(extract(mod1)$beta0,extract(mod1)$beta1,xlab = 'Beta 0', ylab = 'Beta 1',main='Intercept vs. Slope')

# Construct a grid of widths
widths <- seq(20,35,.1)

predth <- matrix(numeric(4000*length(widths)),nrow=4000)
prprob <- function(x,b0,b1){
  return(invlogit(b0+b1*x))  
}
for(i in 1:4000){
  # Obtain the predicted probabilities 
  # at each width based on b0[i]
  # and b1[i]
  predth[i,] <- apply(as.matrix(widths),1,prprob,b0=extract(mod1)$beta0[i],b1=extract(mod1)$beta1[i])
}

with(crab,plot(ylab='Probability of Satellite Presence',main='95% CB',cex=1,col='green',jitter(pa,amount=.05)~width,bty='l'))

meds <- apply(predth,2,median)
low <- apply(predth,2,quantile,.025)
high <- apply(predth,2,quantile,.975)

lines(widths,meds,col='blue')
lines(widths,low,lty=2, col ='red')
lines(widths,high,lty=2,col ='red')
```

##2b.

```{r}
# load("C:/Users/Sapana/Desktop/R/crab.RData")
crab$pa <-ifelse(crab$satell>1,1,0)

cloglog <- function(x) { log(-log(1-x))}
logit <- function(x) { log(x/(1-x))}

library(rstan)

mod2code <- '
  data{
    int<lower=0> n;             // sample size
    int<lower=0,upper=1> y[n];  // presence or absence of satellite
    real x[n];                  // width
  }
  parameters{
    real beta0;
    real beta1;
  }
  model{
    for(i in 1:n){
      y[i] ~ bernoulli(1-exp(-exp(beta0+beta1*x[i])));
    }
    beta0 ~ normal(0,10);
    beta1 ~ normal(0,10);
  }
'

mod2dat <- list(n=nrow(crab),y=crab$pa,x=crab$width)
mod2 <- stan(model_code=mod2code,data=mod2dat)
print(mod2,pars = c('beta0','beta1'),digits=3)


```

```{r}
par(mfrow=c(1,2))
plot(extract(mod2)$beta0,extract(mod2)$beta1,xlab = 'Beta 0', ylab = 'Beta 1',main='Intercept vs. Slope')

# Construct a grid of widths
widths <- seq(20,35,.1)

predth <- matrix(numeric(4000*length(widths)),nrow=4000)
prprob <- function(x,b0,b1){
  return(invlogit(b0+b1*x))  
}
for(i in 1:4000){
  # Obtain the predicted probabilities 
  # at each width based on b0[i]
  # and b1[i]
  predth[i,] <- apply(as.matrix(widths),1,prprob,b0=extract(mod2)$beta0[i],b1=extract(mod2)$beta1[i])
}

with(crab,plot(ylab='Probability of Satellite Presence',main='95% CB',cex=1,col='green',jitter(pa,amount=.05)~width,bty='l'))

meds <- apply(predth,2,median)
low <- apply(predth,2,quantile,.025)
high <- apply(predth,2,quantile,.975)

lines(widths,meds,col='blue')
lines(widths,low,lty=2, col ='red')
lines(widths,high,lty=2,col ='red')

```


##There are two parts..
##2c. Submitted on time with an error NaN
```{r}
theta1 <- rnorm(1000,0,0.1)
theta2 <- rnorm(1000,0,0.1)
b0c <- extract(mod2)$beta0
b1c <- extract(mod2)$beta1

y <- crab$pa
x <- crab$width

nsim <- 4000
r <- numeric(nsim)
logr <- numeric(nsim)

for(i in 1:nsim){
  beta0 <- b0c[i]
  beta1 <- b1c[i]
  
  theta <- invlogit(beta0+beta1*x)
  thetac <- 1-exp(-exp(beta0+beta1*x))
  
  f0 <- dbeta(invlogit(beta0+beta1*22),0,0.1)*
        dbeta(invlogit(beta0+beta1*28),0,0.1)
  logf0 <- log(f0)
  
  f1 <- dbeta(1-exp(-exp(beta0+beta1*22)),0,0.1)*
        dbeta(1-exp(-exp(beta0+beta1*28)),0,0.1)
  logf1 <- log(f1)
 
  logfyM0 <- sum(y*log(theta)+(1-y)*log(1-theta))
  logfyM1 <- sum(y*log(thetac)+(1-y)*log(1-thetac))
  
  logr[i] <- logfyM0+logf0-logfyM1-logf1
}

BFr <- mean(exp(logr))
BFr


```

##2c. This code was written after the deadline of 11:30am in the class..

```{r}
set.seed(889)
theta1 <- rnorm(1000,0,10)
theta2 <- rnorm(1000,0,10)
b0c <- extract(mod2)$beta0
b1c <- extract(mod2)$beta1

y <- crab$pa
x <- crab$width

nsim <- 4000
r <- numeric(nsim)
logr <- numeric(nsim)

for(i in 1:nsim){
  beta0 <- b0c[i]
  beta1 <- b1c[i]
  
  theta <- invlogit(beta0+beta1*x)
  thetac <- 1-exp(-exp(beta0+beta1*x))
  
  f0 <- dnorm(invlogit(beta0+beta1*22),0,10)*
        dnorm(invlogit(beta0+beta1*28),0,10)
  logf0 <- log(f0)
  
  f1 <- dnorm(1-exp(-exp(beta0+beta1*22)),0,10)*
        dnorm(1-exp(-exp(beta0+beta1*28)),0,10)
  logf1 <- log(f1)
 
  logfyM0 <- sum(y*log(theta)+(1-y)*log(1-theta))
  logfyM1 <- sum(y*log(thetac)+(1-y)*log(1-thetac))
  
  logr[i] <- logfyM0+logf0-logfyM1-logf1
}

BFr <- mean(exp(logr))
BFr


```
I get a Bayes Factor of .176. Since it is less than 1. This favors CLL model. 





