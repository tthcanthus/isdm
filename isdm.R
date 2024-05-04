library(spdep)
library(Matrix)
library(spatialreg)

##################################----EM algorithm for ISDM----########################################
estimate.ISDM <- function(size,mat,rho1,rho2,v)
# size is the sample size
# mat is the weight matrix  w
# rho1 is the spatial parameter rho_1 of the center model
# rho2 is the spatial parameter rho_2 of the range model
# when v=0, the error follows normal distribution; when v=2, the error follows t-distribution t(2)
{
    data.f = data.frame(y_c,x_c,x_r,mat%*%x_c,mat%*%x_r) 
    X <- cbind(x_c,x_r,mat%*%x_c,mat%*%x_r)
    formula <- y_c~.
    srk12w <- mat2listw(mat)
    
    begin.p <- lagsarlm(formula,data.f,srk12w)
    rho.h <- begin.p$rho
    beta.h <- begin.p$coefficients[-1]
    inter.h <- begin.p$coefficients[1]
    s2.h <- begin.p$s2
    nu.h <- 3
    
    rho.h0 <- begin.p$rho
    beta.h0 <- begin.p$coefficients[-1]
    
    
    if(v==0)
      nu.h <- 15
    abs <- 1
     while(abs>10^(-7)){
    ##############-------E step--------###################
          epsilon <- y_c-rho.h*mat%*%y_c-X%*%beta.h-inter.h
          delta <- (epsilon)^2/s2.h
          logu <- digamma(0.5*(nu.h+1))-log(0.5*(nu.h+delta))
          u <- as.vector((nu.h+1)/(nu.h+delta))
          x.u <- diag(sqrt(u))%*%X
          data.f <- data.frame(y_c,x.u)
          formula <- y_c~x.u
          rho.b <- rho.h
          beta.b <- beta.h
         s2.b <- s2.h
    #############--------M step--------##################
         log.gamma <- function(nu){
            f <- sum(0.5*nu*logu-0.5*nu*u)+size*(-log(gamma(0.5*nu))+0.5*nu*log(0.5*nu))
            f
          }
          v.loglike <- apply(v.n,2,log.gamma)
          v.index <- which(v.loglike==max(v.loglike))
          nu.h <- v.n[v.index]
          a <- lagsarlm(formula,data.f,srk12w)
          rho.h <- a$rho
          beta.h <- a$coefficients[-1]
          s2.h <- a$s2
          inter.h <- a$coefficients[1]
          abs <- max(c(abs(rho.h-rho.b),abs(beta.h-beta.b),abs(s2.h-s2.b)))
       }
       result1 <- c(rho.h,beta.h,nu.h)
       result10 <- c(rho.h0,beta.h0)
    
    data.fr = data.frame(y_r,x_c,x_r,mat%*%x_c,mat%*%x_r) 
    formula <- y_r~.
    X <- cbind(x_c,x_r,mat%*%x_c,mat%*%x_r)
    begin.pr <- lagsarlm(formula,data.fr,srk12w)
    rho.hr <- begin.pr$rho
    beta.hr <- begin.pr$coefficients[-1]
    inter.hr <- begin.pr$coefficients[1]
    s2.hr <- begin.pr$s2
    nu.h <- 3
    
    rho.hr0 <- begin.pr$rho
    beta.hr0 <- begin.pr$coefficients[-1]
    
    if(v==0)
     nu.h <- 15
    abs <- 1
    while(abs>10^(-7)){
    #      ##############-------E step--------###################
          epsilon <- y_r-rho.hr*mat%*%y_r-X%*%beta.hr-inter.hr
          delta <- (epsilon)^2/s2.hr
          logu <- digamma(0.5*(nu.h+1))-log(0.5*(nu.h+delta))
          u <- as.vector((nu.h+1)/(nu.h+delta))
          x.u <- diag(sqrt(u))%*%X
          data.f <- data.frame(y_r,x.u)
          formula <- y_r~x.u
          rho.b <- rho.hr
          beta.b <- beta.hr
          s2.b <- s2.hr
    #############--------M step--------##################
          log.gamma <- function(nu){
            f <- sum(0.5*nu*logu-0.5*nu*u)+size*(-log(gamma(0.5*nu))+0.5*nu*log(0.5*nu))
            f
          }
          v.loglike <- apply(v.n,2,log.gamma)
          v.index <- which(v.loglike==max(v.loglike))
          nu.h <- v.n[v.index]
          a <- lagsarlm(formula,data.f,srk12w)
          rho.hr <- a$rho
          beta.hr <- a$coefficients[-1]
          inter.hr <- a$coefficients[1]
          s2.hr <- a$s2
         abs <- max(c(abs(rho.hr-rho.b),abs(beta.hr-beta.b),abs(s2.hr-s2.b)))
        }
}


