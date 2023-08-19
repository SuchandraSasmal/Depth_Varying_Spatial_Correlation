#Task 1 : Phi = 0.05


library(mvtnorm)

ms1 <- numeric()
ms2 <- numeric()
ms3 <- numeric()
ms4 <- numeric()
ms5 <- numeric()
phi_opt <- numeric()

for(k in 1:100)
{
  loc <- cbind(runif(100),runif(100))  #simulating 100 locations
  true_phi <- rep(0.05,100)  #true value of phi
  true_sig <- 1     #true value of variance 
  cov1 <- matrix(0,100,100)     #original covariance kernel
  
  for(i in 1:100)
  {
    for(j in 1:100)
    {
      lo <- norm(loc[i,]-loc[j,],"2")
      cov1[i,j] <- exp(-lo/true_phi[i])
    }
  }
  
  X <- cbind(1,loc)   #design matrix
  
  y <- matrix(rmvnorm(n=1,rep(0,100),cov1),100)   #simulating response
  
  
  
  obj <- function(par)   #function giving the covariance matrix for a given phi
  {
    phi <- rep(par,100)
    cov <- matrix(0,100,100)   #covariance matrix for a particular phi
    
    for(i in 1:100)
    {
      for(j in 1:100)
      {
        lo <- norm(loc[i,]-loc[j,],"2")
        cov[i,j] <- exp(-lo/phi[i])
      }
    }
    
    
    #beta estimate for a particular phi
    beta <- solve(t(X)%*%solve(cov)%*%X)%*%t(X)%*%solve(cov)%*%y 
    
    #variance estimate for a particular phi
    sig <- t(y-(X%*%beta))%*%solve(cov)%*%(y-(X%*%beta))/100
    
    #likelihood function for multivariate normal with mean X*beta and dispersion=sig*cov
    cov2 <- as.vector(sig) * cov
    
    like <- exp(-(t(y-(X%*%beta))%*%solve(cov2)%*%(y-(X%*%beta)))/2)/((2*pi)^50 %*% sqrt(det(cov2)))
    like
  }
  
  opt <- optimize(f=obj,interval = c(0.002,0.5),maximum=TRUE)
  
  #optimum value of phi
  phi_opt[k] <- opt$maximum
  
  cov3 <- matrix(0,100,100)       #optimum value of covariance kernel
  
  for(i in 1:100)
  {
    for(j in 1:100)
    {
      lo <- norm(loc[i,]-loc[j,],"2")
      cov3[i,j] <- exp(-lo/phi_opt[k])
    }
  }
  
  #optimum value of beta
  beta_opt <- solve(t(X)%*%solve(cov3)%*%X)%*%t(X)%*%solve(cov3)%*%y
  
  #optimum value of variance
  sig_opt <- t(y-(X%*%beta_opt))%*%solve(cov3)%*%(y-(X%*%beta_opt))/100
  
  
  #value of beta with true phi = 0.05
  beta_hat <- solve(t(X)%*%solve(cov1)%*%X)%*%t(X)%*%solve(cov1)%*%y
  
  #value of variance with true phi = 0.05
  sig_hat <- t(y-(X%*%beta_hat))%*%solve(cov1)%*%(y-(X%*%beta_hat))/100
  
  ms1[k]<- beta_hat[1]-beta_opt[1]
  ms2[k]<- beta_hat[2]-beta_opt[2]
  ms3[k]<- beta_hat[3]-beta_opt[3]
  ms4[k] <- sig_hat - sig_opt
  ms5[k] <- phi_opt[k] - true_phi[1]
}

mse1 <- mean(ms1^2)  #mse of beta_0

mse2 <- mean(ms2^2)  #mse of beta_1

mse3 <- mean(ms3^2)  #mse of beta_2

mse4 <- mean(ms4^2)  #mse of variance

mse5 <- mean(ms5^2)  #mse of phi

mad1 <- mean(abs(ms1))  #mad of beta_0

mad2 <- mean(abs(ms2))  #mad of beta_1

mad3 <- mean(abs(ms3))  #mad of beta_2

mad4 <- mean(abs(ms4))  #mad of variance

mad5 <- mean(abs(ms5))  #mad of phi




#Task2: phi = 0.05 on left and phi = 0.1 on right

ms11 <- numeric()
ms21 <- numeric()
ms31 <- numeric()
ms41 <- numeric()
ms51 <- numeric()
ms61 <- numeric()
phi_opt1 <- numeric()


for(k in 1:100)
{
  loc <- cbind(runif(100),runif(100))  #simulating 100 locations
  true_phi1 <- numeric() #true value of phi
  lon <- loc[,1]
  lat <- loc[,2]
  for(i in 1:100)
  {
    if(lon[i] < 0.5)
    { 
      true_phi1[i] <- 0.05
    }
    else
    {
      true_phi1[i] <- 0.1
    }
  }
  
  true_sig <- rep(1,100)     #true value of variance 
  cov1 <- matrix(0,100,100)     #original covariance kernel
  
  for(i in 1:100)
  {
    for(j in 1:100)
    {
      s1 <- loc[i,]
      s2 <- loc[j,]
      dis <- t(s1-s2)%*%diag(1 / c((true_phi1[i]^2+true_phi1[j]^2)/2,(true_phi1[i]^2+true_phi1[j]^2)/2))%*%(s1-s2)
      rns <- (true_phi1[i]*true_phi1[j]*2*exp(-sqrt(dis)))/(true_phi1[i]^2+true_phi1[j]^2)
      cov1[i,j] <- true_sig[i]*true_sig[j]*rns
    }
  }
  
  X <- cbind(1,loc)   #design matrix
  
  y <- matrix(rmvnorm(n=1,rep(0,100),cov1),100)   #simulating response
  
  
  
  obj <- function(par)   #function giving the covariance matrix for a given phi
  {
    phi <- rep(par,100)
    cov <- matrix(0,100,100)   #covariance matrix for a particular phi
    
    for(i in 1:100)
    {
      for(j in 1:100)
      {
        lo <- norm(loc[i,]-loc[j,],"2")
        cov[i,j] <- exp(-lo/phi[i])
      }
    }
    
    
    #beta estimate for a particular phi
    beta <- solve(t(X)%*%solve(cov)%*%X)%*%t(X)%*%solve(cov)%*%y 
    
    #variance estimate for a particular phi
    sig <- t(y-(X%*%beta))%*%solve(cov)%*%(y-(X%*%beta))/100
    
    #likelihood function for multivariate normal with mean X*beta and dispersion=sig*cov
    cov2 <- as.vector(sig) * cov
    
    like <- exp(-(t(y-(X%*%beta))%*%solve(cov2)%*%(y-(X%*%beta)))/2)/((2*pi)^50 %*% sqrt(det(cov2)))
    like
  }
  
  opt <- optimize(f=obj,interval = c(0.002,0.5),maximum=TRUE)
  
  #optimum value of phi
  phi_opt1[k] <- opt$maximum
  
  cov3 <- matrix(0,100,100)
  
  for(i in 1:100)
  {
    for(j in 1:100)
    {
      lo <- norm(loc[i,]-loc[j,],"2")
      cov3[i,j] <- exp(-lo/phi_opt1[k])
    }
  }
  
  #optimum value of beta
  beta_opt <- solve(t(X)%*%solve(cov3)%*%X)%*%t(X)%*%solve(cov3)%*%y
  
  #optimum value of variance
  sig_opt <- t(y-(X%*%beta_opt))%*%solve(cov3)%*%(y-(X%*%beta_opt))/100
  
  
  #value of beta with true phi = 0.05
  beta_hat <- solve(t(X)%*%solve(cov1)%*%X)%*%t(X)%*%solve(cov1)%*%y
  
  #value of variance with true phi = 0.05
  sig_hat <- t(y-(X%*%beta_hat))%*%solve(cov1)%*%(y-(X%*%beta_hat))/100
  
  ms11[k]<- beta_hat[1]-beta_opt[1]
  ms21[k]<- beta_hat[2]-beta_opt[2]
  ms31[k]<- beta_hat[3]-beta_opt[3]
  ms41[k] <- sig_hat - sig_opt
  ms51[k] <- sum((rep(phi_opt1[k],100) - true_phi)^2)
  ms61[k] <- sum(abs(rep(phi_opt1[k],100) - true_phi))
}


mse11 <- mean(ms1^2)  #mse of beta_0

mse21 <- mean(ms2^2)  #mse of beta_1

mse31 <- mean(ms3^2)  #mse of beta_2

mse41 <- mean(ms4^2)  #mse of variance

mse51 <- ms51/100  #mse of phi

mad11 <- mean(abs(ms1))  #mad of beta_0

mad21 <- mean(abs(ms2))  #mad of beta_1

mad31 <- mean(abs(ms3))  #mad of beta_2

mad41 <- mean(abs(ms4))  #mad of variance

mad51 <- ms61/100  #mad of phi




















