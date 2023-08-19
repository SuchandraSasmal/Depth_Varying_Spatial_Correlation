library(fda)

set.seed(100)

d1 <- expand.grid(seq(0.1, 1, 0.1), seq(0.1, 1, 0.1)) # let's keep over unit square
lon <- d1$Var1  #longitude
lat <- d1$Var2  #latitude
depth <- runif(100,0,200)   #depth

nbasis <- 8
b1 <- create.bspline.basis(c(0,200), nbasis = nbasis)
plot(b1)
B <- eval.basis(depth, b1)
rowSums(B) # should be 1's
alpha <- rnorm(nbasis, mean = -2, sd = 0.1)
phi <- exp(c(B %*% alpha))
phi # values should be around 0.1 * diameter (max dist between two locations) of domain

library(ggmap)
library(readxl)
library(ggplot2)
library(plot3D)
library(viridis)
library(fda)
library(mvtnorm)

dat <- read_xlsx("observed_Bangladesh.xlsx")
dat <- as.data.frame(dat)

bd <- make_bbox(lon=c(min(dat$LON),max(dat$LON)),lat=c(min(dat$LAT),max(dat$LAT)),f=0.1)
bd1 <- get_stamenmap(bbox=bd,zoom=7,maptype="terrain")
as <- ggmap(bd1)+geom_point(data=dat,aes(x=dat$LON,y=dat$LAT,col=log(dat$`As (ug/L)`)))+scale_color_viridis(name="Log-As",option="F")+ggtitle("Scatter PLot of Log of Arsenic Content")

bd <- make_bbox(lon=c(min(dat$LON),max(dat$LON)),lat=c(min(dat$LAT),max(dat$LAT)),f=0.1)
bd1 <- get_stamenmap(bbox=bd,zoom=7,maptype="terrain")
well <- ggmap(bd1)+geom_point(data=dat,aes(x=dat$LON,y=dat$LAT,col=dat$`WELL_DEPTH (m)`))+scale_color_viridis(name="Well_Depth",option="A")+ggtitle("Scatter PLot of Well Depth")

bd <- make_bbox(lon=c(min(dat$LON),max(dat$LON)),lat=c(min(dat$LAT),max(dat$LAT)),f=0.1)
bd1 <- get_stamenmap(bbox=bd,zoom=7,maptype="terrain")
elev <- ggmap(bd1)+geom_point(data=dat,aes(x=dat$LON,y=dat$LAT,col=dat$`ELEV (m)`))+scale_color_viridis(name="Elevation",option="F")+ggtitle("Scatter PLot of Elevation")

set.seed(100)

d1 <- expand.grid(seq(0.1, 1, 0.1), seq(0.1, 1, 0.1)) # let's keep over unit square
lon <- d1$Var1  #longitude
lat <- d1$Var2  #latitude
depth <- runif(100,0,200) #depth
height <- runif(100,0,200)

nbasis <- 8
b1 <- create.bspline.basis(c(0,200), nbasis = nbasis)
plot(b1)
B <- eval.basis(depth, b1)
rowSums(B) # should be 1's
alpha <- rnorm(nbasis, mean = -2, sd = 0.1)
phi <- exp(c(B %*% alpha))
phi # values should be around 0.1 * diameter (max dist between two locations) of domain
sig <- rep(1,nrow(d1))

#Covariance matrix

df <- data.frame(lon,lat,depth,phi,sig)
cov <- matrix(0,nrow(df),nrow(df))

for(i in 1:nrow(df))
{
  for(j in 1:nrow(df))
  {
    s1 <- c(df$lon[i],df$lat[i])
    s2 <- c(df$lon[j],df$lat[j])
    dis <- t(s1-s2)%*%diag(c((df$phi[i]+df$phi[j])/2,(df$phi[i]+df$phi[j])/2),2)%*%(s1-s2)
    rns <- (df$phi[i]*df$phi[j]*2*exp(-sqrt(dis)))/(df$phi[i]^2+df$phi[j]^2)
    cov[i,j] <- sig[i]*sig[j]*rns
  }
}


# mean function using B-Spline
basis <- 2
b11 <- create.bspline.basis(c(0,200), nbasis = 2, norder =2 )
plot(b11)
B1 <- eval.basis(height, b11)
alpha1 <- rnorm(2, mean = 0, sd = 1)
mu <- c(B1 %*% alpha1)
epsi <- mvrnorm(n=1,rep(0,100),cov)
y <- mvrnorm(n=1,mu,cov)
x <- cbind(rep(1,100),height)
beta <- solve(t(x)%*%x)%*%t(x)%*%y
rmse <- sqrt(sum((beta-alpha1)^2)/100)
rmse



#Task 1
library(mvtnorm)

#set.seed(1)
loc <- cbind(runif(100),runif(100))
true_phi <- 0.05
true_sig <- 1
cov1 <- matrix(0,100,100)

for(i in 1:100)
{
  for(j in 1:100)
  {
    lo <- norm(loc[i,]-loc[j,],"2")
    cov1[i,j] <- exp(-lo/true_phi)
  }
}

X <- cbind(1,loc)

y <- matrix(rmvnorm(n=1,rep(0,100),cov1),100)



obj <- function(par)
{
  phi <- par
  cov <- matrix(0,100,100)
  
  for(i in 1:100)
  {
    for(j in 1:100)
    {
      lo <- norm(loc[i,]-loc[j,],"2")
      cov[i,j] <- exp(-lo/phi)
    }
  }
  
  
  
  beta_hat <- solve(t(X)%*%solve(cov)%*%X)%*%t(X)%*%solve(cov)%*%y
  sig_hat <- t(y-(X%*%beta_hat))%*%solve(cov)%*%(y-(X%*%beta_hat))/100
  cov2 <- as.vector(sig_hat) * cov
  
  like <- exp(-(t(y-(X%*%beta_hat))%*%solve(cov2)%*%(y-(X%*%beta_hat)))/2)/((2*pi)^50 %*% sqrt(det(cov2)))
  like
}

optimize(f=obj,interval = c(0.002,0.5),maximum=TRUE)
phi_opt <- opt$maximum
cov <- matrix(0,100,100)

for(i in 1:100)
{
  for(j in 1:100)
  {
    lo <- norm(loc[i,]-loc[j,],"2")
    cov[i,j] <- exp(-lo/phi_opt)
  }
}
beta_opt <- solve(t(X)%*%solve(cov)%*%X)%*%t(X)%*%solve(cov)%*%y
sig_opt <- t(y-(X%*%beta_hat))%*%solve(cov)%*%(y-(X%*%beta_hat))/100
