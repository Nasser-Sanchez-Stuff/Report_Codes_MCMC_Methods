#load packages and set ggplot theme
library(pracma)
library(tidyverse)
library(gridExtra)
ggplot2::theme_set(theme_bw())



# used to create our prior, likelihood and posterior for our model.
# We use the scaled values for plotting
get_bayesian_fns=function(theta,y,alpha,beta,n){
  prior_fn=function(theta,alpha,beta){
     return(theta^(alpha-1)*exp(-beta*theta))
  }
  prior=prior_fn(theta,alpha,beta)
  prior_scaled=prior/trapz(theta,prior)
  # 
  lh_fn=function(theta,n,y){
    return(theta^(n) * exp(-theta*sum(y)) )
  }
  likelihood=lh_fn(theta,n,y)
  likelihood_scaled=likelihood/trapz(theta,likelihood)
  
  posterior_fn=function(theta,y,alpha,beta,n){
    alpha_1=alpha+n
    beta_1=beta+sum(y)
    return(theta^(alpha_1-1)*exp(-beta_1*theta))
  }
  
  
  posterior=posterior_fn(theta,y,alpha,beta,n)
  posterior_scaled=posterior/trapz(theta,posterior)
  
  df=data.frame(x=theta,likelihood=likelihood,likelihood_scaled=likelihood_scaled,
                prior=prior,prior_scaled=prior_scaled,
                posterior=posterior,posterior_scaled=posterior_scaled)
  
  df=df[-c(2,4,6)]
  names(df)=gsub("_scaled","",names(df))
  df$obs=as.factor(n)
  return(df)
}

# initialise variables and pass into function defined above
y_2=c(1.5,3.5)
y_5=c(2.5,2.5,2.5,2.5,2.5)
# we want to create plots for values of theta between 0 and 3
theta=seq(0,3,length.out=1000)
alpha=2
beta=3

# create the datasets of our prior, likelihood and posterior for n=2 and n=5
# observations so we can compare
model_2_samples=get_bayesian_fns(theta,y_2,alpha,beta,2)
model_5_samples=get_bayesian_fns(theta,y_5,alpha,beta,5)
df=rbind(model_2_samples,model_5_samples)
 df_long=df%>%
   pivot_longer(likelihood:posterior,names_to = "type",values_to = "value")

 
# Likelihood function plot for n=2 observations 
ggplot(model_2_samples,aes(x,likelihood))+geom_line()+labs(x="θ",title="Likelihood Function Plot")+coord_cartesian(xlim=c(0,1.5))


# plot prior, likelihood and posterior for n=2 and n=5 observations
ggplot(data=df_long,aes(x,value))+geom_line(aes(color=type),linewidth=1)+
  theme_bw()+labs(
  x="θ",y="Density",title="Plot of Scaled Components of Bayesian Model
Separated by Number of Observations")+facet_wrap(~obs,ncol=2)#+coord_cartesian(ylim=c(0,2.5))

# define gibbs sampler 
gibbs_sampler = function(N, alpha, alpha_1, beta_1, y, init_beta,init_theta) {
  theta=rep(NA,N)
  beta=rep(NA,N)
  
  theta[1]=init_theta
  beta[1]=init_beta
  for (i in 2:N) {
    
    beta[i]=rgamma(1,alpha_1,theta[i-1]+beta_1)
    
    theta[i]=rgamma(1,length(y)+alpha,sum(y)+beta[i])
    
    
  }
  return(data.frame(t=1:N,theta = theta, beta = beta))
  
  
}
# define sampler with inferior performance
gibbs_sampler2 = function(N, alpha, alpha_1, beta_1, y, init_beta,init_theta) {
  theta=rep(NA,N)
  beta=rep(NA,N)
  
  theta[1]=init_theta
  beta[1]=init_beta
  for (i in 2:N) {
    
    beta[i]=rgamma(1,alpha_1,theta[i-1]+beta_1+rpois(1,3/sqrt(i)))
    
    theta[i]=rgamma(1,length(y)+alpha+rpois(1,3/sqrt(i)),sum(y)+beta[i])+
      rnorm(1,0,10/sqrt(i))+
      rpois(1,8/sqrt(i))+
      rnorm(1,0.99^i,5/sqrt(i))
    
    
  }
  return(data.frame(t=1:N,theta = theta, beta = beta))
}


# initialise values and call the samplers

alpha=2
alpha_1=40
beta_1=15
init_beta=1
init_theta=1
y=c(1.5,3.5)

set.seed(123)
samples_100=gibbs_sampler(100,alpha,alpha_1,beta_1,y,init_beta,init_theta)
set.seed(123)
samples_10000=gibbs_sampler(10000,alpha,alpha_1,beta_1,y,init_beta,init_theta)

set.seed(123)
samples_other = gibbs_sampler2(1000,alpha,alpha_1,beta_1,y,init_beta,init_theta)

samples_100$chain=as.factor("100")
samples_10000$chain=as.factor("10000")
samples_other$chain=as.factor("other")

samples=rbind(samples_100,samples_10000,samples_other)

#trace plots of MCMC chains
trace_model_A=ggplot(samples_100,aes(t,theta))+geom_line()+geom_smooth(method="loess")+
  labs(x="Time",y="θ",title="Model A")
trace_model_B=ggplot(samples_10000,aes(t,theta))+geom_line()+geom_smooth(method="loess")+
  labs(x="Time",y="θ",title="Model B") 
trace_model_C=ggplot(samples_other,aes(t,theta))+geom_line()+geom_smooth(method="loess")+
  labs(x="Time",y="θ",title="Model C")
grid.arrange(trace_model_A,trace_model_B,trace_model_C,top=
               "Trace Plots of All Models' Samples of θ Over Time",ncol=3)









# Drop first 10% from each chain as 'burn-in' and use the rest in calculations
# also split remaining data into first 20% and last 50% for comparing means of 
# subsequences
x1=samples_100[10:100,2]
x11=x1[1:18]
x21=x1[45:90]

x2=samples_10000[1001:10000,2]
x12=x2[1:1800]
x22=x2[4501:9000]


x3=samples_other[101:1000,2]
x13=x3[1:180]
x23=x3[451:900]

# Z scores for Geweke's convergence diagnostic
z_1=(mean(x11)-mean(x21))/sqrt((var(x11)/length(x11))+(var(x21)/length(x21)) )

z_2=(mean(x12)-mean(x22))/sqrt((var(x12)/length(x12))+(var(x22)/length(x22)) )

z_3=(mean(x13)-mean(x23))/sqrt((var(x13)/length(x13))+(var(x23)/length(x23)) )

# Plot ACF up to 50 lags for each series
acf1=acf(x1,lag.max=50)
acf1=data.frame(lag=acf1$lag,acf=acf1$acf)
acf2=acf(x2,lag.max=50)
acf2=data.frame(lag=acf2$lag,acf=acf2$acf)
acf3=acf(x3,lag.max=50)
acf3=data.frame(lag=acf3$lag,acf=acf3$acf)

p1=ggplot(acf1,aes(lag,acf))+geom_bar(stat="identity")+coord_cartesian(ylim=c(-0.2,.2))+labs(title="Model A")
p2=ggplot(acf2,aes(lag,acf))+geom_bar(stat="identity")+coord_cartesian(ylim=c(-0.2,.2))+labs(title="Model B")
p3=ggplot(acf3,aes(lag,acf))+geom_bar(stat="identity")+coord_cartesian(ylim=c(-0.2,.2))+labs(title="Model C")


gridExtra::grid.arrange(p1,p2,p3,ncol=3,top="Autocorrelation Functions of Model's A, B and C up to 50 Lags")


# Now calculate ACF series for each chain, using n/10 as number of lags to
# be calculated
acf1=acf(x1,lag.max=9)$acf
acf2=acf(x2,lag.max=900)$acf
acf3=acf(x3,lag.max=90)$acf

# Calculate ESS using the ACF values

ess1=length(x1)/(1+2*sum(acf1))
ess2=length(x2)/(1+2*sum(acf2))
ess3=length(x3)/(1+2*sum(acf3))


# Calculate standard error for each chain
mcmcse1=sqrt(var(x1)/(100*ess1/length(x1)))
mcmcse2=sqrt(var(x2)/(100*ess2/length(x2)))
mcmcse3=sqrt(var(x3)/(100*ess3/length(x3)))
