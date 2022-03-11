# ________________________________________________
# MCMC sampling Theta and Tau normalized posterior. 
# ________________________________________________


# A function to calculate the log unnormalized posterior given Theta, Tau & ti:
Logposterior <- function(theta,tau,mu.theta,mu.tau,t,x,n){
  value <- -1/mu.tau*tau-1/mu.theta*theta 
  w <- 0
  for(i in 1:length(t)){
    val <- log(2/theta)-2/theta*t[i]+x[i]*log(3/4-3/4*exp((-8/3)*(tau+t[i])))+(n[i]-x[i])*log(1/4+3/4*exp((-8/3)*(tau+t[i])))
    w <- w+val
  }
  log.p <- value + w
  return(log.p)
}


# A function to calculate ln alpha for theta:
ln.alpha.theta <- function(theta,theta.new,mu.theta,mu.tau,t,L){
  alpha <- (-1/mu.theta*theta.new+(log(2/theta.new)*L-sum(2/theta.new*t)))-(-1/mu.theta*theta+(log(2/theta)*L-sum(2/theta*t)))
  return(alpha)
}

# A function to calculate the log acceptance probability for each coalescent times:
Logposterior.ti <- function(theta,tau,mu.theta,mu.tau,t,x,n){
  val <- -2/theta*t+x*log(3/4-3/4*exp((-8/3)*(tau+t)))+(n-x)*log(1/4+3/4*exp((-8/3)*(tau+t)))
  log.p  <- val
  return(log.p)
}

mcmc.msc = function(N=1000, theta.init=0.001, tau.init=0.01, t.init=0.001, L=1000, w_theta=0.00045, w_tau=0.0004, w_t=0.0085,x=xi,n=ni,mu.tau=0.005,mu.theta=0.001){
  start_time <- Sys.time() # Time at the start of the program.
  sample.tau <- numeric(N+1) # Vectors to hold visited Tau an Theta for N steps.
  sample.theta <- numeric(N+1) 
  sample.tau[1] <- tau.init # Starting positions are given. 
  sample.theta[1] <- theta.init
  accept.theta <- 0  # Counts the number of acceptances
  accept.tau <- 0 
  accept.t.rate <- c() # Holds acceptance rates of ti for each iteration. 
  tau <- tau.init # Parameters are given initial values.
  theta <- theta.init
  t <- rep(t.init,L)
  lnp.t <- c() # Empty vector to hold lnp at each loci.
  for(i in 1:L){
    lnp.t[i] <- Logposterior.ti(theta,tau,mu.theta,mu.tau,t=t[i],x=x[i],n=n[i]) # compute lnp at each loci. 
  }
  
  for(k in 1:N){ # Starting the chain with N iterations.
    lnp <-  Logposterior(theta,tau,mu.theta,mu.tau,t,x,n) # Keep the lnp up to date.
    # Changes tau, using a sliding window with size w_tau
    tau.new <- tau + w_tau*(runif(1)-0.5)
    if(tau.new<0){
      tau.new <- -tau.new 
    } 
    lnew <- Logposterior(theta,tau.new,mu.theta,mu.tau,t,x,n) #lnp at the proposed position.
    lnalpha <- lnew - lnp # log ratio
    if(lnalpha>0||exp(lnalpha)>runif(1)){ # accept depending on log ratio
      # if accepted replace tau value & add to the accepted count.
      tau <- tau.new 
      accept.tau = accept.tau+1 
    }
    
    # Change theta, using sliding window with size w_theta
    theta.new <- theta + w_theta*(runif(1)-0.5)
    if(theta.new<0){
      theta.new <- -theta.new 
    } 
    lnalpha <- ln.alpha.theta(theta,theta.new,mu.theta,mu.tau,t,L) # log ratio
    if(lnalpha>0||exp(lnalpha)>runif(1)){ # accept depending on log ratio
      # if accepted replace theta value & add to the accepted count.
      theta <- theta.new 
      accept.theta = accept.theta+1
    }
    
    accept.t <- 0 # count acceptance of t moves.
    # Change the coalescent times, using sliding windows with size w_t
    for(i in 1:L){
      tnew = t[i] + w_t*(0.5-runif(1))
      if(tnew<0){
        tnew <- -tnew 
      } 
      lnew <- Logposterior.ti(theta,tau,mu.theta,mu.tau,tnew,x[i],n[i]) # lnp at the proposed position.
      lnalpha <- lnew - lnp.t[i] # log ratio
      if(lnalpha>0||exp(lnalpha)>runif(1)){ # accept depending on log ratio
        lnp.t[i] <- lnew # if accepted replace lnp & ti values & add to the accepted count.
        t[i] <- tnew
        accept.t = accept.t+1
      }
    }
    accept.t.rate <- c(accept.t.rate,accept.t/L) # Store the acceptance rate over L steps.
    
    sample.theta[k+1] <- theta # Store the parameter positions after each interation.
    sample.tau[k+1] <- tau
  }
  ti.accept.rate <- mean(accept.t.rate) # Average acceptance rate over all iterations.
  theta.accept.rate <- accept.theta/N 
  tau.accept.rate <- accept.tau/N 
  end_time <- Sys.time() # Time at the end of the program.
  run_time <- end_time - start_time # Time difference is the running time. 
  results<- list(sample.tau,sample.theta,tau.accept.rate,theta.accept.rate,ti.accept.rate,run_time) # Create a list of all the results.
  return(results)
}
