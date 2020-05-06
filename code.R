library("ape")
library("expm")

#NENTROPY
#returns the node entropies by calculating sum of the state entropies
#prob: matrix of state probabilities

nentropy <- function(prob) {

  k              <- ncol(prob)                       #number of states
  prob[prob>1/k] <- prob[prob>1/k]/(1-k) - 1/(1-k)   #state entropies
  tent           <- apply(prob,1,sum)                #node entropy
 
  #correct absolute 0/1
  tent[tent == 0] <- tent[tent == 0] + runif(1,0,1)/10000
  tent[tent == 1] <- tent[tent == 1] - runif(1,0,1)/10000

  return(tent)
}

#FUNCTION FOR BAYESIAN INFERENCES
#bayesian inferences on the node entropies 
#l0: rate parameter of the exponential prior distribution
#se: standard deviation of the proposal distribution 
#a:  alpha parameter (beta likelihood)
#b:  beta paramter (beta likelihood)
#x:  node entropies

lpalpha <- function(a,b,x,l0) {          #log posterior alpha
  N  <- length(x)
  lp <- N*(lgamma(a+b)-lgamma(a)) - a*(l0-sum(log(x)))
  return(lp)
}

lpbeta  <- function(a,b,x,l0) {          #log posterior beta
  N  <- length(x)
  lp <- N*(lgamma(a+b)-lgamma(b)) - b*(l0-sum(log(1-x)))
  return(lp)
}

mhalpha <- function(a,b,x,l0,se) {       #metropolis hastings alpha
  a0 <- a
  a1 <- exp(rnorm(1,log(a0),se))

  r  <- min(1, exp(lpalpha(a1,b,x,l0) - lpalpha(a0,b,x,l0) ) )

  while (is.na(r) == T) {
    a1 <- exp(rnorm(1,log(a0),se))
    r  <- min(1, exp(lpalpha(a1,b,x,l0) - lpalpha(a0,b,x,l0) ) )
  }

  if (runif(1) < r) {
    return(a1) 
  } else {
    return(a0)
  }
}

mhbeta  <- function(a,b,x,l0,se) {      #metropolis hastings beta
  b0 <- b
  b1 <- exp(rnorm(1,log(b0),se))

  r  <- min(1, exp(lpbeta(a,b1,x,l0) - lpbeta(a,b0,x,l0) ) )

  while (is.na(r) == T) {
    b1 <- exp(rnorm(1,log(b0),se))
    r  <- min(1, exp(lpbeta(a,b1,x,l0) - lpbeta(a,b0,x,l0) ) )
  }  

  if (runif(1) < r) {
    return(b1)
  } else {
    return(b0)
  }
}

#MCMC
#Markov chain monte carlo scheme using the conditional posteriors of alpha and beta
#alpha: initial value of alpha
#beta: initial values of beta
#x: node entropies
#sim: number of iterations
#thin: controles the number of saved iterations = sim/thin
#burn: number of iterates to burn

emcmc <- function(alpha,beta,x,l0,se,sim,thin,burn) {

  usim <- seq(burn,sim,thin)
  gibbs <- matrix(NA,ncol=2,nrow=length(usim))
  p <- 1

  for (i in 1:sim) {
    alpha <- mhalpha(alpha,beta,x,l0,se)
    beta  <- mhbeta(alpha,beta,x,l0,se)

    if (i == usim[p]) {
       gibbs[p,] <- c(alpha,beta)
       p <- p+1
    }
  }  
  return(gibbs)
}

#RATE MATRIX FOR TRAIT EVOLUTION. K=2 TO 5
ratematrix <- function(pi,rho){
  
  k <- length(pi)
  
  if (k==2){
    r <- c(pi[1]*0     ,pi[2]*rho[1],
           pi[1]*rho[1],pi[2]*0)
  }
  
  if (k==3){
    r <- c(pi[1]*0     ,pi[2]*rho[1],pi[3]*rho[2],
           pi[1]*rho[1],pi[2]*0     ,pi[3]*rho[3],
           pi[1]*rho[2],pi[2]*rho[3],pi[3]*0 )
  }
  
  if (k==4){
    r <- c(pi[1]*0     ,pi[2]*rho[1],pi[3]*rho[2],pi[4]*rho[3],
           pi[1]*rho[1],pi[2]*0     ,pi[3]*rho[4],pi[4]*rho[5],
           pi[1]*rho[2],pi[2]*rho[4],pi[3]*0     ,pi[4]*rho[6],
           pi[1]*rho[3],pi[2]*rho[5],pi[3]*rho[6],pi[4]*0 )
  }  
  
  if (k==5){
    r <- c(pi[1]*0     ,pi[2]*rho[1],pi[3]*rho[2],pi[4]*rho[3] ,pi[5]*rho[4],
           pi[1]*rho[1],pi[2]*0     ,pi[3]*rho[5],pi[4]*rho[6] ,pi[5]*rho[7],
           pi[1]*rho[2],pi[2]*rho[5],pi[3]*0     ,pi[4]*rho[8] ,pi[5]*rho[9],
           pi[1]*rho[3],pi[2]*rho[6],pi[3]*rho[8],pi[4]*0      ,pi[5]*rho[10],
           pi[1]*rho[4],pi[2]*rho[7],pi[3]*rho[9],pi[4]*rho[10],pi[5]*0)
  }
  
  R <- matrix(r,ncol=k,nrow=k) 
  diag(R) <- -rowSums(R)
  
  return(R)
}

#RTRAIT
#simulates the evolution of a trait in a given tree
# tree: metric-tree
# R: rate matrix
# nstates: number of states

rtrait <- function(tree,R,nstates) {
  
  nspecis <- length(tree$tip.label)
  
  #tree
  edge <- cbind(tree$edge,tree$edge.length)
  
  ancestral <- rep(NA,2*nspecies-1) 
  ancestral[nspecies+1] <- sample(1:nstates,1,prob=pi) 
  
  #rate change
  inode <- nspecies+1
  while (sum(is.na(ancestral)) > 0) {
    
    inode1 <-  edge[which(edge[,1]==inode)[1],2]
    inode2 <-  edge[which(edge[,1]==inode)[2],2]
    bl1 <- edge[which(edge[,1]==inode)[1],3]
    bl2 <- edge[which(edge[,1]==inode)[2],3]
    
    astate <- rep(0,nstates)
    astate[ancestral[inode]] <- 1 
    
    ancestral[inode1] <- sample(1:nstates,1,prob=astate%*%expm(R*bl1))
    ancestral[inode2] <- sample(1:nstates,1,prob=astate%*%expm(R*bl2))
    
    inode <- inode+1
  }
  return(ancestral[1:nspecies])
  
}

#DELTA
#calculate delta statistic
#trait: trait vector 
delta <- function(trait, tree,lambda0,se,sim,thin,burn) {
  
  ar <- ace(trait,tree,type="discret",method="ML",model="ARD")$lik.anc
  
  # deletes the complex part whenever it pops up
  if (class(ar[1,1]) == "complex"){
    ar <- Re(ar)
  }
  
  x  <- nentropy(ar)
  mc1    <- emcmc(rexp(1),rexp(1),x,lambda0,se,sim,thin,burn)
  mc2    <- emcmc(rexp(1),rexp(1),x,lambda0,se,sim,thin,burn)
  mchain <- rbind(mc1,mc2)
  deltaA <- mean(mchain[,2]/mchain[,1])
  
  return(deltaA)
}

