# Before starting:
#1. make sure you have installed the "ape" package.
#2. the file code.R is in the working directory.

library("ape")
source("code.R")

#SOME PARAMETERS... 
lambda0 <- 0.1   #rate parameter of the proposal 
se      <- 0.5   #standard deviation of the proposal
sim     <- 10000 #number of iterations
thin    <- 10    #we kept only each 10th iterate 
burn    <- 100   #100 iterates are burned-in

#RANDOM TREE: 
#same for both examples, only the trait vector varies.
set.seed(25)
ns   <- 20        #20 species
tree <- rtree(ns)

##########
# CASE A # : with phylogenetic signal
##########
trait <- c(2,1,3,1,1,3,1,3,2,1,1,2,2,2,2,1,1,3,1,1)

#CALCULATE DELTA A
deltaA <- delta(trait,tree,lambda0,se,sim,thin,burn)
print(deltaA)

#DRAW THE TREE...
par(mfrow=c(1,2))
tree$tip.label <- rep("",ns)
plot(tree,main="SCENARIO A")
ar <- ace(trait,tree,type="discret",method="ML",model="ARD")$lik.anc
nodelabels(pie = ar, cex = 1,frame="n") 
mtrait <- matrix(0,ncol=3,nrow=ns)
for ( i in 1:ns) {
  mtrait[i,trait[i]] <- 1
}
tiplabels(pie=mtrait,cex=0.5)



##########
# CASE B # : no phylogenetic signal
##########
trait <- c(2,3,1,3,3,3,3,2,2,3,1,2,1,2,3,1,2,3,1,2) 

#CALCULATE DELTA B
deltaB <-  delta(trait,tree,lambda0,se,sim,thin,burn)
print(deltaA)

#DRAW THE TREE...
ar <- ace(trait,tree,type="discret",method="ML",model="ARD")$lik.anc
plot(tree,main="SCENARIO B")
nodelabels(pie = Re(ar), cex = 1) 
mtrait <- matrix(0,ncol=3,nrow=ns)
for ( i in 1:ns) {  mtrait[i,trait[i]] <- 1 }
tiplabels(pie=mtrait,cex=0.5)


