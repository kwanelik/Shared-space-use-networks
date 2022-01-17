# Functions required for Simulations.R

## Function for loading in data from a given location

loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}

#-------------------

## Function for calculating probability of detecting an individual at a distance, x, from its centroid
## Assuming a negative sigmoidal curve (defined by parameters a & b) describes the change in probability 
## of detection with increasing distance from the centroid (hereafter referred to as the 'home range')

prob_all <- function(x, a, b) {
  return(1/(1 + exp(-(a + b*x))))
}

#-------------------

## Function for generating network based on overlapping home ranges, from known individual centroids (xs and ys) and home range parameters (as & bs)
## Inputs: 	xs = vector of x location
##		ys = vector of y location 
## 		as = vector of a parameter values 
##		bs = vector of b parameter values
##    n = buffer area beyond home ranges
##    res = number of points along both x and y axes (higher = slower)
## Output: 	1. Matrix of overlapping home ranges 

get_network_2D <- function(xs, ys, as, bs, n=2, res=10) {
 
  edges <- matrix(NA,length(xs),length(xs))
  for (i in 1:(nrow(edges)-1)) {
	grid <- expand.grid(x=seq(min(xs)-n,max(xs)+n,length.out=res), ys=seq(min(ys)-n,max(ys)+n,length.out=res))
	d1 <- as.matrix(dist(cbind(c(xs[i],grid$x),c(ys[i],grid$y))))[1,-1]
	probs.1 <- prob_all(log(d1+1), as[i], bs[i])

    for (j in (i+1):nrow(edges)) {

	d2 <- as.matrix(dist(cbind(c(xs[j],grid$x),c(ys[j],grid$y))))[1,-1]
	probs.2 <- prob_all(log(d2+1), as[j], bs[j])

	top <- apply(cbind(probs.1,probs.2),1,min)
	bottom <- apply(cbind(probs.1,probs.2),1,max)

	edges[i,j] <- sum(top)/sum(bottom)
	edges[j,i] <- edges[i,j]
    }
  }
  return(edges)
}

#-------------------

## Function for simulating real network of N nodes with some attribute 
## Example attribute: Sex
## Inputs: 	N = number of nodes / individuals in the network
##		noise = noise drawn from a normal distribution to simulate individual-level variation in home range size & shape
## 		N.classes = Number of classes of individuals (e.g. 2 for "M" and "F")
##		HR.params = Home range parameters: a matrix with 2 columns (parameters a and b) and one row for each class
## Output: 	1. Simulated real network matrix, 
##         	2. Matrix including for each individual: a & b parameters, distances from centroid to each trap, and probability of being detected in each trap
##         	3. Data frame of real centroids [note that Classes are output as letters of the alphabet, starting with 'a' for the first row of HR.params, 'b' for the second row, etc.]
## Note: defaults are the parameters estimated from our own empirical dataset for Male (a = 2.0800; b = -4.8220) and Female (a = 2.8250; b = -6.2072) field voles (Microtus agrestis)
  # Variants: in the MS we simulate real networks with varying effect sizes. We vary the effect size by varying the b parameter for Males and Females.
  # For 1/2 effect size observed in our empirical dataset: Male (a = 2.0800; b = -5.0720) and Female (a = 2.8250; b = -5.9572)
  # For 2 x effect size observed in our empirical dataset: Male (a = 2.0800; b = -4.3220) and Female (a = 2.8250; b = -6.7072)

generate_real_network_data <- function(nodes, noise=0.05, N.classes=2, HR.params=matrix(c(2.0800,-4.8220,2.8250,-6.2072),ncol=2,nrow=2,byrow=TRUE)) {
  
  # Generating real centroids with Class attributes 
  x<-runif(nodes, 1, 10)
  y<-runif(nodes, 1, 10)
  classes<-letters[1:N.classes]
  real_centroids<-data.frame(Tag_ID=1:nodes, x=x, y=y, a=NA, b=NA, Class=sample(classes, nodes, replace=T), stringsAsFactors=FALSE)
  
  # Generate home range parameters for each individual
  for (i in 1:length(classes)) {
    real_centroids$a[real_centroids$Class==classes[i]] <- rnorm(sum(real_centroids$Class==classes[i]),mean=HR.params[i,1],sd=abs(HR.params[i,1])*noise)  
    real_centroids$b[real_centroids$Class==classes[i]] <- rnorm(sum(real_centroids$Class==classes[i]),mean=HR.params[i,2],sd=abs(HR.params[i,2])*noise)  
  }  
  
  # Generating real network using above information
  real_network <- get_network_2D(real_centroids$x, real_centroids$y, real_centroids$a, real_centroids$b,res=10)
  
  # Renaming rows/cols as unique individual (Tag) IDs 
  rownames(real_network)<-real_centroids$Tag_ID
  colnames(real_network)<-real_centroids$Tag_ID
  
  # Creating a matrix of distances between each individual's real centroid and each trap 
  matrix_dists_real <- expand.grid(Tag_ID=real_centroids$Tag_ID,Location=as.character(traps$Location))
  matrix_dists_real$Class <- real_centroids$Class[match(matrix_dists_real$Tag_ID,real_centroids$Tag_ID)]
  matrix_dists_real$x.ind <- real_centroids$x[match(matrix_dists_real$Tag_ID,real_centroids$Tag_ID)]
  matrix_dists_real$y.ind <- real_centroids$y[match(matrix_dists_real$Tag_ID,real_centroids$Tag_ID)]
  matrix_dists_real$x.trap <- traps$x[match(matrix_dists_real$Location,traps$Location)]
  matrix_dists_real$y.trap <- traps$y[match(matrix_dists_real$Location,traps$Location)]
  matrix_dists_real$a.ind <- real_centroids$a[match(matrix_dists_real$Tag_ID,real_centroids$Tag_ID)]
  matrix_dists_real$b.ind <- real_centroids$b[match(matrix_dists_real$Tag_ID,real_centroids$Tag_ID)]
  matrix_dists_real$Dist <- sqrt((matrix_dists_real$x.ind-matrix_dists_real$x.trap)^2 + (matrix_dists_real$y.ind - matrix_dists_real$y.trap)^2)
  
  # Using home range profiles to calculate probabilities of each individual being detected in each trap
  matrix_dists_real$Dist.log<-log(matrix_dists_real$Dist + 1)
  matrix_dists_real$Trap.prob <- NA
  matrix_dists_real$Trap.prob <- prob_all(matrix_dists_real$Dist.log, matrix_dists_real$a.ind, matrix_dists_real$b.ind)
  # NB: Home range profiles and resulting trapping probabilities were calculated given an individual was trapped
  # Reducing probabilities to account for the fact that an individual will not be trapped every time you put a trap out
  matrix_dists_real$Trap.prob <- matrix_dists_real$Trap.prob/(matrix_dists_real$Dist*500)
  matrix_dists_real$Trap.prob[matrix_dists_real$Trap.prob > 0.95] <- 0.95

  return(list(real_network=real_network, matrix_dists_real=matrix_dists_real, real_centroids=real_centroids))
  
}

#-------------------

## Function for simulating (observed) shared trap network, from unique individual (Tag) IDs and locations at which each ID was observed

get_cotrap_network <- function(Tag_ID, Location) {

  a <- t(table(Tag_ID,Location))
	Tags <- colnames(a)
	network <- matrix(NA, nrow=length(Tags), ncol=length(Tags))
	rownames(network) <- Tags
	colnames(network) <- Tags
	for (i in 1:length(Tags)) {
		b <- a+a[,i]
		network[i,] <- colSums(b==2)/colSums(b>0)
	}
	return(network)
}

#-------------------

## Function to compare an observed network to the real (generative) network
## Inputs: 	real_network = original simulated network
##		new_network = network from the simulated dataset 
## 		n.rans = number of randomisations 
##		class = vector of classes for each individual
## Output: 1. Correlation: Mantel correlation between real and observed network, 
##         2. Power: Testing whether significant effect (present in real network) also detected in observed network, by node permutation,
##         3. Accuracy: Sum of absolute differences between real and observed networks.
## Note: This test compares class 2 ("M") to class 1 ("F"). This function should be changed if more classes are included.

get_results_from_sim <- function(real_network, new_network, n.rans=100, class) {
  
	require(vegan)
	require(sna)

  # unique classes
  classes <- c("M", "F") 
  
  # Correlation: Mantel correlation between real and observed network 
  order <- rownames(new_network)
  if (length(order) > 4) {
      m_test<-mantel(real_network[order,order],new_network,na.rm=TRUE, permutations=0)
      corr<-m_test$statistic
  } else {
      corr <- NA
  }
        
  # Power: Testing whether significant effect detected, by node permutation
  # Example effect: Mean degree (sum of edge weights) greater for Males than Females
  deg <- rowSums(new_network,na.rm=T)
  ids <- rownames(new_network)
  obs <- mean(deg[which(class[as.numeric(ids)] == classes[1])]) -  
  mean(deg[which(class[as.numeric(ids)] == classes[2])])
  ran <- rep(NA, n.rans)
  for (k in 1:n.rans) {
      n.r <- rmperm(new_network)
      deg.r <- rowSums(n.r,na.rm=T)
      ran[k] <- mean(deg.r[which(class[as.numeric(ids)] == classes[1])]) -
      mean(deg.r[which(class[as.numeric(ids)] == classes[2])])
  }
  true.pos <- sum(ran >= obs) < 2.5*n.rans/100  # Two-tailed test 
        
  # Accuracy: Calculating sum of absolute differences between real and observed network 
  # Re-including those individuals who were not observed in observed network -> assigning them a zero
  order <- rownames(new_network)
	full.obs.network <- matrix(0,nrow=nrow(real_network),ncol=ncol(real_network))
	diag(full.obs.network) <- NA
	rownames(full.obs.network) <- rownames(real_network)
	colnames(full.obs.network) <- colnames(real_network)
	full.obs.network[order,order] <- new_network
  # Taking the absolute difference
  acc<-sum(abs(real_network-full.obs.network),na.rm=T)

	return(list(corr=corr, true.pos=true.pos, acc=acc))

}

