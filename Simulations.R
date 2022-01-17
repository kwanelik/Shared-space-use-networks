# Simulations

# What this code does:
# 1. Simulates a real network for N individuals, defined as the amount of overlap in the home ranges across all combinations of individuals,
# 2. Simulates observed datasets (given a range of sampling efforts, r) that contain detections of individuals at traps, where the detection probability 
#    for a given individual in a given trap is determined by the position of the trap relative to the home range of the individual, 
# 3. For each simulated observed dataset, generates an observed shared trap network, by calculating the number of traps in which two individuals 
#    were observed divided by the total number of traps that either individual was observed in,
# 4. And an observed overlap network, by estimating individual centroids and home ranges from the observed dataset, then calculating overlap in these 
#    estimated home ranges across all combinations of individuals (as in 1),
# 5. Repeats N times.

# Requires several functions in Functions.R
# Also requires a file with information about the location of traps (unique trap ID, x_coord and y-coord)
# Example: trap location file for 10x10 trapping grid (File: traps).
# Variant: in the MS we repeat this procedure using a denser 19x19 trapping grid (File: ext_traps).

# NB/ Nodes/individuals have attributes - here as an example, this is Sex (Male/"M" or Female/"F")

# Where example values are used these are specified after 'Example:' and refer to those values used to generate Figs 3-5 in the MS (but these can be varied)

#--------------------------------------------------------------------------
### Preparation
#--------------------------------------------------------------------------

# Clear environment
rm(list=ls())

# Reading in required functions in Functions.R
#source(<insert naming structure>)

# Reading in trap location file
#setwd(<insert naming structure>)
traps<-loadRData(file="traps")
traps$x_coord <- as.vector(traps$x_coord)
traps$Location <- as.character(traps$Location)

# Specifying number of simulated real networks/number of repeats of steps 1-4 above  
# Example: 5000
sim<-5000

# Specifying number of individuals required in each real network
# Example: 100
individuals <- 100

# Specifying range of sampling efforts required for simulation of observed datasets (number of draws from the probability of observing an individual in a trap)
# Example: 50-4500 (equivalent to mean number of captures per individual of approximately 1-30)
# For 10 x 10 grid
r_range<-c(50,100,150,200,300,400,500,800,1200,1800,2500,3500) 

#-------------------

# Creating objects to save into
# Check of whether significant effect is indeed present in real network 
true.pos<-matrix(data=NA, sim, 1)

# Mean number of captures per individual in observed dataset
capt<-matrix(data=NA, sim, length(r_range))
colnames(capt)<-r_range

# Number of individuals in observed networks 
no_indivs<-matrix(data=NA, sim, length(r_range))
colnames(no_indivs)<-r_range

# Mean number of captures per trap in observed dataset
xs<-matrix(data=NA, sim, length(r_range))
colnames(xs)<-r_range

# Correlation between real network and observed network (shared trap and overlap)
corr.shared<-corr.overlap.sep<-matrix(data=NA, sim, length(r_range))
colnames(corr.shared)<-colnames(corr.overlap.sep)<-r_range

# Accuracy of observed network (shared trap and overlap)
acc.shared<-acc.overlap.sep<-matrix(data=NA, sim, length(r_range))
colnames(acc.shared)<-colnames(acc.overlap.sep)<-r_range

# Power of observed network (shared trap and overlap)
true.pos.shared<-true.pos.overlap.sep<-matrix(data=NA, sim, length(r_range))
colnames(true.pos.shared)<-colnames(true.pos.overlap.sep)<-r_range
true.pos <- rep(NA, sim)
false.pos.shared<-false.pos.overlap.sep<-matrix(data=NA, sim, length(r_range))
colnames(false.pos.shared)<-colnames(false.pos.overlap.sep)<-r_range
false.pos <- rep(NA, sim)

#--------------------------------------------------------------------------
### Running simulations
#--------------------------------------------------------------------------

for (j in 1:sim){

#----------------------------------------------- 1. Simulating real network ---------------------------------------------#

 data <- generate_real_network_data(individuals,noise=0.05,N.classes=2,HR.params=matrix(c(2.0800,-4.8220,2.8250,-6.2072),ncol=2,nrow=2,byrow=TRUE))
 real_network <- data$real_network
 matrix_dists_real <- data$matrix_dists_real
 real_centroids <- data$real_centroids

 # Update classes to sexes (in our case)
 real_centroids$Sex <- real_centroids$Class
 real_centroids$Sex[which(real_centroids$Sex == "a")] <- "M"
 real_centroids$Sex[which(real_centroids$Sex == "b")] <- "F"
 matrix_dists_real$Sex <- matrix_dists_real$Class
 matrix_dists_real$Sex[which(matrix_dists_real$Sex == "a")] <- "M"
 matrix_dists_real$Sex[which(matrix_dists_real$Sex == "b")] <- "F"

 # Double checking that significant effect present in real network 
 # Example effect: Mean degree (sum of edge weights) greater for Males than Females
 results <- get_results_from_sim(real_network=real_network, new_network=real_network, n.rans=100, class=real_centroids$Sex)
 true.pos[j] <- results$true.pos
 
 #-------------------------------------------- 2. Simulating observed dataset ---------------------------------------------#

 for (i in 1:length(r_range)){ 
  
  r<-r_range[i]
        
  # Adding columns later on, so removing before starting next round of loop
  matrix_dists_real<-matrix_dists_real[,1:13]
        
  # Repeat until all attributes represented in dataset, then move on
  # Example attribute: Sex 
	repeat {
   catches <- rbinom(n=nrow(matrix_dists_real),size=r,prob=(1-(1-matrix_dists_real$Trap.prob)))
   matrix_dists_real$Det.obs <- as.numeric(catches>0) # Whether or not individual was caught in a trap
   matrix_dists_real$Det.count <- catches # Number of times an individual was caught in a trap 
	 if (length(unique(matrix_dists_real$Sex[matrix_dists_real$Det.obs==1])) > 1) { break; }  
	       }
        
  # Calculating mean number of captures per individual & saving
  capt[j,i]<-mean(tapply(matrix_dists_real$Det.count[matrix_dists_real$Det.count>0], list(matrix_dists_real$Tag_ID[matrix_dists_real$Det.count>0]),sum)) 
        
  # Calculating mean number of captures per trap & saving
  xs[j,i]<-mean(tapply(matrix_dists_real$Det.count, list(matrix_dists_real$Location),sum))
        
  #----------------------------------------------- 2. Generating shared trap network ---------------------------------------------#
      
  shared_network <- get_cotrap_network(matrix_dists_real$Tag_ID[which(matrix_dists_real$Det.obs==1)], matrix_dists_real$Location[which(matrix_dists_real$Det.obs==1)])
  diag(shared_network) <- NA 
        
  # Calculating: Correlation, power and accuracy of shared trap network
  results <- get_results_from_sim(real_network=real_network, new_network=shared_network, n.rans=100, class=real_centroids$Sex)
	true.pos.shared[j,i] <- results$true.pos
	corr.shared[j,i] <- results$corr
	acc.shared[j,i] <- results$acc
        
	#----------------------------------------------------- 3. Generating overlap network  ---------------------------------------------#
       
  # Recalculating (weighted) centroids 
	# Weighted because a trap in which an individual was caught multiple times will have a greater influence on its centroid than a trap in which it was caught only once
	matrix_dists_real2<-matrix_dists_real[,c(1,3,6,7,14,15)]
	matrix_dists_real2$Det.count[matrix_dists_real2$Det.count==0]<-1 
	matrix_dists_real2 <- as.data.frame(lapply(matrix_dists_real2, rep, matrix_dists_real2$Det.count)) 
	matrix_dists_real2<-matrix_dists_real2[,-6]
	x <- sapply(by(matrix_dists_real2$x.trap[matrix_dists_real2$Det.obs==1], matrix_dists_real2$Tag_ID[matrix_dists_real2$Det.obs==1],mean),identity)
  y <- sapply(by(matrix_dists_real2$y.trap[matrix_dists_real2$Det.obs==1], matrix_dists_real2$Tag_ID[matrix_dists_real2$Det.obs==1],mean),identity)
  obs_centroids <- data.frame(Tag_ID=real_centroids$Tag_ID,Sex=real_centroids$Sex)
  obs_centroids$x <- x[match(obs_centroids$Tag_ID,names(x))]
  obs_centroids$y <- y[match(obs_centroids$Tag_ID,names(y))]
        
  # Creating a matrix of distances between each individual's observed centroid and each trap
  matrix_dists_obs <- matrix_dists_real
  matrix_dists_obs$x.ind <- obs_centroids$x[match(matrix_dists_obs$Tag_ID,obs_centroids$Tag_ID)]
  matrix_dists_obs$y.ind <- obs_centroids$y[match(matrix_dists_obs$Tag_ID,obs_centroids$Tag_ID)]
  matrix_dists_obs$Dist <- sqrt((matrix_dists_obs$x.ind-matrix_dists_obs$x.trap)^2 + (matrix_dists_obs$y.ind - matrix_dists_obs$y.trap)^2)
  matrix_dists_obs$Dist.log<-log(matrix_dists_obs$Dist + 1)
  # Removing unobserved individuals
  matrix_dists_obs <- matrix_dists_obs[which(!is.na(matrix_dists_obs$Dist)),] 
        
  # Calculating total number of detections per individual
  a<-tapply(matrix_dists_obs$Det.count, list(matrix_dists_obs$Tag_ID),sum)
  matrix_dists_obs$Det.count.total <- a[match(matrix_dists_obs$Tag_ID,names(a))]
        
  # Estimating home range profiles (negative sigmoidal curves) using this data 
  # Example: Sex-specific curves (male home range profile = fit.males; female home range profile = fit.females)
  fit.males <- glm(Det.obs ~ Dist.log, data=matrix_dists_obs[which(matrix_dists_obs$Sex=="M"),], family=binomial,control = list(maxit = 50))
  fit.females <- glm(Det.obs ~ Dist.log, data=matrix_dists_obs[which(matrix_dists_obs$Sex=="F"),], family=binomial,control = list(maxit = 50))
  
  # Generating overlap network 
  params <- data.frame(sex=c("M","F"), a=c(coef(fit.males)[1],coef(fit.females)[1]), b=c(coef(fit.males)[2],coef(fit.females)[2]))
  as <- params$a[match(obs_centroids$Sex,params$sex)]
  bs <- params$b[match(obs_centroids$Sex,params$sex)]

  overlap_network_sep <- get_network_2D(obs_centroids$x[which(!is.na(obs_centroids$x))], obs_centroids$y[which(!is.na(obs_centroids$y))],as,bs)
  rownames(overlap_network_sep)=colnames(overlap_network_sep)=na.omit(obs_centroids)$Tag_ID
        
  # Calculating: Correlation, power and accuracy of overlap network
  results <- get_results_from_sim(real_network=real_network, new_network=overlap_network_sep, n.rans=100, class=real_centroids$Sex)
	true.pos.overlap.sep[j,i] <- results$true.pos
	corr.overlap.sep[j,i] <- results$corr
	acc.overlap.sep[j,i] <- results$acc
        
	#--------------------------------------------------------------------------------------------------------------------------------------#
        
  # Calculating number of individuals that are present in both observed networks (shared trap and overlap)
  no_indivs[j,i]<-nrow(shared_network)
        
 }}

# NB: Check true.pos: if FALSE (i.e. significant effect not present in real network) then turn TRUES to NA in true.pos.shared and true.pos.overlap.sep


