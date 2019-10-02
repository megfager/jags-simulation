# Meghan Fager
# Generate and Estimate GRM bifactor models with interactions via JAGS (terminal)

onHPC = FALSE

### Load required packages
if (!require(truncnorm)) install.packages("truncnorm")
library(truncnorm)
if (!require(mvtnorm)) install.packages("mvtnorm")
library(mvtnorm)
if (!require(MASS)) install.packages("MASS")
library(MASS)
if (!require(MCMCpack)) install.packages("MCMCpack")
library(MCMCpack)
if (!require(mcmcplots)) install.packages("mcmcplots")
library(mcmcplots)
if (!require(gtools)) install.packages("gtools")
library(gtools)

# Define directories
odir <- "../output-GRM"
dir.create(odir)
modelFile <- "BGRM-Interactions.bug"
initsFile1 <- "../inits1.R" # Starting values for chain 1
initsFile2 <- "../inits2.R" # Starting values for chain 2

# Load script building function
source("../R/build_cmd_file.R")

### Data-generation conditions
reps <- 1
set.seed(1738)
K <- 6 # Number of Response Options
nBreaks <- K-1
P <- c(0, 0.25, 0.5) # Proportion of items with Interactions
estC <- c("None", "Match", "All")
nPersons <- c(500, 1000) # Number of Respondents
nItems <- 32 # Number of Items
# Number of grouping dimensions - nItems/nGroups, nItems/2, AND nItems/4 must be an integer
# Note: This does not guarantee that the number of interacting items per grouping dimension is the same
nGroups <- 4

# Create sample size directories
for(N in 1:length(nPersons)){
  dir.create(paste0(odir, "/N", nPersons[N]))
}

# MCMC Estimation Specifications
nBurn <- 5000
nSamples <- 10000
nThin <- max(1, floor((nSamples - nBurn) / 1000))

# Define seeds
chain1 <- seq(from = 289457, to = 289756, by = 3)
chain2 <- seq(from = 92834, to = 93133, by = 3)

# Data-generating Conditions, per model type
conditions <- data.frame("P" = sort(rep(P, 2)), "N" = rep(sort(rep(nPersons, 1)), 1))

### Create indicator by factor pattern matrix
perms <- unique(permutations(n = nGroups, r = nGroups, v = c(1, rep(0, nGroups-1)), set = FALSE, repeats.allowed = FALSE))
save <- NULL
for(d in 1:(nItems/nGroups)){
  save <- rbind(save, perms)
}
colnames(save) <- paste0("q", seq(1:nGroups))
qmat <- data.frame("qg" = rep(1, nItems), save)

### Theta covariance matrix
D <- nGroups + 1
mu <- rep(0, D)
SIG <- diag(1, nrow = D)

g <- 1
for(g in 1:nrow(conditions)){  # Start loop for data generation
  
  ### Model conditions
  p <- conditions[g,"P"]
  n <- conditions[g,"N"]
  
  ### Data-gen Output directory
  gen <- paste0("GRM", "-", p, "-", n)
  dir.create(file.path(odir, paste0("N", n)))
  odir2 <- file.path(odir, paste0("N", n), gen)
  dir.create(odir2)
  
  r <- 1
  for(r in 1:reps){  # Start loop for replications
    
    ### Replication Output directory
    odir3 <- file.path(odir2, r)
    dir.create(odir3)
    
    ### Generate Thetas
    true_thetas <- data.frame(mvrnorm(n, mu, SIG))
    
    ### Update Q-matrix to also indicate interactions
    first <- 1:(nItems*p)
    second <- 1+(nItems*p):(nItems-1)
    qmat$beta2 <- NULL
    qmat[first,"beta2"] <- 1
    qmat[second,"beta2"] <- 0
    
    ### Generate main effects and interactions
    amat <- data.frame("beta1" = matrix(runif(nItems*(D), 0, 2), ncol = D))    ## Main Effects
    amat <- cbind(amat, "beta2" = matrix(runif(nItems, -0.5, 0.5)))     ## Interactions
    amat <- amat * qmat
    
    
    ### Generate thresholds
    kmat <- data.frame(t(matrix(runif(nItems*D, -3, 3), ncol = 5)))
    for(k in 1:nItems){
      kmat[,k] <- sort(kmat[,k])
    }
    kmat <- t(kmat)
    colnames(kmat) <- paste0("beta0_", seq(1:nBreaks))
    
    # Item parameters
    item_params <- cbind(amat, kmat)
    
    ### Generate data
    dat <- matrix(NA, nrow = n, ncol = nItems)
    prob <- array(NA, c(n, nItems, K))
    
    for (item in 1:nItems){
      for (person in 1:n){ # Add Main effects and interactions to probit
        probit <- item_params[item, 1]*true_thetas[person, 1]  # General dimension main effect
        for(group in 2:(D)){  # Grouping dimensions main effects and interactions
          probit <- probit + item_params[item, group]*true_thetas[person, group]
          probit <- probit + item_params[item, "beta2"]*true_thetas[person, 1]*true_thetas[person, group]
        }
        probit1 <- item_params[item,"beta0_1"] + probit  # Lowest category
        probit2 <- item_params[item,"beta0_2"] + probit
        probit3 <- item_params[item,"beta0_3"] + probit
        probit4 <- item_params[item,"beta0_4"] + probit
        probit5 <- item_params[item,"beta0_5"] + probit
        probit6 <- 1
        prob[person,item,6] <- probit6 - pnorm(probit5) # Probability that person responds to category 6
        prob[person,item,5] <- pnorm(probit5) - pnorm(probit4)
        prob[person,item,4] <- pnorm(probit4) - pnorm(probit3)
        prob[person,item,3] <- pnorm(probit3) - pnorm(probit2)
        prob[person,item,2] <- pnorm(probit2) - pnorm(probit1)
        prob[person,item,1] <- pnorm(probit1) - 0
        
        
        dat[person,item] <- which(rmultinom(n = 1, size = 1,
                                            prob = prob[person,item,])==1)
      }
    }
    
    item_params <- cbind(amat, kmat)
    item_params$beta2 <- NULL
    betas <- qmat[,1:(D)] * amat$beta2
    names(betas) <- paste0("beta2.", seq(from = 1, to = D))
    betas$beta2.1 <- 0
    item_params <- cbind(item_params, betas)
    
    # Export files
    write.csv(dat, file = file.path(odir3, "data.csv"))
    write.csv(item_params, file = file.path(odir3, "item_params.csv"))
    write.csv(true_thetas, file = file.path(odir3, "thetas.csv"))
    
    #### End data-generation ===================================

    ### Initial Values ==========================================
    ### Write out replications seeds
    .RNG.seed <- chain1[r]
    set.seed(chain1[r])

    k_star <- as.matrix(cbind(rnorm(nItems, mean = -3), 
                              rnorm(nItems, mean = -1), 
                              rnorm(nItems, mean = 0), 
                              rnorm(nItems, mean = 1), 
                              rnorm(nItems, mean = 3)))
    k_star <- t(apply(t(k_star), 2, sort))
    .RNG.name <- "base::Super-Duper"
    dump(list = c(".RNG.seed", ".RNG.name", "k_star"), file = file.path(odir3, "inits1.R"), control = c("showAttributes"), append = FALSE)


    .RNG.seed <- chain2[r]
    set.seed(chain2[r])
    k_star <- as.matrix(cbind(rnorm(nItems, mean = -3), 
                              rnorm(nItems, mean = -1), 
                              rnorm(nItems, mean = 0), 
                              rnorm(nItems, mean = 1), 
                              rnorm(nItems, mean = 3)))
    k_star <- t(apply(t(k_star), 2, sort))
    .RNG.name <- "base::Wichmann-Hill"
    dump(list = c(".RNG.seed", ".RNG.name", "k_star"), file = file.path(odir3, "inits2.R"), control = c("showAttributes"), append = FALSE)
    
    #### ESTIMATION - DATA CONDITIONS ==========================
    cond = 1
    for(cond in 1:length(estC)){
      designMain <- qmat[,1:D] # Main effects
      est <- estC[cond]

      if((p==0 & est=="Match")==TRUE){ next } else { # skip iteration and go to next iteration
        ### Estimation Output directory
        odir4 <- file.path(odir3, est)
        dir.create(odir4)

        # JAGS conditions
        if(est == "All"){
          designInt <- designMain
          designInt[,1] <- 0
        }
        if(est == "Match"){
          designInt <- designMain*qmat[,"beta2"]
          designInt[,1] <- 0
        }
        if(est == "None"){
          designInt <- designMain*0
        }
        
        # Create data file that is readable for jags (via terminal)
        designMain <- unname(as.matrix(designMain))
        designInt <- unname(as.matrix(designInt))
        beta0 <- matrix(NA, nrow = nItems, ncol = nBreaks)
        dump(list = c("nItems", "n", "dat", "designMain", "designInt", "D", "mu", "SIG", "beta0", "K"),
             file = file.path(odir4, "data.R"), control = c("showAttributes"))

        # Create jags script
        script <- build_cmd_file(est = est, r = r, nThin = nThin, nBurn = nBurn, nSamples = nSamples, initsFile1 = initsFile1,
                                 initsFile2 = initsFile2, modelFile = modelFile)
        write.table(script, file = file.path(odir4, "script.jags"), quote = FALSE, row.names = FALSE, col.names = FALSE)

      }

    } # Estimation
  } # Replications
} # Conditions
