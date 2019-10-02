# Meghan Fager
# UPDATE JAGS SCRIPT WITH REVISED MCMC ESTIMATION SPECS

# Define directories
odir <- "../../output"
modelFile <- "B2PL-Interactions.bug"

# Data Estimation Conditions
nBurn <- 3500
nSamples <- 7000
nThin <- max(1, floor((nSamples - nBurn) / 1000))

### Data-generation conditions
m <- 2 # Number of response options
reps <- 100
set.seed(1738)
P <- c(0, 0.25, 0.5) # Proportion of items with Interactions
estC <- c("None", "Match", "All")
nPersons <- c(500, 3000) # Number of Respondents

# Load script building function
source("build_cmd_file.R")

# Data-generating Conditions, per model type
conditions <- data.frame("P" = sort(rep(P, 2)), "N" = rep(sort(rep(nPersons, 1)), 1))

g <- 1
for(g in 1:nrow(conditions)){  # Start loop
  
  ### Model conditions
  p <- conditions[g,"P"]
  n <- conditions[g,"N"]
  
  ### Data-gen Output directory
  mod <- ifelse(m==2, "C2PL", "GRM")
  gen <- paste0(mod, "-", p, "-", n)
  odir2 <- file.path(odir, paste0("N", n), gen)
  
  r <- 1
  for(r in 1:reps){  # Start loop for replications
    
    ### Replication Output directory
    odir3 <- file.path(odir2, r)
    
    initsFile1 <- file.path(odir3, "inits1.R")
    initsFile2 <- file.path(odir3, "inits2.R")
    
    #### ESTIMATION - DATA CONDITIONS ==========================
    cond = 1
    for(cond in 1:length(estC)){
      est <- estC[cond]
      
      if((p==0 & est=="Match")==TRUE){ next } else { # skip iteration and go to next iteration
        ### Estimation Output directory
        odir4 <- file.path(odir3, est)
        
        script <- build_cmd_file(est = est, r = r, nThin = nThin, nBurn = nBurn, nSamples = nSamples, initsFile1 = initsFile1,
                                 initsFile2 = initsFile2, modelFile = modelFile)
        write.table(script, file = file.path(odir4, "script.jags"), quote = FALSE, row.names = FALSE, col.names = FALSE)
      }
      
    } # Estimation
  } # Replications
} # Conditions
