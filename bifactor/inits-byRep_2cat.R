# Meghan Fager
### Write out init files in replication directories
## C2PL conditions

# Define directories
odir <- "../output"

## Conditions
P <- c(0, 0.25, 0.5) # Proportion of items with Interactions
estC <- c("None", "Match", "All")
nPersons <- c(500, 3000) # Number of Respondents
reps <- 1

# Data-generating Conditions, per model type
conditions <- data.frame("P" = sort(rep(P, 2)), "N" = rep(sort(rep(nPersons, 1)), 1))

# Define seeds
chain1 <- seq(from = 289457, to = 289756, by = 3)
chain2 <- seq(from = 92834, to = 93133, by = 3)

g <- 1
for(g in 1:nrow(conditions)){  # Start loop for data generation
  
  ### Model conditions
  p <- conditions[g,"P"]
  n <- conditions[g,"N"]
  
  ### Data-gen Output directory
  gen <- paste0("C2PL", "-", p, "-", n)
  odir2 <- file.path(odir, paste0("N", n), gen)
  
  for(r in 1:reps){
    
    ### Replication Output directory
    odir3 <- file.path(odir2, r)
    
    ### Write out replications seeds
    .RNG.seed <- chain1[r]
    .RNG.name <- "base::Super-Duper"
    dump(list = c(".RNG.seed", ".RNG.name"), file = file.path(odir3, "inits1.R"), control = c("showAttributes"), append = FALSE)
    
    .RNG.seed <- chain2[r]
    .RNG.name <- "base::Wichmann-Hill"
    dump(list = c(".RNG.seed", ".RNG.name"), file = file.path(odir3, "inits2.R"), control = c("showAttributes"), append = FALSE)
  }
}
