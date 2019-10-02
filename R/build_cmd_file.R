# Meghan Fager
# Build jags file to estimate bifactor models with interaction effects

build_cmd_file = function(est, r, nThin, nBurn, nSamples, initsFile1, initsFile2, modelFile){
  cmdfile = NULL
  cmdfile = c(cmdfile, paste0("model in ", modelFile))
  cmdfile = c(cmdfile, paste0("load dic"))
  cmdfile = c(cmdfile, paste0("data in data.R"))
  cmdfile = c(cmdfile, paste0("compile, nchains(2)"))
  cmdfile = c(cmdfile, paste0("parameters in ", initsFile1, ", chain(1)"))
  cmdfile = c(cmdfile, paste0("parameters in ", initsFile2, ", chain(2)"))
  cmdfile = c(cmdfile, paste0("initialize"))
  cmdfile = c(cmdfile, paste0("samplers to samplers.txt"))
  cmdfile = c(cmdfile, paste0("update ", nBurn))
  cmdfile = c(cmdfile, paste0("monitor theta, thin(", nThin, ")"))
  cmdfile = c(cmdfile, paste0("monitor beta0, thin(", nThin, ")"))
  cmdfile = c(cmdfile, paste0("monitor beta1, thin(", nThin, ")"))
  cmdfile = c(cmdfile, paste0("monitor beta2, thin(", nThin, ")"))
  cmdfile = c(cmdfile, paste0("monitor deviance, thin(", nThin, ")"))
  cmdfile = c(cmdfile, paste0("monitor pD, thin(", nThin, ")"))
  cmdfile = c(cmdfile, paste0("update ", nSamples))
  cmdfile = c(cmdfile, paste0("coda *, stem(MCMC-)"))
  cmdfile = c(cmdfile, paste0("coda pD, stem(pD-)"))
  cmdfile = c(cmdfile, paste0("coda deviance, stem(DIC-)"))
  cmdfile = c(cmdfile, paste0("exit"))
}

#test <- build_cmd_file(est = est, r = r, nThin = 3, nBurn = 1000, nSamples = 3000, initsFile1 = "init1.R",
#                       initsFile2 = "init2.R", modelFile = "B2PL-Interactions.bug")
#write.table(test, file = "test.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
