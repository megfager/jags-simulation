# Meghan Fager
# Build text file with directories to sim replications

top = "/home/megfager/cat2-sim"
reps = 100
model = "C2PL"
N = c("500", "3000")
P = c(0, 0.25, 0.5)
est = c("None", "All", "Match")

build_array_file = function(reps, top, model, N, P, est){
  array_file = NULL
  all = NULL
  for(r in 1:reps){
    for(n in 1:length(N)){
      for(e in 1:length(est)){
        for(p in 1:length(P)){
          if((P[p]==0 & est[e]=="Match")==TRUE){ next } else {
            array_file = c(file.path(paste0(top, "/", paste0("N", N[n]), "/", model, "-", P[p], "-", N[n], "/", r, "/", est[e])))
            all = rbind(all, array_file)
          }
        }
      }
    }
  }
  return(all)
}

table.500 <- build_array_file(reps, top, model, N = "500", P, est)
table.3000 <- build_array_file(reps, top, model, N = "3000", P, est)
write.table(table.500, "table-500.txt", row.names = FALSE, quote = FALSE, col.names = FALSE)
write.table(table.3000, "table-3000.txt", row.names = FALSE, quote = FALSE, col.names = FALSE)
