# code for fitting Jolly Seber model to North Atlantic right whale sightings (1990-)

# Load the data/inputs----
#source("NARW_JS_1data.R")
load("NARW_JS_inputs.Rdata")
library(parallel)
library(MCMCvis)
library(mcmcplots)
library(nimble)

# Create CSV of sighting histories
obs <- whale.data$y[,-1]
obs[obs==2] <- 0 #turn to 0/1s
obs <- obs[which(apply(obs,1,sum) > 0),]  #limit to known individuals
write.csv(obs,paste0("sightings_NARW_1990-",Year.Max,".csv"),row.names = F)

# Parameters monitored
params <- c(
  #"phiaf","phiam","phi0", "phi1","phi2","phi3","phi4", 
  "b0","b.age", "b.female","b.regime","a0","a.female",
  #"pcap1", "pcap2", 
  "pi.female", "psi",             
  "sigma.phi", "sigma.p.i", "sigma.p.t",
  "epsilon.phi.t","epsilon.p.t",
  "gamma", "entry", "N", "NF","NM", "B", "Nd","aliveT"
  )

# multiple chain inits
info <- list(
  list(seed = 1,
       inits = inits),
  list(seed = 2,
       inits = inits),
  list(seed = 3,
       inits = inits)
)

# NIMBLE model fitting code ----
run_MCMC_allcode <- function(info, whale.data, whale.constants, params,
                             n.iter, n.burnin, n.thin){
  library(nimble)
  
  source("Pace2017_NARW_model.txt")
  
  # testing only, not when running cluster
  # n.iter=2000; n.burnin=1000; n.thin=1
  # info <- info[[1]]
  
  print(Sys.time())
  whale.model <- nimbleModel(code=whale.code.pop, name="whale", data=whale.data,
                             constants = whale.constants, inits = info$inits)
  whale.model$calculate()
  
  print(Sys.time())
  whaleConf <- configureMCMC(whale.model, print=F)
  
  whaleConf$addMonitors(params)
  whaleConf$removeSamplers(c('a0','a.female','b0','b.female','b.regime','b.age'))
  whaleConf$addSampler(target=c('a0','a.female','b0','b.female','b.regime','b.age'),type="AF_slice")
  
  print(Sys.time())
  
  whaleMCMC <- buildMCMC(whaleConf)
  CwhaleMCMC <- compileNimble(whale.model, whaleMCMC, resetFunctions = TRUE)
  print(Sys.time())
  
  out <- runMCMC(CwhaleMCMC$whaleMCMC, niter=n.iter, nburnin=n.burnin, thin=n.thin, 
                 setSeed = FALSE) #info$seed)
  print(Sys.time())
  
  return(out)
}

# Fit multiple chains to model with cluster ----
this_cluster <- makeCluster(3)

## Basic MCMC parameters (~7.4hrs for 25k)
n.iterX   <-  25000  #600000
n.burninX <-  5000  #20000
n.thinX <- 5

(start <- Sys.time())
outL <- parLapply(cl = this_cluster, X = info[1:length(this_cluster)],
                  fun = run_MCMC_allcode,
                  params = params,
                  whale.data = whale.data,
                  whale.constants = whale.constants,
                  n.iter=n.iterX, n.burnin=n.burninX, n.thin=n.thinX)
(end <- Sys.time()-start)

stopCluster(this_cluster)

save(outL, run_MCMC_allcode, whale.data, whale.constants, #whale.tbl,
     file = paste0("./out/outL_NARW_1990-",Year.Max,".Rdata"))

# Examine output----
load(file = paste0("./out/outL_NARW_1990-",Year.Max,".Rdata"))

param.rm <- paste(c("N","B","entry","epsilon","gamma","aliveT"),collapse="|")

MCMCsummary(outL, round = 3, 
            excl = colnames(outL[[1]])[grep(pat = param.rm,colnames(outL[[1]]))])

mcmcplot(as.mcmc.list(lapply(outL,as.mcmc)), dir = paste0("./out"), filename="whale_popsize_parms",
         parms = colnames(outL[[1]])[-grep(pat = param.rm,colnames(outL[[1]]))])

mcmcplot(as.mcmc.list(lapply(outL,as.mcmc)), dir = paste0("./out"), filename="whale_popsize_N",
         parms = colnames(outL[[1]])[grep(pat = "N\\[", colnames(outL[[1]]))])
