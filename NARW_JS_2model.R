# code for fitting Jolly Seber model to North Atlantic right whale sightings (1990-)

# Load the data/inputs----
Year.Max <- 2024
#source("NARW_JS_1data.R")
load(paste0("./data/NARW_JS_inputs_",Year.Max,".Rdata"))
library(parallel)
library(MCMCvis)
library(mcmcplots)
library(nimble)
library(dplyr)

yrs <- 1990:Year.Max
calves <- read.csv("./data/Calves_Observed_1990-present.csv")

whale.data$calves <- (calves %>% filter(Year %in% yrs) %>%
                        #only include calves not known to be lost
                        mutate(Recruits = Tot.Calves-Known.Losses)
                      )$Recruits

# for calf-integrated model
inits$gamma <- c(0.25,rep(NA,whale.constants$n.occasions-2))
inits$phi_calf <- 0.95

# Create CSV of sighting histories
obs <- whale.data$y[,-1]
obs[obs==2] <- 0 #turn to 0/1s
obs <- obs[which(apply(obs,1,sum) > 0),]  #limit to known individuals
write.csv(obs,paste0("./data/sightings_NARW_1990-",Year.Max,".csv"),row.names = F)

# Parameters monitored
params <- c(
  #"phiaf","phiam","phi0", "phi1","phi2","phi3","phi4", 
  "b0","b.age", "b.female","b.regime","a0","a.female",
  #"pcap1", "pcap2", "psi",
  "pi.female", "phi_calf",           
  "sigma.phi", "sigma.p.i", "sigma.p.t",
  "epsilon.phi.t","epsilon.p.t",
  "N", "NF", "NM", "NadF", "Nd",
  "gamma", "B", "aliveT" #,"z"
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
  
  source("NARW_model_w_calfs.txt")
  #source("Pace2017_NARW_model.txt")
  
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
n.iterX   <-  25000
n.burninX <-  5000
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
