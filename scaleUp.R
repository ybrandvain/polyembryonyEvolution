# Simulating the evolution of polyembryony across parameter space
# October 11th 2022
# Yaniv Brandvain, Alex Harkness & Tanja Pyhajarvi
# Accompanies Brandvain et al. "Reproductive compensation and selection among viable embryos drive the evolution of polyembryony"

# Run a bunch of invasions across parameter values
# Be sure this code is in the same folder as the files below
source("polyembryonyEvolution.R")
source("processPoly.R")


intro.poly.run <- function(simID, n.introductions, n.gen.after.fix, n.times.to.track.fix = 10, n.initial.poly.copies = 1,which.model = 1:5){
  load(simID)
  poly.models       <- data.frame( equalizedW = c(T,T,F,F),  compete = c(T,F,T,F), e1_survive_ok =F)
  poly.models       <- rbind(poly.models ,c(T,F,T))
  n.inds            <- z$params["n.inds"][1,1]
  selfing.rate      <- z$params["selfing.rate"][1,1]
  polyemb.p0        <- as.numeric(n.initial.poly.copies / (2 * z$params["n.inds"]))
  fitness.effects   <- ifelse(  as.character(z$params["fitness.effects"][1,1]) == "uniform", -1,as.numeric(as.character(z$params["fitness.effects"][1,1])))
  dom.effects       <- ifelse(  as.character(z$params["dom.effects"][1,1]) == "uniform", -1,as.numeric(as.character(z$params["dom.effects"][1,1])))
  U       <- z$params["U"][1,1]
  poly.models.list  <- list(equalizedW_compete = 1, equalizedW_nocompete = 2, nonequalizedW_compete = 3,nonequalizedW_nocompete = 4, resurrection = 5)
  poly.models.list  <- poly.models.list[which.model] 
  dist.timing       <- as.numeric(strsplit(as.character(z$params["dist.timing"][1,1]),":")[[1]])
  names(dist.timing)<- c("E","B","L")
  hard.embryo.selection <- ifelse(length(which(names(z$params) == "hard.embryo.selection")) == 0, TRUE, z$params["hard.embryo.selection"][1,1])
  p.poly.mono.geno      <- ifelse(length(which(names(z$params) == "p.poly.mono.geno")) == 0, 0, z$params["p.poly.mono.geno"][1,1])
  intro.results <- lapply(poly.models.list,function(X){
    #poly.models[X,]
    fixed.outcome <- list()
    times.fixed   <- 0
    #    sim.summary   <- data.frame(matrix(nrow = n.introductions, ncol = length(z$params)))
    #    colnames(sim.summary) <- colnames(z$params)
    for(i in 1:n.introductions) { # should be n.introductions
      print(i)
      this.sim <- runSim(n.inds                = n.inds,  
                         n.gen                 = 0, 
                         selfing.rate          = selfing.rate,
                         U                     = U,
                         fitness.effects       = fitness.effects, 
                         dom.effects           = dom.effects ,  
                         gen.after.loss        = 1,  
                         gen.after.fix         = ifelse(times.fixed < n.times.to.track.fix, n.gen.after.fix,1), 
                         polyemb.p0            = polyemb.p0, 
                         introduce.polyem      = 0,
                         equalizedW            = poly.models[X,"equalizedW"] ,
                         compete               = poly.models[X,"compete"],
                         genomes               = z$genome,
                         genome.id             = simID,
                         hard.embryo.selection = hard.embryo.selection,
                         p.poly.mono.geno      = p.poly.mono.geno,
                         just.return.genomes   = times.fixed >= n.times.to.track.fix,
                         dist.timing           = dist.timing,
                         e1_survive_ok         = poly.models[X,"e1_survive_ok"]
      )
      print(":)")
      if(class(this.sim$params$fitness.effects) == "factor"){this.sim$params$fitness.effects <- as.numeric(fitness.effects)}
      if(class(this.sim$params$dom.effects) == "factor"){this.sim$params$dom.effects <- as.numeric(dom.effects)}
      if(i == 1){sim.summary <- as_tibble(this.sim$params)}
      if(i > 1){sim.summary[i,] <- this.sim$params}
      if(this.sim$params[1,"fixed"]){     
        times.fixed <- times.fixed + 1
        if( times.fixed <= n.times.to.track.fix){fixed.outcome[[times.fixed]] <- this.sim}
      }
      print(sprintf("SCALEUPDATE %s, completedIntroduction %s, of %s, timesfixed %s, model %s, genomeID %s", paste(rbind(strsplit( date(), " ",)[[1]][2:4]),collapse="_"),  i, n.introductions, times.fixed ,names(poly.models.list)[X], simID ))
    }
    return(list(sim.summary = sim.summary,fixed.outcome = fixed.outcome))
  })
  sim.summary   <- as_tibble(do.call(rbind,lapply(intro.results, function(X){X$sim.summary})))
  fixed.details <- bind_rows(lapply(intro.results, function(TMp){
    bind_rows(lapply(TMp$fixed.outcome, function(TMP){
      bind_cols(as_tibble(nest(TMP$gen.summary, data = everything())) %>% rename(gen.summary = data),
                as_tibble(nest(TMP$genome, data = everything())  %>% rename(genome = data) ),
                as_tibble(TMP$params, data = everything()))
    }))
  }))
  return(list(sim.summary = sim.summary, fixed.details = fixed.details))
}