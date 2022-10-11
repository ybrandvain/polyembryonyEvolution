# Simulating the evolution of polyembryony
# October 11th 2022
# Yaniv Brandvain, Alex Harkness & Tanja Pyhajarvi
# Accompanies Brandvain et al. "Reproductive compensation and selection among viable embryos drive the evolution of polyembryony"
library(tidyverse)

##########################
# How to live 
##########################
initializeGenomes <- function(n.inds, genomes = NULL){
  if(!is.null(genomes)){return(genomes)}
  tibble(ind    = rep(1:n.inds,2), 
         id     = 10,    # 10 means monoembryony. 11 means poly
         s      =  0,
         h      =  1,    # might change for convinience.. currently meaningless
         timing = "D")
}
introducePoly     <- function(genomes, polyemb.p0){
  genomes %>% 
    dplyr::mutate(id = case_when(id %in% 10:11 ~ 
                                   sample(x = as.double(10:11),
                                          size = n(),
                                          replace = TRUE,
                                          prob = c(1 - polyemb.p0, polyemb.p0)), #10 means monoembryony. 11 means poly
                                 !id %in% 10:11 ~ id))
}
addMutations      <- function(tmp.genomes, U, fitness.effects, dom.effects, dist.timing, n.inds, smallest.s){
  getVal <- function(thing, num, this.min = 0, this.max = 1, prelim.vals = NULL){
    if(is.null(prelim.vals)){ prelim.vals <- runif(n = num, min = this.min, this.max) }
    if(thing == "uniform") {   return(prelim.vals)  }
    if(is.numeric(thing))  {   return(rep(thing, num))}
    recover() # we can transform normal to an dist here
  }
  inds   <-   unique(tmp.genomes$ind)
  n.muts <-   rpois(1, U * length(inds) )
  bind_rows(tmp.genomes,                                                     # old genome
            tibble(ind    = sample(inds, size = n.muts, replace = TRUE),     # assigning muts to inds
                   timing = sample(x       = names(dist.timing),             # timing of selection
                                   size    = n.muts,                         # number of muts defined
                                   replace = TRUE,                           # obviously
                                   prob    = dist.timing),                   # right now all muts equi-probable. can change this
                   id = getVal(thing = "uniform", num = n.muts, this.min = smallest.s/n.inds),
                   s  = getVal(thing = fitness.effects, num = n.muts, this.min = smallest.s/n.inds, prelim.vals = id),
                   h  = getVal(thing = dom.effects, num = n.muts))  %>%    # s from uniform as described in ms. can change
              dplyr::mutate(s = ifelse(timing == "B", (1-sqrt(1-s)),s)))   
}
getFitness        <- function(tmp.genomes, dev.to.exclude, adult = TRUE){
  ind.genomes <- tmp.genomes  %>% ungroup()                           %>%
    dplyr::mutate(dup = as.numeric(duplicated(tmp.genomes)))          %>%
    dplyr::filter(!duplicated(tmp.genomes, fromLast = TRUE) & 
                    !timing %in% dev.to.exclude)                               %>%
    dplyr::mutate(w.loc = 1 - (1-dup) * h * s - dup *s)       
  if(adult) {ind.genomes  <- ind.genomes %>% mutate(mono = NA, selfed = NA)%>% group_by(ind)}
  if(!adult){ind.genomes  <- ind.genomes %>% group_by(mating, embryo)}
  ind.genomes                                                         %>%   
    dplyr::summarise(w = prod(w.loc), mono = mean(mono), selfed = mean(as.numeric(selfed)))    
}
findMates         <- function(adult.fitness, selfing.rate, n.inds, epsilon = 1e-16){
  outbred.parents <- replicate(3,with(adult.fitness, sample(ind, size = n.inds, replace = TRUE, prob = (w + epsilon)))) 
  self <- data.frame(matrix(rbinom(n = n.inds * 2, size = 1, prob = selfing.rate),ncol = 2))
  colnames(outbred.parents) <- c("mom","dad1","dad2")
  as_tibble(outbred.parents) %>%
    mutate(dad1 = ifelse(self$X1 == 1,mom,dad1), # selfing
           dad2 = ifelse(self$X2 == 1,mom,dad2), # selfing
           mating = 1:n.inds)                    # indexing
}
embryoFitness     <- function(tmp.embryos, p.poly.mono.geno, d2exclude = "L"){
  tmp.embryos                                                 %>% 
    dplyr::group_by(mating)                                   %>%
    dplyr::mutate(mono_geno    = sum(parent == "mat" & id == 10)/2,
                  mono_chance  = rbinom(n = 1, size = 1, prob = 1 - p.poly.mono.geno),
                  mono = mono_geno * mono_chance)   %>%    
    #dplyr::filter( !(mono == 1 & embryo == "e2") )           %>% 
    dplyr::group_by(embryo, add = TRUE)                       %>% ungroup() %>%
    select(- parent)                                          %>%
    getFitness(dev.to.exclude = d2exclude, adult = FALSE)           %>% ungroup() 
}

favoriteChild     <- function(temp.kidsW, equalizedW = TRUE, compete = TRUE, epsilon = 1e-16, hard.embryo.selection,e1_survive_ok){
  if( hard.embryo.selection){
    temp.kidsW <- temp.kidsW                                    %>%
      mutate(alive = rbinom(n = n(),size = 1,prob = w ))    
    if(equalizedW &  hard.embryo.selection){
      temp.kidsW <- temp.kidsW                                   %>%
        group_by(mating)                                         %>%
        mutate(alive2 = sum(alive * as.numeric(embryo == "e1"))) %>% ungroup() %>%
        filter(alive2 == 1) %>% 
        select(-alive2)
      if(e1_survive_ok){temp.kidsW <- temp.kidsW %>% mutate(alive = 1)} #weird flag to allow dead embryos to surive if e1 survived under experiental form of polyembryony 
    }
    temp.kidsW <- temp.kidsW                                       %>% 
      dplyr::filter(alive == 1)                                    
  }
  temp.kidsW <- temp.kidsW                          %>% 
    dplyr::filter( !( mono == 1 & embryo == "e2") ) %>%  
    dplyr::group_by(mating) 
  if(!compete){
    temp.kidsW <- temp.kidsW  %>% 
      mutate(w = max(w))
    # just flip
  }
  temp.kidsW <- temp.kidsW   %>%  
    sample_n(1,weight = w + epsilon )  %>% ungroup()  
  return(temp.kidsW)
}

grabInds          <- function(selectedEmbryos, embryos){
  embryoId  <- dplyr::mutate(selectedEmbryos, winners = paste(mating,embryo)) %>% 
    dplyr::select(winners) %>% 
    pull()
  embryos %>% 
    dplyr::mutate(combo = paste(mating, embryo)) %>% 
    dplyr::filter(combo %in% embryoId) %>%
    dplyr::select(ind = mating, id = id, s = s, h = h, timing = timing) 
}
doMeiosis         <- function(tmp.genomes, parents){
  #  to.meios <- nest(tmp.genomes, data = c(id, s, h, timing))         %>%    # to new
  to.meios <- nest(tmp.genomes, -ind)         %>% 
    arrange(ind)                              %>%  # major bug fix
    dplyr::slice(parents)                     %>% 
    dplyr::mutate(mating = 1:length(parents)) %>%
    dplyr::select(mating, data)               %>%
    #unnest(cols = c(data) )           # to new
    unnest()   
  poly.allele     <- to.meios$id %in% 10:11 
  to.transmit.hom <- duplicated(to.meios, fromLast = FALSE) & !poly.allele
  pick.one        <- duplicated(to.meios, fromLast = TRUE) | poly.allele  | to.transmit.hom
  # first tranmit het alleels at random
  random.alleles <- to.meios                                    %>% 
    dplyr::filter(!pick.one)                                    %>%
    dplyr::filter(rbinom(n = sum(!pick.one), size = 1, prob = .5) == 1)
  # next be sure to transmit exactly one hom allele
  hom.transmit   <-  to.meios                                   %>% 
    dplyr::filter(to.transmit.hom)
  # finally, pick one of the alelles at the developement locus
  dev.transmit <- to.meios                                      %>% 
    dplyr::filter(poly.allele)                                  %>%
    dplyr::group_by(mating)                                     %>%
    sample_n(size = 1, replace = FALSE)                         %>%
    ungroup()
  # now shove them all together
  bind_rows(random.alleles, hom.transmit, dev.transmit )        %>% 
    nest(-mating)                                               %>%
    #   nest(data = c(id, s, h, timing))                            %>%  # to new
    dplyr::arrange(mating) 
}
makeBabies        <- function(tmp.genomes, mates){   
  selfed <- mates%>%
    mutate(e1 = mom == dad1, e2 = mom ==dad2)%>% select(mating, e1, e2) %>% 
    gather(key = embryo, value = selfed, - mating) 
  bind_cols(
    doMeiosis(tmp.genomes, mates$mom)  %>% select(mating = mating, mat = data),
    doMeiosis(tmp.genomes, mates$dad1) %>% select(e1 = data),
    doMeiosis(tmp.genomes, mates$dad2) %>% select(e2 = data))                       %>%
    gather(key = embryo, value = pat, -mating, - mat)                               %>%
    full_join(selfed, by = c("mating", "embryo"))                                   %>%
    gather(key = parent, value = haplo, -mating, - embryo,-selfed)                  %>% # syngamy
    unnest(haplo)    
  #   unnest(cols = c(haplo))                                                         # to new
}
summarizeGen      <- function(tmp.genomes, mates, embryos, selectedEmbryos){
  n.intial  <- embryos$mating%>%unique()%>%length()
  n.survive <- tmp.genomes$ind %>%unique()%>%length()
  muts_e1      <- embryos %>% filter(embryo == "e1")%>% 
    mutate(s = ifelse(id %in%  c(10,11), id,s)) %>% 
    group_by(id,s,h,timing) %>% 
    tally() %>% 
    ungroup()
  #
  muts_survive <- tmp.genomes %>%   # these are allele counts after embryo selection 
    mutate(s = ifelse(id %in%  c(10,11), id,s)) %>% 
    group_by(id,s,h,timing) %>% 
    tally() %>% 
    ungroup()
  selfing.info <- tibble(   realized_selfing = selectedEmbryos %>% summarise(mean(selfed)) %>% pull(), 
                            primary_selfing = embryos %>% summarise(mean(selfed))%>% pull())
  realized_selfing_by_mono <- selectedEmbryos  %>% 
    group_by(mono) %>% 
    summarize(reself = mean(selfed)) %>% 
    summarize(mono_realized_selfing =sum(as.numeric(mono == 1) *reself ) / sum(as.numeric(mono == 1)),
              poly_realized_selfing =sum(as.numeric(mono == 0) *reself )/ sum(as.numeric(mono == 0)))
  pop.stats    <- bind_cols(
    muts_survive %>% 
      filter(timing == "D") %>%
      rename(n_alleles = n) %>%
      summarise(two_n = sum(n_alleles), 
                freq_poly_after_sel = sum(as.numeric(s == 11) * n_alleles/ two_n )),
    muts_e1%>% 
      summarise(freq_poly_e1 = sum(as.numeric(s == 11))/ (n.intial*2)  )
  )
  # muts / ind by tming survivors
  mut.per.ind.survive <- muts_survive %>% 
    filter(timing !="D") %>% 
    group_by(timing) %>% 
    rename(n_alleles = n) %>%
    tally(wt = n_alleles ) %>%
    mutate(timing = paste(timing,"survive",sep = "_"))
  names(mut.per.ind.survive)[2] <- "n" 
  mut.per.ind.survive <- mut.per.ind.survive %>%
    mutate(muts_per_ind_survive = n / pop.stats$two_n)%>% 
    select(-n) %>%
    spread(key = timing, value = muts_per_ind_survive) 
  if(nrow(mut.per.ind.survive) >0){
    names(mut.per.ind.survive) <- paste(names(mut.per.ind.survive),"perhaploidgenome",sep="_")
  }
  # muts / ind by tming all 
  mut.per.ind.e1 <- muts_e1 %>% 
    filter(timing !="D") %>% 
    group_by(timing) %>% 
    rename(n_alleles = n) %>%
    tally(wt = n_alleles ) %>%
    mutate(timing = paste(timing,"e1",sep = "_"))
  names(mut.per.ind.e1)[2] <- "n" 
  mut.per.ind.e1 <- mut.per.ind.e1 %>%
    mutate(muts_per_ind_e1 = n / pop.stats$two_n)%>% 
    select(-n) %>%
    spread(key = timing, value = muts_per_ind_e1) 
  if(nrow(mut.per.ind.e1) >0){
    names(mut.per.ind.e1) <- paste(names(mut.per.ind.e1),"perhaploidgenome",sep="_")
  }
  
  w.all.stats <- full_join(
    embryoFitness(embryos,p.poly.mono.geno=0,d2exclude = "L") %>% 
      select(mating, embryo, w_early_all = w) %>%ungroup(),
    embryoFitness(embryos,p.poly.mono.geno=0,d2exclude = "E") %>% 
      select(mating, embryo, w_late_all = w)%>%ungroup(), 
    by = c("mating" , "embryo") ) %>% 
    summarize(mean_w_early_all = mean(w_early_all), 
              mean_w_late_all  = mean(w_late_all),
              var_w_early_all  = var(w_early_all), 
              var_w_late_all   = var(w_late_all),
              cor_w_early_late_all = cor(w_early_all, w_late_all))
  
  
  w.summary <- bind_cols(tibble(
    early_w_selected = getFitness(tmp.genomes,dev.to.exclude = "L", adult = TRUE) %>% select(w) %>%pull(),
    late_w           = getFitness(tmp.genomes,dev.to.exclude = "E", adult = TRUE) %>% select(w) %>%pull()) %>%
      summarise(mean_w_late_survivors = mean(late_w), 
                mean_w_early_survivors = mean(early_w_selected), 
                var_w_late_survivors  = var(late_w), 
                var_w_early_survivors = var(early_w_selected), 
                cor_w_early_late_survivors = cor(early_w_selected,late_w)) ,
    w.all.stats)
  
  w_early_by_self_mono_survive <- selectedEmbryos %>% 
    group_by(mono,selfed)%>% 
    summarize(w = mean(w)) %>% 
    ungroup() %>% 
    summarize(w_early_mono_selfed_survive = sum(w * as.numeric(mono==1 & selfed == 1)) /sum(as.numeric(mono==1 & selfed == 1)),
              w_early_mono_out_survive    = sum(w * as.numeric(mono==1 & selfed == 0)) /sum(as.numeric(mono==1 & selfed == 0)),
              w_early_poly_selfed_survive = sum(w * as.numeric(mono==0 & selfed == 1)) /sum(as.numeric(mono==0 & selfed == 1)),
              w_early_poly_out_survive    = sum(w * as.numeric(mono==0 & selfed == 0)) /sum(as.numeric(mono==0 & selfed == 0)))
  list(genome    = tmp.genomes, 
       #summaries = bind_cols(nest(muts, data= everything() ) %>% rename(muts = data),  # to new
       summaries = bind_cols(nest(muts_survive)  %>% rename(muts_survive = data),
                             nest(muts_e1)       %>% rename(mutse1 = data),
                             pop.stats,  
                             selfing.info , 
                             realized_selfing_by_mono,
                             mut.per.ind.e1, 
                             mut.per.ind.survive, 
                             w.summary,
                             w_early_by_self_mono_survive ))
}

##########################
##########################
##########################

# running one generation
oneGen <- function(tmp.genomes, n.inds, selfing.rate, U, fitness.effects, dom.effects, dist.timing, equalizedW, compete, just.return.genomes,  p.poly.mono.geno,  hard.embryo.selection, smallest.s, e1_survive_ok){
  tmp.genomes     <- addMutations(tmp.genomes, U, fitness.effects, dom.effects, dist.timing, n.inds, smallest.s)
  adult.fitness   <- getFitness(tmp.genomes, dev.to.exclude = "E")                       # Adult Fitness
  mates           <- findMates(adult.fitness, selfing.rate, n.inds = n.inds)             # Mating / Selection
  embryos         <- makeBabies(tmp.genomes, mates)                                      # meiosis is in here too
  kidsW           <- embryoFitness(embryos, p.poly.mono.geno)
  selectedEmbryos <- favoriteChild(kidsW, equalizedW = equalizedW, compete = compete,   hard.embryo.selection =  hard.embryo.selection, e1_survive_ok = e1_survive_ok)                                                # pick your child !
  tmp.genomes     <- grabInds(selectedEmbryos = selectedEmbryos, embryos = embryos) %>%  # extract the genomes of selected embryos from our chosen children
    mutate(ind = as.numeric(factor(rank(ind, ties.method = "min"))))                     # ugh.. this last line is kinda gross. but necessary. in means inds are numbered 1:n... this is importnat for  other bits above
  if(just.return.genomes){return(list(genome = tmp.genomes))}
  summarizeGen(tmp.genomes, mates, embryos, selectedEmbryos)
}
# running for a bunch of generations
runSim <- function(n.inds = 1000, selfing.rate = 0, U = .5, fitness.effects  = "uniform", 
                   dom.effects = "uniform", n.gen  = 1000, dist.timing  = c(E = 1/2, B = 0, L = 1/2), 
                   equalizedW = TRUE, compete = TRUE , p.poly.mono.geno = 0,  hard.embryo.selection = TRUE,
                   introduce.polyem = Inf, polyemb.p0  = .01, genomes = NULL, genome.id = NA,
                   gen.after.loss   = 1,gen.after.fix    = 1, just.return.genomes = FALSE, smallest.s=20, e1_survive_ok = FALSE){
  # n.inds           =      1000, 
  # selfing.rate     =         0, # recall selfing = 0 is RANDOM MATING and DOES NOT PRECLUDE SELFING
  # U                =         1, 
  # fitness.effects  = "uniform" or -1, mean uniform, any other number is a fixed value  # need to implement more options. currently takes a fixed val or "uniform"
  # dom.effects      = "uniform",or -1, mean uniform # need to implement more options. currently takes a fixed val or "uniform"
  # n.gen            =      1000, # prob should add an option like "until lost/fixed"
  # dist.timing      = c(E = 1/3, B = 1/3, L = 1/3), add what you like! 
  #                  = numeric options: 1 = c(E = 1/2, B = 0, L = 1/2)
  #                                     2 = c(E = 1, B = 0, L = 0) 
  #                                     3 = c(E = 0, B = 0, L = 1)
  #                                     4 = c(E = 0, B = 1, L = 0)
  # equalizedW       = TRU. Eshould the expected number of embryos produced by mono and poplyembryonic genos be equivalent? Achieved by group sel at level of mom
  # compete          = compete = TRUE , should embryos be chosen at random or with respect to their fitnesses? 
  # introduce.polyem = Inf       , # gen at which we introduce polyembryony allele
  # polyemb.p0       = .01       , # freq of polyembryony allele once introduced
  # genomes          = NULL        # An option to hand genomes from a previous run
  # genome.id        = NULL 
  # gen.after.loss   = 1
  # gen.after.fix    = 1
  # p.poly.mono.geno = 0. This is the proportion of "monoembryonic" genotypes that are polyrmbryonic
  # recommend a value of approx 0.1 for invasion of / bunn in with soft selection
  # note: change the range of s here
  # hard.embryo.selection = TRUE. is selection on embryos hard or soft? 
  # ignoring just.return.genomes
  # smallest.s the lowest selection coeficeint times N (so min.fitness effect is  smallest.s/n)
  if(fitness.effects == -1){fitness.effects <- "uniform"}
  if(dom.effects == -1){dom.effects <- "uniform"}
  if(length(dist.timing) == 1){
    dist.timing <- list(c(E = 1/2, B = 0, L = 1/2), c(E = 1, B = 0, L = 0), c(E = 0, B = 0, L = 1),c(E = 0, B = 1, L = 0))[[dist.timing]]
  }
  g             <- 0
  g.after.fix   <- 0 
  g.after.loss  <- 0 
  #g.since.fixed <- 0
  ans           <- list(genome = initializeGenomes(n.inds = n.inds, genomes)) # Make genomes   # will need to keep track of things... but what?
  keep.going = TRUE
  while(keep.going){  # or stopping rule tbd   # i realize this should be a for loop, but sense that a while loop will give me flexibility for broader stopping rules
    if((g %% 100) == 0){print(sprintf("SIMULATION progress, gen %s, genomeID %s",g,genome.id))}
    if(g == introduce.polyem){ans$genome <- introducePoly(ans$genome, polyemb.p0)} # introduce polyembryony allele
    g                 <- g + 1
    print(g)
    ans               <- oneGen(ans$genome, n.inds, selfing.rate, U, fitness.effects, 
                                dom.effects, dist.timing, equalizedW = equalizedW, compete = compete, 
                                just.return.genomes = just.return.genomes,  
                                p.poly.mono.geno = p.poly.mono.geno,  
                                hard.embryo.selection =  hard.embryo.selection, 
                                smallest.s = smallest.s, e1_survive_ok =e1_survive_ok)
    if(!just.return.genomes){
      if(g == 1){ gen.summary    <- as_tibble(ans$summaries %>% mutate(gen = g))}
      if(g > 1){ gen.summary[g,] <- ans$summaries %>% mutate(gen = g)}
    }
    status <- t(ans$genome  %>% 
                  filter(timing == "D")  %>% 
                  summarise(loss = sum(id == 11) == 0 , fix = sum(id == 10) == 0) )
    g.after.fix    <- g.after.fix   + as.numeric(status["fix",])
    g.after.loss   <- g.after.loss  + as.numeric(status["loss",])
    if(g >= n.gen   &  (g.after.loss >=  gen.after.loss)   |   (g.after.fix >=  gen.after.fix)){keep.going = FALSE} 
  }
  fixed <- ifelse(introduce.polyem == Inf, NA, ifelse(g.after.fix >0, TRUE, FALSE))
  lost <- ifelse(introduce.polyem == Inf, NA, ifelse(g.after.fix >0, TRUE, FALSE))
  gen.after.fixed.or.lost <- ifelse(introduce.polyem == Inf, NA, ifelse(g.after.fix >0, g.after.fix, g.after.loss ))
  final.status <- ifelse(polyemb.p0 == 0 | (introduce.polyem > n.gen), "burnin", ifelse(status["loss",1], "loss", ifelse(status["fix",1], "fix","stopped")))
  params <- data.frame(n.inds = n.inds, 
                       selfing.rate = selfing.rate, 
                       U = U, 
                       fitness.effects = fitness.effects,
                       dom.effects = dom.effects, 
                       n.gen = n.gen, 
                       g = g, 
                       dist.timing = paste(round(dist.timing, digits = 2), collapse = ":"),
                       introduce.polyem = introduce.polyem, 
                       polyemb.p0  = polyemb.p0 , 
                       existing.genome = !is.null(genomes), 
                       genom.id = genome.id,  
                       last.gen = g, 
                       gen.after.fixed.or.lost  = gen.after.fixed.or.lost, 
                       fixed = fixed, 
                       equalizedW = equalizedW, 
                       compete = compete,
                       e1_survive_ok = e1_survive_ok,
                       hard.embryo.selection = hard.embryo.selection,
                       p.poly.mono.geno=p.poly.mono.geno, 
                       smallest.s= smallest.s,
                       e1_survive_ok =e1_survive_ok)
  if(just.return.genomes){gen.summary<-NULL}
  if(!just.return.genomes){
    ans$summaries <- oneGen(ans$genome, n.inds, selfing.rate, U, fitness.effects, 
                            dom.effects, dist.timing, equalizedW = equalizedW, 
                            compete = compete, just.return.genomes = just.return.genomes,  
                            p.poly.mono.geno = p.poly.mono.geno,  hard.embryo.selection =  hard.embryo.selection, 
                            smallest.s = smallest.s,e1_survive_ok =e1_survive_ok)$summaries 
    final.mean_w_early_all <- round(ans$summaries$mean_w_early_all,digits = 3)
    final.mean_w_late_all  <- round(ans$summaries$mean_w_late_all,digits = 3)
    final.n <- ans$summaries$two_n/2
    final.L.perhaploidgenome <- ifelse("L_perhaploidgenome" %in% names(ans$summaries),round(ans$summaries$L_perhaploidgenome,digits = 3),0)
    final.E.perhaploidgenome <- ifelse("E_perhaploidgenome" %in% names(ans$summaries),round(ans$summaries$E_perhaploidgenome,digits = 3),0)
    print(sprintf("SIMULATION done, stoppedGen %s, status %s, finalMeanWlateAll %s, finalMeanWearlyAll %s, Lperdiploidgenome %s, Eperdiploidgenome %s, genomeID %s", g, final.status, final.mean_w_late_all, final.mean_w_early_all, final.L.perhaploidgenome, final.E.perhaploidgenome, genome.id))
  }
  return(list(genome = ans$genome, gen.summary = gen.summary,params = params))
}


