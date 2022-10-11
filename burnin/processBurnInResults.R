library(tidyverse)
library(ggthemes)
library(ggrepel)
library(RColorBrewer)
library(wesanderson)
library(gridExtra)
library(infer)

getGenSummaryNoMuts <- function(file){
  load(file)
  z$gen.summary   %>% 
    select(-muts) %>% 
    mutate(dist.timing      = z$params$dist.timing,
           U                = z$params$U,
           selfing.rate     = z$params$selfing.rate,
           fitness.effects  = z$params$fitness.effects,
           dom.effects      = z$params$dom.effects,
           equalizedW       = z$params$equalizedW,
           compete          = z$params$compete,   
           run              = str_remove(file, ".*_i"),
           file             = file)
}
# EXAMPLE: getGenSummaryNoMuts(file = ("BurnIn-1_-1_0.5_1/BurnIn-1_-1_0.5_1_0/BurninGenome_w-1_h-1_U0.5_t1_S0_i1"))


getGenSummaryNoMutsOverReps <- function(path.to.params){
  files <- paste(path.to.params,list.files(path.to.params), sep = "/")
  param.summary <- do.call(rbind,lapply(files, getGenSummaryNoMuts))
  return(param.summary)
} 

#getGenSummaryNoMutsOverReps(path.to.params = "LongBurnIns/BurnIn2000-1_-1_0.5_1_0/")

burnInFinal <- function(){
  lapply(paste("LongBurnIns", list.files("LongBurnIns/"),sep = "/"), 
         function(PARAMZ){
           print(PARAMZ)
           recover()
           getGenSummaryNoMutsOverReps(PARAMZ) %>%
             filter(gen == max(gen)) %>%
             mutate(fitness.effects = as.character(fitness.effects ), 
                    dom.effects     = as.character(dom.effects))
         }) %>%
    bind_rows()
}


# For muts over time 
mutsOverTime <- function(){
  lapply(paste("LongBurnIns", list.files("LongBurnIns/"),sep = "/"), 
         function(PARAMZ){
           print(PARAMZ)
           getGenSummaryNoMutsOverReps(PARAMZ) %>%
             select(two_n,E_perhaploidgenome,L_perhaploidgenome, gen, U, selfing.rate, fitness.effects, dom.effects, run, file) %>%
             mutate(fitness.effects = as.character(fitness.effects ), 
                    dom.effects     = as.character(dom.effects))
         }) %>%
    bind_rows()
}

# Generate a high level summary of burn in results for each generation
all_burn_ins <- mutsOverTime()   
write_csv(all_burn_ins, path = "allburninsgen.csv")

# Generate a high level summary of burn in results for the final generation 
burn_in_final <- burnInFinal()
write_csv(burn_in_final, path = "allburninsfinal.csv")
