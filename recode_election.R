#Recoding Election Results: (Cleaned March 2021 - no changes however)
rm(list = ls())
library(tidyverse)

################################################################################
#### ELECTION RESULTS ##########################################################
################################################################################
load_to_env <- function (rdata, env = new.env()) {
    load(rdata, env)
    return(env)
}
path_dat = "/Users/Ciara/Dropbox/kpop/application/data/"

pres_env <-  load_to_env(paste0(path_dat, "1976-2016-president.RData"))

election_2016 <- pres_env$x %>%
    filter(year == 2016) %>%
    select(state, party, candidatevotes, totalvotes) %>% 
    group_by(state) %>%
    summarize(demtotal = sum(candidatevotes[party == "democrat"], na.rm = TRUE),
              reptotal = sum(candidatevotes[party == "republican"],
                             na.rm = TRUE),
              totalvotes = unique(totalvotes)) %>%
    ungroup()

#not updated march 2021 (no changes)
saveRDS(election_2016, paste0(path_dat, "election.rds"))
