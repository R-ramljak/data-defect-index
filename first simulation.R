library(tidyverse)

set.seed(3)


pop <- tibble(ID = paste0("A", 1:100000), 
              y = sample(1:100, 100000, replace = T)) %>% 
  mutate(prob.help = case_when(between(y, 90, 100) ~ 1,
                               TRUE ~ 0)) %>% 
  mutate(prob = case_when(prob.help == "1" ~ (1 / sum(prob.help)) / length(prob.help),
                          prob.help == "0" ~ (1 / (length(ID) - sum(prob.help))) / length(prob.help))) # possible sampling error, inclusion probs dont add up to 1 perfectly


# 100 SRS without replacement
n.SRS <- map(1:100, ~sample(pop$ID, 1000, replace = F)) %>%
  set_names(paste("Sample", 1:100)) %>% 
  map(as_tibble) %>% 
  map(~select(., ID = value)) %>% 
  map(~mutate(., selected = 1)) %>% 
  map(~right_join(., pop, by = "ID")) %>% 
  map(~mutate(., selected = case_when(is.na(selected) ~ 0,
                                      TRUE ~ selected)))
hist.SRS <- n.SRS %>% 
  map_df(~summarise(., mean = mean(y)))


rho.SRS2 <- n.SRS %>% 
  map_df(~summarise(., rho = cor(selected, y))) %>% 
  summarise(mean.rho = mean(rho), var.rho = var(rho)) # mean and variance is nearly 0 for SRS

# 100 Non probability samples without replacement
n.NP <- map(1:100, ~sample(pop$ID, 1000, replace = F, prob = pop$prob)) %>%
  set_names(paste("Sample", 1:100)) %>% 
  map(as_tibble) %>% 
  map(~select(., ID = value)) %>% 
  map(~mutate(., selected = 1)) %>% 
  map(~right_join(., pop, by = "ID")) %>% 
  map(~mutate(., selected = case_when(is.na(selected) ~ 0,
                                      TRUE ~ selected)))

hist.NP <- n.NP %>% 
  map_df(~summarise(., mean = mean(y)))


rho.NP2 <- n.NP %>% 
  map_df(~summarise(., rho = cor(selected, y))) %>% 
  summarise(mean.rho = mean(rho), var.rho = var(rho)) # mean and variance is nearly 0 for SRS


# expected value of rho is dependent on the difference in inlcusion probabilities (selection bias) and the size of range of y 
