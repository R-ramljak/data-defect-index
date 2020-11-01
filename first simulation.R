library(tidyverse)
library(Hmisc)

set.seed(3)


pop <- tibble(ID = paste0("A", 1:100000), 
              y = rexp(100000, rate = .2)) %>% 
  mutate(gen.inclusion.prob = 1 / (length(ID))) %>% 
  arrange(y) %>% 
  mutate(rank = percent_rank(y)) %>% 
  mutate(prob.help = case_when(between(percent_rank(y), 0, 0.2) | # min(y) < y < 1
                                 between(percent_rank(y), 0.9, 1) ~ 1, # 10 < y < max(y)
                               TRUE ~ 0)) %>% 
  mutate(adj.inclusion.prob = case_when(prob.help == 1 ~ 1 / sum(prob.help),
                                        prob.help == 0 ~ 1 / sum(length(ID) - prob.help)))

pop <- tibble(ID = paste0("A", 1:100000), 
              y = rnorm(100000, mean = 40, 10)) %>% 
  mutate(gen.inclusion.prob = 1 / (length(ID))) %>% 
  arrange(y) %>% 
  mutate(rank = percent_rank(y)) %>% 
  mutate(prob.help = case_when(between(percent_rank(y), 0, 0.5) ~ 1, #| # min(y) < y < 1
                                 # between(percent_rank(y), 0.9, 1) ~ 1, # 10 < y < max(y)
                               TRUE ~ 0)) %>% 
  mutate(adj.inclusion.prob = case_when(prob.help == 1 ~ 1 / sum(prob.help),
                                        prob.help == 0 ~ 1 / sum(length(ID) - prob.help)))


# 100 SRS without replacement
n.SRS <- map(1:100, ~sample(pop$ID, 1000, replace = T)) %>%
  set_names(paste("Sample", 1:100)) %>% 
  map(as_tibble) %>% 
  map(~select(., ID = value)) %>% 
  map(~mutate(., selected = 1)) %>% 
  map(~right_join(., pop, by = "ID")) %>% 
  map(~mutate(., selected = case_when(is.na(selected) ~ 0,
                                      TRUE ~ selected)))


# 100 Non probability samples without replacement
n.NP <- map(1:100, ~sample(pop$ID, 10000, replace = T, prob = pop$adj.inclusion.prob)) %>%
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
  map_df(~summarise(., target.mean = weighted.mean(y, w = selected), target.sd = sqrt(wtd.var(y, weights = selected)),
                    rho = cor(selected, y))) %>% 
  summarise(mean(target.mean), mean(target.sd),
            mean.rho = mean(rho), var.rho = var(rho)) # mean and variance is nearly 0 for SRS



hist.NP <- n.NP %>% 
  map_df(~summarise(., mean = mean(y)))


rho.NP2 <- n.NP %>% 
  map_df(~summarise(., target.mean = weighted.mean(y, w = selected), target.sd = sqrt(wtd.var(y, weights = selected)),
                    rho = cor(selected, y))) %>%  
  summarise(mean(target.mean), mean(target.sd),
            mean.rho = mean(rho), var.rho = var(rho)) # mean and variance is nearly 0 for SRS


rho.SRS2 %>% 
  ggplot() +
  geom_density(aes(target.var))

# expected value of rho is dependent on the difference in inlcusion probabilities (selection bias) and the size of range of y 


# the problem difficulty always needs to be the true value of the sd and not the one from the biased sample. this means one needs auxiliary info to correctly estimate the sd in order to correct for the bias
