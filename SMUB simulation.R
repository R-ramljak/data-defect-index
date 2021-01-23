library(tidyverse)
library(svglite)


r <- seq(0, 1, 0.1)
phi <- seq(0, 1, 0.1)

f <- crossing(r, phi)

g <- function(phi, r) {
  num <- phi + (1 - phi) * r
  denom <- phi * r + (1 - phi)
  
  final <- num / denom
  return(final)
}

e <- f %>% 
  mutate(num = phi + (1 - phi) * r) %>% 
  mutate(denom = phi * r + (1 - phi)) %>% 
  mutate(new = num / denom) 


e %>% 
  ggplot() +
  geom_point(aes(x = phi, y = r)) 

e %>% 
  ggplot() +
  geom_point(aes(x = r, y = new)) +
  facet_wrap(~phi)

plot <- e %>% 
  filter(!is.infinite(new)) %>% 
  ggplot() +
  geom_point(aes(x = phi, y = new)) +
  geom_line(aes(x = phi, y  = new)) +
  facet_wrap(~r, labeller = label_bquote(r[ZY]^(1) ~"="~  .(r))) +
  labs(x = expression(varphi), y = "Weight adjustment factor"
       #title = "Range of Weight adjustment factor given r and phi"
       )

ggsave("Plots/weight.adjustment.range.svg", plot, device = "svg")


N.pop <- 10000

pop <- tibble(ID = paste0("A", 1:N.pop), 
              y = rnorm(N.pop, mean = 40, 10)) %>% 
  mutate(gen.inclusion.prob = 1 / (length(ID))) %>% 
  arrange(y) %>% 
  mutate(rank = percent_rank(y)) %>% 
  mutate(prob.help = case_when(between(percent_rank(y), 0, 0.2) ~ 1, #| # min(y) < y < 1
                               # between(percent_rank(y), 0.9, 1) ~ 1, # 10 < y < max(y)
                               TRUE ~ 0)) %>% 
  mutate(adj.inclusion.prob = case_when(prob.help == 1 ~ 1 / sum(prob.help),
                                        prob.help == 0 ~ 1 / sum(length(ID) - prob.help)))


f <- sample(pop$ID, 1000, replace = F, prob = pop$adj.inclusion.prob) %>%
  as_tibble() %>% 
  select(ID = value) %>% 
  mutate(selected = 1) %>% 
  right_join(pop, by = "ID") %>% 
  mutate(selected = case_when(is.na(selected) ~ 0,
                              TRUE ~ selected)) %>% 
  group_by(selected) %>% 
  mutate(number = n()) %>% 
  mutate(variance.x = n() / N.pop * (1 - n() / N.pop)) %>% 
  mutate(variance.y = var(y)) %>% 
  ungroup() %>% 
  mutate(x.star = selected * sqrt(variance.y / variance.x))

