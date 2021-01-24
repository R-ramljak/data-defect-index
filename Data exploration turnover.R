library(tidyverse)


raw <- readRDS("Data/sample.rds")


data <- raw %>% 
  filter(!omzet_excl_btw <= 0 ,
         !wp <= 0) %>% 
  mutate(ID = row_number()) %>% 
  select(ID, everything()) %>% 
  mutate(wp = as.numeric(wp)) %>% 
  mutate(year.quarter = paste0("Y:", Year, " ", "Q:0", Quarter)) %>% 
  arrange(desc(year.quarter, daysSinceEOQ)) %>%
  group_by(year.quarter) %>% 
  mutate(day.rank = 1 / n()) %>% 
  mutate(day.cum = cumsum(day.rank)) %>% 
  mutate(selective.sam = case_when(day.cum <= 0.1 ~ 1,
                                   TRUE ~ 0)) %>% 
  mutate(log.turnover = log10(omzet_excl_btw)) %>% 
  ungroup()


## Descriptives

# n per quarter and year
data %>% 
  group_by(Year, Quarter) %>% 
  summarise(count = n())

# n per NACE code
data %>% 
  group_by(kerncel) %>% 
  summarise(count = n())

# density turnover
data %>% 
  ggplot() +
  geom_density(aes(log10(omzet_excl_btw), group = selective.sam, color = selective.sam)) +
  facet_wrap(~year.quarter)

# density days
data %>% 
  ggplot() +
  geom_density(aes(daysSinceEOQ)) +
  facet_wrap(~year.quarter)

# density turnover
data %>% 
  ggplot() +
  geom_density(aes(log10(wp))) +
  facet_wrap(~year.quarter)

# turnover total per period
e <- data %>% 
  group_by(year.quarter) %>% 
  summarise(tot.turnover = sum(omzet_excl_btw))

# correlation plot
data %>% 
  ggplot() +
  geom_point(aes(x = log10(wp), y = log.turnover)) +
  facet_wrap(~year.quarter)
  


## Modelling
mod1 <- lm(log.turnover ~ wp, data = data)
summary(mod1)

y.name = "Year"
z.name = "Quarter"


#### SMUB

## input is a data frame with 3 variables: y, z, sample. sample is a dummy, indicating if the unit belongs to the sample or not
## This mean the whole data frame without missings will be implemented into the function
## phi is another input which is limited between 0 and 1 and can be a constant or a numeric vector indicating the values of interest
## the output will be tibble with as many rows as the phi vector contains and the following variables:
## phi, r, sdm (standardized mean difference in terms of z), waf (weight adjustment factor), SMUB.estimate 

SMUB <- function(data, phi) {
  
  ## helpers
  
  # datasets
  sample.data <- data %>% 
    filter(sample == 1)
  
  population.data <- data
  
  # sample descriptives
  
  # rescale z
  
  sample.mean.y <- mean(sample.data$y)
  sample.var.y <- var(sample.data$y)
  
  sample.mean.z <- mean(sample.data$z)
  sample.var.z <- var(sample.data$z)
  
  r <- cor(sample.data$y, sample.data$z)
  
  # population descriptives
  
  pop.mean.y <- mean(population.data$y)
  pop.sd.y <- sd(population.data$y)
  
  
  pop.mean.z <- mean(population.data$z)
  
  rho <- cor(population.data$y, as.numeric(population.data$sample))
  
  # SMUB formula helpers
  sdm <- ( mean(sample.data$z) - mean(population.data$z) ) / sd(sample.data$z)
  
  
  ## calculation of output
  
  # tibble creation with final SMUB.estimate
  sens.data <- tibble(phi) %>% 
    mutate(r = r, sdm = sdm) %>% 
    mutate(num = phi + (1 - phi) * r, denom = phi * r + (1 - phi)) %>% 
    mutate(waf = num / denom) %>% 
    mutate(SMUB.estimate = waf * sdm) %>% 
    mutate(estimated.mean = sample.mean.y + waf * sqrt(sample.var.y / sample.var.z) * (pop.mean.z - sample.mean.z))
  
  # true bias and true phi
  bias <- sample.mean.y - pop.mean.y
  std.bias <- (sample.mean.y - pop.mean.y) / pop.sd.y
  var.yz <- sqrt(sample.var.y / sample.var.z)
  mean.diff.z <- sample.mean.z - pop.mean.z

  true.waf <- bias / (var.yz * mean.diff.z)
  true.phi <- (bias - r * var.yz * mean.diff.z) / ( (1 - r) * (bias + var.yz * mean.diff.z) )
  true.estimated.mean <- sample.mean.y + true.waf * var.yz * (pop.mean.z - sample.mean.z) # different than mean.diff.z
  
  
  ## output organization
  # Phi range tibble
  estimate <- sens.data %>% 
    select(-c(num, denom))
  
  # Global Statistics
  global.characteristics <- list(n.prop = sum(as.numeric(population.data$sample)) / length(population.data$sample),
                                 sample.mean.y = sample.mean.y,
                                 pop.mean.y = pop.mean.y,
                                 bias = bias,
                                 std.bias = std.bias,
                                 r = r,
                                 true.waf = true.waf,
                                 true.phi = true.phi,
                                 rho = rho,
                                 sample.mean.z = sample.mean.z,
                                 sample.var.y = sample.var.y,
                                 sample.var.z = sample.var.z,
                                 var.yz = var.yz,
                                 mean.diff.z = mean.diff.z,
                                 pop.mean.z = pop.mean.z,
                                 pop.sd.y = pop.sd.y,
                                 true.estimated.mean1 = true.estimated.mean)

  return(list(sensitivity.analysis = estimate,
              global.characteristics = global.characteristics))
}


selected.data <- data %>% 
  select(y = log.turnover, z = wp, sample = selective.sam, year.quarter) %>%
  filter(sample == 1) %>% 
  group_by(year.quarter) %>% 
  summarise(std.var.sel.yz = sqrt( var(y) * var(z) ), .groups = "drop") %>% 
  group_split(year.quarter) %>% 
  map_df(~select(., std.var.sel.yz)) %>% 
  deframe()

data.final <- data %>% 
  split(.$year.quarter) %>%
  map(~select(., y = log.turnover, z = wp, sample = selective.sam)) %>% 
  map2(., selected.data, ~mutate(.x, z = z * .y))
  
smub1 <- data.final %>% 
  map(~SMUB(data = ., phi = seq(0, 5, 0.1)))

# Comp. table global stats all quarters
smub.ind.global <- smub1 %>% 
  map_depth(1, ~as_tibble(.$global.characteristics)) %>% 
  map_dfr(~select(., n.prop, bias, std.bias, mean.diff.z, r, true.phi, mean.diff.z), .id = "year.quarter")




# e.helper <- data %>%
#   filter(year.quarter == "Y:2015 Q:04") %>%
#   
# 
# e <- data %>% 
#   filter(year.quarter == "Y:2015 Q:04") %>% 
#   select(y = log.turnover, z = wp, sample = selective.sam) %>% 
#   mutate(z = z * sqrt( var(e.helper$y) / var(e.helper$z)))
# 
# d1 <- SMUB(data = e, phi = seq(0, 5, 0.1))


