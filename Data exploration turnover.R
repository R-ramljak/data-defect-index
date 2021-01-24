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


