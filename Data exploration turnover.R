library(tidyverse)
library(broom)
library(modelr)


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
  geom_density(aes(log10(omzet_excl_btw), 
                   color = as.factor(selective.sam))) +
  geom_density(aes(log10(omzet_excl_btw)), color = "black") +
  theme(legend.position="bottom") +
  facet_wrap(~year.quarter) +
  labs(x = "log10(turnover)", y = "density", color = "Selection") +
  ggtitle("Density: Log10 Turnover by Year/Quarter", subtitle = "Grouped by selection, black line population")

# density days
data %>% 
  ggplot() +
  geom_density(aes(daysSinceEOQ)) +
  facet_wrap(~year.quarter)

# density working people
data %>% 
  ggplot() +
  geom_density(aes(log10(wp), 
                   color = as.factor(selective.sam))) +
  geom_density(aes(log10(wp)), color = "black") +
  theme(legend.position="bottom") +
  facet_wrap(~year.quarter) +
  labs(x = "log10(working people)", y = "density", color = "Selection") +
  ggtitle("Density: Log10 Working people by Year/Quarter", subtitle = "Grouped by selection, black line population")


# turnover total per period
e <- data %>% 
  group_by(year.quarter) %>% 
  summarise(tot.turnover = sum(omzet_excl_btw))


## Modelling wp and log10(turnover)
quarter_model <- function(df) {
  lm(log.turnover ~ wp, data = df)
}

by.quarter <- data %>% 
  group_by(year.quarter) %>% 
  nest() %>% 
  mutate(model = map(data, quarter_model)) %>% 
  mutate(resids = map2(data, model, add_residuals))

resids <- unnest(by.quarter, resids)
  
resids %>% 
  ggplot(aes(sample = resid)) +
  stat_qq() +
  stat_qq_line() +
  facet_wrap(~year.quarter) +
  ggtitle("Normal Q-Q")
