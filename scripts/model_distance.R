#' ##############################
#' load libraries and set seed
#' ##############################
library(tidyverse)
library(gt)
library(tidybayes)
library(brms)

set.seed(16)


#' ####################################################
#' read data
#' ####################################################
read_csv("./data/dat_long.csv") %>%
  identity() -> dat_long
dat_long


dat_long %>%
  group_by(subject_room_id) %>%
  filter(subject_room_day == min(subject_room_day)) %>%
  ungroup() %>%
  identity() -> dat_first
dat_first


#' ####################################################
#' generative model for SARS-CoV-2 contamination
#' - initial sample only
#' - elevated vs floor (no accounting for high-touch)
#' - binomial probability of detection
#' ####################################################

dat_first %>%
  # exclude bathroom sites and nursing sites (different distance scheme)
  filter(grepl("bathroom|sink|toilet|nurse",site_descriptor) == FALSE) %>%
  mutate(scv2_detected = !is.na(copies_max),
         site_category = case_when(grepl("floor",site_descriptor) == TRUE ~ "floor",
                                   grepl("floor",site_descriptor) == FALSE ~ "elevated"),
         high_touch = touch == "High") %>%
  #count(site_descriptor, site_category, high_touch) %>% gt()
  select(redcap_id, subject_covid_day, subject_hosp_day, subject_covid_day, swab_site, unit, site_category, touch, high_touch, distance, copies_max, scv2_detected) %>%
  mutate(site_category = stringr::str_to_title(gsub("_"," ",site_category))) %>%
  distinct() %>%
  identity() -> dat

dat


#' get prior
dat %>%
  brms::get_prior(data = ., family = bernoulli,
                  scv2_detected ~ (1 + distance | site_category)
  ) %>%
  gt::gt()


#' run binomial model
# dat %>%
#   select(scv2_detected, distance, site_category) %>%
#   brm(formula = scv2_detected ~ 0 + (1 + distance | site_category),
#       data = .,
#       family = bernoulli,
#       chains = 4,
#       cores = 4,
#       control = list("adapt_delta" = 0.9999, max_treedepth = 10),
#       backend = "cmdstanr",
#       seed = 16) -> m_binom_scv2_distance_mix_category
# 
# m_binom_scv2_distance_mix_category %>% write_rds(file = "./models/binomial/m_binom_scv2_distance_mix_category.rds.gz", compress = "gz")
m_binom_scv2_distance_mix_category <- read_rds(file = "./models/binomial/m_binom_scv2_distance_mix_category.rds.gz")

m_binom_scv2_distance_mix_category
rstan::check_hmc_diagnostics(m_binom_scv2_distance_mix_category$fit)
m_binom_scv2_distance_mix_category %>% pp_check()

m_binom_scv2_distance_mix_category %>%
  posterior_summary() %>%
  as_tibble(rownames = "param") %>%
  gt::gt() %>%
  gt::fmt_number(columns = 2:5, n_sigfig = 3)


#' fitted
m_binom_scv2_distance_mix_category$data %>%
  as_tibble() %>%
  expand(distance = modelr::seq_range(distance, n = 100),
         site_category = unique(site_category)#,
         #high_touch = unique(high_touch)
         ) %>%
  #filter(!(high_touch == TRUE & site_category == "Floor")) %>%
  add_epred_draws(m_binom_scv2_distance_mix_category) %>%
  identity() -> m_binom_scv2_distance_mix_category_fitted
m_binom_scv2_distance_mix_category_fitted


m_binom_scv2_distance_mix_category_fitted %>%
  #mutate(high_touch = case_when(high_touch == TRUE ~ "High Touch",
  #                              high_touch == FALSE ~ "Low Touch")) %>%
  ggplot(aes(x = distance, y = .epred)) +
  #geom_point(data = m_binom_scv2_distance_mix_category$data, aes(x = distance, y = scv2_detected), color = "grey", alpha = 0.5) +
  stat_lineribbon() +
  #facet_wrap(facets = ~ site_category + high_touch, scales = "fixed") +
  facet_wrap(facets = ~ site_category) +
  scale_fill_brewer(palette = "Reds") +
  labs(x = "Distance from Patient (meters)",
       y = "Probability of SARS-CoV-2 Detection by RT-PCR",
       fill = "Posterior Credible Interval") +
  theme_bw() +
  theme(strip.text = ggtext::element_markdown(color = "black", size = 8),
        axis.text.x = ggtext::element_markdown(color = "black"),
        axis.text.y = ggtext::element_markdown(color = "black"),
        legend.position = "top",
        #legend.background = element_rect(fill = "white", color = "black", size = 0.25),
        strip.background = element_blank()) -> p_scv2_binomial_distance_site_category
p_scv2_binomial_distance_site_category


# p_scv2_binomial_distance_site_category %>%
#   ggsave(filename = "./figs/p_scv2_binomial_distance_site_category.pdf", height = 5, width = 8, units = "in")
# p_scv2_binomial_distance_site_category %>%
#   ggsave(filename = "./figs/p_scv2_binomial_distance_site_category.png", height = 5, width = 8, units = "in", dpi = 600)
# p_scv2_binomial_distance_site_category %>%
#   ggsave(filename = "./figs/p_scv2_binomial_distance_site_category.svg", height = 5, width = 8, units = "in")





#' ####################################################
#' generative model for SARS-CoV-2 contamination
#' - initial sample only
#' - elevated vs floor -- accounting for high-touch
#' - binomial probability of detection
#' ####################################################

#' get prior
dat %>%
  brms::get_prior(data = ., family = bernoulli,
                  scv2_detected ~ (1 + distance | site_category / high_touch)
  ) %>%
  gt::gt()


#' run binomial model
# dat %>%
#   select(scv2_detected, distance, site_category, high_touch) %>%
#   brm(formula = scv2_detected ~ 0 + (1 + distance | site_category / high_touch),
#       data = .,
#       family = bernoulli,
#       chains = 4,
#       cores = 4,
#       control = list("adapt_delta" = 0.999, max_treedepth = 10),
#       backend = "cmdstanr",
#       seed = 16) -> m_binom_scv2_distance_mix_category_touch
# 
# m_binom_scv2_distance_mix_category_touch %>% write_rds(file = "./models/binomial/m_binom_scv2_distance_mix_category_touch.rds.gz", compress = "gz")
m_binom_scv2_distance_mix_category_touch <- read_rds(file = "./models/binomial/m_binom_scv2_distance_mix_category_touch.rds.gz")

m_binom_scv2_distance_mix_category_touch
rstan::check_hmc_diagnostics(m_binom_scv2_distance_mix_category_touch$fit)
m_binom_scv2_distance_mix_category_touch %>% pp_check()

m_binom_scv2_distance_mix_category_touch %>%
  posterior_summary() %>%
  as_tibble(rownames = "param") %>%
  gt::gt() %>%
  gt::fmt_number(columns = 2:5, n_sigfig = 3)


#' fitted
m_binom_scv2_distance_mix_category_touch$data %>%
  as_tibble() %>%
  expand(distance = modelr::seq_range(distance, n = 100),
         site_category = unique(site_category),
         high_touch = unique(high_touch)
  ) %>%
  filter(!(high_touch == TRUE & site_category == "Floor")) %>%
  add_epred_draws(m_binom_scv2_distance_mix_category_touch) %>%
  identity() -> m_binom_scv2_distance_mix_category_touch_fitted
m_binom_scv2_distance_mix_category_touch_fitted


m_binom_scv2_distance_mix_category_touch_fitted %>%
  mutate(high_touch = case_when(high_touch == TRUE ~ "High Touch",
                                high_touch == FALSE ~ "Low Touch")) %>%
  ggplot(aes(x = distance, y = .epred)) +
  #geom_point(data = m_binom_scv2_distance_mix_category_touch$data, aes(x = distance, y = scv2_detected), color = "grey", alpha = 0.5) +
  stat_lineribbon() +
  facet_wrap(facets = ~ site_category + high_touch, scales = "fixed") +
  #facet_wrap(facets = ~ site_category) +
  scale_fill_brewer(palette = "Reds") +
  labs(x = "Distance from Patient (meters)",
       y = "Probability of SARS-CoV-2 Detection by RT-PCR",
       fill = "Posterior Credible Interval") +
  theme_bw() +
  theme(strip.text = ggtext::element_markdown(color = "black", size = 8),
        axis.text.x = ggtext::element_markdown(color = "black"),
        axis.text.y = ggtext::element_markdown(color = "black"),
        legend.position = "top",
        #legend.background = element_rect(fill = "white", color = "black", size = 0.25),
        strip.background = element_blank()) -> p_scv2_binomial_distance_site_category_touch
p_scv2_binomial_distance_site_category_touch


# p_scv2_binomial_distance_site_category_touch %>%
#   ggsave(filename = "./figs/p_scv2_binomial_distance_site_category_touch.pdf", height = 5, width = 8, units = "in")
# p_scv2_binomial_distance_site_category_touch %>%
#   ggsave(filename = "./figs/p_scv2_binomial_distance_site_category_touch.png", height = 5, width = 8, units = "in", dpi = 600)
# p_scv2_binomial_distance_site_category_touch %>%
#   ggsave(filename = "./figs/p_scv2_binomial_distance_site_category_touch.svg", height = 5, width = 8, units = "in")







#' ####################################################
#' generative model for SARS-CoV-2 contamination
#' - initial sample only
#' - elevated vs floor (no accounting for high-touch)
#' - linear quantity detected
#' ####################################################

dat_first %>%
  # exclude bathroom sites and nursing sites (different distance scheme)
  filter(grepl("bathroom|sink|toilet|nurse",site_descriptor) == FALSE) %>%
  mutate(scv2_detected = !is.na(copies_max),
         site_category = case_when(grepl("floor",site_descriptor) == TRUE ~ "floor",
                                   grepl("floor",site_descriptor) == FALSE ~ "elevated"),
         high_touch = touch == "High") %>%
  #count(site_descriptor, site_category, high_touch)
  select(redcap_id, subject_covid_day, subject_hosp_day, subject_covid_day, swab_site, unit, site_category, touch, high_touch, distance, copies_max, scv2_detected) %>%
  mutate(site_category = stringr::str_to_title(gsub("_"," ",site_category))) %>%
  distinct() %>%
  identity() -> dat

dat



#' get prior
dat %>%
  brms::get_prior(data = ., family = bernoulli,
                  log10(copies_max) ~ (1 + distance | site_category)
  ) %>%
  gt::gt()


#' run linear model: does NOT account for specimens with unmeasured copy numbers
# dat %>%
#   select(copies_max, distance, site_category) %>%
#   brm(formula = log10(copies_max) ~ 0 + (1 + distance | site_category),
#       data = .,
#       family = gaussian,
#       chains = 4,
#       cores = 4,
#       control = list("adapt_delta" = 0.9999, max_treedepth = 10),
#       backend = "cmdstanr",
#       seed = 16) -> m_linear_scv2_distance_mix_category
# 
# m_linear_scv2_distance_mix_category %>% write_rds(file = "./models/linear/m_linear_scv2_distance_mix_category.rds.gz", compress = "gz")
m_linear_scv2_distance_mix_category <- read_rds(file = "./models/linear/m_linear_scv2_distance_mix_category.rds.gz")

m_linear_scv2_distance_mix_category
rstan::check_hmc_diagnostics(m_linear_scv2_distance_mix_category$fit)
m_linear_scv2_distance_mix_category %>% pp_check()

m_linear_scv2_distance_mix_category %>%
  posterior_summary() %>%
  as_tibble(rownames = "param") %>%
  gt::gt() %>%
  gt::fmt_number(columns = 2:5, n_sigfig = 3)


#' fitted
m_linear_scv2_distance_mix_category$data %>%
  as_tibble() %>%
  expand(distance = modelr::seq_range(distance, n = 100),
         site_category = unique(site_category)#,
         #high_touch = unique(high_touch)
  ) %>%
  #filter(!(high_touch == TRUE & site_category == "Floor")) %>%
  add_epred_draws(m_linear_scv2_distance_mix_category) %>%
  identity() -> m_linear_scv2_distance_mix_category_fitted
m_linear_scv2_distance_mix_category_fitted


m_linear_scv2_distance_mix_category_fitted %>%
  #mutate(high_touch = case_when(high_touch == TRUE ~ "High Touch",
  #                              high_touch == FALSE ~ "Low Touch")) %>%
  ggplot(aes(x = distance, y = .epred)) +
  stat_lineribbon() +
  geom_point(data = m_linear_scv2_distance_mix_category$data, aes(x = distance, y = log10(copies_max)), color = "grey", alpha = 0.5) +
  #facet_wrap(facets = ~ site_category + high_touch, scales = "fixed") +
  facet_wrap(facets = ~ site_category) +
  scale_fill_brewer(palette = "Reds") +
  labs(x = "Distance from Patient (meters)",
       y = "log<sub>10</sub> copies SARS-CoV-2 by RT-PCR",
       fill = "Posterior Credible Interval") +
  theme_bw() +
  theme(strip.text = ggtext::element_markdown(color = "black", size = 8),
        axis.text.x = ggtext::element_markdown(color = "black"),
        axis.text.y = ggtext::element_markdown(color = "black"),
        axis.title.x = ggtext::element_markdown(color = "black"),
        axis.title.y = ggtext::element_markdown(color = "black"),
        legend.position = "top",
        #legend.background = element_rect(fill = "white", color = "black", size = 0.25),
        strip.background = element_blank()) -> p_scv2_linear_distance_site_category
p_scv2_linear_distance_site_category


# p_scv2_linear_distance_site_category %>%
#   ggsave(filename = "./figs/p_scv2_linear_distance_site_category.pdf", height = 5, width = 8, units = "in")
# p_scv2_linear_distance_site_category %>%
#   ggsave(filename = "./figs/p_scv2_linear_distance_site_category.png", height = 5, width = 8, units = "in", dpi = 600)
# p_scv2_linear_distance_site_category %>%
#   ggsave(filename = "./figs/p_scv2_linear_distance_site_category.svg", height = 5, width = 8, units = "in")





#' ####################################################
#' generative model for SARS-CoV-2 contamination
#' - initial sample only
#' - elevated vs floor -- accounting for high-touch
#' - linear quantity detected
#' ####################################################

dat_first %>%
  # exclude bathroom sites and nursing sites (different distance scheme)
  filter(grepl("bathroom|sink|toilet|nurse",site_descriptor) == FALSE) %>%
  mutate(scv2_detected = !is.na(copies_max),
         site_category = case_when(grepl("floor",site_descriptor) == TRUE ~ "floor",
                                   grepl("floor",site_descriptor) == FALSE ~ "elevated"),
         high_touch = touch == "High") %>%
  #count(site_descriptor, site_category, high_touch)
  select(redcap_id, subject_covid_day, subject_hosp_day, subject_covid_day, swab_site, unit, site_category, touch, high_touch, distance, copies_max, scv2_detected) %>%
  mutate(site_category = stringr::str_to_title(gsub("_"," ",site_category))) %>%
  distinct() %>%
  identity() -> dat

dat



#' get prior
dat %>%
  brms::get_prior(data = ., family = bernoulli,
                  log10(copies_max) ~ (1 + distance | site_category / high_touch)
  ) %>%
  gt::gt()


#' run linear model
# dat %>%
#   select(copies_max, distance, site_category, high_touch) %>%
#   brm(formula = log10(copies_max) ~ 0 + (1 + distance | site_category / high_touch),
#       data = .,
#       family = gaussian,
#       chains = 4,
#       cores = 4,
#       control = list("adapt_delta" = 0.999999, max_treedepth = 10),
#       backend = "cmdstanr",
#       seed = 16) -> m_linear_scv2_distance_mix_category_touch
# 
# m_linear_scv2_distance_mix_category_touch %>% write_rds(file = "./models/linear/m_linear_scv2_distance_mix_category_touch.rds.gz", compress = "gz")
m_linear_scv2_distance_mix_category_touch <- read_rds(file = "./models/linear/m_linear_scv2_distance_mix_category_touch.rds.gz")

m_linear_scv2_distance_mix_category_touch
rstan::check_hmc_diagnostics(m_linear_scv2_distance_mix_category_touch$fit)
m_linear_scv2_distance_mix_category_touch %>% pp_check()

m_linear_scv2_distance_mix_category_touch %>%
  posterior_summary() %>%
  as_tibble(rownames = "param") %>%
  gt::gt() %>%
  gt::fmt_number(columns = 2:5, n_sigfig = 3)


#' fitted
m_linear_scv2_distance_mix_category_touch$data %>%
  as_tibble() %>%
  expand(distance = modelr::seq_range(distance, n = 100),
         site_category = unique(site_category),
         high_touch = unique(high_touch)
  ) %>%
  filter(!(high_touch == TRUE & site_category == "Floor")) %>%
  add_epred_draws(m_linear_scv2_distance_mix_category_touch) %>%
  identity() -> m_linear_scv2_distance_mix_category_touch_fitted
m_linear_scv2_distance_mix_category_touch_fitted


m_linear_scv2_distance_mix_category_touch_fitted %>%
  mutate(high_touch = case_when(high_touch == TRUE ~ "High Touch",
                                high_touch == FALSE ~ "Low Touch")) %>%
  ggplot(aes(x = distance, y = .epred)) +
  stat_lineribbon() +
  geom_point(data = mutate(as_tibble(m_linear_scv2_distance_mix_category_touch$data),
                           high_touch = case_when(high_touch == TRUE ~ "High Touch",
                                                  high_touch == FALSE ~ "Low Touch")
  ),
  aes(x = distance, y = log10(copies_max)), color = "grey", alpha = 0.5) +
  facet_wrap(facets = ~ site_category + high_touch, scales = "fixed") +
  #facet_wrap(facets = ~ site_category) +
  scale_fill_brewer(palette = "Reds") +
  labs(x = "Distance from Patient (meters)",
       y = "log<sub>10</sub> copies SARS-CoV-2 by RT-PCR",
       fill = "Posterior Credible Interval") +
  theme_bw() +
  theme(strip.text = ggtext::element_markdown(color = "black", size = 8),
        axis.text.x = ggtext::element_markdown(color = "black"),
        axis.text.y = ggtext::element_markdown(color = "black"),
        axis.title.x = ggtext::element_markdown(color = "black"),
        axis.title.y = ggtext::element_markdown(color = "black"),
        legend.position = "top",
        #legend.background = element_rect(fill = "white", color = "black", size = 0.25),
        strip.background = element_blank()) -> p_scv2_linear_distance_site_category_touch
p_scv2_linear_distance_site_category_touch


# p_scv2_linear_distance_site_category_touch %>%
#   ggsave(filename = "./figs/p_scv2_linear_distance_site_category_touch.pdf", height = 5, width = 8, units = "in")
# p_scv2_linear_distance_site_category_touch %>%
#   ggsave(filename = "./figs/p_scv2_linear_distance_site_category_touch.png", height = 5, width = 8, units = "in", dpi = 600)
# p_scv2_linear_distance_site_category_touch %>%
#   ggsave(filename = "./figs/p_scv2_linear_distance_site_category_touch.svg", height = 5, width = 8, units = "in")




#' #############################################
#' HURDLE MODELS
#' #############################################

#' ####################################################
#' generative model for SARS-CoV-2 contamination
#' - initial sample only
#' - elevated vs floor (no accounting for high-touch)
#' - hurdle quantity detected
#' ####################################################

dat_first %>%
  # exclude bathroom sites and nursing sites (different distance scheme)
  filter(grepl("bathroom|sink|toilet|nurse",site_descriptor) == FALSE) %>%
  mutate(scv2_detected = !is.na(copies_max),
         site_category = case_when(grepl("floor",site_descriptor) == TRUE ~ "floor",
                                   grepl("floor",site_descriptor) == FALSE ~ "elevated"),
         high_touch = touch == "High") %>%
  #count(site_descriptor, site_category, high_touch)
  select(redcap_id, subject_covid_day, subject_hosp_day, subject_covid_day, swab_site, unit, site_category, touch, high_touch, distance, copies_max, scv2_detected) %>%
  mutate(site_category = stringr::str_to_title(gsub("_"," ",site_category))) %>%
  distinct() %>%
  identity() -> dat

dat



#' get prior
dat %>%
  mutate(copies_max = replace(copies_max, is.na(copies_max), 0)) %>%
  brms::get_prior(data = ., family = hurdle_lognormal(),
                  copies_max ~ (1 + distance | site_category)
  ) %>%
  gt::gt()


#' run hurdle model: does NOT account for specimens with unmeasured copy numbers
# dat %>%
#   mutate(copies_max = replace(copies_max, is.na(copies_max), 0)) %>%
#   select(copies_max, distance, site_category) %>%
#   brm(formula = copies_max ~ 1 + (1 + distance | site_category),
#       data = .,
#       family = hurdle_lognormal(),
#       chains = 4,
#       cores = 4,
#       control = list("adapt_delta" = 0.9999, max_treedepth = 10),
#       backend = "cmdstanr",
#       seed = 16) -> m_hurdle_scv2_distance_mix_category
# 
# m_hurdle_scv2_distance_mix_category %>% write_rds(file = "./models/hurdle/m_hurdle_scv2_distance_mix_category.rds.gz", compress = "gz")
m_hurdle_scv2_distance_mix_category <- read_rds(file = "./models/hurdle/m_hurdle_scv2_distance_mix_category.rds.gz")

m_hurdle_scv2_distance_mix_category
rstan::check_hmc_diagnostics(m_hurdle_scv2_distance_mix_category$fit)
m_hurdle_scv2_distance_mix_category %>% pp_check() + scale_x_log10()

m_hurdle_scv2_distance_mix_category %>%
  posterior_summary() %>%
  as_tibble(rownames = "param") %>%
  gt::gt() %>%
  gt::fmt_number(columns = 2:5, n_sigfig = 3)


#' fitted
m_hurdle_scv2_distance_mix_category$data %>%
  as_tibble() %>%
  expand(distance = modelr::seq_range(distance, n = 100),
         site_category = unique(site_category)#,
         #high_touch = unique(high_touch)
  ) %>%
  #filter(!(high_touch == TRUE & site_category == "Floor")) %>%
  add_epred_draws(m_hurdle_scv2_distance_mix_category) %>%
  identity() -> m_hurdle_scv2_distance_mix_category_fitted
m_hurdle_scv2_distance_mix_category_fitted


m_hurdle_scv2_distance_mix_category_fitted %>%
  #mutate(high_touch = case_when(high_touch == TRUE ~ "High Touch",
  #                              high_touch == FALSE ~ "Low Touch")) %>%
  ggplot(aes(x = distance, y = .epred)) +
  stat_lineribbon() +
  geom_point(data = m_hurdle_scv2_distance_mix_category$data, aes(x = distance, y = copies_max), color = "grey", alpha = 0.5) +
  #facet_wrap(facets = ~ site_category + high_touch, scales = "fixed") +
  facet_wrap(facets = ~ site_category) +
  scale_y_log10() +
  scale_fill_brewer(palette = "Reds") +
  labs(x = "Distance from Patient (meters)",
       y = "Copies SARS-CoV-2 by RT-PCR",
       fill = "Posterior Credible Interval") +
  theme_bw() +
  theme(strip.text = ggtext::element_markdown(color = "black", size = 8),
        axis.text.x = ggtext::element_markdown(color = "black"),
        axis.text.y = ggtext::element_markdown(color = "black"),
        axis.title.x = ggtext::element_markdown(color = "black"),
        axis.title.y = ggtext::element_markdown(color = "black"),
        legend.position = "top",
        #legend.background = element_rect(fill = "white", color = "black", size = 0.25),
        strip.background = element_blank()) -> p_scv2_hurdle_distance_site_category
p_scv2_hurdle_distance_site_category


# p_scv2_hurdle_distance_site_category %>%
#   ggsave(filename = "./figs/p_scv2_hurdle_distance_site_category.pdf", height = 5, width = 8, units = "in")
# p_scv2_hurdle_distance_site_category %>%
#   ggsave(filename = "./figs/p_scv2_hurdle_distance_site_category.png", height = 5, width = 8, units = "in", dpi = 600)
# p_scv2_hurdle_distance_site_category %>%
#   ggsave(filename = "./figs/p_scv2_hurdle_distance_site_category.svg", height = 5, width = 8, units = "in")





#' ####################################################
#' generative model for SARS-CoV-2 contamination
#' - initial sample only
#' - elevated vs floor -- accounting for high-touch
#' - hurdle quantity detected
#' ####################################################

dat_first %>%
  # exclude bathroom sites and nursing sites (different distance scheme)
  filter(grepl("bathroom|sink|toilet|nurse",site_descriptor) == FALSE) %>%
  mutate(scv2_detected = !is.na(copies_max),
         site_category = case_when(grepl("floor",site_descriptor) == TRUE ~ "floor",
                                   grepl("floor",site_descriptor) == FALSE ~ "elevated"),
         high_touch = touch == "High") %>%
  #count(site_descriptor, site_category, high_touch)
  select(redcap_id, subject_covid_day, subject_hosp_day, subject_covid_day, swab_site, unit, site_category, touch, high_touch, distance, copies_max, scv2_detected) %>%
  mutate(site_category = stringr::str_to_title(gsub("_"," ",site_category))) %>%
  distinct() %>%
  identity() -> dat

dat



#' get prior
dat %>%
  brms::get_prior(data = ., family = bernoulli,
                  log10(copies_max) ~ (1 + distance | site_category / high_touch)
  ) %>%
  gt::gt()


#' run hurdle model
# dat %>%
#   mutate(copies_max = replace(copies_max, is.na(copies_max), 0)) %>%
#   select(copies_max, distance, site_category, high_touch) %>%
#   brm(formula = copies_max ~ 1 + (1 + distance | site_category / high_touch),
#       data = .,
#       family = hurdle_lognormal(),
#       chains = 4,
#       cores = 4,
#       control = list("adapt_delta" = 0.99999999, max_treedepth = 10),
#       backend = "cmdstanr",
#       seed = 16) -> m_hurdle_scv2_distance_mix_category_touch
# 
# m_hurdle_scv2_distance_mix_category_touch %>% write_rds(file = "./models/hurdle/m_hurdle_scv2_distance_mix_category_touch.rds.gz", compress = "gz")
m_hurdle_scv2_distance_mix_category_touch <- read_rds(file = "./models/hurdle/m_hurdle_scv2_distance_mix_category_touch.rds.gz")

m_hurdle_scv2_distance_mix_category_touch
rstan::check_hmc_diagnostics(m_hurdle_scv2_distance_mix_category_touch$fit)
m_hurdle_scv2_distance_mix_category_touch %>% pp_check() + scale_x_log10()

m_hurdle_scv2_distance_mix_category_touch %>%
  posterior_summary() %>%
  as_tibble(rownames = "param") %>%
  gt::gt() %>%
  gt::fmt_number(columns = 2:5, n_sigfig = 3)


#' fitted
m_hurdle_scv2_distance_mix_category_touch$data %>%
  as_tibble() %>%
  expand(distance = modelr::seq_range(distance, n = 100),
         site_category = unique(site_category),
         high_touch = unique(high_touch)
  ) %>%
  filter(!(high_touch == TRUE & site_category == "Floor")) %>%
  add_epred_draws(m_hurdle_scv2_distance_mix_category_touch) %>%
  identity() -> m_hurdle_scv2_distance_mix_category_touch_fitted
m_hurdle_scv2_distance_mix_category_touch_fitted


m_hurdle_scv2_distance_mix_category_touch_fitted %>%
  mutate(high_touch = case_when(high_touch == TRUE ~ "High Touch",
                                high_touch == FALSE ~ "Low Touch")) %>%
  ggplot(aes(x = distance, y = .epred)) +
  stat_lineribbon() +
  geom_point(data = mutate(as_tibble(m_hurdle_scv2_distance_mix_category_touch$data),
                           high_touch = case_when(high_touch == TRUE ~ "High Touch",
                                                  high_touch == FALSE ~ "Low Touch")
  ),
  aes(x = distance, y = copies_max), color = "grey", alpha = 0.5) +
  facet_wrap(facets = ~ site_category + high_touch, scales = "fixed") +
  #facet_wrap(facets = ~ site_category) +
  scale_fill_brewer(palette = "Reds") +
  scale_y_log10() +
  labs(x = "Distance from Patient (meters)",
       y = "Copies SARS-CoV-2 by RT-PCR",
       fill = "Posterior Credible Interval") +
  theme_bw() +
  theme(strip.text = ggtext::element_markdown(color = "black", size = 8),
        axis.text.x = ggtext::element_markdown(color = "black"),
        axis.text.y = ggtext::element_markdown(color = "black"),
        axis.title.x = ggtext::element_markdown(color = "black"),
        axis.title.y = ggtext::element_markdown(color = "black"),
        legend.position = "top",
        #legend.background = element_rect(fill = "white", color = "black", size = 0.25),
        strip.background = element_blank()) -> p_scv2_hurdle_distance_site_category_touch
p_scv2_hurdle_distance_site_category_touch


# p_scv2_hurdle_distance_site_category_touch %>%
#   ggsave(filename = "./figs/p_scv2_hurdle_distance_site_category_touch.pdf", height = 5, width = 8, units = "in")
# p_scv2_hurdle_distance_site_category_touch %>%
#   ggsave(filename = "./figs/p_scv2_hurdle_distance_site_category_touch.png", height = 5, width = 8, units = "in", dpi = 600)
# p_scv2_hurdle_distance_site_category_touch %>%
#   ggsave(filename = "./figs/p_scv2_hurdle_distance_site_category_touch.svg", height = 5, width = 8, units = "in")




#' ############################
#' MULTIVARIABLE MODELS to COMPARE with HEIRARCHICAL MODELS
#' ############################

#' ####################################################
#' generative model for SARS-CoV-2 contamination
#' - initial sample only
#' - elevated vs floor -- accounting for high-touch
#' - binomial probability of detection
#' ####################################################

dat_first %>%
  # exclude bathroom sites and nursing sites (different distance scheme)
  filter(grepl("bathroom|sink|toilet|nurse",site_descriptor) == FALSE) %>%
  mutate(scv2_detected = !is.na(copies_max),
         site_category = case_when(grepl("floor",site_descriptor) == TRUE ~ "floor",
                                   grepl("floor",site_descriptor) == FALSE ~ "elevated"),
         high_touch = touch == "High") %>%
  #count(site_descriptor, site_category, high_touch) %>% gt()
  select(redcap_id, subject_covid_day, subject_hosp_day, subject_covid_day, swab_site, unit, site_category, touch, high_touch, distance, copies_max, scv2_detected) %>%
  mutate(site_category = stringr::str_to_title(gsub("_"," ",site_category))) %>%
  distinct() %>%
  identity() -> dat

dat

#' get prior
dat %>%
  brms::get_prior(data = ., family = bernoulli,
                  scv2_detected ~ 1 + distance + site_category + high_touch
  ) %>%
  gt::gt()


#' run binomial model
# dat %>%
#   select(scv2_detected, distance, site_category, high_touch) %>%
#   brm(formula = scv2_detected ~ 1 + distance + site_category + high_touch,
#       data = .,
#       family = bernoulli,
#       chains = 4,
#       cores = 4,
#       control = list("adapt_delta" = 0.999, max_treedepth = 10),
#       backend = "cmdstanr",
#       seed = 16) -> m_mvbinom_scv2_distance_fix_category_touch
# 
# m_mvbinom_scv2_distance_fix_category_touch %>% write_rds(file = "./models/binomial/m_mvbinom_scv2_distance_fix_category_touch.rds.gz", compress = "gz")
m_mvbinom_scv2_distance_fix_category_touch <- read_rds(file = "./models/binomial/m_mvbinom_scv2_distance_fix_category_touch.rds.gz")

m_mvbinom_scv2_distance_fix_category_touch
rstan::check_hmc_diagnostics(m_mvbinom_scv2_distance_fix_category_touch$fit)
m_mvbinom_scv2_distance_fix_category_touch %>% pp_check()

m_mvbinom_scv2_distance_fix_category_touch %>%
  posterior_summary() %>%
  as_tibble(rownames = "param") %>%
  gt::gt() %>%
  gt::fmt_number(columns = 2:5, n_sigfig = 3)


#' fitted
m_mvbinom_scv2_distance_fix_category_touch$data %>%
  as_tibble() %>%
  expand(distance = modelr::seq_range(distance, n = 100),
         site_category = unique(site_category),
         high_touch = unique(high_touch)
  ) %>%
  filter(!(high_touch == TRUE & site_category == "Floor")) %>%
  add_epred_draws(m_mvbinom_scv2_distance_fix_category_touch) %>%
  identity() -> m_mvbinom_scv2_distance_fix_category_touch_fitted
m_mvbinom_scv2_distance_fix_category_touch_fitted


m_mvbinom_scv2_distance_fix_category_touch_fitted %>%
  mutate(high_touch = case_when(high_touch == TRUE ~ "High Touch",
                                high_touch == FALSE ~ "Low Touch")) %>%
  ggplot(aes(x = distance, y = .epred)) +
  #geom_point(data = m_mvbinom_scv2_distance_fix_category_touch$data, aes(x = distance, y = scv2_detected), color = "grey", alpha = 0.5) +
  stat_lineribbon() +
  facet_wrap(facets = ~ site_category + high_touch, scales = "fixed") +
  #facet_wrap(facets = ~ site_category) +
  scale_fill_brewer(palette = "Reds") +
  labs(x = "Distance from Patient (meters)",
       y = "Probability of SARS-CoV-2 Detection by RT-PCR",
       fill = "Posterior Credible Interval") +
  theme_bw() +
  theme(strip.text = ggtext::element_markdown(color = "black", size = 8),
        axis.text.x = ggtext::element_markdown(color = "black"),
        axis.text.y = ggtext::element_markdown(color = "black"),
        legend.position = "top",
        #legend.background = element_rect(fill = "white", color = "black", size = 0.25),
        strip.background = element_blank()) -> p_scv2_mvbinomial_distance_site_category_touch
p_scv2_mvbinomial_distance_site_category_touch


# p_scv2_mvbinomial_distance_site_category_touch %>%
#   ggsave(filename = "./figs/p_scv2_mvbinomial_distance_site_category_touch.pdf", height = 5, width = 8, units = "in")
# p_scv2_mvbinomial_distance_site_category_touch %>%
#   ggsave(filename = "./figs/p_scv2_mvbinomial_distance_site_category_touch.png", height = 5, width = 8, units = "in", dpi = 600)
# p_scv2_mvbinomial_distance_site_category_touch %>%
#   ggsave(filename = "./figs/p_scv2_mvbinomial_distance_site_category_touch.svg", height = 5, width = 8, units = "in")


#' ######################
#' model comparison
#' ######################
loo::loo_compare(loo(m_mvbinom_scv2_distance_fix_category_touch), loo(m_binom_scv2_distance_mix_category_touch), loo(m_binom_scv2_distance_mix_category)) %>%
  identity() -> binom_distance_category_touch_loo_comparison
binom_distance_category_touch_loo_comparison

binom_distance_category_touch_loo_comparison %>%
  as_tibble(rownames = "model") %>%
  gt()



#' ####################################################
#' generative model for SARS-CoV-2 contamination
#' - initial sample only
#' - elevated vs floor -- accounting for high-touch
#' - hurdle quantity detected
#' ####################################################

dat_first %>%
  # exclude bathroom sites and nursing sites (different distance scheme)
  filter(grepl("bathroom|sink|toilet|nurse",site_descriptor) == FALSE) %>%
  mutate(scv2_detected = !is.na(copies_max),
         site_category = case_when(grepl("floor",site_descriptor) == TRUE ~ "floor",
                                   grepl("floor",site_descriptor) == FALSE ~ "elevated"),
         high_touch = touch == "High") %>%
  #count(site_descriptor, site_category, high_touch)
  select(redcap_id, subject_covid_day, subject_hosp_day, subject_covid_day, swab_site, unit, site_category, touch, high_touch, distance, copies_max, scv2_detected) %>%
  mutate(site_category = stringr::str_to_title(gsub("_"," ",site_category))) %>%
  distinct() %>%
  identity() -> dat

dat



#' get prior
dat %>%
  brms::get_prior(data = ., family = bernoulli,
                  log10(copies_max) ~ 1 + distance + site_category + high_touch
  ) %>%
  gt::gt()


#' run hurdle model
# dat %>%
#   mutate(copies_max = replace(copies_max, is.na(copies_max), 0)) %>%
#   select(copies_max, distance, site_category, high_touch) %>%
#   brm(formula = copies_max ~ 1 + distance + site_category + high_touch,
#       data = .,
#       family = hurdle_lognormal(),
#       chains = 4,
#       cores = 4,
#       control = list("adapt_delta" = 0.999, max_treedepth = 10),
#       backend = "cmdstanr",
#       seed = 16) -> m_mvhurdle_scv2_distance_fix_category_touch
# 
# m_mvhurdle_scv2_distance_fix_category_touch %>% write_rds(file = "./models/hurdle/m_mvhurdle_scv2_distance_fix_category_touch.rds.gz", compress = "gz")
m_mvhurdle_scv2_distance_fix_category_touch <- read_rds(file = "./models/hurdle/m_mvhurdle_scv2_distance_fix_category_touch.rds.gz")

m_mvhurdle_scv2_distance_fix_category_touch
rstan::check_hmc_diagnostics(m_mvhurdle_scv2_distance_fix_category_touch$fit)
m_mvhurdle_scv2_distance_fix_category_touch %>% pp_check() + scale_x_log10()

m_mvhurdle_scv2_distance_fix_category_touch %>%
  posterior_summary() %>%
  as_tibble(rownames = "param") %>%
  gt::gt() %>%
  gt::fmt_number(columns = 2:5, n_sigfig = 3)


#' fitted
m_mvhurdle_scv2_distance_fix_category_touch$data %>%
  as_tibble() %>%
  expand(distance = modelr::seq_range(distance, n = 100),
         site_category = unique(site_category),
         high_touch = unique(high_touch)
  ) %>%
  filter(!(high_touch == TRUE & site_category == "Floor")) %>%
  add_epred_draws(m_mvhurdle_scv2_distance_fix_category_touch) %>%
  identity() -> m_mvhurdle_scv2_distance_fix_category_touch_fitted
m_mvhurdle_scv2_distance_fix_category_touch_fitted


m_mvhurdle_scv2_distance_fix_category_touch_fitted %>%
  mutate(high_touch = case_when(high_touch == TRUE ~ "High Touch",
                                high_touch == FALSE ~ "Low Touch")) %>%
  ggplot(aes(x = distance, y = .epred)) +
  stat_lineribbon() +
  geom_point(data = mutate(as_tibble(m_mvhurdle_scv2_distance_fix_category_touch$data),
                           high_touch = case_when(high_touch == TRUE ~ "High Touch",
                                                  high_touch == FALSE ~ "Low Touch")
  ),
  aes(x = distance, y = copies_max), color = "grey", alpha = 0.5) +
  facet_wrap(facets = ~ site_category + high_touch, scales = "fixed") +
  #facet_wrap(facets = ~ site_category) +
  scale_fill_brewer(palette = "Reds") +
  scale_y_log10() +
  labs(x = "Distance from Patient (meters)",
       y = "Copies SARS-CoV-2 by RT-PCR",
       fill = "Posterior Credible Interval") +
  theme_bw() +
  theme(strip.text = ggtext::element_markdown(color = "black", size = 8),
        axis.text.x = ggtext::element_markdown(color = "black"),
        axis.text.y = ggtext::element_markdown(color = "black"),
        axis.title.x = ggtext::element_markdown(color = "black"),
        axis.title.y = ggtext::element_markdown(color = "black"),
        legend.position = "top",
        #legend.background = element_rect(fill = "white", color = "black", size = 0.25),
        strip.background = element_blank()) -> p_scv2_mvhurdle_distance_site_category_touch
p_scv2_mvhurdle_distance_site_category_touch


# p_scv2_mvhurdle_distance_site_category_touch %>%
#   ggsave(filename = "./figs/p_scv2_mvhurdle_distance_site_category_touch.pdf", height = 5, width = 8, units = "in")
# p_scv2_mvhurdle_distance_site_category_touch %>%
#   ggsave(filename = "./figs/p_scv2_mvhurdle_distance_site_category_touch.png", height = 5, width = 8, units = "in", dpi = 600)
# p_scv2_mvhurdle_distance_site_category_touch %>%
#   ggsave(filename = "./figs/p_scv2_mvhurdle_distance_site_category_touch.svg", height = 5, width = 8, units = "in")



#' ######################
#' model comparison
#' ######################
loo::loo_compare(loo(m_mvhurdle_scv2_distance_fix_category_touch), loo(m_hurdle_scv2_distance_mix_category_touch), loo(m_hurdle_scv2_distance_mix_category)) %>%
  identity() -> hurdle_distance_category_touch_loo_comparison
hurdle_distance_category_touch_loo_comparison

hurdle_distance_category_touch_loo_comparison %>%
  as_tibble(rownames = "model") %>%
  gt()





#' ############################
#' MODELS with SUBJECT-LEVEL RANDOM EFFECTS
#' PSIS diagnositics: plot(loo(m_binom_scv2_distance_mix_category_touch))
#' outliers analysis via random effects
#' ############################


#' ####################################################
#' generative model for SARS-CoV-2 contamination ~ days from COVID diagnosis
#' - elevated vs floor -- accounting for high-touch
#' - binomial probability of detection
#' ####################################################

dat_long %>%
  # exclude bathroom sites and nursing sites (different distance scheme)
  filter(grepl("bathroom|sink|toilet|nurse",site_descriptor) == FALSE) %>%
  mutate(scv2_detected = !is.na(copies_max),
         site_category = case_when(grepl("floor",site_descriptor) == TRUE ~ "floor",
                                   grepl("floor",site_descriptor) == FALSE ~ "elevated"),
         high_touch = touch == "High") %>%
  #count(site_descriptor, site_category, high_touch)
  select(subject_room_id, redcap_id, subject_room_day, distance, subject_hosp_day, subject_covid_day, swab_site, unit, site_category, touch, high_touch, distance, copies_max, scv2_detected) %>%
  mutate(redcap_id = paste0("decon_subject_", redcap_id)) %>%
  mutate(site_category = stringr::str_to_title(gsub("_"," ",site_category))) %>%
  distinct() %>%
  identity() -> dat

dat


#' get prior
dat %>%
  brms::get_prior(data = ., family = bernoulli,
                  scv2_detected ~ (1 + distance + site_category + high_touch | subject_room_id)
  ) %>%
  gt::gt()


#' run binomial model
# dat %>%
#   select(scv2_detected, distance, site_category, high_touch, subject_room_id) %>%
#   brm(formula = scv2_detected ~ (1 + distance + site_category + high_touch | subject_room_id),
#       data = .,
#       family = bernoulli,
#       prior = c(set_prior(prior = "normal(0, 1)",class = "Intercept"), set_prior(prior = "normal(0,1)", class = "sd")),
#       chains = 4,
#       cores = 4,
#       iter = 3000,
#       warmup = 1000,
#       control = list("adapt_delta" = 0.999, max_treedepth = 10),
#       backend = "cmdstanr",
#       seed = 16) -> m_mvbinom_scv2_distance_category_touch_mix_subject
# 
# m_mvbinom_scv2_distance_category_touch_mix_subject %>% write_rds(file = "./models/binomial/m_mvbinom_scv2_distance_category_touch_mix_subject.rds.bz2", compress = "bz2")
m_mvbinom_scv2_distance_category_touch_mix_subject <- read_rds(file = "./models/binomial/m_mvbinom_scv2_distance_category_touch_mix_subject.rds.bz2")

m_mvbinom_scv2_distance_category_touch_mix_subject
rstan::check_hmc_diagnostics(m_mvbinom_scv2_distance_category_touch_mix_subject$fit)
m_mvbinom_scv2_distance_category_touch_mix_subject %>% pp_check()

m_mvbinom_scv2_distance_category_touch_mix_subject %>%
  posterior_summary() %>%
  as_tibble(rownames = "param") %>%
  gt::gt() %>%
  gt::fmt_number(columns = 2:5, n_sigfig = 3)


#' fitted
m_mvbinom_scv2_distance_category_touch_mix_subject$data %>%
  as_tibble() %>%
  expand(distance = modelr::seq_range(distance, n = 30),
         site_category = unique(site_category),
         high_touch = unique(high_touch),
         subject_room_id = unique(subject_room_id)
  ) %>%
  filter(!(high_touch == TRUE & site_category == "Floor")) %>%
  add_epred_draws(m_mvbinom_scv2_distance_category_touch_mix_subject, ndraws = 1000) %>%
  identity() -> m_mvbinom_scv2_distance_category_touch_mix_subject_fitted
m_mvbinom_scv2_distance_category_touch_mix_subject_fitted


m_mvbinom_scv2_distance_category_touch_mix_subject_fitted %>%
  mutate(high_touch = case_when(high_touch == TRUE ~ "High Touch",
                                high_touch == FALSE ~ "Low Touch")) %>%
  ggplot(aes(x = distance, y = .epred)) +
  #geom_point(data = m_mvbinom_scv2_distance_category_touch_mix_subject$data, aes(x = distance, y = scv2_detected), color = "grey", alpha = 0.5) +
  stat_lineribbon() +
  facet_wrap(facets = ~ site_category + high_touch, scales = "fixed") +
  #facet_wrap(facets = ~ site_category) +
  scale_fill_brewer(palette = "Reds") +
  labs(x = "Distance from Patient (meters)",
       y = "Probability of SARS-CoV-2 Detection by RT-PCR",
       fill = "Posterior Credible Interval") +
  theme_bw() +
  theme(strip.text = ggtext::element_markdown(color = "black", size = 8),
        axis.text.x = ggtext::element_markdown(color = "black"),
        axis.text.y = ggtext::element_markdown(color = "black"),
        axis.title.x = ggtext::element_markdown(color = "black"),
        axis.title.y = ggtext::element_markdown(color = "black"),
        legend.position = "top",
        #legend.background = element_rect(fill = "white", color = "black", size = 0.25),
        strip.background = element_blank()) -> p_scv2_mvbinomial_distance_site_category_touch_mix_subject
p_scv2_mvbinomial_distance_site_category_touch_mix_subject


# p_scv2_mvbinomial_distance_site_category_touch_mix_subject %>%
#   ggsave(filename = "./figs/p_scv2_mvbinomial_distance_site_category_touch_mix_subject.pdf", height = 5, width = 8, units = "in")
# p_scv2_mvbinomial_distance_site_category_touch_mix_subject %>%
#   ggsave(filename = "./figs/p_scv2_mvbinomial_distance_site_category_touch_mix_subject.png", height = 5, width = 8, units = "in", dpi = 600)
# p_scv2_mvbinomial_distance_site_category_touch_mix_subject %>%
#   ggsave(filename = "./figs/p_scv2_mvbinomial_distance_site_category_touch_mix_subject.svg", height = 5, width = 8, units = "in")


#' compare subjects
m_mvbinom_scv2_distance_category_touch_mix_subject_fitted %>%
  group_by(subject_room_id, distance, site_category, high_touch) %>%
  summarise(.epred = median(.epred, na.rm = TRUE)) %>%
  mutate(high_touch = case_when(high_touch == TRUE ~ "High Touch",
                                high_touch == FALSE ~ "Low Touch")) %>%
  ggplot(aes(x = distance, y = .epred, group = subject_room_id, color = subject_room_id)) +
  #geom_line(color = "black", alpha = 0.5) +
  geom_line() +
  facet_wrap(facets = ~ site_category + high_touch, scales = "fixed") +
  colorspace::scale_color_discrete_sequential(palette = "YlGnBu") +
  labs(x = "Distance from Patient (meters)",
       y = "Probability of SARS-CoV-2 Detection by RT-PCR",
       fill = "Posterior Credible Interval") +
  theme_bw() +
  theme(strip.text = ggtext::element_markdown(color = "black", size = 8),
        axis.text.x = ggtext::element_markdown(color = "black"),
        axis.text.y = ggtext::element_markdown(color = "black"),
        axis.title.x = ggtext::element_markdown(color = "black"),
        axis.title.y = ggtext::element_markdown(color = "black"),
        legend.position = "none",
        #legend.background = element_rect(fill = "white", color = "black", size = 0.25),
        strip.background = element_blank()) -> p_scv2_mvbinomial_distance_site_category_touch_mix_subject_lines
p_scv2_mvbinomial_distance_site_category_touch_mix_subject_lines


# p_scv2_mvbinomial_distance_site_category_touch_mix_subject_lines %>%
#   ggsave(filename = "./figs/p_scv2_mvbinomial_distance_site_category_touch_mix_subject_lines.pdf", height = 5, width = 8, units = "in")
# p_scv2_mvbinomial_distance_site_category_touch_mix_subject_lines %>%
#   ggsave(filename = "./figs/p_scv2_mvbinomial_distance_site_category_touch_mix_subject_lines.png", height = 5, width = 8, units = "in", dpi = 600)
# p_scv2_mvbinomial_distance_site_category_touch_mix_subject_lines %>%
#   ggsave(filename = "./figs/p_scv2_mvbinomial_distance_site_category_touch_mix_subject_lines.svg", height = 5, width = 8, units = "in")



#' patchwork plot joining fixed-effects and random-subject-effects models

((p_scv2_mvbinomial_distance_site_category_touch + labs(y = "Probability of SARS-CoV-2<br>Detection by RT-PCR")) / 
    (p_scv2_mvbinomial_distance_site_category_touch_mix_subject_lines + labs(y = "Probability of SARS-CoV-2<br>Detection by RT-PCR"))) +
  patchwork::plot_annotation(tag_levels = "A") -> p_distance_models_combined
p_distance_models_combined


# p_distance_models_combined %>%
#   ggsave(filename = "./figs/p_distance_models_combined.pdf", height = 7, width = 8, units = "in")
# p_distance_models_combined %>%
#   ggsave(filename = "./figs/p_distance_models_combined.png", height = 7, width = 8, units = "in", dpi = 600)
# p_distance_models_combined %>%
#   ggsave(filename = "./figs/p_distance_models_combined.svg", height = 7, width = 8, units = "in")


#' examine relationship between subjects' floors and high-touch elevated
m_mvbinom_scv2_distance_category_touch_mix_subject_fitted %>%
  ungroup() %>%
  filter(distance == min(distance, na.rm = TRUE)) %>%
  filter((site_category == "Elevated" & high_touch == TRUE) | site_category == "Floor") %>%
  group_by(subject_room_id, distance, site_category, high_touch) %>%
  summarise(.epred = median(.epred, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(numeric_id = gsub("decon_","",subject_room_id)) %>%
  select(subject_room_id, numeric_id, site_category, .epred) %>%
  pivot_wider(id_cols = c(subject_room_id, numeric_id), names_from = site_category, values_from = .epred) %>%
  left_join(summarise(group_by(dat_long, subject_room_id), total_study_day = min(total_study_day, na.rm = TRUE), subject_room_day = min(subject_room_day, na.rm = TRUE), unit = unit[1]), by = "subject_room_id") %>%
  #mutate(numeric_id = replace(numeric_id, Floor <= 0.75 | Elevated <= 0.75, NA)) %>%
  ggplot(data = ., aes(x = Floor, y = Elevated, color = subject_room_id, label = numeric_id)) +
  #geom_point() +
  geom_point(shape = 21, color = "black") +
  #ggrepel::geom_text_repel(color = "black", max.overlaps = 50) +
  coord_equal(xlim = c(0,1), ylim = c(0,1)) +
  #colorspace::scale_color_discrete_sequential() +
  labs(x = "Probability of SARS-CoV-2<br>Detection on Floor",
       y = "Probability of SARS-CoV-2<br>Detection on High-Touch Elevated",
  ) +
  theme_bw() +
  theme(strip.text = ggtext::element_markdown(color = "black", size = 8),
        axis.text.x = ggtext::element_markdown(color = "black"),
        axis.text.y = ggtext::element_markdown(color = "black"),
        axis.title.x = ggtext::element_markdown(color = "black"),
        axis.title.y = ggtext::element_markdown(color = "black"),
        legend.position = "none",
        #legend.background = element_rect(fill = "white", color = "black", size = 0.25),
        strip.background = element_blank())




#' ####################################################
#' generative model for SARS-CoV-2 contamination ~ days from COVID diagnosis
#' - elevated vs floor -- accounting for high-touch
#' - hurdle quantity detected
#' ####################################################

dat_long %>%
  # exclude bathroom sites and nursing sites (different distance scheme)
  filter(grepl("bathroom|sink|toilet|nurse",site_descriptor) == FALSE) %>%
  mutate(scv2_detected = !is.na(copies_max),
         site_category = case_when(grepl("floor",site_descriptor) == TRUE ~ "floor",
                                   grepl("floor",site_descriptor) == FALSE ~ "elevated"),
         high_touch = touch == "High") %>%
  #count(site_descriptor, site_category, high_touch)
  select(subject_room_id, redcap_id, subject_room_day, distance, subject_hosp_day, distance, swab_site, unit, site_category, touch, high_touch, distance, copies_max, scv2_detected) %>%
  mutate(redcap_id = paste0("decon_subject_", redcap_id)) %>%
  mutate(site_category = stringr::str_to_title(gsub("_"," ",site_category))) %>%
  distinct() %>%
  identity() -> dat

dat


#' get prior
dat %>%
  brms::get_prior(formula = copies_max ~ (1 + distance + site_category + high_touch | subject_room_id),
                  data = .,
                  family = hurdle_lognormal()
  ) %>%
  gt::gt()


#' run hurdle model
# dat %>%
#   select(copies_max, distance, site_category, high_touch, subject_room_id) %>%
#   brm(formula = copies_max ~ (1 + distance + site_category + high_touch | subject_room_id),
#       data = .,
#       family = hurdle_lognormal(),
#       prior = c(set_prior(prior = "normal(0, 1)",class = "Intercept"), set_prior(prior = "normal(0,1)", class = "sd")),
#       chains = 4,
#       cores = 4,
#       iter = 2000,
#       warmup = 1000,
#       control = list("adapt_delta" = 0.999, max_treedepth = 10),
#       backend = "cmdstanr",
#       seed = 16) -> m_mvhurdle_scv2_distance_category_touch_mix_subject
# 
# m_mvhurdle_scv2_distance_category_touch_mix_subject %>% write_rds(file = "./models/hurdle/m_mvhurdle_scv2_distance_category_touch_mix_subject.rds.bz2", compress = "bz2")
m_mvhurdle_scv2_distance_category_touch_mix_subject <- read_rds(file = "./models/hurdle/m_mvhurdle_scv2_distance_category_touch_mix_subject.rds.bz2")

m_mvhurdle_scv2_distance_category_touch_mix_subject
rstan::check_hmc_diagnostics(m_mvhurdle_scv2_distance_category_touch_mix_subject$fit)
m_mvhurdle_scv2_distance_category_touch_mix_subject %>% pp_check() + scale_x_log10()

m_mvhurdle_scv2_distance_category_touch_mix_subject %>%
  posterior_summary() %>%
  as_tibble(rownames = "param") %>%
  gt::gt() %>%
  gt::fmt_number(columns = 2:5, n_sigfig = 3)


#' fitted
m_mvhurdle_scv2_distance_category_touch_mix_subject$data %>%
  as_tibble() %>%
  expand(distance = modelr::seq_range(distance, n = 30),
         site_category = unique(site_category),
         high_touch = unique(high_touch),
         subject_room_id = unique(subject_room_id)
  ) %>%
  filter(!(high_touch == TRUE & site_category == "Floor")) %>%
  add_epred_draws(m_mvhurdle_scv2_distance_category_touch_mix_subject, ndraws = 1000) %>%
  identity() -> m_mvhurdle_scv2_distance_category_touch_mix_subject_fitted
m_mvhurdle_scv2_distance_category_touch_mix_subject_fitted


m_mvhurdle_scv2_distance_category_touch_mix_subject_fitted %>%
  mutate(high_touch = case_when(high_touch == TRUE ~ "High Touch",
                                high_touch == FALSE ~ "Low Touch")) %>%
  ggplot(aes(x = distance, y = .epred)) +
  #geom_point(data = m_mvhurdle_scv2_distance_category_touch_mix_subject$data, aes(x = distance, y = scv2_detected), color = "grey", alpha = 0.5) +
  stat_lineribbon() +
  facet_wrap(facets = ~ site_category + high_touch, scales = "fixed") +
  #facet_wrap(facets = ~ site_category) +
  scale_fill_brewer(palette = "Reds") +
  scale_y_log10() +
  labs(x = "Distance from Patient (meters)",
       y = "Copies SARS-CoV-2 by RT-PCR",
       fill = "Posterior Credible Interval") +
  theme_bw() +
  theme(strip.text = ggtext::element_markdown(color = "black", size = 8),
        axis.text.x = ggtext::element_markdown(color = "black"),
        axis.text.y = ggtext::element_markdown(color = "black"),
        axis.title.x = ggtext::element_markdown(color = "black"),
        axis.title.y = ggtext::element_markdown(color = "black"),
        legend.position = "top",
        #legend.background = element_rect(fill = "white", color = "black", size = 0.25),
        strip.background = element_blank()) -> p_scv2_mvhurdle_distance_site_category_touch_mix_subject
p_scv2_mvhurdle_distance_site_category_touch_mix_subject


#' compare subjects
m_mvhurdle_scv2_distance_category_touch_mix_subject_fitted %>%
  group_by(subject_room_id, distance, site_category, high_touch) %>%
  summarise(.epred = median(.epred, na.rm = TRUE)) %>%
  mutate(high_touch = case_when(high_touch == TRUE ~ "High Touch",
                                high_touch == FALSE ~ "Low Touch")) %>%
  ggplot(aes(x = distance, y = .epred, group = subject_room_id, color = subject_room_id)) +
  #geom_line(color = "black", alpha = 0.5) +
  geom_line() +
  facet_wrap(facets = ~ site_category + high_touch, scales = "fixed") +
  colorspace::scale_color_discrete_sequential(palette = "YlGnBu") +
  scale_y_log10() +
  labs(x = "Distance from Patient (meters)",
       y = "Copies SARS-CoV-2 by RT-PCR",
       fill = "Posterior Credible Interval") +
  theme_bw() +
  theme(strip.text = ggtext::element_markdown(color = "black", size = 8),
        axis.text.x = ggtext::element_markdown(color = "black"),
        axis.text.y = ggtext::element_markdown(color = "black"),
        axis.title.x = ggtext::element_markdown(color = "black"),
        axis.title.y = ggtext::element_markdown(color = "black"),
        legend.position = "none",
        #legend.background = element_rect(fill = "white", color = "black", size = 0.25),
        strip.background = element_blank()) -> p_scv2_mvhurdle_distance_site_category_touch_mix_subject_lines
p_scv2_mvhurdle_distance_site_category_touch_mix_subject_lines







