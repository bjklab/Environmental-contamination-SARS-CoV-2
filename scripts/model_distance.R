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
  #mutate(unit = replace(unit, unit=="M1", "R1")) %>%
  identity() -> dat_long
dat_long

dat_long %>%
  group_by(redcap_id) %>%
  filter(subject_study_day == min(subject_study_day)) %>%
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
#       control = list("adapt_delta" = 0.999, max_treedepth = 10),
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
  add_fitted_draws(m_binom_scv2_distance_mix_category) %>%
  identity() -> m_binom_scv2_distance_mix_category_fitted
m_binom_scv2_distance_mix_category_fitted


m_binom_scv2_distance_mix_category_fitted %>%
  #mutate(high_touch = case_when(high_touch == TRUE ~ "High Touch",
  #                              high_touch == FALSE ~ "Low Touch")) %>%
  ggplot(aes(x = distance, y = .value)) +
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
  add_fitted_draws(m_binom_scv2_distance_mix_category_touch) %>%
  identity() -> m_binom_scv2_distance_mix_category_touch_fitted
m_binom_scv2_distance_mix_category_touch_fitted


m_binom_scv2_distance_mix_category_touch_fitted %>%
  mutate(high_touch = case_when(high_touch == TRUE ~ "High Touch",
                                high_touch == FALSE ~ "Low Touch")) %>%
  ggplot(aes(x = distance, y = .value)) +
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
#       control = list("adapt_delta" = 0.999, max_treedepth = 10),
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
  add_fitted_draws(m_linear_scv2_distance_mix_category) %>%
  identity() -> m_linear_scv2_distance_mix_category_fitted
m_linear_scv2_distance_mix_category_fitted


m_linear_scv2_distance_mix_category_fitted %>%
  #mutate(high_touch = case_when(high_touch == TRUE ~ "High Touch",
  #                              high_touch == FALSE ~ "Low Touch")) %>%
  ggplot(aes(x = distance, y = .value)) +
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
#       control = list("adapt_delta" = 0.99999, max_treedepth = 10),
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
  add_fitted_draws(m_linear_scv2_distance_mix_category_touch) %>%
  identity() -> m_linear_scv2_distance_mix_category_touch_fitted
m_linear_scv2_distance_mix_category_touch_fitted


m_linear_scv2_distance_mix_category_touch_fitted %>%
  mutate(high_touch = case_when(high_touch == TRUE ~ "High Touch",
                                high_touch == FALSE ~ "Low Touch")) %>%
  ggplot(aes(x = distance, y = .value)) +
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
#       control = list("adapt_delta" = 0.999, max_treedepth = 10),
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
  add_fitted_draws(m_hurdle_scv2_distance_mix_category) %>%
  identity() -> m_hurdle_scv2_distance_mix_category_fitted
m_hurdle_scv2_distance_mix_category_fitted


m_hurdle_scv2_distance_mix_category_fitted %>%
  #mutate(high_touch = case_when(high_touch == TRUE ~ "High Touch",
  #                              high_touch == FALSE ~ "Low Touch")) %>%
  ggplot(aes(x = distance, y = .value)) +
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
#       control = list("adapt_delta" = 0.99999, max_treedepth = 10),
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
  add_fitted_draws(m_hurdle_scv2_distance_mix_category_touch) %>%
  identity() -> m_hurdle_scv2_distance_mix_category_touch_fitted
m_hurdle_scv2_distance_mix_category_touch_fitted


m_hurdle_scv2_distance_mix_category_touch_fitted %>%
  mutate(high_touch = case_when(high_touch == TRUE ~ "High Touch",
                                high_touch == FALSE ~ "Low Touch")) %>%
  ggplot(aes(x = distance, y = .value)) +
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






