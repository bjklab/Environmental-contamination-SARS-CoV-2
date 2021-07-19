#' ##############################
#' load libraries and set seed
#' ##############################
library(tidyverse)
library(gt)
library(tidybayes)
library(brms)

set.seed(16)


#' ####################################################
#' read data and generate descriptive tables and plots
#' ####################################################
read_csv("./data/dat_long.csv") %>%
  #mutate(unit = replace(unit, unit=="M1", "R1")) %>%
  identity() -> dat_long

dat_long %>%
  group_by(redcap_id) %>%
  filter(subject_study_day == min(subject_study_day)) %>%
  ungroup() %>%
  identity() -> dat_first



#' ####################################################
#' generative model for SARS-CoV-2 contamination
#' - linear abundance data
#' - binomial probability of detection
#' ####################################################

dat_long %>%
  mutate(scv2_detected = !is.na(copies_max),
         site_category = case_when(grepl("wall", site_descriptor) ~ "wall",
                                   grepl("floor", site_descriptor) ~ "floor",
                                   grepl("sink|toilet", site_descriptor) ~ "bathroom",
                                   grepl("doorknob|bed|keyboard|mouse", site_descriptor) ~ "high_touch")) %>%
  #count(site_descriptor, site_category)
  select(redcap_id, subject_covid_day, subject_hosp_day, subject_covid_day, swab_site, unit, site_category, touch, distance, copies_max, scv2_detected) %>%
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
dat %>%
  filter(site_category != "Bathroom") %>%
  brm(formula = scv2_detected ~ 0 + (1 + distance | site_category),
      data = .,
      family = bernoulli,
      chains = 4,
      cores = 4,
      control = list("adapt_delta" = 0.999, max_treedepth = 18),
      backend = "cmdstanr",
      seed = 16) -> m_binom_scv2_distance_mix_category

m_binom_scv2_distance_mix_category %>% write_rds(file = "./models/binomial/m_binom_scv2_distance_mix_category.rds.gz", compress = "gz")
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
         site_category = unique(site_category)) %>%
  add_fitted_draws(m_binom_scv2_distance_mix_category) %>%
  identity() -> m_binom_scv2_distance_mix_category_fitted
m_binom_scv2_distance_mix_category_fitted


m_binom_scv2_distance_mix_category_fitted %>%
  ggplot(aes(x = distance, y = .value)) +
  #geom_point(data = m_binom_scv2_distance_mix_category$data, aes(x = distance, y = scv2_detected), color = "grey", alpha = 0.5) +
  stat_lineribbon() +
  facet_wrap(facets = ~ site_category, scales = "fixed") +
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



#' run linear model
# dat %>%
#   filter(site_category != "Bathroom") %>%
#   brm(formula = log10(copies_max) ~ 0 + (1 + distance | site_category),
#       data = .,
#       family = gaussian,
#       chains = 4,
#       cores = 4,
#       control = list("adapt_delta" = 0.999, max_treedepth = 18),
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
         site_category = unique(site_category)) %>%
  add_fitted_draws(m_linear_scv2_distance_mix_category) %>%
  identity() -> m_linear_scv2_distance_mix_category_fitted
m_linear_scv2_distance_mix_category_fitted


m_linear_scv2_distance_mix_category_fitted %>%
  ggplot(aes(x = distance, y = .value)) +
  geom_point(data = m_linear_scv2_distance_mix_category$data, aes(x = distance, y = log10(copies_max)), color = "grey", alpha = 0.5) +
  stat_lineribbon() +
  facet_wrap(facets = ~ site_category, scales = "fixed") +
  scale_fill_brewer(palette = "Reds") +
  labs(x = "Distance from Patient (meters)",
       y = "Expected Log<sub>10</sub> Copies of SARS-CoV-2 by RT-PCR",
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



#' ###################################
#' DISTANCE MODELS: FLOOR vs ELEVATED
#' - nested random effects? for high-touch
#' ###################################




