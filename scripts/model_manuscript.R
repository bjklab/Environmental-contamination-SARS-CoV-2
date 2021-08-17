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
#' ####################################################
#' SPATIAL MODELS
#' ####################################################
#' ####################################################

dat_first %>%
  # exclude bathroom sites and nursing sites (different distance scheme)
  filter(grepl("bathroom|sink|toilet|nurse",site_descriptor) == FALSE) %>%
  mutate(scv2_detected = !is.na(copies_max),
         site_category = case_when(grepl("floor",site_descriptor) == TRUE ~ "floor",
                                   grepl("floor",site_descriptor) == FALSE ~ "elevated"),
         high_touch = touch == "High") %>%
  #count(site_descriptor, site_category, high_touch) %>% gt()
  select(redcap_id, subject_covid_day, subject_hosp_day, subject_room_day, total_study_day, swab_site, unit, site_category, touch, high_touch, distance, copies_max, scv2_detected) %>%
  mutate(site_category = stringr::str_to_title(gsub("_"," ",site_category))) %>%
  distinct() %>%
  identity() -> dat

dat


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
  select(redcap_id, subject_covid_day, subject_hosp_day, subject_room_day, total_study_day, swab_site, unit, site_category, touch, high_touch, distance, copies_max, scv2_detected) %>%
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
  scale_y_continuous(limits = c(0,1)) +
  facet_wrap(facets = ~ site_category + high_touch, scales = "fixed") +
  #facet_wrap(facets = ~ site_category) +
  scale_fill_brewer(palette = "Reds") +
  labs(x = "Distance from Patient (meters)",
       y = "Probability of SARS-CoV-2<br>Detection by RT-PCR",
       fill = "Posterior<br>Credible<br>Interval") +
  theme_bw() +
  theme(strip.text = ggtext::element_markdown(color = "black", size = 10),
        axis.text.x = ggtext::element_markdown(color = "black"),
        axis.text.y = ggtext::element_markdown(color = "black"),
        axis.title.x = ggtext::element_markdown(color = "black"),
        axis.title.y = ggtext::element_markdown(color = "black"),
        legend.title = ggtext::element_markdown(color = "black"),
        legend.position = "right",
        #legend.background = element_rect(fill = "white", color = "black", size = 0.25),
        strip.background = element_blank()) -> p_scv2_mvbinomial_distance_site_category_touch
p_scv2_mvbinomial_distance_site_category_touch


# p_scv2_mvbinomial_distance_site_category_touch %>%
#   ggsave(filename = "./figs/p_scv2_mvbinomial_distance_site_category_touch.pdf", height = 5, width = 8, units = "in")
# p_scv2_mvbinomial_distance_site_category_touch %>%
#   ggsave(filename = "./figs/p_scv2_mvbinomial_distance_site_category_touch.png", height = 5, width = 8, units = "in", dpi = 600)
# p_scv2_mvbinomial_distance_site_category_touch %>%
#   ggsave(filename = "./figs/p_scv2_mvbinomial_distance_site_category_touch.svg", height = 5, width = 8, units = "in")







#' ####################################################
#' ####################################################
#' TEMPORAL MODELS
#' ####################################################
#' ####################################################

dat_long %>%
  # exclude bathroom sites and nursing sites (different distance scheme)
  filter(grepl("bathroom|sink|toilet|nurse",site_descriptor) == FALSE) %>%
  mutate(scv2_detected = !is.na(copies_max),
         site_category = case_when(grepl("floor",site_descriptor) == TRUE ~ "floor",
                                   grepl("floor",site_descriptor) == FALSE ~ "elevated"),
         high_touch = touch == "High") %>%
  #count(site_descriptor, site_category, high_touch)
  select(subject_room_id, redcap_id, subject_room_day, subject_covid_day, subject_hosp_day, total_study_day, swab_site, unit, site_category, touch, high_touch, distance, copies_max, scv2_detected) %>%
  mutate(redcap_id = paste0("decon_subject_", redcap_id)) %>%
  mutate(site_category = stringr::str_to_title(gsub("_"," ",site_category))) %>%
  distinct() %>%
  identity() -> dat

dat




#' ####################################################
#' generative model for SARS-CoV-2 contamination ~ days from COVID diagnosis
#' - elevated vs floor -- accounting for high-touch
#' - binomial probability of detection
#' ####################################################

#' get prior
dat %>%
  brms::get_prior(data = ., family = bernoulli,
                  scv2_detected ~ 1 + subject_covid_day + site_category + high_touch
  ) %>%
  gt::gt()


#' run binomial model
# dat %>%
#   select(scv2_detected, subject_covid_day, site_category, high_touch) %>%
#   brm(formula = scv2_detected ~ 1 + subject_covid_day + site_category + high_touch,
#       data = .,
#       family = bernoulli,
#       chains = 4,
#       cores = 4,
#       control = list("adapt_delta" = 0.999, max_treedepth = 10),
#       backend = "cmdstanr",
#       seed = 16) -> m_mvbinom_scv2_time_fix_category_touch
# 
# m_mvbinom_scv2_time_fix_category_touch %>% write_rds(file = "./models/binomial/m_mvbinom_scv2_time_fix_category_touch.rds.gz", compress = "gz")
m_mvbinom_scv2_time_fix_category_touch <- read_rds(file = "./models/binomial/m_mvbinom_scv2_time_fix_category_touch.rds.gz")

m_mvbinom_scv2_time_fix_category_touch
rstan::check_hmc_diagnostics(m_mvbinom_scv2_time_fix_category_touch$fit)
m_mvbinom_scv2_time_fix_category_touch %>% pp_check()

m_mvbinom_scv2_time_fix_category_touch %>%
  posterior_summary() %>%
  as_tibble(rownames = "param") %>%
  gt::gt() %>%
  gt::fmt_number(columns = 2:5, n_sigfig = 3)


#' fitted
m_mvbinom_scv2_time_fix_category_touch$data %>%
  as_tibble() %>%
  expand(subject_covid_day = modelr::seq_range(subject_covid_day, n = 100),
         site_category = unique(site_category),
         high_touch = unique(high_touch)
  ) %>%
  filter(!(high_touch == TRUE & site_category == "Floor")) %>%
  add_epred_draws(m_mvbinom_scv2_time_fix_category_touch, ndraws = 1000) %>%
  identity() -> m_mvbinom_scv2_time_fix_category_touch_fitted
m_mvbinom_scv2_time_fix_category_touch_fitted


m_mvbinom_scv2_time_fix_category_touch_fitted %>%
  mutate(high_touch = case_when(high_touch == TRUE ~ "High Touch",
                                high_touch == FALSE ~ "Low Touch")) %>%
  ggplot(aes(x = subject_covid_day, y = .epred)) +
  #geom_point(data = m_mvbinom_scv2_time_fix_category_touch$data, aes(x = subject_covid_day, y = scv2_detected), color = "grey", alpha = 0.5) +
  stat_lineribbon() +
  facet_wrap(facets = ~ site_category + high_touch, scales = "fixed") +
  #facet_wrap(facets = ~ site_category) +
  scale_fill_brewer(palette = "Reds") +
  scale_y_continuous(limits = c(0,1)) +
  labs(x = "Days after COVID-19 Diagnosis",
       y = "Probability of SARS-CoV-2<br>Detection by RT-PCR",
       fill = "Posterior<br>Credible<br>Interval") +
  theme_bw() +
  theme(strip.text = ggtext::element_markdown(color = "black", size = 10),
        axis.text.x = ggtext::element_markdown(color = "black"),
        axis.text.y = ggtext::element_markdown(color = "black"),
        axis.title.x = ggtext::element_markdown(color = "black"),
        axis.title.y = ggtext::element_markdown(color = "black"),
        legend.title = ggtext::element_markdown(color = "black"),
        legend.position = "right",
        #legend.background = element_rect(fill = "white", color = "black", size = 0.25),
        strip.background = element_blank()) -> p_scv2_mvbinomial_subject_covid_day_site_category_touch
p_scv2_mvbinomial_subject_covid_day_site_category_touch


# p_scv2_mvbinomial_subject_covid_day_site_category_touch %>%
#   ggsave(filename = "./figs/p_scv2_mvbinomial_subject_covid_day_site_category_touch.pdf", height = 5, width = 8, units = "in")
# p_scv2_mvbinomial_subject_covid_day_site_category_touch %>%
#   ggsave(filename = "./figs/p_scv2_mvbinomial_subject_covid_day_site_category_touch.png", height = 5, width = 8, units = "in", dpi = 600)
# p_scv2_mvbinomial_subject_covid_day_site_category_touch %>%
#   ggsave(filename = "./figs/p_scv2_mvbinomial_subject_covid_day_site_category_touch.svg", height = 5, width = 8, units = "in")




#' ####################################################
#' ####################################################
#' COMBINE UNADJUSTED MODEL PLOTS
#' ####################################################
#' ####################################################

#' patchwork plot joining fixed-effects and random-subject-effects models

((p_scv2_mvbinomial_distance_site_category_touch) / 
    (p_scv2_mvbinomial_subject_covid_day_site_category_touch)) + #theme(strip.text = element_blank()))) +
  patchwork::plot_layout(guides = 'collect') +
  patchwork::plot_annotation(tag_levels = "A") -> p_combined_distance_time_unadjusted
p_combined_distance_time_unadjusted


# p_combined_distance_time_unadjusted %>%
#   ggsave(filename = "./figs/p_combined_distance_time_unadjusted.pdf", height = 7, width = 8, units = "in")
# p_combined_distance_time_unadjusted %>%
#   ggsave(filename = "./figs/p_combined_distance_time_unadjusted.png", height = 7, width = 8, units = "in", dpi = 600)
# p_combined_distance_time_unadjusted %>%
#   ggsave(filename = "./figs/p_combined_distance_time_unadjusted.svg", height = 7, width = 8, units = "in")








#' ####################################################
#' ####################################################
#' SUBJECT-LEVEL RANDOM EFFECTS MODELS
#' ####################################################
#' ####################################################

dat_long %>%
  # exclude bathroom sites and nursing sites (different distance scheme)
  filter(grepl("bathroom|sink|toilet|nurse",site_descriptor) == FALSE) %>%
  mutate(scv2_detected = !is.na(copies_max),
         site_category = case_when(grepl("floor",site_descriptor) == TRUE ~ "floor",
                                   grepl("floor",site_descriptor) == FALSE ~ "elevated"),
         high_touch = touch == "High") %>%
  #count(site_descriptor, site_category, high_touch)
  select(subject_room_id, redcap_id, subject_room_day, subject_covid_day, subject_hosp_day, total_study_day, swab_site, unit, site_category, touch, high_touch, distance, copies_max, scv2_detected) %>%
  mutate(redcap_id = paste0("decon_subject_", redcap_id)) %>%
  mutate(site_category = stringr::str_to_title(gsub("_"," ",site_category))) %>%
  distinct() %>%
  identity() -> dat

dat



#' ####################################################
#' generative model for SARS-CoV-2 contamination ~ days from COVID diagnosis
#' - elevated vs floor -- accounting for high-touch
#' - binomial probability of detection
#' ####################################################

#' get prior
dat %>%
  brms::get_prior(data = ., family = bernoulli,
                  scv2_detected ~ (1 + subject_covid_day + site_category + high_touch | subject_room_id)
  ) %>%
  gt::gt()


#' run binomial model
# dat %>%
#   select(scv2_detected, subject_covid_day, site_category, high_touch, subject_room_id) %>%
#   brm(formula = scv2_detected ~ (1 + subject_covid_day + site_category + high_touch | subject_room_id),
#       data = .,
#       family = bernoulli,
#       prior = c(set_prior(prior = "normal(0, 1)",class = "Intercept"), set_prior(prior = "normal(0,1)", class = "sd")),
#       chains = 4,
#       cores = 4,
#       iter = 3000,
#       warmup = 1000,
#       control = list("adapt_delta" = 0.999, max_treedepth = 10),
#       backend = "cmdstanr",
#       seed = 16) -> m_mvbinom_scv2_time_category_touch_mix_subject
# 
# m_mvbinom_scv2_time_category_touch_mix_subject %>% write_rds(file = "./models/binomial/m_mvbinom_scv2_time_category_touch_mix_subject.rds.bz2", compress = "bz2")
m_mvbinom_scv2_time_category_touch_mix_subject <- read_rds(file = "./models/binomial/m_mvbinom_scv2_time_category_touch_mix_subject.rds.bz2")

m_mvbinom_scv2_time_category_touch_mix_subject
rstan::check_hmc_diagnostics(m_mvbinom_scv2_time_category_touch_mix_subject$fit)
m_mvbinom_scv2_time_category_touch_mix_subject %>% pp_check()

m_mvbinom_scv2_time_category_touch_mix_subject %>%
  posterior_summary() %>%
  as_tibble(rownames = "param") %>%
  gt::gt() %>%
  gt::fmt_number(columns = 2:5, n_sigfig = 3)


#' fitted
m_mvbinom_scv2_time_category_touch_mix_subject$data %>%
  as_tibble() %>%
  expand(subject_covid_day = modelr::seq_range(subject_covid_day, n = 40),
         site_category = unique(site_category),
         high_touch = unique(high_touch),
         subject_room_id = unique(subject_room_id)
  ) %>%
  filter(!(high_touch == TRUE & site_category == "Floor")) %>%
  add_epred_draws(m_mvbinom_scv2_time_category_touch_mix_subject, ndraws = 1000) %>%
  identity() -> m_mvbinom_scv2_time_category_touch_mix_subject_fitted
m_mvbinom_scv2_time_category_touch_mix_subject_fitted


# m_mvbinom_scv2_time_category_touch_mix_subject_fitted %>%
#   mutate(high_touch = case_when(high_touch == TRUE ~ "High Touch",
#                                 high_touch == FALSE ~ "Low Touch")) %>%
#   ggplot(aes(x = subject_covid_day, y = .epred)) +
#   #geom_point(data = m_mvbinom_scv2_time_category_touch_mix_subject$data, aes(x = subject_covid_day, y = scv2_detected), color = "grey", alpha = 0.5) +
#   stat_lineribbon() +
#   facet_wrap(facets = ~ site_category + high_touch, scales = "fixed") +
#   #facet_wrap(facets = ~ site_category) +
#   scale_fill_brewer(palette = "Reds") +
#   labs(x = "Days after COVID-19 Diagnosis",
#        y = "Probability of SARS-CoV-2 Detection by RT-PCR",
#        fill = "Posterior Credible Interval") +
#   theme_bw() +
#   theme(strip.text = ggtext::element_markdown(color = "black", size = 8),
#         axis.text.x = ggtext::element_markdown(color = "black"),
#         axis.text.y = ggtext::element_markdown(color = "black"),
#         axis.title.x = ggtext::element_markdown(color = "black"),
#         axis.title.y = ggtext::element_markdown(color = "black"),
#         legend.position = "top",
#         #legend.background = element_rect(fill = "white", color = "black", size = 0.25),
#         strip.background = element_blank())


#' compare subjects
m_mvbinom_scv2_time_category_touch_mix_subject_fitted %>%
  group_by(subject_room_id, subject_covid_day, site_category, high_touch) %>%
  summarise(.epred = median(.epred, na.rm = TRUE)) %>%
  mutate(high_touch = case_when(high_touch == TRUE ~ "High Touch",
                                high_touch == FALSE ~ "Low Touch")) %>%
  ungroup() %>%
  left_join(select(ungroup(summarise(group_by(dat,subject_room_id), pandemic_index = min(total_study_day,na.rm=TRUE))), subject_room_id, pandemic_index), by = "subject_room_id") %>%
  ggplot(aes(x = subject_covid_day, y = .epred, group = subject_room_id, color = pandemic_index)) +
  #geom_line(color = "black", alpha = 0.5) +
  geom_line() +
  scale_y_continuous(limits = c(0,1)) +
  facet_wrap(facets = ~ site_category + high_touch, scales = "fixed") +
  colorspace::scale_color_continuous_sequential(palette = "YlGnBu") +
  labs(x = "Days after COVID-19 Diagnosis",
       y = "Probability of SARS-CoV-2<br>Detection by RT-PCR",
       color = "Days from<br>Start of<br>Local<br>COVID-19<br>Wave") +
  theme_bw() +
  theme(strip.text = ggtext::element_markdown(color = "black", size = 10),
        axis.text.x = ggtext::element_markdown(color = "black"),
        axis.text.y = ggtext::element_markdown(color = "black"),
        axis.title.x = ggtext::element_markdown(color = "black"),
        axis.title.y = ggtext::element_markdown(color = "black"),
        legend.title = ggtext::element_markdown(color = "black"),
        legend.position = "right",
        #legend.background = element_rect(fill = "white", color = "black", size = 0.25),
        strip.background = element_blank()) -> p_scv2_mvbinomial_time_site_category_touch_mix_subject_lines
p_scv2_mvbinomial_time_site_category_touch_mix_subject_lines


# p_scv2_mvbinomial_time_site_category_touch_mix_subject_lines %>%
#   ggsave(filename = "./figs/p_scv2_mvbinomial_time_site_category_touch_mix_subject_lines.pdf", height = 5, width = 8, units = "in")
# p_scv2_mvbinomial_time_site_category_touch_mix_subject_lines %>%
#   ggsave(filename = "./figs/p_scv2_mvbinomial_time_site_category_touch_mix_subject_lines.png", height = 5, width = 8, units = "in", dpi = 600)
# p_scv2_mvbinomial_time_site_category_touch_mix_subject_lines %>%
#   ggsave(filename = "./figs/p_scv2_mvbinomial_time_site_category_touch_mix_subject_lines.svg", height = 5, width = 8, units = "in")



#' ####################################################
#' ####################################################
#' ADJUST FOR PANDEMIC WAVE INDEX
#' ####################################################
#' ####################################################


dat_long %>%
  # exclude bathroom sites and nursing sites (different distance scheme)
  filter(grepl("bathroom|sink|toilet|nurse",site_descriptor) == FALSE) %>%
  mutate(scv2_detected = !is.na(copies_max),
         site_category = case_when(grepl("floor",site_descriptor) == TRUE ~ "floor",
                                   grepl("floor",site_descriptor) == FALSE ~ "elevated"),
         high_touch = touch == "High") %>%
  #count(site_descriptor, site_category, high_touch)
  select(subject_room_id, redcap_id, subject_room_day, subject_covid_day, subject_hosp_day, total_study_day, swab_site, unit, site_category, touch, high_touch, distance, copies_max, scv2_detected) %>%
  mutate(redcap_id = paste0("decon_subject_", redcap_id)) %>%
  mutate(site_category = stringr::str_to_title(gsub("_"," ",site_category))) %>%
  distinct() %>%
  identity() -> dat

dat




#' ####################################################
#' generative model for SARS-CoV-2 contamination ~ days from COVID diagnosis, adjusted for total days since second wave
#' - elevated vs floor -- accounting for high-touch
#' - binomial probability of detection
#' ####################################################

#' get prior
dat %>%
  brms::get_prior(data = ., family = bernoulli,
                  scv2_detected ~ 1 + subject_covid_day + site_category * total_study_day + high_touch
  ) %>%
  gt::gt()


#' run binomial model
# dat %>%
#   select(scv2_detected, subject_covid_day, total_study_day, site_category, high_touch) %>%
#   brm(formula = scv2_detected ~ 1 + subject_covid_day + site_category * total_study_day + high_touch,
#       data = .,
#       family = bernoulli,
#       chains = 4,
#       cores = 4,
#       control = list("adapt_delta" = 0.999, max_treedepth = 10),
#       backend = "cmdstanr",
#       seed = 16) -> m_mvbinom_scv2_time_fix_category_touch_adjust_wave
# 
# m_mvbinom_scv2_time_fix_category_touch_adjust_wave %>% write_rds(file = "./models/binomial/m_mvbinom_scv2_time_fix_category_touch_adjust_wave.rds.gz", compress = "gz")
m_mvbinom_scv2_time_fix_category_touch_adjust_wave <- read_rds(file = "./models/binomial/m_mvbinom_scv2_time_fix_category_touch_adjust_wave.rds.gz")

m_mvbinom_scv2_time_fix_category_touch_adjust_wave
rstan::check_hmc_diagnostics(m_mvbinom_scv2_time_fix_category_touch_adjust_wave$fit)
m_mvbinom_scv2_time_fix_category_touch_adjust_wave %>% pp_check()

m_mvbinom_scv2_time_fix_category_touch_adjust_wave %>%
  posterior_summary() %>%
  as_tibble(rownames = "param") %>%
  gt::gt() %>%
  gt::fmt_number(columns = 2:5, n_sigfig = 3)


#' fitted
m_mvbinom_scv2_time_fix_category_touch_adjust_wave$data %>%
  as_tibble() %>%
  expand(subject_covid_day = modelr::seq_range(subject_covid_day, n = 10),
         site_category = unique(site_category),
         high_touch = unique(high_touch),
         total_study_day = modelr::seq_range(total_study_day, n = 10)
  ) %>%
  filter(!(high_touch == TRUE & site_category == "Floor")) %>%
  add_epred_draws(m_mvbinom_scv2_time_fix_category_touch_adjust_wave, ndraws = 1000) %>%
  identity() -> m_mvbinom_scv2_time_fix_category_touch_adjust_wave_fitted
m_mvbinom_scv2_time_fix_category_touch_adjust_wave_fitted


# m_mvbinom_scv2_time_fix_category_touch_adjust_wave_fitted %>%
#   mutate(high_touch = case_when(high_touch == TRUE ~ "High Touch",
#                                 high_touch == FALSE ~ "Low Touch")) %>%
#   ggplot(aes(x = subject_covid_day, y = .epred)) +
#   #geom_point(data = m_mvbinom_scv2_time_fix_category_touch_adjust_wave$data, aes(x = subject_covid_day, y = scv2_detected), color = "grey", alpha = 0.5) +
#   stat_lineribbon() +
#   facet_wrap(facets = ~ site_category + high_touch, scales = "fixed") +
#   #facet_wrap(facets = ~ site_category) +
#   scale_fill_brewer(palette = "Purples") +
#   scale_y_continuous(limits = c(0,1)) +
#   labs(x = "Days after COVID-19 Diagnosis",
#        y = "Probability of SARS-CoV-2<br>Detection by RT-PCR",
#        fill = "Posterior<br>Credible<br>Interval") +
#   theme_bw() +
#   theme(strip.text = ggtext::element_markdown(color = "black", size = 10),
#         axis.text.x = ggtext::element_markdown(color = "black"),
#         axis.text.y = ggtext::element_markdown(color = "black"),
#         axis.title.x = ggtext::element_markdown(color = "black"),
#         axis.title.y = ggtext::element_markdown(color = "black"),
#         legend.title = ggtext::element_markdown(color = "black"),
#         legend.position = "right",
#         #legend.background = element_rect(fill = "white", color = "black", size = 0.25),
#         strip.background = element_blank())


m_mvbinom_scv2_time_fix_category_touch_adjust_wave_fitted %>%
  mutate(high_touch = case_when(high_touch == TRUE ~ "High Touch",
                                high_touch == FALSE ~ "Low Touch")) %>%
  ggplot(aes(x = total_study_day, y = .epred)) +
  #geom_point(data = m_mvbinom_scv2_time_fix_category_touch_adjust_wave$data, aes(x = subject_covid_day, y = scv2_detected), color = "grey", alpha = 0.5) +
  stat_lineribbon() +
  facet_wrap(facets = ~ site_category + high_touch, scales = "fixed") +
  #facet_wrap(facets = ~ site_category) +
  scale_fill_brewer(palette = "Blues") +
  scale_y_continuous(limits = c(0,1)) +
  labs(x = "Days from Start of Second COVID-19 Wave",
       y = "Probability of SARS-CoV-2<br>Detection by RT-PCR",
       fill = "Posterior<br>Credible<br>Interval") +
  theme_bw() +
  theme(strip.text = ggtext::element_markdown(color = "black", size = 10),
        axis.text.x = ggtext::element_markdown(color = "black"),
        axis.text.y = ggtext::element_markdown(color = "black"),
        axis.title.x = ggtext::element_markdown(color = "black"),
        axis.title.y = ggtext::element_markdown(color = "black"),
        legend.title = ggtext::element_markdown(color = "black"),
        legend.position = "right",
        #legend.background = element_rect(fill = "white", color = "black", size = 0.25),
        strip.background = element_blank()) -> p_scv2_mvbinomial_subject_covid_day_site_category_touch_adjust_wave
p_scv2_mvbinomial_subject_covid_day_site_category_touch_adjust_wave



# p_scv2_mvbinomial_subject_covid_day_site_category_touch_adjust_wave %>%
#   ggsave(filename = "./figs/p_scv2_mvbinomial_subject_covid_day_site_category_touch_adjust_wave.pdf", height = 5, width = 8, units = "in")
# p_scv2_mvbinomial_subject_covid_day_site_category_touch_adjust_wave %>%
#   ggsave(filename = "./figs/p_scv2_mvbinomial_subject_covid_day_site_category_touch_adjust_wave.png", height = 5, width = 8, units = "in", dpi = 600)
# p_scv2_mvbinomial_subject_covid_day_site_category_touch_adjust_wave %>%
#   ggsave(filename = "./figs/p_scv2_mvbinomial_subject_covid_day_site_category_touch_adjust_wave.svg", height = 5, width = 8, units = "in")



#' ####################################################
#' ####################################################
#' COMBINE SUBJECT-LEVEL & ADJUSTED MODEL PLOTS
#' ####################################################
#' ####################################################

#' patchwork plot joining fixed-effects and random-subject-effects models

((p_scv2_mvbinomial_time_site_category_touch_mix_subject_lines) / 
    (p_scv2_mvbinomial_subject_covid_day_site_category_touch_adjust_wave)) + #theme(strip.text = element_blank()))) +
  #patchwork::plot_layout(guides = 'collect') +
  patchwork::plot_annotation(tag_levels = "A") -> p_combined_subject_wave_time_adjusted
p_combined_subject_wave_time_adjusted


# p_combined_subject_wave_time_adjusted %>%
#   ggsave(filename = "./figs/p_combined_subject_wave_time_adjusted.pdf", height = 7, width = 8, units = "in")
# p_combined_subject_wave_time_adjusted %>%
#   ggsave(filename = "./figs/p_combined_subject_wave_time_adjusted.png", height = 7, width = 8, units = "in", dpi = 600)
# p_combined_subject_wave_time_adjusted %>%
#   ggsave(filename = "./figs/p_combined_subject_wave_time_adjusted.svg", height = 7, width = 8, units = "in")



















