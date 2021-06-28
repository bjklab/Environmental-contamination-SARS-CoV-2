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
  mutate(unit = replace(unit, unit=="M1", "R1")) %>%
  identity() -> dat_long

dat_long %>%
  group_by(redcap_id) %>%
  filter(subject_study_day == min(subject_study_day)) %>%
  ungroup() %>%
  identity() -> dat_first

# sample site summary
dat_first %>%
  group_by(swab_site, site_descriptor) %>%
  summarise(prop_positive = sum(!is.na(copies_max)) / n()) %>%
  ungroup() %>%
  arrange(desc(prop_positive)) %>%
  gt::gt() %>%
  gt::fmt_percent(columns = prop_positive, decimals = 1) %>%
  gt::tab_header(title = "Summary of Sampling Sites", )


# unit summary
dat_first %>%
  group_by(unit) %>%
  summarise(prop_positive = sum(!is.na(copies_max)) / n()) %>%
  ungroup() %>%
  arrange(desc(prop_positive)) %>%
  gt::gt() %>%
  gt::fmt_percent(columns = prop_positive, decimals = 1) %>%
  gt::tab_header(title = "Summary of Hospital Units")



# subject summary
dat_first %>%
  group_by(redcap_id) %>%
  summarise(prop_positive = sum(!is.na(copies_max)) / n()) %>%
  ungroup() %>%
  arrange(desc(prop_positive)) %>%
  gt::gt() %>%
  gt::fmt_percent(columns = prop_positive, decimals = 1) %>%
  gt::tab_header(title = "Summary of Hospital Units")



# time summary
dat_long %>%
  group_by(subject_study_day, swab_site) %>%
  summarise(prop_positive = sum(!is.na(copies_max)) / n()) %>%
  ungroup() %>%
  qplot(data = ., x = subject_study_day, y = prop_positive, facets = ~ swab_site, geom = "point") +
  labs(title = "SARS-CoV-2 Detection vs Days Post Enrollment") +
  theme(plot.title.position = "plot")


dat_long %>%
  group_by(subject_hosp_day, swab_site) %>%
  summarise(prop_positive = sum(!is.na(copies_max)) / n()) %>%
  ungroup() %>%
  qplot(data = ., x = subject_hosp_day, y = prop_positive, facets = ~ swab_site, geom = "point") +
  labs(title = "SARS-CoV-2 Detection vs Days Post Hospitalization") +
  theme(plot.title.position = "plot")


dat_long %>%
  group_by(subject_covid_day, swab_site) %>%
  summarise(prop_positive = sum(!is.na(copies_max)) / n()) %>%
  ungroup() %>%
  qplot(data = ., x = subject_covid_day, y = prop_positive, facets = ~ swab_site, geom = "point") +
  labs(title = "SARS-CoV-2 Detection vs Days Post COVID-19 Diagnosis") +
  theme(plot.title.position = "plot")
