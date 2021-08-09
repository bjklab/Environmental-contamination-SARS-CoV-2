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
  # arrange(redcap_id, subject_study_day) %>%
  # mutate(subject_room_id = as.numeric(factor(paste0("decon_", redcap_id, "_", room))),
  #        subject_room_id = paste0("decon_",stringr::str_pad(subject_room_id, width = 3, side = "left", pad = "0"))) %>%
  # select(subject_room_id, everything()) %>%
  # group_by(subject_room_id) %>%
  # arrange(subject_study_day) %>%
  # mutate(subject_room_day = subject_study_day - min(subject_study_day, na.rm = TRUE)) %>%
  # ungroup() %>%
  # select(subject_room_id, redcap_id, unit, room, unit_room, contains("day"), everything()) %>%
  identity() -> dat_long
dat_long


dat_long %>%
  group_by(subject_room_id) %>%
  filter(subject_room_day == min(subject_room_day)) %>%
  ungroup() %>%
  identity() -> dat_first
dat_first


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



# subject - room summary
dat_first %>%
  group_by(subject_room_id) %>%
  summarise(prop_positive = sum(!is.na(copies_max)) / n()) %>%
  ungroup() %>%
  arrange(desc(prop_positive)) %>%
  gt::gt() %>%
  gt::fmt_percent(columns = prop_positive, decimals = 1) %>%
  gt::tab_header(title = "Summary of Subjects")


# subject & unit summary
dat_first %>%
  group_by(subject_room_id, unit) %>%
  summarise(prop_positive = sum(!is.na(copies_max)) / n()) %>%
  ungroup() %>%
  group_by(unit) %>%
  arrange(desc(prop_positive)) %>%
  gt::gt() %>%
  gt::fmt_percent(columns = prop_positive, decimals = 1) %>%
  gt::tab_header(title = "Summary of Subjects & Units")



# time summary
dat_long %>%
  group_by(subject_room_day, swab_site) %>%
  summarise(prop_positive = sum(!is.na(copies_max)) / n()) %>%
  ungroup() %>%
  qplot(data = ., x = subject_room_day, y = prop_positive, facets = ~ swab_site, geom = "point") +
  labs(title = "SARS-CoV-2 Detection vs Days Post Enrollment") +
  theme(plot.title.position = "plot")


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



# distance summary
dat_first %>%
  qplot(data = ., x = distance, y = copies_max, geom = "point") +
  facet_wrap(facets = ~ swab_site, scales = "free_x") +
  scale_y_log10() +
  labs(title = "SARS-CoV-2 Abundance vs Distance from Patient's Bed") +
  theme(plot.title.position = "plot")
