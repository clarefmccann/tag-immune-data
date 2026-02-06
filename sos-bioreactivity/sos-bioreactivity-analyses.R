pacman::p_load("dplyr", "ggplot2", "tidyr", "stringr", "purrr", "here", "lme4", "lmerTest", "broom", "broom.mixed")

proj_path <- here()
data_path <- Sys.getenv("CAS_DIR")

cytokines <- read.csv(paste0(data_path, "/TAG/data/SOS/saliva/Cleaned_datasets/SOS_cytokines_cleaned.csv"))
crp <- read.csv(paste0(data_path, "/TAG/data/SOS/saliva/Cleaned_datasets/SOS_CRP_cleaned.csv"))
hormones <- read.csv(paste0(data_path, "/TAG/data/SOS/saliva/Cleaned_datasets/SOS_hormones_cleaned.csv"))

str(cytokines)
str(crp)
str(hormones)

full_dataset <- left_join(cytokines, left_join(crp, hormones, by = c("SampleID", "Time"))) %>% 
  select(-ends_with(c(".x", ".y"))) 

selected_dataset <- full_dataset %>% 
  select(SampleID, Time, ends_with("log"))

paired_t_one_tailed <- function(x, y, alternative = c("greater", "less")) {
  alternative <- match.arg(alternative)
  tt <- t.test(x, y, paired = TRUE, alternative = alternative)
  tibble(
    estimate = unname(tt$estimate),
    statistic = unname(tt$statistic),
    df = unname(tt$parameter),
    p_value = tt$p.value,
    conf_low = tt$conf.int[1],
    conf_high = tt$conf.int[2]
  )
}

paired_wilcox <- function(x, y, alternative = c("two.sided", "greater", "less")) {
  alternative <- match.arg(alternative)
  wt <- suppressWarnings(wilcox.test(x, y, paired = TRUE, exact = FALSE, alternative = alternative))
  tibble(
    statistic = unname(wt$statistic),
    p_value = wt$p.value
  )
}

pairwise_cor_fdr <- function(df, method = "pearson") {
  vars <- names(df)
  pairs <- combn(vars, 2, simplify = FALSE)
  
  out <- map_dfr(pairs, function(p) {
    a <- p[1]; b <- p[2]
    ct <- cor.test(df[[a]], df[[b]], method = method)
    tibble(var1 = a, var2 = b,
           r = unname(ct$estimate),
           p_value = ct$p.value)
  })
  
  out %>% mutate(p_fdr = p.adjust(p_value, method = "fdr"))
}


# -----------------------------
# hypothesis 1a
# -----------------------------
# 1a(i) linear + quadratic effect of timepoint on each biomarker
# random intercept + random slope of time (linear)
biomarkers_long <- selected_dataset %>%
  mutate(
    Time_2 = Time^2
  ) %>% 
  pivot_longer(cols = 3:10, names_to = "biomarker", values_to = "value", )

fit_biomarker_models <- function(df_one_bio) {

  m <- lmer(value ~ Time + Time_2 + (1 + Time | SampleID), data = biomarkers_long, REML = FALSE)
  sm <- summary(m)
  
  fe <- broom.mixed::tidy(m, effects = "fixed") %>%
    select(term, estimate, std.error, statistic, p.value)
  
  list(model = m, fixed = fe)
}

biomarker_model_results <- biomarkers_long %>%
  group_by(biomarker) %>%
  group_split() %>%
  set_names(map_chr(., ~ unique(.x$biomarker))) %>%
  map(fit_biomarker_models)

biomarker_fixed_table <- imap_dfr(biomarker_model_results, function(res, bio) {
  res$fixed %>% mutate(biomarker = bio, .before = 1)
}) %>%
  group_by(term) %>% 
  mutate(p_fdr = p.adjust(p.value, method = "fdr")) %>%
  ungroup()

# 1a(ii) generate peak levels (max of T3/T4) and compare to T1 with one-tailed paired t-tests
# make peak among T3,T4 per subject per biomarker; keep T1
bio_peaks <- biomarkers_long %>%
  filter(Time %in% c(1, 3, 4)) %>%
  group_by(SampleID, biomarker) %>%
  summarise(
    t1 = value[Time == 1][1],
    t3 = value[Time == 3][1],
    t4 = value[Time == 4][1],
    peak = max(c(t3, t4), na.rm = TRUE),
    peak_tp = if_else(!is.na(t3) & !is.na(t4) & t3 >= t4, 3,
                      if_else(is.na(t3) & !is.na(t4), 4,
                              if_else(!is.na(t3) & is.na(t4), 3, 4))),
    .groups = "drop"
  ) %>%
  mutate(delta_peak_t1 = peak - t1)

# one-tailed paired t-tests: peak > t1
bio_peak_tests <- bio_peaks %>%
  filter(!is.na(t1) & !is.na(peak)) %>%
  group_by(biomarker) %>%
  summarise(
    paired_t_one_tailed(peak, t1, alternative = "greater"),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(p_fdr = p.adjust(p_value, method = "fdr"))

# outputs:
# - bio_peaks has peak levels + magnitude of change (delta_peak_t1)
# - bio_peak_tests has one-tailed p + FDR


# -----------------------------
# hypothesis 1b
# -----------------------------

# rest average across two resting scans
hr_wide <- hr_wide %>%
  mutate(
    hr_rest_avg = rowMeans(across(c(hr_rest1, hr_rest2)), na.rm = TRUE),
    delta_task_rest = hr_task - hr_rest_avg
  )

# 1b(i) task HR > rest HR (one-tailed paired t-test)
hr_task_vs_rest <- paired_t_one_tailed(hr_wide$hr_task, hr_wide$hr_rest_avg, alternative = "greater") %>%
  mutate(test = "task_vs_rest")

# 1b(ii) self-connected vs average of other three conditions (one-tailed paired t-test)
hr_wide <- hr_wide %>%
  mutate(
    hr_other_avg = rowMeans(across(c(hr_disconnected, hr_other1, hr_other2)), na.rm = TRUE),
    delta_connected_vs_disconnected = hr_connected - hr_disconnected,
    delta_connected_vs_otheravg = hr_connected - hr_other_avg,
    delta_connected_vs_rest = hr_connected - hr_rest_avg
  )

hr_connected_vs_others <- paired_t_one_tailed(hr_wide$hr_connected, hr_wide$hr_other_avg, alternative = "greater") %>%
  mutate(test = "connected_vs_otheravg")

hr_connected_vs_disconnected <- paired_t_one_tailed(hr_wide$hr_connected, hr_wide$hr_disconnected, alternative = "greater") %>%
  mutate(test = "connected_vs_disconnected")

hr_tests <- bind_rows(hr_task_vs_rest, hr_connected_vs_others, hr_connected_vs_disconnected) %>%
  mutate(p_fdr = p.adjust(p_value, method = "fdr"))


# -----------------------------
# hypothesis 1c
# -----------------------------

# wilcoxon signed-rank: connected vs disconnected and connected vs other

wilcox_alt <- "greater"  

ratings_wide <- ratings_wide %>%
  mutate(
    delta_rating_conn_vs_disc = rating_connected - rating_disconnected,
    delta_rating_conn_vs_other = rating_connected - rating_other
  )

rating_conn_vs_disc <- paired_wilcox(ratings_wide$rating_connected, ratings_wide$rating_disconnected, alternative = wilcox_alt) %>%
  mutate(test = "rating_connected_vs_disconnected")

rating_conn_vs_other <- paired_wilcox(ratings_wide$rating_connected, ratings_wide$rating_other, alternative = wilcox_alt) %>%
  mutate(test = "rating_connected_vs_other")

rating_tests <- bind_rows(rating_conn_vs_disc, rating_conn_vs_other) %>%
  mutate(p_fdr = p.adjust(p_value, method = "fdr"))


# -----------------------------
# hypothesis 1d
# -----------------------------

# build one dataframe with ALL magnitude-of-change (difference score) variables
# include biomarker deltas + HR deltas + rating deltas
diff_scores <- bio_peaks %>%
  select(id, biomarker, delta_peak_t1) %>%
  pivot_wider(names_from = biomarker, values_from = delta_peak_t1, names_prefix = "delta_bio_") %>%
  left_join(
    hr_wide %>%
      select(
        id,
        delta_task_rest,
        delta_connected_vs_disconnected,
        delta_connected_vs_otheravg,
        delta_connected_vs_rest
      ),
    by = "id"
  ) %>%
  left_join(
    ratings_wide %>%
      select(id, delta_rating_conn_vs_disc, delta_rating_conn_vs_other),
    by = "id"
  )

# pearson correlations among all diff scores (pairwise complete obs)
cor_results <- pairwise_cor_fdr(diff_scores %>% select(-id), method = "pearson")

# -----------------------------
# to export
# -----------------------------
# biomarker_fixed_table
# bio_peaks
# bio_peak_tests
# hr_tests
# rating_tests
# diff_scores
# cor_results
