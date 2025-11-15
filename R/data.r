library(dplyr)
library(here)
library(janitor)
library(readr)
library(stringr)
library(tidyr)

data_dir <- here("data")
alsfrs_data_path <- file.path(data_dir, "02_sen2025_alsfrs.csv")
cognitive_data_path <- file.path(data_dir, "06_sen2025_cognitive.csv")
demographics_data_path <- file.path(data_dir, "01_sen2025_basico.csv")
genetics_data_path <- file.path(data_dir, "11_sen2025_genet.csv")
kings_data_path <- file.path(data_dir, "04_sen2025_kings.csv")
niv_data_path <- file.path(data_dir, "12_sen2025_ventilation.csv")
treatment_data_path <- file.path(data_dir, "08_sen2025_trmnt.csv")

as_als_diagnosis <- function(x) {
    x |>
        str_remove(" \\([^)]*\\)") |>
        factor(levels = c(
                   "Classical ALS",
                   "Progressive Muscular Atrophy",
                   "Primary Lateral Sclerosis",
                   "FOSMN"
               ))
}

as_eer_category <- function(x) {
    x |>
        str_remove(" \\([^)]*\\)") |>
        str_replace("laboratory results supported", "lab supported") |>
        factor(
            levels = c(
                "Definite ALS",
                "Probable ALS",
                "Probable ALS, lab supported",
                "Possible ALS",
                "Not meeting criteria"
            )
        )
}

as_family_history_category <- function(x) {
  factor(x, levels = c("Sporadic ALS", "Familial ALS"))
}

as_genetic_result <- function(x) {
  factor(x, levels = c("Negative", "Positive"))
}

as_sex <- function(x) {
  factor(x, levels = c("Male", "Female"))
}

as_onset_category <- function(x) {
  factor(x, levels = c("Spinal", "Bulbar", "Respiratory", "Generalized"))
}

as_vital_status <- function(x) {
  factor(x, levels = c("Alive", "Deceased"))
}

as_yesno_category <- function(x) {
  case_when(
    x %in% c("Sí", "Yes", TRUE) ~ "Sí",
    x %in% c("No", FALSE) ~ "No"
  ) |>
    factor(levels = c("No", "Sí"))
}

demographics_data <- read_csv(demographics_data_path) |>
    clean_names() |>
    mutate(
        across(sex, as_sex),
        across(status, as_vital_status),
        across(als_fam_hist_type, as_family_history_category),
        across(eer_category, as_eer_category),
        across(final_diagnosis, as_als_diagnosis),
        across(c(ends_with("_yn"), ends_with("_dic")), as_yesno_category),
        site_of_onset = as_onset_category(case_when(
        (
            !is.na(als_symp_bulbar) +
            !is.na(als_symp_spinal) +
            !is.na(als_symp_resp)
        ) > 1 ~ "Generalized",
            !is.na(als_symp_bulbar) ~ "Bulbar",
            !is.na(als_symp_spinal) ~ "Spinal",
            !is.na(als_symp_resp) ~ "Respiratory"
        )),
        diagnostic_delay_yrs = age_diagnosis - age_onset
    )

alsfrs_data <- read_csv(alsfrs_data_path) |>
    clean_names() |>
    left_join(
        demographics_data |> select(record_id, age_onset),
        by = "record_id"
    ) |>
    mutate(
        years_since_onset = round(alsfrs_age - age_onset, 1),
        delta_fs = if_else(
            is.finite((48-alsfrs_total)/(years_since_onset*12)),
            (48-alsfrs_total)/(years_since_onset*12), NA
        )
    )

cognitive_data <- read_csv(cognitive_data_path) |>
    clean_names()

genetics_data <- read_csv(genetics_data_path) |>
    clean_names()

kings_data <- read_csv(kings_data_path) |>
    clean_names()

niv_data <- read_csv(niv_data_path) |>
    clean_names()

treatment_data <- read_csv(treatment_data_path) |>
    clean_names()

baseline_alsfrs <- alsfrs_data |>
    slice_min(alsfrs_age, by = record_id, n = 1, with_ties = FALSE) |>
    transmute(
        record_id,
        baseline_deltafs = delta_fs,
        progression_category = case_when(
            baseline_deltafs < 0.5 ~ "SP",
            baseline_deltafs |> between(0.5, 1) ~ "NP",
            baseline_deltafs > 1 ~ "FP"
        ) |> factor(levels = c("SP", "NP", "FP"))
    )

baseline_ecas <- cognitive_data |>
  drop_na(ecas_als_specific, ecas_als_non_specific, ecas_total) |>
  slice_min(cognitive_age, by = record_id, n = 1, with_ties = FALSE) |>
  select(record_id, ecas_als_specific, ecas_als_non_specific, ecas_total)

baseline_genetics <- genetics_data |>
    summarize(
        c9_status = if_else(
            any(gen_gene == "C9orf72", na.rm = TRUE),
            "Positive", "Negative"
        ),
        sod1_status = if_else(
            any(gen_gene == "SOD1", na.rm = TRUE),
            "Positive", "Negative"
        ),
        tardbp_status = if_else(
            any(gen_gene == "TARDBP", na.rm = TRUE),
            "Positive", "Negative"
        ),
        fus_status = if_else(
            any(gen_gene == "FUS", na.rm = TRUE),
            "Positive", "Negative"
        ),
        other_status = if_else(
            any(!is.na(gen_gene) & !(gen_gene %in% c("C9orf72", "SOD1", "TARDBP", "FUS"))),
            "Positive", "Negative"
        ),
        .by = record_id
    ) |>
    mutate(
        across(ends_with("_status"), as_genetic_result),
    )

baseline_treatment <- treatment_data |>
  summarize(
        treated_riluzole = any(treat_type == "Riluzole", na.rm = TRUE),
        treated_tofersen = any(treat_type == "Tofersen", na.rm = TRUE),
        riluzole_start = suppressWarnings(if_else(
          treated_riluzole,
          min(med_start_age[which(treat_type == "Riluzole")]),
          NA
        )),
        tofersen_start = suppressWarnings(if_else(
          treated_tofersen,
          min(med_start_age[which(treat_type == "Tofersen")]),
          NA
        )),
        .by = record_id
    )

baseline_data <- demographics_data |>
  left_join(baseline_alsfrs, by = "record_id") |>
  left_join(baseline_ecas, by = "record_id") |>
  left_join(baseline_genetics, by = "record_id") |>
  left_join(baseline_treatment, by = "record_id") |>
  mutate(
    riluzole_delay = if_else(
      is.finite(riluzole_start-age_onset),
      riluzole_start - age_onset, NA
    ),
    tofersen_delay = if_else(
      is.finite(tofersen_start-age_onset),
      tofersen_start - age_onset, NA
    )
  )
