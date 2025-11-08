library(dplyr)
library(here)
library(janitor)
library(readr)
library(stringr)

alsfrs_data_path <- here("data", "02_sen2025_alsfrs.csv")
cognitive_data_path <- here("data", "06_sen2025_cognitive.csv")
demographics_data_path <- here("data", "01_sen2025_basico_mod.csv")
genetics_data_path <- here("data", "10_sen2025_genet.csv")
kings_data_path <- here("data", "04_sen2025_kings.csv")
niv_data_path <- here("data", "11_sen2025_ventilation.csv")
treatment_data_path <- here("data", "08_sen2025_trmnt.csv")

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
        factor(levels = c(
                   "Definite ALS",
                   "Probable ALS",
                   "Probable ALS, lab supported",
                   "Possible ALS",
                   "Not meeting criteria"
               ))
}

as_genetic_result <- \(x) factor(x, levels = c("Positive", "Negative"))
as_sex <- \(x) factor(x, levels = c("Male", "Female"))
as_vital_status <- \(x) factor(x, levels = c("Alive", "Deceased"))
as_yesno_category <- \(x) factor(x, levels = c("Sí", "No"))

demographics_data <- read_csv(demographics_data_path) |>
    clean_names() |>
    mutate(
        across(sex, as_sex),
        across(status, as_vital_status),
        across(eer_category, as_eer_category),
        across(final_diagnosis, as_als_diagnosis),
        across(gold_coast_crit_yn, as_yesno_category),
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
        ) |> factor(levels = c("NP", "SP", "FP"))
    )

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
        .by = record_id
    )

baseline_data <- demographics_data |>
    left_join(baseline_alsfrs, by = "record_id") |>
    left_join(baseline_genetics, by = "record_id") |>
    left_join(baseline_treatment, by = "record_id")
