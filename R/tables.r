library(gt)
library(gtsummary)
library(here)
library(rstpm2)
library(survival)

source("R/data.r")

output_dir <- here("output", "tables")
dir.create(output_dir, showWarnings = FALSE)

baseline_data |>
    select(
        sex, age_onset, diagnostic_delay_yrs,
        site_of_onset, final_diagnosis,
        eer_category, gold_coast_crit_yn,
        baseline_deltafs, progression_category,
        treated_riluzole
    ) |>
    tbl_summary(missing = "no") |>
    as_gt() |>
    gtsave(file.path(output_dir, "summary.docx"))

baseline_genetics |>
    select(ends_with("_status")) |>
    tbl_summary() |>
    as_gt() |>
    gtsave(file.path(output_dir, "genetics.docx"))

baseline_data |>
    select(
        c9_status, age_onset, diagnostic_delay_yrs,
        site_of_onset, final_diagnosis,
        eer_category, gold_coast_crit_yn,
        baseline_deltafs, progression_category
    ) |>
    tbl_summary(by = "c9_status", missing = "no") |>
    add_p() |>
    as_gt() |>
    gtsave(file.path(output_dir, "summary-by-c9_status.docx"))

baseline_data |>
    select(
        sod1_status, age_onset, diagnostic_delay_yrs,
        site_of_onset, final_diagnosis,
        eer_category, gold_coast_crit_yn,
        baseline_deltafs, progression_category,
        treated_tofersen
    ) |>
    tbl_summary(by = "sod1_status", missing = "no") |>
    add_p() |>
    as_gt() |>
    gtsave(file.path(output_dir, "summary-by-sod1_status.docx"))

baseline_data |>
    select(
        tardbp_status, age_onset, diagnostic_delay_yrs,
        site_of_onset, final_diagnosis,
        eer_category, gold_coast_crit_yn,
        baseline_deltafs, progression_category
    ) |>
    tbl_summary(by = "tardbp_status", missing = "no") |>
    add_p() |>
    as_gt() |>
    gtsave(file.path(output_dir, "summary-by-tardbp_status.docx"))

baseline_data |>
    select(
        fus_status, age_onset, diagnostic_delay_yrs,
        site_of_onset, final_diagnosis,
        eer_category, gold_coast_crit_yn,
        baseline_deltafs, progression_category
    ) |>
    tbl_summary(by = "fus_status", missing = "no") |>
    add_p() |>
    as_gt() |>
    gtsave(file.path(output_dir, "summary-by-fus_status.docx"))

survival_cox.fit <- coxph(
  Surv(time, status == "Deceased") ~
    sex + age_onset + log(diagnostic_delay_yrs) +
    I(baseline_deltafs^(1/3)) +
    site_of_onset + c9_status,
  data = baseline_data |>
    filter(age > age_onset, diagnostic_delay_yrs > 0) |>
    mutate(time = age - age_onset),
)

tbl_regression(survival_cox.fit, exponentiate = TRUE) |>
  as_gt() |>
  gtsave(file.path(output_dir, "survival_coxph.docx"))

tofersen_cox.fit <- coxph(
  Surv(age-age_onset, status == "Deceased") ~
    age_onset + progression_category + treated_tofersen,
  data = baseline_data |>
    filter(c9_status == "Negative", sod1_status == "Positive")
)

tbl_regression(tofersen_cox.fit, exponentiate = TRUE) |>
  as_gt() |>
  gtsave(file.path(output_dir, "tofersen_coxph.docx"))
