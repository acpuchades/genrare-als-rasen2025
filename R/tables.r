library(gt)
library(gtsummary)
library(here)

source("R/data.r")

output_dir <- here("output")
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
