library(dplyr)
library(forcats)
library(ggplot2)
library(ggsurvfit)
library(stringr)
library(survival)
library(tidyr)

source("R/data.r")

theme_color_primary <- "#143771"
theme_color_secondary <- "#4E5F68"

output_dir <- here("output", "figures")
dir.create(output_dir, showWarnings = FALSE)

baseline_data |>
    select(record_id) |>
    mutate(
        alsfrs_available = record_id %in% alsfrs_data$record_id,
        cognitive_available = record_id %in% cognitive_data$record_id,
        genetics_available = record_id %in% genetics_data$record_id,
        kings_available = record_id %in% kings_data$record_id,
        niv_available = record_id %in% niv_data$record_id
    ) |>
    summarise(
        across(ends_with("_available"), \(x) sum(x, na.rm = TRUE)),
        patients_available = n(),
        total_available = n()
    ) |>
    pivot_longer(ends_with("_available") & -total_available) |>
    mutate(
        across(name, ~.x |>
                         fct_reorder(desc(value)) |>
                         fct_relevel("patients_available") |>
                         fct_rev()
    )) |>
    ggplot(aes(name, value)) +
    geom_col(fill = theme_color_primary) +
    geom_text(
        aes(y = value/2, label = str_glue("{round(value / total_available * 100)} %")),
        color = "white", hjust = 0.5, size = 4
    ) +
    coord_flip() +
    labs(y = "Nº de registros", x = NULL) +
    theme_minimal()
ggsave(
    file.path(output_dir, "data-availability.png"),
    bg = "white", width = 10, height = 5, dpi = 300
)

baseline_data |>
    filter(age >= age_onset) |>
    with(survfit2(Surv(age-age_onset, status == "Deceased") ~ 1)) |>
    ggsurvfit(color = theme_color_primary) +
    add_confidence_interval(fill = theme_color_primary) +
    add_quantile(0.5, color = theme_color_secondary) +
    scale_ggsurvfit()
ggsave(
    file.path(output_dir, "survival.png"),
    width = 10, height = 7, dpi = 300
)

baseline_data |>
    filter(age >= age_onset) |>
    with(survfit2(Surv(age-age_onset, status == "Deceased") ~ site_of_onset)) |>
    ggsurvfit() +
    scale_ggsurvfit() +
    add_confidence_interval() +
    add_legend_title("Categoría de inicio") +
    add_pvalue("annotation") +
    theme_minimal() +
    theme(legend.position = "bottom")
ggsave(
    file.path(output_dir, "survival-by-site_of_onset.png"),
    bg = "white", width = 10, height = 7, dpi = 300
)

baseline_data |>
    filter(age >= age_onset) |>
    with(survfit2(Surv(age-age_onset, status == "Deceased") ~ eer_category)) |>
    ggsurvfit() +
    scale_ggsurvfit() +
    add_confidence_interval() +
    add_legend_title("Categoría EER") +
    add_pvalue("annotation") +
    theme_minimal() +
    theme(legend.position = "bottom")
ggsave(
    file.path(output_dir, "survival-by-eer_category.png"),
    bg = "white", width = 10, height = 7, dpi = 300
)

baseline_data |>
    filter(final_diagnosis != "FOSMN", age >= age_onset) |>
    with(survfit2(Surv(age-age_onset, status == "Deceased") ~ final_diagnosis)) |>
    ggsurvfit() +
    scale_ggsurvfit() +
    add_confidence_interval() +
    add_legend_title("Fenotipo") +
    add_pvalue("annotation") +
    theme_minimal() +
    theme(legend.position = "bottom")
ggsave(
    file.path(output_dir, "survival-by-final_diagnosis.png"),
    bg = "white", width = 10, height = 7, dpi = 300
)

baseline_data |>
    filter(age >= age_onset) |>
    with(survfit2(Surv(age-age_onset, status == "Deceased") ~ progression_category)) |>
    ggsurvfit() +
    scale_ggsurvfit() +
    add_confidence_interval() +
    add_legend_title("Categoría de progresión") +
    add_pvalue("annotation") +
    theme_minimal() +
    theme(legend.position = "bottom")
ggsave(
    file.path(output_dir, "survival-by-progression_category.png"),
    bg = "white", width = 10, height = 7, dpi = 300
)

baseline_data |>
    filter(age >= age_onset) |>
    with(survfit2(Surv(age-age_onset, status == "Deceased") ~ c9_status)) |>
    ggsurvfit() +
    add_confidence_interval() +
    add_legend_title("C9orf72") +
    add_pvalue("annotation") +
    scale_ggsurvfit() +
    scale_color_manual(values = c(theme_color_primary, theme_color_secondary)) +
    scale_fill_manual(values = c(theme_color_primary, theme_color_secondary)) +
    theme_minimal() +
    theme(legend.position = "bottom")
ggsave(
    file.path(output_dir, "survival-by-c9_status.png"),
    bg = "white", width = 10, height = 7, dpi = 300
)

baseline_data |>
    filter(age >= age_onset) |>
    with(survfit2(Surv(age-age_onset, status == "Deceased") ~ sod1_status)) |>
    ggsurvfit() +
    add_confidence_interval() +
    add_legend_title("SOD1") +
    add_pvalue("annotation") +
    scale_ggsurvfit() +
    scale_color_manual(values = c(theme_color_primary, theme_color_secondary)) +
    scale_fill_manual(values = c(theme_color_primary, theme_color_secondary)) +
    theme_minimal() +
    theme(legend.position = "bottom")
ggsave(
    file.path(output_dir, "survival-by-sod1_status.png"),
    bg = "white", width = 10, height = 7, dpi = 300
)

genetics_data |>
    drop_na(gen_gene) |>
    mutate(across(gen_gene,
                  ~.x |>
                      fct_infreq() |>
                      fct_lump_min(5) |>
                      fct_rev()
    )) |>
    ggplot(aes(gen_gene)) +
    geom_bar(fill = theme_color_primary) +
    geom_text(
        aes(y = after_stat(count/2), label = after_stat(count)),
        hjust = 0.5, color = "white", stat = "count"
    ) +
    coord_flip() +
    labs(x = "Gen", y = "Pacientes") +
    theme_minimal()
ggsave(
    file.path(output_dir, "altered-genes.png"),
    bg = "white", width = 10, height = 7, dpi = 300
)

alsfrs_data |>
    filter(years_since_onset <= 10) |>
    ggplot(aes(years_since_onset, alsfrs_total, group = record_id)) +
    geom_line(linewidth = 0.1, color = theme_color_secondary) +
    geom_jitter(color = theme_color_secondary, size = 0.1, alpha = 0.5) +
    scale_x_continuous(breaks = seq(0, 10, 2), minor_breaks = seq(0, 10, 1)) +
    scale_y_continuous(breaks = seq(0, 48, 6), minor_breaks = seq(0, 48, 3)) +
    labs(x = "Años desde el inicio", y = "ALSFRS-R") +
    theme_minimal()
ggsave(
    file.path(output_dir, "alsfrs.png"),
    bg = "white", width = 10, height = 6, dpi = 300
)

alsfrs_data |>
    filter(years_since_onset <= 10) |>
    left_join(baseline_data, by = "record_id") |>
    ggplot(aes(
        years_since_onset, alsfrs_total, group = record_id,
        color = progression_category, fill = progression_category
    )) +
    geom_line(linewidth = 0.1) + geom_jitter(size = 0.1, alpha = 0.5) +
    geom_smooth(aes(group = NULL, fill = progression_category)) +
    scale_x_continuous(breaks = seq(0, 10, 2), minor_breaks = seq(0, 10, 1)) +
    scale_y_continuous(breaks = seq(0, 48, 6), minor_breaks = seq(0, 48, 3)) +
    labs(
        x = "Años desde el inicio", y = "ALSFRS-R",
        color = "Tipo de progresión", fill = "Tipo de progresión"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
ggsave(
    file.path(output_dir, "alsfrs-by-progression_category.png"),
    bg = "white", width = 10, height = 6, dpi = 300
)

alsfrs_data |>
    filter(years_since_onset <= 10) |>
    left_join(baseline_data, by = "record_id") |>
    ggplot(aes(
        years_since_onset, alsfrs_total, group = record_id,
        color = c9_status, fill = c9_status
    )) +
    geom_line(linewidth = 0.1) + geom_jitter(size = 0.1, alpha = 0.5) +
    geom_smooth(aes(group = NULL, fill = c9_status)) +
    scale_x_continuous(breaks = seq(0, 10, 2), minor_breaks = seq(0, 10, 1)) +
    scale_y_continuous(breaks = seq(0, 48, 6), minor_breaks = seq(0, 48, 3)) +
    scale_color_manual(
        values = c(
            Positive = alpha(theme_color_primary, 1),
            Negative = alpha(theme_color_secondary, 0.25)
        )
    ) +
    scale_fill_manual(
        values = c(
            Positive = alpha(theme_color_primary, 1),
            Negative = alpha(theme_color_secondary, 0.25)
        )
    ) +
    labs(
        x = "Años desde el inicio", y = "ALSFRS-R",
        color = "C9orf72", fill = "C9orf72"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
ggsave(
    file.path(output_dir, "alsfrs-by-c9_status.png"),
    bg = "white", width = 10, height = 6, dpi = 300
)

alsfrs_data |>
    filter(years_since_onset <= 10) |>
    left_join(baseline_data, by = "record_id") |>
    ggplot(aes(
        years_since_onset, alsfrs_total, group = record_id,
        color = sod1_status, fill = sod1_status
    )) +
    geom_line(linewidth = 0.1) + geom_jitter(size = 0.1, alpha = 0.5) +
    geom_smooth(aes(group = NULL, fill = sod1_status)) +
    scale_x_continuous(breaks = seq(0, 10, 2), minor_breaks = seq(0, 10, 1)) +
    scale_y_continuous(breaks = seq(0, 48, 6), minor_breaks = seq(0, 48, 3)) +
    scale_color_manual(
        values = c(
            Positive = alpha(theme_color_primary, 1),
            Negative = alpha(theme_color_secondary, 0.25)
        )
    ) +
    scale_fill_manual(
        values = c(
            Positive = alpha(theme_color_primary, 1),
            Negative = alpha(theme_color_secondary, 0.25)
        )
    ) +
    labs(
        x = "Años desde el inicio", y = "ALSFRS-R",
        color = "SOD1", fill = "SOD1"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
ggsave(
    file.path(output_dir, "alsfrs-by-sod1_status.png"),
    bg = "white", width = 10, height = 6, dpi = 300
)
