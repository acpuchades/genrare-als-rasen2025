library(dplyr)
library(forcats)
library(ggplot2)
library(ggsci)
library(ggsignif)
library(ggsurvfit)
library(stringr)
library(survival)
library(tidyr)
library(visdat)

source("R/data.r")

theme_color_primary <- "#143771"
theme_color_secondary <- "#4E5F68"

custom_scale_fill_brewer <- scale_fill_lancet
custom_scale_color_brewer <- scale_color_lancet

output_dir <- here("output", "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

vis_miss(baseline_data) +
  theme(
    plot.margin = margin(10, 65, 10, 10),
  )
ggsave(
  file.path(output_dir, "data-missing.png"),
  bg = "white", dpi = 300
)

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

genetics_data |>
  summarize(
    gen_result = as_genetic_result(if_else(
      any(!is.na(gen_gene)), "Positive", "Negative"
    )),
    .by = record_id
  ) |>
  left_join(baseline_data, by = "record_id") |>
  ggplot(aes(gen_result, age_onset, fill = gen_result)) +
  geom_boxplot() +
  geom_signif(
    comparisons = list(c("Positive", "Negative")),
    map_signif_level = TRUE
  ) +
  scale_fill_manual(values = c(
    Positive = theme_color_primary,
    Negative = theme_color_secondary
  )) +
  labs(
    x = NULL, y = "Edad de inicio",
    fill = "Resultado genético"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave(
  file.path(output_dir, "age_onset-by-gen_result.png"),
  bg = "white", width = 10, height = 8, dpi = 300
)

genetics_data |>
  summarize(
    gen_result = as_genetic_result(if_else(
      any(!is.na(gen_gene)), "Positive", "Negative"
    )),
    .by = record_id
  ) |>
  left_join(baseline_data, by = "record_id") |>
  mutate(age_group = age_onset |> cut(
    breaks = c(0, 18, 30, 40, 50, 60, Inf),
    labels = c(
      "0 - 18", "19 - 30", "31 - 40",
      "41 - 50", "51 - 60", "> 60"
    )
  )) |>
  drop_na(age_group, gen_result) |>
  ggplot(aes(age_group, fill = gen_result)) +
  geom_bar(position = "fill") +
  geom_text(
    aes(
      label = after_stat(count / tapply(count, x, sum)[x]) |>
        scales::percent(accuracy = 1)
    ),
    stat = "count", position = position_fill(vjust = 0.5),
    color = "white", subset = ~gen_result == "Positive"
  ) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c(
    Positive = theme_color_primary,
    Negative = theme_color_secondary
  )) +
  labs(
    x = "Edad de inicio",
    fill = "Resultado genético"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave(
  file.path(output_dir, "gen_positive-by-age_group.png"),
  bg = "white", width = 10, height = 8, dpi = 300
)

bind_rows(
  baseline_data |> mutate(gene = "C9orf72", gene_status = c9_status),
  baseline_data |> mutate(gene = "SOD1", gene_status = sod1_status),
  baseline_data |> mutate(gene = "TARDBP", gene_status = tardbp_status),
  baseline_data |> mutate(gene = "FUS", gene_status = fus_status)
) |>
  ggplot(aes(gene, age_onset, fill = gene_status)) +
  geom_boxplot() +
  scale_fill_manual(values = c(
    Positive = theme_color_primary,
    Negative = theme_color_secondary
  )) +
  labs(x = "Gen", y = "Edad de inicio", fill = "Resultado") +
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave(
  file.path(output_dir, "age_at_onset-by-gene.png"),
  bg = "white", width = 10, height = 8, dpi = 300
)

ggplot(baseline_data, aes(diagnostic_delay_yrs)) +
  geom_histogram(color = "black", fill = theme_color_primary) +
  geom_vline(
    xintercept = median(baseline_data$diagnostic_delay_yrs, na.rm = TRUE),
    color = "red", linetype = 2
  ) +
  scale_x_log10() +
  labs(title = "Retraso diagnóstico", x = "Años", y = "N") +
  theme_minimal()
ggsave(
  file.path(output_dir, "distribution-diagnostic_delay.png"),
  bg = "white", dpi = 300, width = 10, height = 10
)


baseline_data |>
  drop_na(sex) |>
  ggplot(aes(x = "", fill = sex)) +
  geom_bar(color = "white", width = 0.5) +
  geom_text(
    aes(x = 1, label = after_stat(str_glue("{fill}\n{count} ({round(count/sum(count) * 100)}%)"))),
    color = "white", size = 6, stat = "count", hjust = 0.5,
    position = position_stack(vjust = 0.5)
  ) +
  coord_polar("y", start = 0) +
  custom_scale_fill_brewer() +
  labs(fill = "Sexo") +
  theme_void() +
  theme(legend.position = "none")
ggsave(
  file.path(output_dir, "distribution-sex.png"),
  bg = "white", width = 10, height = 10, dpi = 300
)

baseline_data |>
  drop_na(als_fam_hist_type) |>
  ggplot(aes(x = "", fill = als_fam_hist_type)) +
  geom_bar(color = "white", width = 0.5) +
  geom_text(
    aes(x = 1, label = after_stat(str_glue("{fill}\n{count} ({round(count/sum(count) * 100)}%)"))),
    color = "white", size = 6, stat = "count", hjust = 0.5,
    position = position_stack(vjust = 0.5), vjust = 0.5
  ) +
  coord_polar("y", start = 0, clip = "off") +
  custom_scale_fill_brewer() +
  labs(fill = "Historia familiar") +
  theme_void() +
  theme(
    legend.position = "none",
    plot.margin = unit(c(1, 3, 1, 3), "cm")
  )
ggsave(
  file.path(output_dir, "distribution-als_fam_hist_type.png"),
  bg = "white", width = 10, height = 8, dpi = 300
)

baseline_data |>
  drop_na(site_of_onset) |>
  ggplot(aes(x = "", fill = site_of_onset)) +
  geom_bar(color = "white", width = 0.5) +
  geom_text(
    aes(x = 1, label = after_stat(str_glue("{fill}\n{count} ({round(count/sum(count) * 100)}%)"))),
    color = "white", size = 6, stat = "count", hjust = 0.5,
    position = position_stack(vjust = 0.5)
  ) +
  coord_polar("y", start = 0) +
  custom_scale_fill_brewer() +
  labs(fill = "Categoría de inicio") +
  theme_void() +
  theme(legend.position = "none")
ggsave(
  file.path(output_dir, "distribution-site_of_onset.png"),
  bg = "white", width = 10, height = 10, dpi = 300
)

baseline_data |>
  filter(final_diagnosis != "FOSMN") |>
  drop_na(final_diagnosis) |>
  ggplot(aes(x = "", fill = final_diagnosis)) +
  geom_bar(color = "white", width = 0.5) +
  geom_text(
    aes(x = 1.4, label = after_stat(str_glue("{fill}\n{count} ({round(count/sum(count) * 100)}%)"))),
    color = "black", size = 6, stat = "count", hjust = 0,
    position = position_stack(vjust = 0.5),
  ) +
  coord_polar("y", start = 0, clip = "off") +
  custom_scale_fill_brewer() +
  labs(fill = "Fenotipo") +
  theme_void() +
  theme(
    legend.position = "none",
    plot.margin = unit(c(1, 3, 1, 3), "cm")
  )
ggsave(
  file.path(output_dir, "distribution-final_diagnosis.png"),
  bg = "white", width = 10, height = 8, dpi = 300
)

baseline_data |>
  semi_join(
    genetics_data |> drop_na(gen_gene),
    by = "record_id"
  ) |>
  drop_na(als_fam_hist_type) |>
  ggplot(aes(x = "", fill = als_fam_hist_type)) +
  geom_bar(color = "white", width = 1) +
  geom_text(
    aes(x = 1, label = after_stat(str_glue("{fill}\n{count} ({round(count/sum(count) * 100)}%)"))),
    color = "white", size = 5, stat = "count", hjust = 0.5,
    position = position_stack(vjust = 0.5)
  ) +
  coord_polar("y", start = 0) +
  custom_scale_fill_brewer() +
  labs(fill = "Historia familiar") +
  theme_void() +
  theme(legend.position = "none")
ggsave(
  file.path(output_dir, "distribution-gen_gene-als_fam_hist_type.png"),
  bg = "white", width = 10, height = 8, dpi = 300
)

genetics_data |>
  drop_na(gen_gene) |>
  summarize(
    gen_gene = case_when(
      n_distinct(gen_gene) > 1 ~ "Multiple genes",
      first(gen_gene) %in% c("C9orf72", "SOD1", "TARDBP", "FUS") ~ first(gen_gene),
      TRUE ~ "Other genes"
    ) |>
      factor(levels = c(
        "C9orf72", "SOD1", "TARDBP", "FUS", "Other genes", "Multiple genes"
      )),
    .by = record_id
  ) |>
  ggplot(aes(x = "", fill = gen_gene)) +
  geom_bar(color = "white", width = 1) +
  geom_text(
    aes(x = 1.7, label = after_stat(str_glue("{fill}\n{count} ({round(count/sum(count) * 100)}%)"))),
    color = "black", size = 3, stat = "count", hjust = 0.5,
    position = position_stack(vjust = 0.5)
  ) +
  coord_polar("y", start = 0) +
  custom_scale_fill_brewer() +
  labs(fill = "Estudio genético") +
  theme_void() +
  theme(legend.position = "none")
ggsave(
  file.path(output_dir, "distribution-gen_gene.png"),
  bg = "white", width = 10, height = 8, dpi = 300
)

baseline_data |>
  filter(als_fam_hist_type == "Familial ALS") |>
  inner_join(
    genetics_data |> transmute(
      record_id,
      gen_gene = case_when(
        gen_gene %in% c("C9orf72", "SOD1", "TARDBP", "FUS") ~ gen_gene,
        !is.na(gen_gene) ~ "Other genes",
        TRUE ~ "None found"
      ) |> factor(levels = c(
        "C9orf72", "SOD1", "TARDBP", "FUS",
        "Other genes", "None found"
      ))
    ),
    by = "record_id"
  ) |>
  ggplot(aes(x = "", fill = gen_gene)) +
  geom_bar(color = "white", width = 1) +
  geom_text(
    aes(x = 1.7, label = after_stat(str_glue("{fill}\n{count} ({round(count/sum(count) * 100)}%)"))),
    color = "black", size = 4, stat = "count", hjust = 0.5,
    position = position_stack(vjust = 0.5)
  ) +
  coord_polar("y", start = 0) +
  custom_scale_fill_brewer() +
  labs(fill = "Estudio genético") +
  theme_void() +
  theme(legend.position = "none")
ggsave(
  file.path(output_dir, "distribution-fALS-gen_gene.png"),
  bg = "white", width = 10, height = 8, dpi = 300
)

demographics_data |>
  filter(status == "Alive") |>
  drop_na(gast_carrier_dic) |>
  ggplot(aes(x = "", fill = gast_carrier_dic)) +
  geom_bar(color = "white", width = 1) +
  geom_text(
    aes(x = 1.2, label = after_stat(str_glue("{fill}\n{count} ({round(count/sum(count) * 100)}%)"))),
    color = "white", size = 5, stat = "count", hjust = 0.5,
    position = position_stack(vjust = 0.5)
  ) +
  coord_polar("y", start = 0) +
  custom_scale_fill_brewer() +
  labs(fill = "Gastrostomía") +
  theme_void() +
  theme(legend.position = "none")
ggsave(
  file.path(output_dir, "distribution-gast_carrier.png"),
  bg = "white", width = 10, height = 8, dpi = 300
)

kings_data |>
  slice_max(kings_age, n = 1, by = "record_id", with_ties = FALSE) |>
  semi_join(baseline_data |> filter(status == "Alive"), by = "record_id") |>
  mutate(across(kings_total, ~.x |>
                               case_match(
                                 1 ~ "I",
                                 2 ~ "II",
                                 3 ~ "III",
                                 4 ~ "IV"
                               ) |> fct_drop())) |>
  drop_na(kings_total) |>
  ggplot(aes(x = "", fill = kings_total)) +
  geom_bar(color = "white", size = 0.5) +
  geom_text(
    aes(x = 1, label = after_stat(str_glue("{fill}\n{count} ({round(count/sum(count) * 100)}%)"))),
    color = "white", size = 5, stat = "count", hjust = 0.5, position = position_stack(vjust = 0.5),
  ) +
  coord_polar("y", start = 0) +
  labs(fill = "King's") +
  custom_scale_fill_brewer() +
  theme_void() +
  theme(legend.position = "none")
ggsave(
  file.path(output_dir, "estado-actual-by-kings.png"),
  bg = "white", width = 10, height = 10, dpi = 300
)

overall_survival.fit <- baseline_data |>
  filter(age >= age_onset) |>
  with(survfit2(Surv(age-age_onset, status == "Deceased") ~ 1))

overall_survival.med <- quantile(overall_survival.fit, probs = 0.5)$quantile

ggsurvfit(overall_survival.fit, color = theme_color_primary) +
  add_confidence_interval(fill = theme_color_primary) +
  add_quantile(0.5, color = theme_color_secondary) +
  annotate("text", x = overall_survival.med + 0.25, y = 0.5,
           label = str_c(round(overall_survival.med, 1), " años"),
           hjust = 0, vjust = -0.5) +
  scale_ggsurvfit() +
  theme_minimal()
ggsave(
    file.path(output_dir, "survival.png"),
    bg = "white", width = 10, height = 7, dpi = 300
)

baseline_data |>
  filter(age >= age_onset) |>
  with(survfit2(Surv(age-age_onset, status == "Deceased") ~ site_of_onset)) |>
  ggsurvfit() +
  scale_ggsurvfit() +
  custom_scale_color_brewer() +
  custom_scale_fill_brewer() +
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
  custom_scale_color_brewer() +
  custom_scale_fill_brewer() +
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
  custom_scale_color_brewer() +
  custom_scale_fill_brewer() +
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
  filter(final_diagnosis != "FOSMN", (age - age_onset) >= 4) |>
  with(survfit2(Surv(age-age_onset, status == "Deceased") ~ final_diagnosis)) |>
  ggsurvfit() +
  scale_ggsurvfit() +
  scale_x_continuous(limits = c(4, 51)) +
  custom_scale_color_brewer() +
  custom_scale_fill_brewer() +
  add_confidence_interval() +
  add_legend_title("Fenotipo") +
  add_pvalue("annotation") +
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave(
    file.path(output_dir, "survival-by-final_diagnosis-since_4yrs.png"),
    bg = "white", width = 10, height = 7, dpi = 300
)

baseline_data |>
  filter(age >= age_onset) |>
  mutate(across(progression_category, ~factor(.x, levels = c("SP", "FP", "NP")))) |>
  with(survfit2(Surv(age-age_onset, status == "Deceased") ~ progression_category)) |>
  ggsurvfit() +
  scale_ggsurvfit() +
  custom_scale_color_brewer() +
  custom_scale_fill_brewer() +
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
  filter(
    age >= age_onset,
    final_diagnosis == "Classical ALS",
    is.na(sod1_status) | sod1_status == "Negative"
  ) |>
  with(survfit2(Surv(age-age_onset, status == "Deceased") ~ c9_status)) |>
  ggsurvfit() +
  add_confidence_interval() +
  add_legend_title("C9orf72") +
  scale_ggsurvfit() +
  scale_x_continuous(limits = c(0, 10)) +
  scale_color_manual(values = c(theme_color_primary, theme_color_secondary)) +
  scale_fill_manual(values = c(theme_color_primary, theme_color_secondary)) +
  add_pvalue("annotation") +
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave(
    file.path(output_dir, "survival-by-c9_status.png"),
    bg = "white", width = 10, height = 7, dpi = 300
)

baseline_data |>
  filter(
    age >= age_onset,
    final_diagnosis == "Classical ALS",
    is.na(c9_status) | c9_status == "Negative"
  ) |>
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

baseline_data |>
  filter(
    age >= age_onset,
    final_diagnosis == "Classical ALS",
    is.na(c9_status) | c9_status == "Negative"
  ) |>
  mutate(
    group = factor(case_when(
      sod1_status == "Negative" ~ "Non-SOD1",
      sod1_status == "Positive" & !treated_tofersen ~ "SOD1 Positive / No treatment",
      sod1_status == "Positive" &  treated_tofersen ~ "SOD1 Positive / Tofersen"
    ), levels = c("Non-SOD1", "SOD1 Positive / No treatment", "SOD1 Positive / Tofersen"))
  ) |>
  with(survfit2(Surv(age-age_onset, status == "Deceased") ~ group)) |>
  ggsurvfit() +
  add_confidence_interval() +
  add_legend_title("SOD1") +
  scale_ggsurvfit() +
  scale_x_continuous(limits = c(0, 10)) +
  add_pvalue("annotation") +
  custom_scale_color_brewer() +
  custom_scale_fill_brewer() +
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave(
    file.path(output_dir, "survival-by-sod1_status-tofersen.png"),
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
  drop_na(progression_category) |>
  mutate(across(progression_category, ~factor(.x, levels = c("SP", "FP", "NP")))) |>
  ggplot(aes(
      years_since_onset, alsfrs_total, group = record_id,
      color = progression_category, fill = progression_category
  )) +
  geom_line(linewidth = 0.1) + geom_jitter(size = 0.1, alpha = 0.5) +
  geom_smooth(aes(group = NULL, fill = progression_category)) +
  scale_x_continuous(breaks = seq(0, 10, 2), minor_breaks = seq(0, 10, 1)) +
  scale_y_continuous(breaks = seq(0, 48, 6), minor_breaks = seq(0, 48, 3)) +
  custom_scale_color_brewer() +
  custom_scale_fill_brewer() +
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
