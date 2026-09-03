# Step 4: Merge AMR and Socioeconomic Data (script version)

suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
  library(writexl)
  library(ggplot2)
  library(cowplot)
  library(sf)
  library(geodata)
})

# ------------------------------------------------------------------
# Load inputs
# ------------------------------------------------------------------
step2_path <- file.path("data", "step2_imputation_gam_output.xlsx")
step3_path <- file.path("data", "step3_cleanedAMRdata_output.xlsx")

if (!file.exists(step2_path)) {
  stop("Expected file not found: ", step2_path, ". Run step2_imputation_gam first.")
}

if (!file.exists(step3_path)) {
  stop("Expected file not found: ", step3_path, ". Run step3_cleanAMRdata first.")
}

socio_data <- read_xlsx(step2_path)
amr_data <- read_xlsx(step3_path)

# ------------------------------------------------------------------
# Prepare AMR summary by district and year
# ------------------------------------------------------------------
vars <- c(
  "abaumanii_carbapenem",
  "ecoli_3gc",
  "saureus_penicillinase",
  "kpneumoniae_carbapenem",
  "paeruginosa_aminoglycoside",
  "enterobacter_3gc",
  "efaecium_glycopeptide",
  "efaecalis_glycopeptide"
)

cleaned_amr_data <- amr_data %>%
  mutate(
    district_code = as.character(district_code),
    year_taken = as.integer(year_taken)
  ) %>%
  filter(!is.na(district_code), trimws(district_code) != "", !is.na(year_taken))

district_year_summary <- cleaned_amr_data %>%
  group_by(district_code, year_taken) %>%
  summarise(
    across(
      all_of(vars),
      list(
        resistant = ~ sum(. == "Resistant", na.rm = TRUE),
        sensitive = ~ sum(. == "Sensitive", na.rm = TRUE),
        tested = ~ sum(!is.na(.) & . != "", na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )

district_year_summary_L1 <- cleaned_amr_data %>%
  filter(hospital_tier == "Level1") %>%
  group_by(district_code, year_taken) %>%
  summarise(
    across(
      all_of(vars),
      list(
        resistant = ~ sum(. == "Resistant", na.rm = TRUE),
        sensitive = ~ sum(. == "Sensitive", na.rm = TRUE),
        tested = ~ sum(!is.na(.) & . != "", na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )

district_year_summary_without_L1 <- cleaned_amr_data %>%
  filter(is.na(hospital_tier) | hospital_tier != "Level1") %>%
  group_by(district_code, year_taken) %>%
  summarise(
    across(
      all_of(vars),
      list(
        resistant = ~ sum(. == "Resistant", na.rm = TRUE),
        sensitive = ~ sum(. == "Sensitive", na.rm = TRUE),
        tested = ~ sum(!is.na(.) & . != "", na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )

write.csv(district_year_summary, file.path("data", "step4_amr_district_year_summary_alllevels.csv"), row.names = FALSE)
write.csv(district_year_summary_L1, file.path("data", "step4_amr_district_year_summary_level1.csv"), row.names = FALSE)
write.csv(
  district_year_summary_without_L1,
  file.path("data", "step4_amr_district_year_summary_without_level1.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# Split summaries by bug–drug and save per file
# ------------------------------------------------------------------
output_dir <- file.path("data", "amr_data_bug_drug")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

make_subset <- function(df, bug_drug) {
  tested_col <- paste0(bug_drug, "_tested")
  resistant_col <- paste0(bug_drug, "_resistant")
  resistance_per_col <- paste0(bug_drug, "_resistance_per")
  sub <- dplyr::select(
    df,
    district_code,
    year_taken,
    all_of(resistant_col),
    all_of(paste0(bug_drug, "_sensitive")),
    all_of(tested_col)
  )
  sub <- dplyr::filter(
    sub,
    !is.na(.data[[tested_col]]),
    .data[[tested_col]] != 0,
    .data[[tested_col]] != ""
  )
  dplyr::mutate(
    sub,
    `above 30` = dplyr::if_else(.data[[tested_col]] >= 30, "yes", "no"),
    !!resistance_per_col := .data[[resistant_col]] / .data[[tested_col]]
  )
}

bug_drugs <- c(
  "abaumanii_carbapenem",
  "ecoli_3gc",
  "saureus_penicillinase",
  "kpneumoniae_carbapenem",
  "paeruginosa_aminoglycoside",
  "enterobacter_3gc",
  "efaecium_glycopeptide",
  "efaecalis_glycopeptide"
)

bug_drug_files_all <- lapply(bug_drugs, function(bd) {
  sub_df <- make_subset(district_year_summary, bd)
  write.csv(sub_df, file.path(output_dir, paste0(bd, "_alllevels.csv")), row.names = FALSE)
  sub_df
})

bug_drug_files_L1 <- lapply(bug_drugs, function(bd) {
  sub_df <- make_subset(district_year_summary_L1, bd)
  write.csv(sub_df, file.path(output_dir, paste0(bd, "_level1.csv")), row.names = FALSE)
  sub_df
})

bug_drug_files_without_L1 <- lapply(bug_drugs, function(bd) {
  sub_df <- make_subset(district_year_summary_without_L1, bd)
  write.csv(sub_df, file.path(output_dir, paste0(bd, "_without_level1.csv")), row.names = FALSE)
  sub_df
})

# ------------------------------------------------------------------
# Maps: resistance by district/year for all bug–drug combos (all levels and Level1)
# ------------------------------------------------------------------
zaf_gadm2 <- geodata::gadm("ZAF", level = 2, path = "data") |>
  sf::st_as_sf()

map_dir_base <- file.path("maps", "available data maps")
if (!dir.exists(map_dir_base)) dir.create(map_dir_base, recursive = TRUE)

map_inputs <- list(
  alllevels = district_year_summary,
  level1 = district_year_summary_L1,
  without_level1 = district_year_summary_without_L1
)

run_panel_only <- identical(Sys.getenv("RUN_PANEL_ONLY", "0"), "1")

if (!run_panel_only) {
invisible(lapply(names(map_inputs), function(level_name) {
  level_df <- map_inputs[[level_name]]

  lapply(bug_drugs, function(bd) {
    map_dir <- file.path(map_dir_base, level_name, bd, "maps")
    if (!dir.exists(map_dir)) dir.create(map_dir, recursive = TRUE)

    resistance_per_col <- paste0(bd, "_resistance_per")
    tested_col <- paste0(bd, "_tested")

    bd_df <- make_subset(level_df, bd) |>
      dplyr::mutate(district_code = as.character(district_code))

    years <- sort(unique(bd_df$year_taken))
    plot_list <- list()

    for (yr in years) {
      df_year <- dplyr::filter(bd_df, year_taken == yr)

      map_data <- zaf_gadm2 |>
        dplyr::left_join(df_year, by = c("CC_2" = "district_code"))

      message(sprintf("Creating map: %s | %s | %s | year %d", level_name, bd, "all data", yr))
      map_pts <- map_data %>%
        dplyr::filter(!is.na(.data[[tested_col]])) %>%
        sf::st_point_on_surface()

      p <- ggplot(map_data) +
        geom_sf(aes(fill = .data[[resistance_per_col]]), color = "grey50", size = 0.1) +
        geom_sf(
          data = map_pts,
          aes(size = .data[[tested_col]]),
          shape = 21,
          fill = "black",
          color = "white",
          alpha = 0.5,
          stroke = 0.2
        ) +
        scale_fill_gradientn(
          colors = c("#FFE52A", "#F79A19", "#CF0F0F"),
          limits = c(0, 1),
          na.value = "white",
          name = "Resistance proportion"
        ) +
        scale_size_continuous(
          name = "Tested isolates",
          trans = "sqrt",
          range = c(0.5, 4)
        ) +
        labs(
          title = paste(bd, "Resistance —", yr, sprintf("(%s)", level_name)),
          subtitle = "District-level South Africa",
          caption = "Fill = resistance proportion; point size = tested isolates. Source: district_year_summary + GADM (geodata)"
        ) +
        theme_minimal(base_size = 14) +
        theme(
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.box = "horizontal",
          panel.border = element_blank()
        )

      out_file <- file.path(map_dir, sprintf("%s_resistance_%s_%d.png", bd, level_name, yr))
      if (!file.exists(out_file)) {
        ggsave(
          filename = out_file,
          plot = p,
          width = 10, height = 8, dpi = 300
        )
      } else {
        message(sprintf("Skipping existing map: %s", out_file))
      }
      plot_list[[length(plot_list) + 1]] <- p
    }

    if (length(plot_list) > 0) {
      message(sprintf("Creating panel map: %s | %s | %s", level_name, bd, "all data"))
      combined <- cowplot::plot_grid(plotlist = plot_list, ncol = 1)
      per_plot_height <- 6
      panel_height <- min(per_plot_height * length(plot_list), 45)  # stay well under ggsave 50\" limit
      out_panel <- file.path(map_dir, sprintf("%s_resistance_%s_panel.png", bd, level_name))
      if (!file.exists(out_panel)) {
        ggsave(
          filename = out_panel,
          plot = combined,
          width = 10, height = panel_height, dpi = 300
        )
      } else {
        message(sprintf("Skipping existing panel: %s", out_panel))
      }
    }
  })
}))
}

# ------------------------------------------------------------------
# Panel: all levels above threshold, year 2023 only (6 bug–drug combos)
# ------------------------------------------------------------------
panel_year <- 2023
panel_level <- "alllevels"
panel_bug_drugs <- setdiff(
  bug_drugs,
  c("efaecium_glycopeptide", "efaecalis_glycopeptide")
)
panel_level_df <- map_inputs[[panel_level]]

panel_bd_list <- lapply(panel_bug_drugs, function(bd) {
  make_subset(panel_level_df, bd) %>%
    dplyr::filter(`above 30` == "yes", year_taken == panel_year) %>%
    dplyr::mutate(district_code = as.character(district_code))
})
names(panel_bd_list) <- panel_bug_drugs

panel_map_dir <- file.path(map_dir_base, paste0(panel_level, "_abovethreshold"), "panel")
if (!dir.exists(panel_map_dir)) dir.create(panel_map_dir, recursive = TRUE)

panel_out_file <- file.path(panel_map_dir, sprintf("alllevels_abovethreshold_panel_%d.png", panel_year))
if (run_panel_only || !file.exists(panel_out_file)) {
  panel_max_tested <- max(
    unlist(lapply(panel_bug_drugs, function(bd) panel_bd_list[[bd]][[paste0(bd, "_tested")]])),
    na.rm = TRUE
  )

  make_panel_map <- function(bd, with_legend = FALSE, panel_label = NULL) {
    resistance_per_col <- paste0(bd, "_resistance_per")
    tested_col <- paste0(bd, "_tested")

    bd_df <- panel_bd_list[[bd]]

    map_data <- zaf_gadm2 |>
      dplyr::left_join(bd_df, by = c("CC_2" = "district_code"))

    map_pts <- map_data %>%
      dplyr::filter(!is.na(.data[[tested_col]])) %>%
      sf::st_point_on_surface()

    p <- ggplot(map_data) +
      geom_sf(aes(fill = .data[[resistance_per_col]]), color = "grey50", size = 0.1) +
      geom_sf(
        data = map_pts,
        aes(size = .data[[tested_col]]),
        shape = 21,
        fill = "black",
        color = "white",
        alpha = 0.5,
        stroke = 0.2
      ) +
      scale_fill_gradientn(
        colors = c("#FFE52A", "#F79A19", "#CF0F0F"),
        limits = c(0, 1),
        na.value = "white",
        name = "Resistance proportion"
      ) +
      scale_size_continuous(
        name = "Tested isolates",
        trans = "sqrt",
        limits = c(0, panel_max_tested),
        range = c(0.5, 4),
        breaks = seq(250, panel_max_tested, by = 250)
      ) +
      guides(
        fill = guide_colorbar(title.position = "top", barwidth = unit(6, "cm")),
        size = guide_legend(title.position = "top", override.aes = list(alpha = 0.6))
      ) +
      labs(title = NULL, subtitle = NULL, caption = NULL) +
      theme_minimal(base_size = 12) +
      theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        plot.margin = margin(0, 0, 0, 0),
        legend.position = if (with_legend) "bottom" else "none",
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.title.align = 0.5
      )

    if (!is.null(panel_label)) {
      bbox <- sf::st_bbox(map_data)
      label_x <- bbox["xmin"] + 0.02 * (bbox["xmax"] - bbox["xmin"])
      label_y <- bbox["ymax"] - 0.02 * (bbox["ymax"] - bbox["ymin"])
      p <- p + annotate(
        "text",
        x = label_x,
        y = label_y,
        label = panel_label,
        hjust = 0,
        vjust = 1,
        size = 5,
        fontface = "bold"
      )
    }

    p
  }

  panel_plots <- lapply(seq_along(panel_bug_drugs), function(i) {
    make_panel_map(panel_bug_drugs[[i]], panel_label = LETTERS[i])
  })

  # Build a dedicated legend to avoid missing guides when data are sparse
  legend_df <- data.frame(
    lon = c(16, 18),
    lat = c(-33, -31),
    resistance = c(0.1, 0.9),
    tested = c(10, panel_max_tested)
  )
  legend_plot <- ggplot(legend_df, aes(x = lon, y = lat)) +
    geom_point(aes(size = tested, fill = resistance), shape = 21, color = "white", alpha = 0.6, stroke = 0.2) +
    scale_fill_gradientn(
      colors = c("#FFE52A", "#F79A19", "#CF0F0F"),
      limits = c(0, 1),
      name = "Resistance proportion"
    ) +
    scale_size_continuous(
      name = "Tested isolates",
      trans = "sqrt",
      limits = c(0, panel_max_tested),
      range = c(0.5, 4),
      breaks = seq(250, panel_max_tested, by = 250)
    ) +
    guides(
      fill = guide_colorbar(title.position = "top", barwidth = unit(6, "cm")),
      size = guide_legend(
        title.position = "top",
        override.aes = list(alpha = 0.8, shape = 21, fill = "grey50", color = "grey50"),
        nrow = 1,
        label.position = "bottom"
      )
    ) +
    theme_void(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "vertical",
      legend.title.align = 0.5,
      legend.margin = margin(t = 6, r = 6, b = 10, l = 6),
      legend.box.margin = margin(t = 4, r = 4, b = 8, l = 4),
      legend.text = element_text(color = "black"),
      legend.title = element_text(color = "black"),
      legend.key.height = unit(0.8, "cm"),
      legend.key.width = unit(0.8, "cm"),
      plot.background = element_rect(fill = "white", color = NA),
      legend.background = element_rect(fill = "white", color = NA),
      legend.key = element_rect(fill = "white", color = NA)
    )

  panel_legend <- legend_plot

  # Use patchwork for tight grid, then cowplot to place legend on the right
  panel_grid <- patchwork::wrap_plots(panel_plots, ncol = 3) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(theme = theme(plot.margin = margin(0, 0, 0, 0)))

  combined_panel <- cowplot::plot_grid(
    panel_grid,
    panel_legend,
    ncol = 2,
    rel_widths = c(1, 0.30),
    align = "h",
    axis = "tb"
  )

  ggsave(
    filename = panel_out_file,
    plot = combined_panel,
    width = 14,
    height = 9,
    dpi = 300
  )
} else {
  message(sprintf("Skipping existing panel: %s", panel_out_file))
}
if (run_panel_only) {
  message("Panel-only run complete.")
  quit(save = "no")
}

# Maps: only rows with tested >= 30 (above threshold)
if (!run_panel_only) {
invisible(lapply(names(map_inputs), function(level_name) {
  level_df <- map_inputs[[level_name]]

  lapply(bug_drugs, function(bd) {
    map_dir <- file.path(map_dir_base, paste0(level_name, "_abovethreshold"), bd, "maps")
    csv_dir <- file.path("data", "amr_data_bug_drug_abovethreshold", level_name, bd)
    if (!dir.exists(map_dir)) dir.create(map_dir, recursive = TRUE)
    if (!dir.exists(csv_dir)) dir.create(csv_dir, recursive = TRUE)

    resistance_per_col <- paste0(bd, "_resistance_per")
    tested_col <- paste0(bd, "_tested")

    bd_df <- make_subset(level_df, bd) |>
      dplyr::mutate(district_code = as.character(district_code)) |>
      dplyr::filter(`above 30` == "yes")

    write.csv(
      bd_df,
      file.path(csv_dir, sprintf("%s_%s_abovethreshold.csv", bd, level_name)),
      row.names = FALSE
    )

    years <- sort(unique(bd_df$year_taken))
    plot_list <- list()

    for (yr in years) {
      df_year <- dplyr::filter(bd_df, year_taken == yr)

      map_data <- zaf_gadm2 |>
        dplyr::left_join(df_year, by = c("CC_2" = "district_code"))

      message(sprintf("Creating map: %s | %s | %s | year %d", level_name, bd, "above threshold", yr))
      map_pts <- map_data %>%
        dplyr::filter(!is.na(.data[[tested_col]])) %>%
        sf::st_point_on_surface()

      p <- ggplot(map_data) +
        geom_sf(aes(fill = .data[[resistance_per_col]]), color = "grey50", size = 0.1) +
        geom_sf(
          data = map_pts,
          aes(size = .data[[tested_col]]),
          shape = 21,
          fill = "black",
          color = "white",
          alpha = 0.5,
          stroke = 0.2
        ) +
        scale_fill_gradientn(
          colors = c("#FFE52A", "#F79A19", "#CF0F0F"),
          limits = c(0, 1),
          na.value = "white",
          name = "Resistance proportion"
        ) +
        scale_size_continuous(
          name = "Tested isolates",
          trans = "sqrt",
          range = c(0.5, 4)
        ) +
        labs(
          title = paste(bd, "Resistance —", yr, sprintf("(%s above threshold)", level_name)),
          subtitle = "District-level South Africa (tested >= 30)",
          caption = "Fill = resistance proportion; point size = tested isolates. Source: district_year_summary + GADM (geodata)"
        ) +
        theme_minimal(base_size = 14) +
        theme(
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.box = "horizontal",
          panel.border = element_blank()
        )

      out_file <- file.path(map_dir, sprintf("%s_resistance_%s_abovethreshold_%d.png", bd, level_name, yr))
      if (!file.exists(out_file)) {
        ggsave(
          filename = out_file,
          plot = p,
          width = 10, height = 8, dpi = 300
        )
      } else {
        message(sprintf("Skipping existing map: %s", out_file))
      }
      plot_list[[length(plot_list) + 1]] <- p
    }

    if (length(plot_list) > 0) {
      message(sprintf("Creating panel map: %s | %s | %s", level_name, bd, "above threshold"))
      combined <- cowplot::plot_grid(plotlist = plot_list, ncol = 1)
      per_plot_height <- 6
      panel_height <- min(per_plot_height * length(plot_list), 45)  # stay well under ggsave 50\" limit
      out_panel <- file.path(map_dir, sprintf("%s_resistance_%s_abovethreshold_panel.png", bd, level_name))
      if (!file.exists(out_panel)) {
        ggsave(
          filename = out_panel,
          plot = combined,
          width = 10, height = panel_height, dpi = 300
        )
      } else {
        message(sprintf("Skipping existing panel: %s", out_panel))
      }
    }
  })
}))
}
