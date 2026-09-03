suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(scales)
  library(sf)
})

input_file <- "outputs/annualized_rates/stgpr2023_only/all_bugdrug_annualized_rates_stgpr2023_only.csv"
output_dir <- "outputs/annualized_rates/priority_scoring"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
excluded_bug_drugs <- c(
  "efaecalis_glycopeptide_observed_concat_stgpr",
  "efaecium_glycopeptide_observed_concat_stgpr"
)

# Equal-sized decile bins (1-10) across non-missing rows within each bug-drug group.
assign_decile_points <- function(x) {
  if (all(is.na(x))) return(rep(NA_integer_, length(x)))
  dplyr::ntile(x, 10)
}

get_province_lookup <- function() {
  gpkg_path <- "data/spatial/gadm41_ZAF_2_sf.gpkg"
  rds_path <- "data/spatial/gadm41_ZAF_2_sf.rds"
  if (file.exists(gpkg_path)) {
    shp <- sf::st_read(gpkg_path, quiet = TRUE)
  } else if (file.exists(rds_path)) {
    shp <- readRDS(rds_path)
  } else {
    return(tibble(district = character(), province = character()))
  }
  shp %>%
    sf::st_drop_geometry() %>%
    distinct(CC_2, NAME_1) %>%
    transmute(
      district = as.character(CC_2),
      province = as.character(NAME_1)
    )
}

get_sa_district_sf <- function() {
  gpkg_path <- "data/spatial/gadm41_ZAF_2_sf.gpkg"
  rds_path <- "data/spatial/gadm41_ZAF_2_sf.rds"
  if (file.exists(gpkg_path)) {
    return(sf::st_read(gpkg_path, quiet = TRUE))
  }
  if (file.exists(rds_path)) {
    return(readRDS(rds_path))
  }
  NULL
}

get_lesotho_sf <- function() {
  lesotho_sf <- NULL
  try({
    if (requireNamespace("geodata", quietly = TRUE)) {
      lso_v <- geodata::gadm("LSO", level = 0, path = "data/spatial")
      lesotho_sf <- sf::st_as_sf(lso_v)
    }
  }, silent = TRUE)
  if (is.null(lesotho_sf)) {
    try({
      if (requireNamespace("rnaturalearth", quietly = TRUE)) {
        lesotho_sf <- rnaturalearth::ne_countries(
          country = "Lesotho",
          scale = "medium",
          returnclass = "sf"
        )
      }
    }, silent = TRUE)
  }
  if (is.null(lesotho_sf)) {
    # Fallback polygon roughly around Lesotho if external boundary fetch fails.
    coords <- matrix(
      c(
        27.0, -30.7,
        29.5, -30.7,
        29.5, -28.6,
        27.0, -28.6,
        27.0, -30.7
      ),
      ncol = 2,
      byrow = TRUE
    )
    lesotho_sf <- sf::st_sf(
      geometry = sf::st_sfc(sf::st_polygon(list(coords)), crs = 4326)
    )
  }
  lesotho_sf
}

raw_df <- read_csv(input_file, show_col_types = FALSE)

# Map current schema to requested variable names.
mapped_df <- raw_df %>%
  transmute(
    district = district_code,
    bug_drug = bugdrug,
    modelled_resistance_2023 = end_value,
    annual_absolute_change_2013_2023 = annual_absolute_change
  ) %>%
  filter(!bug_drug %in% excluded_bug_drugs)

# Exclude rows with NA district/values per bug-drug before rank/decile assignment.
scored_df <- mapped_df %>%
  filter(
    !is.na(district),
    !is.na(bug_drug),
    !is.na(modelled_resistance_2023),
    !is.na(annual_absolute_change_2013_2023)
  ) %>%
  group_by(bug_drug) %>%
  mutate(
    trend_nonnegative = pmax(annual_absolute_change_2013_2023, 0),
    prevalence_decile_points = assign_decile_points(modelled_resistance_2023),
    trend_decile_points = assign_decile_points(trend_nonnegative),
    priority_score = prevalence_decile_points + trend_decile_points,
    priority_score_scaled = (priority_score - 2) / 18 * 100
  ) %>%
  ungroup() %>%
  select(
    district,
    bug_drug,
    modelled_resistance_2023,
    annual_absolute_change_2013_2023,
    prevalence_decile_points,
    trend_decile_points,
    priority_score,
    priority_score_scaled
  )

district_summary <- scored_df %>%
  group_by(district) %>%
  summarise(
    district_priority_max = max(priority_score, na.rm = TRUE),
    district_priority_sum = sum(priority_score, na.rm = TRUE),
    critical_bugdrug_count = sum(priority_score >= 18, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(district_priority_sum), desc(district_priority_max), district)

# Write outputs
full_out <- file.path(output_dir, "district_bugdrug_priority_scores.csv")
summary_out <- file.path(output_dir, "district_priority_summary.csv")
write_csv(scored_df, full_out)
write_csv(district_summary, summary_out)

# Visual 1: Top districts by priority sum
province_lookup <- get_province_lookup()
province_levels <- c(
  "Western Cape",
  "Northern Cape",
  "Eastern Cape",
  "Free State",
  "KwaZulu-Natal",
  "North West",
  "Gauteng",
  "Mpumalanga",
  "Limpopo"
)

p1_df <- district_summary %>%
  slice_head(n = 25) %>%
  left_join(province_lookup, by = "district") %>%
  mutate(
    province = ifelse(is.na(province), "Unknown", province),
    province = factor(province, levels = c(province_levels, "Unknown"))
  ) %>%
  group_by(province) %>%
  arrange(desc(district_priority_sum), district, .by_group = TRUE) %>%
  ungroup() %>%
  arrange(province, desc(district_priority_sum), district) %>%
  mutate(district = factor(district, levels = rev(unique(district))))

tier_breaks <- as.numeric(
  quantile(
    district_summary$district_priority_sum,
    probs = c(0.25, 0.50, 0.75),
    na.rm = TRUE
  )
)

p1 <- p1_df %>%
  ggplot(aes(x = district, y = district_priority_sum, fill = district_priority_max)) +
  geom_col(width = 0.8) +
  geom_hline(
    yintercept = tier_breaks,
    linetype = "dashed",
    linewidth = 0.35,
    color = "grey35"
  ) +
  coord_flip() +
  scale_fill_gradient(
    low = "#FEE5D9",
    high = "#A50F15"
  ) +
  facet_grid(province ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(
    title = "Top 25 Districts by Total AMR Priority Score",
    x = "District",
    y = "District Priority Sum",
    fill = "Highest single-\ncombination score",
    caption = paste0(
      "Each bar represents the sum of decile-based prioritisation scores across six pathogen-drug combinations ",
      "(range 12-120).\nColour intensity reflects the highest single-combination score."
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.y = element_text(size = 8.7),
    strip.text.y.left = element_text(size = 9, face = "bold", angle = 0, hjust = 0),
    strip.background = element_rect(fill = "grey96", color = "grey82", linewidth = 0.25),
    strip.placement = "outside",
    panel.spacing.y = grid::unit(0.22, "lines"),
    panel.border = element_rect(fill = NA, color = "grey88", linewidth = 0.2),
    plot.caption = element_text(size = 8.5, hjust = 0),
    legend.position = "right",
    plot.margin = margin(10, 16, 10, 30)
  )

ggsave(
  filename = file.path(output_dir, "district_priority_sum_top25.png"),
  plot = p1,
  width = 12.2,
  height = 8.6,
  dpi = 300
)

# Visual 1c: Systemic vs acute risk scatter
p1_scatter_df <- district_summary %>%
  left_join(province_lookup, by = "district") %>%
  mutate(
    province = ifelse(is.na(province), "Unknown", province),
    province = factor(province, levels = c(province_levels, "Unknown"))
  )

p1_scatter <- ggplot(
  p1_scatter_df,
  aes(
    x = district_priority_sum,
    y = district_priority_max,
    size = critical_bugdrug_count,
    color = province
  )
) +
  geom_point(alpha = 0.9) +
  scale_size_continuous(range = c(2.5, 10)) +
  scale_color_viridis_d(option = "plasma", end = 0.92) +
  labs(
    title = "District AMR Priority: Systemic vs Acute Risk",
    subtitle = "X = total district priority sum, Y = highest single-combination score, size = count of critical combinations",
    x = "District priority sum",
    y = "District priority max",
    size = "Critical\ncombinations",
    color = "Province"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = file.path(output_dir, "district_priority_sum_vs_max_scatter.png"),
  plot = p1_scatter,
  width = 11,
  height = 7.5,
  dpi = 300
)

# Visual 1b: All districts by priority sum, grouped by province
p1_all_df <- district_summary %>%
  left_join(province_lookup, by = "district") %>%
  mutate(
    province = ifelse(is.na(province), "Unknown", province),
    province = factor(province, levels = c(province_levels, "Unknown"))
  ) %>%
  arrange(province, desc(district_priority_sum), district) %>%
  mutate(district = factor(district, levels = rev(unique(district))))

p1_all <- p1_all_df %>%
  ggplot(aes(x = district, y = district_priority_sum, fill = district_priority_max)) +
  geom_col(width = 0.8) +
  coord_flip() +
  scale_fill_gradient(low = "#8ecae6", high = "#c1121f") +
  facet_grid(province ~ ., scales = "free_y", space = "free_y") +
  labs(
    title = "Modelled 2023 Districts by Total AMR Priority Score (Grouped by Province)",
    subtitle = paste0("N = ", nrow(p1_all_df), " districts with STGPR predictions in 2023"),
    x = "District",
    y = "District Priority Sum",
    fill = "District\nPriority Max"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  filename = file.path(output_dir, "district_priority_sum_all_districts.png"),
  plot = p1_all,
  width = 10,
  height = 14,
  dpi = 300
)

# Visual 2: Heatmap (top 30 districts by sum x all bug-drugs)
top_districts <- district_summary %>%
  slice_head(n = 30) %>%
  pull(district)

heat_df <- scored_df %>%
  filter(district %in% top_districts) %>%
  mutate(district = factor(district, levels = rev(top_districts)))

p2 <- ggplot(heat_df, aes(x = bug_drug, y = district, fill = priority_score)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_gradientn(colors = c("#e9f5db", "#90be6d", "#f8961e", "#d00000")) +
  labs(
    title = "District-BugDrug Priority Score Heatmap (Top 30 Districts)",
    x = "Bug-Drug",
    y = "District",
    fill = "Priority\nScore"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1)
  )

ggsave(
  filename = file.path(output_dir, "district_bugdrug_priority_heatmap_top30.png"),
  plot = p2,
  width = 14,
  height = 9,
  dpi = 300
)

# Visual 3: Priority score distribution
p3 <- ggplot(scored_df, aes(x = priority_score)) +
  geom_histogram(binwidth = 1, fill = "#457b9d", color = "white", boundary = 1.5) +
  scale_x_continuous(breaks = 2:20) +
  labs(
    title = "Distribution of District x BugDrug Priority Scores",
    x = "Priority Score (2-20)",
    y = "Count"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  filename = file.path(output_dir, "priority_score_distribution.png"),
  plot = p3,
  width = 11,
  height = 6,
  dpi = 300
)

# Visual 4: Map of districts included in current priority scoring set
sa_gadm2 <- get_sa_district_sf()
if (!is.null(sa_gadm2) && "CC_2" %in% names(sa_gadm2)) {
  map_df <- sa_gadm2 %>%
    left_join(district_summary, by = c("CC_2" = "district")) %>%
    mutate(
      district = as.character(CC_2)
    )
  lesotho_sf <- get_lesotho_sf()

  p4 <- ggplot() +
    geom_sf(data = map_df, aes(fill = district_priority_sum), color = "grey55", linewidth = 0.15) +
    geom_sf(data = lesotho_sf, fill = "grey75", color = "grey55", linewidth = 0.2) +
    scale_fill_gradientn(
      colors = c("#e9f5db", "#90be6d", "#f8961e", "#d00000"),
      na.value = "#f2f2f2",
      name = "Priority Sum"
    ) +
    labs(
      title = "District Priority Score Map",
      subtitle = "Color coding by district priority sum (STGPR 2023 districts scored)",
      x = NULL,
      y = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3),
      legend.position = "bottom"
    )

  ggsave(
    filename = file.path(output_dir, "district_priority_sum_map.png"),
    plot = p4,
    width = 9,
    height = 8,
    dpi = 300
  )
}

cat("Saved full scored data:", full_out, "\n")
cat("Saved district summary:", summary_out, "\n")
cat("Rows (scored_df):", nrow(scored_df), "\n")
cat("Rows (district_summary):", nrow(district_summary), "\n")
cat("Distinct bug_drug:", dplyr::n_distinct(scored_df$bug_drug), "\n")
cat("Distinct districts:", dplyr::n_distinct(scored_df$district), "\n")
