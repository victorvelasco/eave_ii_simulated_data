library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

out_dir <- "figures"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

results <- readRDS("results.RDS")

# ── 1. Clean numeric columns ──────────────────────────────────────────────────
num_cols <- c("True value", "Estimate (unweighted)", "SD (unweighted)",
              "Estimate (weighted)", "SD (weighted)")
results <- results |>
  mutate(across(all_of(num_cols), \(x) {
    if (is.numeric(x)) return(x)
    v <- str_replace_all(as.character(x), "−", "-")
    suppressWarnings(as.numeric(v))
  }))

# ── 2. Summarise across simulations ──────────────────────────────────────────
summ <- results |>
  group_by(outcome, parameter, rho) |>
  summarise(
    true       = first(`True value`),
    unw_mean   = mean(`Estimate (unweighted)`, na.rm = TRUE),
    unw_mc_sd  = sd(`Estimate (unweighted)`,   na.rm = TRUE),
    w_mean     = mean(`Estimate (weighted)`,   na.rm = TRUE),
    w_mc_sd    = sd(`Estimate (weighted)`,     na.rm = TRUE),
    .groups = "drop"
  )

# ── 3. Pivot to long form (one row per method) ────────────────────────────────
long <- bind_rows(
  summ |> transmute(
    outcome, parameter, rho, true,
    method = "Unweighted",
    mean   = unw_mean,
    lo     = unw_mean - 1.96 * unw_mc_sd,
    hi     = unw_mean + 1.96 * unw_mc_sd
  ),
  summ |> transmute(
    outcome, parameter, rho, true,
    method = "Weighted (MSM)",
    mean   = w_mean,
    lo     = w_mean - 1.96 * w_mc_sd,
    hi     = w_mean + 1.96 * w_mc_sd
  )
) |>
  mutate(
    method    = factor(method, levels = c("Weighted (MSM)", "Unweighted")),
    rho_label = factor(
      paste0("ρ = ", rho),
      levels = paste0("ρ = ", sort(unique(rho)))
    ),
    # y position: group by rho, offset by method within each rho
    y_id = as.numeric(rho_label) * 3 + as.numeric(method)
  )

# ── 4. True-value lookup (one per outcome × parameter) ────────────────────────
true_vals <- long |>
  distinct(outcome, parameter, rho_label, true)

# ── 5. Nice parameter labels ──────────────────────────────────────────────────
param_labels <- c(
  a_0 = expression(alpha[0]~"(intercept)"),
  a_1 = expression(alpha[1]~"(log OR, AZ vs unvaccinated)")
)

# ── 6. Colour palette ─────────────────────────────────────────────────────────
method_colours <- c("Unweighted" = "#E69F00", "Weighted (MSM)" = "#0072B2")
method_shapes  <- c("Unweighted" = 16,        "Weighted (MSM)" = 17)

# ── 7. Build a single combined figure (outcome × parameter grid) ──────────────

# True value per (outcome, parameter) panel — constant across rho
true_df <- long |>
  distinct(outcome, parameter, true)

# y-axis tick positions and labels (same for every panel)
rho_tick_pos <- long |>
  group_by(rho_label) |>
  summarise(y_mid = mean(y_id), .groups = "drop")

p <- ggplot(long, aes(x = mean, y = y_id, colour = method, shape = method)) +
  geom_errorbar(
    aes(xmin = lo, xmax = hi),
    width = 0.35, linewidth = 0.6,
    orientation = "y"
  ) +
  geom_point(size = 2.5) +
  geom_vline(
    data     = true_df,
    aes(xintercept = true),
    linetype = "dashed", colour = "black", linewidth = 0.6,
    inherit.aes = FALSE
  ) +
  facet_grid(
    outcome ~ parameter,
    scales   = "free_x",
    labeller = labeller(
      parameter = as_labeller(c(
        a_0 = "alpha[0]~(intercept)",
        a_1 = "alpha[1]~(log~OR)"
      ), label_parsed)
    )
  ) +
  scale_y_continuous(
    breaks = rho_tick_pos$y_mid,
    labels = rho_tick_pos$rho_label,
    expand = expansion(add = 1)
  ) +
  scale_colour_manual(values = method_colours, name = "Method") +
  scale_shape_manual( values = method_shapes,  name = "Method") +
  labs(
    x       = "Estimate (mean ± 1.96 × Monte Carlo SD)",
    y       = NULL,
    caption = "Dashed line: true parameter value"
  ) +
  theme_bw(base_size = 22) +
  theme(
    strip.text         = element_text(size = 20),
    strip.text.y       = element_text(size = 19, angle = 0),
    legend.position    = "bottom",
    legend.text        = element_text(size = 20),
    legend.title       = element_text(size = 20),
    axis.text          = element_text(size = 18),
    axis.title.x       = element_text(size = 20),
    plot.caption       = element_text(size = 17),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank()
  )

fname <- file.path(out_dir, "forest_all.png")
ggsave(fname, plot = p, width = 14, height = 18, dpi = 300)
message("Saved: ", fname)
