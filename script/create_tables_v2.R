# --- Packages ---
library(dplyr)
library(stringr)
library(knitr)
library(kableExtra)
library(glue)

# ============= SETTINGS =============
rotate_pages <- FALSE   # TRUE wraps each table in \begin{landscape}...\end{landscape}
options(knitr.kable.NA = "")

# ============= 1) CLEAN NUMERIC COLUMNS =============
num_cols <- c("True value","Estimate (unweighted)","SD (unweighted)",
              "Estimate (weighted)","SD (weighted)","rho")

results_clean <- results %>%
  mutate(
    across(all_of(num_cols), \(x) {
      if (is.numeric(x)) return(x)
      v <- as.character(x)
      v <- str_replace_all(v, fixed("\u2212"), "-")  # Unicode minus → ASCII
      v <- str_replace_all(v, "[,\\s]", "")         # remove commas/spaces
      suppressWarnings(as.numeric(v))
    })
  )

# ============= 2) SUMMARISE ACROSS REPLICATES =============
summ <- results_clean %>%
  mutate(rho_num = as.numeric(rho)) %>%
  group_by(outcome, parameter, rho_num) %>%
  summarise(
    N         = dplyr::n(),
    True      = first(`True value`),
    
    # Unweighted
    unw_mean    = mean(`Estimate (unweighted)`, na.rm = TRUE),
    unw_mc_sd   = sd(  `Estimate (unweighted)`, na.rm = TRUE),
    unw_mean_se = mean(`SD (unweighted)`,       na.rm = TRUE),
    unw_bias    = mean(`Estimate (unweighted)` - True, na.rm = TRUE),
    unw_rmse    = sqrt(mean((`Estimate (unweighted)` - True)^2, na.rm = TRUE)),
    unw_cov95   = mean(
      True >= `Estimate (unweighted)` - 1.96 * `SD (unweighted)` &
        True <= `Estimate (unweighted)` + 1.96 * `SD (unweighted)`,
      na.rm = TRUE
    ),
    
    # Weighted
    w_mean      = mean(`Estimate (weighted)`, na.rm = TRUE),
    w_mc_sd     = sd(  `Estimate (weighted)`, na.rm = TRUE),
    w_mean_se   = mean(`SD (weighted)`,       na.rm = TRUE),
    w_bias      = mean(`Estimate (weighted)` - True, na.rm = TRUE),
    w_rmse      = sqrt(mean((`Estimate (weighted)` - True)^2, na.rm = TRUE)),
    w_cov95     = mean(
      True >= `Estimate (weighted)` - 1.96 * `SD (weighted)` &
        True <= `Estimate (weighted)` + 1.96 * `SD (weighted)`,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  arrange(outcome, parameter, rho_num)

fmt4 <- function(x) ifelse(is.na(x), "", formatC(x, digits = 4, format = "f"))

stack_mean_sd <- function(mean, sd) {
  sprintf("\\makecell[r]{%s\\\\\\scriptsize(%s)}", fmt4(mean), fmt4(sd))
}

tab_fmt <- summ %>%
  mutate(
    parameter = gsub("a_(\\d+)", "$a_{\\1}$", parameter),
    rho_disp  = paste0("$", rho_num, "$"),
    
    unw_mean_sd = stack_mean_sd(unw_mean, unw_mc_sd),
    w_mean_sd   = stack_mean_sd(w_mean,  w_mc_sd),
    
    unw_cov95 = sprintf("%.1f\\%%", 100 * unw_cov95),
    w_cov95   = sprintf("%.1f\\%%", 100 * w_cov95)
  ) %>%
  mutate(across(
    c(True, unw_mean_se, unw_bias, unw_rmse,
      w_mean_se, w_bias, w_rmse),
    ~ round(.x, 4)
  )) %>%
  arrange(outcome, parameter, rho_num)

out_dir <- "tables"
dir.create(out_dir, showWarnings = FALSE)

outcomes <- unique(tab_fmt$outcome)

for (oc in outcomes) {
  
  base <- str_to_lower(str_replace_all(oc, "[^A-Za-z0-9]+", "_"))
  
  df <- tab_fmt %>%
    filter(outcome == oc) %>%
    select(
      parameter, rho_disp, True, N,
      
      unw_mean_sd, unw_mean_se, unw_bias, unw_rmse, unw_cov95,
      w_mean_sd,   w_mean_se,   w_bias,   w_rmse, w_cov95
      
    )
  
  col_names <- c(
    "Parameter", "$\\rho$", "True", "$N$",
    "Mean (MC SD)", "Mean SE", "Bias", "RMSE", "Cov. 95\\%",
    "Mean (MC SD)", "Mean SE", "Bias", "RMSE", "Cov. 95\\%"
  )
  
  align <- c("c","c","r","r", rep("r", 8))
  
  kbl <- kable(
    df,
    format = "latex",
    booktabs = TRUE,
    longtable = TRUE,
    escape = FALSE,
    align = align,
    col.names = col_names,
    caption = glue("{oc}: Unweighted vs Weighted simulation results"),
    label = paste0("tab:", base, "-combined")
  ) |>
    add_header_above(c(
      " " = 4,
      "Unweighted" = 5,
      "Weighted"   = 5
    )) |>
    kable_styling(
      latex_options = c("striped"),
      full_width = FALSE,
      position = "center"
    ) |>
    collapse_rows(columns = 1, latex_hline = "major")
  
  save_kable(
    kbl,
    file = file.path(out_dir, paste0("results_", base, "_combined.tex"))
  )
}

