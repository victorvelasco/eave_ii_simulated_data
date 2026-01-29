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

# ============= 3) FORMAT COMMON PIECES =============
fmt4 <- function(x) ifelse(is.na(x), "", formatC(x, format = "f", digits = 4))
h3  <- function(a,b,c="") {           # 3-line header helper
  if (c == "") sprintf("\\makecell[c]{%s\\\\%s}", a, b)
  else         sprintf("\\makecell[c]{%s\\\\%s\\\\%s}", a, b, c)
}

summ_fmt <- summ %>%
  mutate(
    parameter = gsub("a_(\\d+)", "$a_{\\1}$", parameter),  # a_0 -> $a_{0}$
    rho_disp  = paste0("$", rho_num, "$"),
    # combine mean and MC SD in one cell, MC SD on next line in parentheses
    unw_mean_sd = sprintf("\\makecell[r]{%s\\\\\\scriptsize(%s)}", fmt4(unw_mean), fmt4(unw_mc_sd)),
    w_mean_sd   = sprintf("\\makecell[r]{%s\\\\\\scriptsize(%s)}", fmt4(w_mean),   fmt4(w_mc_sd)),
    # percentages
    unw_cov95 = sprintf("%.1f\\%%", 100 * unw_cov95),
    w_cov95   = sprintf("%.1f\\%%", 100 * w_cov95)
  )

# Round remaining numeric display columns
round_cols <- c("True","unw_mean_se","unw_bias","unw_rmse","w_mean_se","w_bias","w_rmse")
summ_fmt[round_cols] <- lapply(summ_fmt[round_cols], function(x) round(x, 4))

# ============= 4) BUILD TWO TABLES PER OUTCOME (with 3-line headers) =============
out_dir <- "tables"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Column sets
cols_unw <- c("parameter","rho_disp","True","N",
              "unw_mean_sd","unw_mean_se","unw_bias","unw_rmse","unw_cov95")
cols_w   <- c("parameter","rho_disp","True","N",
              "w_mean_sd","w_mean_se","w_bias","w_rmse","w_cov95")

# 3-line headers (short + tall to save width)
col_names_unw <- c(
  h3("","Parameter",""),
  h3("","$\\rho$",""),
  h3("","True",""),
  h3("","$N$",""),
  h3("Unweighted","Mean","(MC SD)"),
  h3("Unweighted","Mean","SE"),
  h3("Unweighted","Bias",""),
  h3("Unweighted","RMSE",""),
  h3("Unweighted","Cov","95\\%")
)
col_names_w <- c(
  h3("","Parameter",""),
  h3("","$\\rho$",""),
  h3("","True",""),
  h3("","$N$",""),
  h3("Weighted","Mean","(MC SD)"),
  h3("Weighted","Mean","SE"),
  h3("Weighted","Bias",""),
  h3("Weighted","RMSE",""),
  h3("Weighted","Cov","95\\%")
)

# Alignment: center first two cols; right-align the rest
align_unw <- c("c","c", rep("r", length(cols_unw) - 2))
align_w   <- c("c","c", rep("r", length(cols_w)   - 2))

outcomes <- sort(unique(summ_fmt$outcome))
files_map <- list()

for (oc in outcomes) {
  base <- str_to_lower(str_replace_all(oc, "[^A-Za-z0-9]+", "_"))
  df_oc <- filter(summ_fmt, outcome == oc)
  
  # --- Unweighted table ---
  df_unw <- select(df_oc, all_of(cols_unw))
  cap_unw <- glue("{oc}: Unweighted summary across simulations (N = {{df_unw$N[1]}} per cell).")
  file_unw <- file.path(out_dir, paste0("results_", base, "_unweighted.tex"))
  
  tbl_unw <- kable(
    df_unw %>% arrange(rho_disp),
    format    = "latex",
    booktabs  = TRUE,
    caption   = cap_unw,
    label     = paste0("tab:", base, "-unweighted"),
    col.names = col_names_unw,
    align     = align_unw,
    escape    = FALSE,
    longtable = TRUE
  ) |>
    kable_styling(latex_options = c("striped"),
                  full_width = FALSE, position = "center") |>
    row_spec(0, extra_css = "white-space:normal;", font_size = 9) |>  # smaller header
    collapse_rows(columns = 1, latex_hline = "major")
  
  save_kable(tbl_unw, file = file_unw)
  
  # --- Weighted table ---
  df_w <- select(df_oc, all_of(cols_w))
  cap_w <- glue("{oc}: Weighted summary across simulations (N = {{df_w$N[1]}} per cell).")
  file_w <- file.path(out_dir, paste0("results_", base, "_weighted.tex"))
  
  tbl_w <- kable(
    df_w %>% arrange(rho_disp),
    format    = "latex",
    booktabs  = TRUE,
    caption   = cap_w,
    label     = paste0("tab:", base, "-weighted"),
    col.names = col_names_w,
    align     = align_w,
    escape    = FALSE,
    longtable = TRUE
  ) |>
    kable_styling(latex_options = c("striped"),
                  full_width = FALSE, position = "center") |>
    row_spec(0, extra_css = "white-space:normal;", font_size = 9) |>  # smaller header
    collapse_rows(columns = 1, latex_hline = "major")
  
  save_kable(tbl_w, file = file_w)
  
  files_map[[oc]] <- list(unweighted = basename(file_unw), weighted = basename(file_w))
}

# ============= 5) WRITE main.tex (includes both tables per outcome) =============
header <- c(
  "\\documentclass[12pt,a4paper]{article}",
  "\\usepackage{geometry}",
  "\\geometry{margin=1in}",
  "\\usepackage{amsmath}",
  "\\usepackage{booktabs}",
  "\\usepackage{array}",
  "\\usepackage{longtable}",
  "\\usepackage{multirow}",
  "\\usepackage{makecell}",     # <-- needed for 3-line headers & stacked cells
  "\\usepackage{threeparttable}",
  "\\usepackage{threeparttablex}",
  "\\usepackage[table]{xcolor}",
  if (rotate_pages) "\\usepackage{pdflscape}" else NULL,
  "\\usepackage{hyperref}",
  "",
  "\\title{Simulation Study Results: Unweighted vs Weighted}",
  "\\author{Your Name}",
  "\\date{\\today}",
  "",
  "\\begin{document}",
  "\\maketitle",
  "",
  "\\section{Results by Outcome}"
)

body <- character(0)
for (oc in outcomes) {
  base <- str_to_lower(str_replace_all(oc, "[^A-Za-z0-9]+", "_"))
  fn_unw <- files_map[[oc]]$unweighted
  fn_w   <- files_map[[oc]]$weighted
  
  body <- c(body,
            paste0("\\subsection{", oc, "}"),
            # Unweighted
            "\\subsubsection{Unweighted}",
            if (rotate_pages) "\\begin{landscape}" else NULL,
            paste0("\\input{tables/", fn_unw, "}"),
            if (rotate_pages) "\\end{landscape}" else NULL,
            "",
            # Weighted
            "\\subsubsection{Weighted}",
            if (rotate_pages) "\\begin{landscape}" else NULL,
            paste0("\\input{tables/", fn_w, "}"),
            if (rotate_pages) "\\end{landscape}" else NULL,
            "")
}

footer <- c("\\end{document}")
writeLines(c(header, body, footer), "main.tex")

message("Done. Split tables with 3-line headers written to 'tables/', and 'main.tex' created.")
