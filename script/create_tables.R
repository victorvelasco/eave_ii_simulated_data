# --- Packages ---
library(dplyr)
library(stringr)
library(knitr)
library(kableExtra)
library(glue)
library(officer)
library(flextable)

# ============= SETTINGS =============
rotate_pages <- FALSE   # TRUE wraps each table in \begin{landscape}...\end{landscape}
options(knitr.kable.NA = "")

# ============= 1) CLEAN NUMERIC COLUMNS =============
num_cols <- c("True value","Estimate (unweighted)","SD (unweighted)",
              "Estimate (weighted)","SD (weighted)","rho")

results <- readRDS("results.RDS")
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

# ============= 4) BUILD ONE COMBINED TABLE PER OUTCOME =============
out_dir <- "tables"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Combined column set: shared cols + unweighted cols + weighted cols
cols_combined <- c(
  "parameter", "rho_disp", "True", "N",
  "unw_mean_sd", "unw_mean_se", "unw_bias", "unw_rmse", "unw_cov95",
  "w_mean_sd",   "w_mean_se",   "w_bias",   "w_rmse",   "w_cov95"
)

# Column names (2-line with makecell)
col_names_combined <- c(
  h3("", "Parameter"),
  h3("", "$\\rho$"),
  h3("", "True"),
  h3("", "$N$"),
  h3("Unweighted", "Mean (MC SD)"),
  h3("Unweighted", "Mean SE"),
  h3("Unweighted", "Bias"),
  h3("Unweighted", "RMSE"),
  h3("Unweighted", "Cov 95\\%"),
  h3("Weighted", "Mean (MC SD)"),
  h3("Weighted", "Mean SE"),
  h3("Weighted", "Bias"),
  h3("Weighted", "RMSE"),
  h3("Weighted", "Cov 95\\%")
)

# Alignment
align_combined <- c("c", "c", rep("r", length(cols_combined) - 2))

outcomes <- sort(unique(summ_fmt$outcome))
files_map <- list()

for (oc in outcomes) {
  base  <- str_to_lower(str_replace_all(oc, "[^A-Za-z0-9]+", "_"))
  df_oc <- filter(summ_fmt, outcome == oc) %>% arrange(parameter, rho_num)
  df_combined <- select(df_oc, all_of(cols_combined))

  cap <- glue("{oc}: Unweighted and weighted summary across {df_oc$N[1]} simulations per cell.")
  file_combined <- file.path(out_dir, paste0("results_", base, ".tex"))

  tbl <- kable(
    df_combined,
    format    = "latex",
    booktabs  = TRUE,
    caption   = cap,
    label     = paste0(base, "-combined"),   # kable prepends tab: automatically
    col.names = col_names_combined,
    align     = align_combined,
    escape    = FALSE,
    longtable = TRUE
  ) |>
    kable_styling(latex_options = c("striped"),
                  full_width = FALSE, position = "center") |>
    add_header_above(
      c(" " = 4, "Unweighted" = 5, "Weighted" = 5),
      bold = TRUE, escape = FALSE
    ) |>
    collapse_rows(columns = 1, latex_hline = "major")

  save_kable(tbl, file = file_combined)
  files_map[[oc]] <- basename(file_combined)
}

# ============= 5) WRITE main.tex (one combined table per outcome) =============
header <- c(
  "\\documentclass[12pt,a4paper]{article}",
  "\\usepackage{geometry}",
  "\\geometry{margin=0.75in}",
  "\\usepackage{amsmath}",
  "\\usepackage{booktabs}",
  "\\usepackage{array}",
  "\\usepackage{longtable}",
  "\\usepackage{multirow}",
  "\\usepackage{makecell}",
  "\\usepackage{threeparttable}",
  "\\usepackage{threeparttablex}",
  "\\usepackage[table]{xcolor}",
  "\\usepackage{pdflscape}",
  "\\usepackage{adjustbox}",
  "\\usepackage{hyperref}",
  "",
  "\\title{Simulation Study Results}",
  "\\author{}",
  "\\date{\\today}",
  "",
  "\\begin{document}",
  "\\maketitle",
  "",
  "\\section{Results by Outcome}"
)

body <- character(0)
for (oc in outcomes) {
  fn <- files_map[[oc]]
  body <- c(body,
            paste0("\\subsection{", oc, "}"),
            "\\begin{landscape}",
            paste0("\\input{tables/", fn, "}"),
            "\\end{landscape}",
            "")
}

footer <- "\\end{document}"
writeLines(c(header, body, footer), "main.tex")

message("Done. One combined table per outcome written to 'tables/', assembled in 'main.tex'.")

# ============= 6) WORD DOCUMENT =============

# Build a clean (no LaTeX) data frame from the raw numeric summary for Word output
make_word_df <- function(df_summ) {
  df_summ %>%
    transmute(
      Parameter = gsub("a_(\\d+)", "a₀", parameter) %>%   # a_0 -> a₀, a_1 -> a₁
                  { ifelse(grepl("0", parameter), "a₀", "a₁") },
      rho       = as.character(rho_num),
      True      = round(True, 4),
      N         = N,
      unw_mean  = sprintf("%.4f (%.4f)", unw_mean, unw_mc_sd),
      unw_se    = round(unw_mean_se, 4),
      unw_bias  = round(unw_bias, 4),
      unw_rmse  = round(unw_rmse, 4),
      unw_cov   = sprintf("%.1f%%", 100 * unw_cov95),
      w_mean    = sprintf("%.4f (%.4f)", w_mean, w_mc_sd),
      w_se      = round(w_mean_se, 4),
      w_bias    = round(w_bias, 4),
      w_rmse    = round(w_rmse, 4),
      w_cov     = sprintf("%.1f%%", 100 * w_cov95)
    )
}

ft_col_keys   <- c("Parameter","rho","True","N",
                   "unw_mean","unw_se","unw_bias","unw_rmse","unw_cov",
                   "w_mean",  "w_se", "w_bias",  "w_rmse",  "w_cov")
ft_col_labels <- c("Parameter", "ρ", "True value", "N",
                   "Mean (MC SD)", "Mean SE", "Bias", "RMSE", "Cov 95%",
                   "Mean (MC SD)", "Mean SE", "Bias", "RMSE", "Cov 95%")

make_flextable <- function(df_word, caption_text) {
  ft <- flextable(df_word, col_keys = ft_col_keys) %>%
    set_header_labels(values = setNames(as.list(ft_col_labels), ft_col_keys)) %>%
    add_header_row(
      values    = c("", "", "", "", "Unweighted", "Weighted"),
      colwidths = c(1,  1,  1,  1,  5,            5)
    ) %>%
    bold(part = "header") %>%
    bg(part = "header", bg = "#D9E1F2") %>%
    bg(i = seq(2, nrow(df_word), by = 2), bg = "#F5F5F5", part = "body") %>%
    align(align = "center", part = "header") %>%
    align(j = 1:2, align = "center", part = "body") %>%
    align(j = 3:14, align = "right", part = "body") %>%
    merge_v(j = "Parameter", part = "body") %>%
    valign(j = "Parameter", valign = "top", part = "body") %>%
    border_outer(part = "all",    border = fp_border(color = "#2E5FA3", width = 1.5)) %>%
    border_inner_h(part = "body", border = fp_border(color = "#BFBFBF", width = 0.5)) %>%
    border_inner_v(part = "body", border = fp_border(color = "#BFBFBF", width = 0.5)) %>%
    hline(part = "header",        border = fp_border(color = "#2E5FA3", width = 1.5)) %>%
    fontsize(size = 9, part = "all") %>%
    font(fontname = "Calibri", part = "all") %>%
    set_table_properties(layout = "autofit", width = 1) %>%
    set_caption(caption = as_paragraph(
      as_chunk(caption_text, props = fp_text(bold = FALSE, font.size = 10,
                                              font.family = "Calibri"))
    ))
  ft
}

doc <- read_docx()

for (i in seq_along(outcomes)) {
  oc     <- outcomes[i]
  df_s   <- filter(summ, outcome == oc) %>% arrange(parameter, rho_num)
  df_w   <- make_word_df(df_s)
  cap    <- paste0("Table ", i, ". Simulation results (unweighted and weighted) for the ",
                   oc, " outcome. Based on ", df_s$N[1], " simulations per cell.")
  ft     <- make_flextable(df_w, cap)

  doc <- body_add_par(doc, "", style = "Normal")
  doc <- body_add_flextable(doc, ft, align = "center")
  # Each table gets its own landscape section (= page break + landscape)
  doc <- body_end_section_landscape(doc)
}

print(doc, target = "simulation_results.docx")
message("Word document written to simulation_results.docx")
