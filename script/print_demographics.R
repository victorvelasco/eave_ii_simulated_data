library(arsenal)

# demographics: data frame with columns sex ("Female"/"Male"), age_group, total

# Expand aggregated counts to one row per individual so tableby can compute
# both N and column-% per cell.
demographics_long <- demographics |>
  select(sex, age_group, total) |>
  uncount(total) |>
  mutate(
    sex       = factor(sex, levels = c("Female", "Male")),
    age_group = factor(age_group)
  )

# Build table: rows = age_group, columns = Female / Male / Total.
# countpct -> "n (xx.x%)" in every cell. pct.col gives column percentages
# (each sex column sums to 100%); switch to "row" if you want row percentages.
tab <- tableby(
  sex ~ age_group,
  data    = demographics_long,
  total   = TRUE,
  control = tableby.control(
    cat.stats  = "countpct",
    cat.simplify = FALSE,
    pct.col    = TRUE,
    digits     = 0,
    digits.pct = 1
  ), 
  test = FALSE
)

# Console / R Markdown
summary(tab,
        labelTranslations = list(age_group = "Age group (years)"),
        text = TRUE)

# Word export
write2word(tab, file.path(getwd(), "tables/table1_age_by_sex.docx"),
           title = "Table 1. Baseline characteristics of the simulated population",
           labelTranslations = list(age_group = "Age group (years)"))
