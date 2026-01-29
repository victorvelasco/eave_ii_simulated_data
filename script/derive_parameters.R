### visits <- c("6d", "13d", "20d", "27d")
###
### tmt_models <- list()
###
### treatments_coded[is.na(treatments_coded)] <- 1
###
### ## Tmt model for 0-6 days
### treatments1 <- treatments_coded |>
###   filter(endpoint_group == "any_throm_haem") |>
###   mutate(cases = `AZ_0-6`, controls = `AZ_7-13`+`AZ_14-20`+`AZ_21-27`+`AZ_27+`+Unvaccinated) |>
###   select(starts_with("AgeGroup"), SexMale, cases, controls)
### tmt_models[[1]] <- glm(cbind(cases, controls) ~ ., family = "binomial", data = treatments1)
###
### ## Tmt model for 7-13 days
### treatments1 <- treatments_coded |>
###   filter(endpoint_group == "any_throm_haem") |>
###   mutate(cases = `AZ_7-13`, controls = `AZ_14-20`+`AZ_21-27`+`AZ_27+`+Unvaccinated) |>
###   select(starts_with("AgeGroup"), SexMale, cases, controls)
### tmt_models[[2]] <- glm(cbind(cases, controls) ~ ., family = "binomial", data = treatments1)
###
### ## Tmt model for 14-20 days
### treatments1 <- treatments_coded |>
###   filter(endpoint_group == "any_throm_haem") |>
###   mutate(cases = `AZ_14-20`, controls = `AZ_21-27`+`AZ_27+`+Unvaccinated) |>
###   select(starts_with("AgeGroup"), SexMale, cases, controls)
### tmt_models[[3]] <- glm(cbind(cases, controls) ~ ., family = "binomial", data = treatments1)
###
### ## Tmt model for 21-27 days
### treatments1 <- treatments_coded |>
###   filter(endpoint_group == "any_throm_haem") |>
###   mutate(cases = `AZ_21-27`, controls = + `AZ_27+`+Unvaccinated) |>
###   select(starts_with("AgeGroup"), SexMale, cases, controls)
### tmt_models[[4]] <- glm(cbind(cases, controls) ~ ., family = "binomial", data = treatments1)
###
### gammas <- lapply(tmt_models, coef)

tmt_model <- glm(cbind(AZ, Unvaccinated) ~ . - PB - total, family = "binomial", data = demographics_coded)

gammas0 <- coef(tmt_model)
gammas0
gammas <- list(gammas0)

betas <- list(
  Thrombocytopenia = list(c(qlogis(0.005), log(1.42))),
  ITP              = list(c(qlogis(0.005), log(5.77))),
  Venous           = list(c(qlogis(0.005), log(1.03))),
  Arterial         = list(c(qlogis(0.005), log(1.22))),
  Hemorrhagic      = list(c(qlogis(0.005), log(1.48)))
)

### betas <- list(
###   c(qlogis(0.005), log(3.43)),
###   c(qlogis(0.005), log(4.60)),
###   c(qlogis(0.005), log(7.81)),
###   c(qlogis(0.005), log(14.07))
### )
### ### gammas <- list(gammas0)
### #   gammas0,
### #   c(gammas0, 0),
### #   c(gammas0, 0, 0),
### #   c(gammas0, 0, 0, 0)
### # )
###
