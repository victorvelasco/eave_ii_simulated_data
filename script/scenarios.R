## sample.size <- function(outcome) {
##   switch(
##     outcome,
##     Thrombocytopenia = { 2343  + 407  + 414      },
##     ITP              = { 702   + 48   + 142      },
##     Venous.          = { 26843 + 3616 + 4110     },
##     Arterial         = { 67599 + 13925  + 13157  },
##     
##   )
## }

betas <- list(
  Thrombocytopenia = list(c(qlogis(0.005), log(1.42))),
  ITP              = list(c(qlogis(0.005), log(5.77))),
  Venous.          = list(c(qlogis(0.005), log(1.03))),
  Arterial.        = list(c(qlogis(0.005), log(1.22))),
  Hemorrhagic      = list(c(qlogis(0.005), log(1.48)))
)