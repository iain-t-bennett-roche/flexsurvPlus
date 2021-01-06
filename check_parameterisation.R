adtte <- sim_adtte(seed = 2020, rho = 0.6)

PFS_data <- adtte %>%
  filter(PARAMCD=="PFS") %>%
  transmute(USUBJID,
            ARMCD,
            PFS_days = AVAL,
            PFS_event = 1- CNSR
  )

checkit <- runPSM(data=PFS_data,
                  time_var="PFS_days",
                  event_var="PFS_event",
                  model.type=c("Common shape"),
                  strata_var = "ARMCD",
                  int_name = "B",
                  ref_name = "A")

#OLD


summary(checkit$models$comshp.exp, type = "mean")
mean_exp(rate = exp(checkit$parameters$exp.rate.ref))
mean_exp(rate = exp(checkit$parameters$exp.rate.int))

summary(checkit$models$comshp.weibull, type = "mean")
mean_weibull(shape = exp(checkit$parameters$weibull.shape.ref), scale = exp(checkit$parameters$weibull.scale.ref))
mean_weibull(shape = exp(checkit$parameters$weibull.shape.int), scale = exp(checkit$parameters$weibull.scale.int))

summary(checkit$models$comshp.gompertz, type = "mean")
mean_gompertz(shape = checkit$parameters$gompertz.shape.ref, rate = exp(checkit$parameters$gompertz.rate.ref))
mean_gompertz(shape = checkit$parameters$gompertz.shape.int, rate = exp(checkit$parameters$gompertz.rate.int))

summary(checkit$models$comshp.lnorm, type = "mean")
mean_lnorm(meanlog = checkit$parameters$lnorm.meanlog.ref, sdlog = exp(checkit$parameters$lnorm.sdlog.ref))
mean_lnorm(meanlog = checkit$parameters$lnorm.meanlog.int, sdlog = exp(checkit$parameters$lnorm.sdlog.int))

summary(checkit$models$comshp.llogis, type = "mean")
mean_llogis(shape = exp(checkit$parameters$llogis.shape.ref), scale = exp(checkit$parameters$llogis.scale.ref))
mean_llogis(shape = exp(checkit$parameters$llogis.shape.int), scale = exp(checkit$parameters$llogis.scale.int))

summary(checkit$models$comshp.gengamma, type = "mean")
mean_gengamma(mu = checkit$parameters$gengamma.mu.ref, sigma = exp(checkit$parameters$gengamma.sigma.ref), Q = checkit$parameters$gengamma.Q.ref)
mean_gengamma(mu = checkit$parameters$gengamma.mu.int, sigma = exp(checkit$parameters$gengamma.sigma.int), Q = checkit$parameters$gengamma.Q.int)

summary(checkit$models$comshp.gamma, type = "mean")
mean_gamma(shape = exp(checkit$parameters$gamma.shape.ref), rate = exp(checkit$parameters$gamma.rate.ref))
mean_gamma(shape = exp(checkit$parameters$gamma.shape.int), rate = exp(checkit$parameters$gamma.rate.int))

#NEW

summary(checkit$models$comshp.exp, type = "mean")
mean_exp(rate = checkit$parameters$exp.rate.ref)
mean_exp(rate = checkit$parameters$exp.rate.int)

summary(checkit$models$comshp.weibull, type = "mean")
mean_weibull(shape = checkit$parameters$weibull.shape.ref, scale = checkit$parameters$weibull.scale.ref)
mean_weibull(shape = checkit$parameters$weibull.shape.int, scale = checkit$parameters$weibull.scale.int)

summary(checkit$models$comshp.gompertz, type = "mean")
mean_gompertz(shape = checkit$parameters$gompertz.shape.ref, rate = checkit$parameters$gompertz.rate.ref)
mean_gompertz(shape = checkit$parameters$gompertz.shape.int, rate = checkit$parameters$gompertz.rate.int)

summary(checkit$models$comshp.lnorm, type = "mean")
mean_lnorm(meanlog = checkit$parameters$lnorm.meanlog.ref, sdlog = checkit$parameters$lnorm.sdlog.ref)
mean_lnorm(meanlog = checkit$parameters$lnorm.meanlog.int, sdlog = checkit$parameters$lnorm.sdlog.int)

summary(checkit$models$comshp.llogis, type = "mean")
mean_llogis(shape = checkit$parameters$llogis.shape.ref, scale = checkit$parameters$llogis.scale.ref)
mean_llogis(shape = checkit$parameters$llogis.shape.int, scale = checkit$parameters$llogis.scale.int)

summary(checkit$models$comshp.gengamma, type = "mean")
mean_gengamma(mu = checkit$parameters$gengamma.mu.ref, sigma = checkit$parameters$gengamma.sigma.ref, Q = checkit$parameters$gengamma.Q.ref)
mean_gengamma(mu = checkit$parameters$gengamma.mu.int, sigma = checkit$parameters$gengamma.sigma.int, Q = checkit$parameters$gengamma.Q.int)

summary(checkit$models$comshp.gamma, type = "mean")
mean_gamma(shape = checkit$parameters$gamma.shape.ref, rate = checkit$parameters$gamma.rate.ref)
mean_gamma(shape = checkit$parameters$gamma.shape.int, rate = checkit$parameters$gamma.rate.int)