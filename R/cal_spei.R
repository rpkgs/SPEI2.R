#' @importFrom lmom cdfnor cdfgam cdfpe3 cdfglo
#' @importFrom lmom pelnor pelgam pelpe3 pelglo
#' @importFrom lmomco parnor pargam parpe3 parglo
#' @importFrom lmomco pwm2lmom are.lmom.valid
NULL

par_nor = setNames(rep(NA, 2), c("mu", "sigma"))
par_gam = setNames(rep(NA, 2), c("alpha", "beta"))
par_pe3 = setNames(rep(NA, 3), c("mu", "sigma", "gamma"))
par_glo = setNames(rep(NA, 3), c("xi", "alpha", "kappa"))
# par_gev = setNames(rep(NA, 3), c("xi", "alpha", "kappa"))

.options <- list(
  "normal"       = list(pel = pelnor, par = parnor, cdf = cdfnor, param = par_nor),
  "Gamma"        = list(pel = pelgam, par = pargam, cdf = cdfgam, param = par_gam),
  "PearsonIII"   = list(pel = pelpe3, par = parpe3, cdf = cdfpe3, param = par_pe3),
  "log-Logistic" = list(pel = pelglo, par = parglo, cdf = cdfglo, param = par_glo)
)

#' cal_spei
#'
#' @param x A monthly time series of certain month
#' @param distribution Distribution name, one of `normal`, `Gamma`, `PearsonIII`, `log-Logistic`.
#' @param fit Fitting method, one of `ub-pwm`, `pp-pwm`, `max-lik`.
#'
#' @example 
#' x = read.table("data-raw/data.txt")$V1
#' cal_spei(x)
#' @export
cal_spei <- function(x, distribution = "log-Logistic", fit = "ub-pwm", ...) {
  r <- list(z = x * NA, coef = .options[[distribution]]$param)

  x.mon <- x[!is.na(x)]
  if (distribution != "log-Logistic") {
    pze <- sum(x.mon == 0) / length(x.mon)
    x.mon <- x.mon[x.mon > 0]
  }

  # Get functions
  dist <- .options[[distribution]]
  pel_lmom <- dist$pel
  par_lmomco <- dist$par
  cdf <- dist$cdf

  # Fit distribution parameters
  if (length(x.mon) < 4) {
    return(r)
  }

  x.mon_sd <- sd(x.mon, na.rm = TRUE)
  if (is.na(x.mon_sd) || (x.mon_sd == 0)) {
    return(r)
  }

  # Calculate probability weighted moments based on `lmomco` or
  # `TLMoments`
  pwm <- switch(fit,
    "pp-pwm" = pwm.pp(x.mon, -0.35, 0, nmom = 3, sort = TRUE),
    "ub-pwm" = PWM(x.mon, order = 0:2)
  )

  # Check L-moments validity
  lmom <- pwm2lmom(pwm)
  if (!are.lmom.valid(lmom) || anyNA(lmom[[1]]) || any(is.nan(lmom[[1]]))) {
    return(r)
  }

  # `lmom` fortran functions need specific inputs L1, L2, T3
  # This is handled internally by `lmomco` with `lmorph`
  fortran_vec <- c(lmom$lambdas[1:2], lmom$ratios[3])

  # Calculate parameters based on distribution with `lmom`, then `lmomco`
  f_params <- tryCatch(
    pel_lmom(fortran_vec),
    error = function(e) {
      par_lmomco(lmom)$para
    }
  )

  # Adjust if user chose `log-Logistic` and `max-lik`
  if (distribution == "log-Logistic" && fit == "max-lik") {
    f_params <- parglo.maxlik(x.mon, f_params)$para
  }

  cdf_res <- cdf(x, f_params) # exclude nans
  z <- qnorm(cdf_res)

  # Adjust for `pze` if distribution is Gamma or PearsonIII
  if (distribution == "Gamma" | distribution == "PearsonIII") {
    z <- qnorm(pze + (1 - pze) * pnorm(z))
  }

  list(z = z, coef = f_params)
}
