#' @importFrom lmom cdfnor cdfgam cdfpe3 cdfglo
#' @importFrom lmom pelnor pelgam pelpe3 pelglo
#' @importFrom lmomco parnor pargam parpe3 parglo
#' @importFrom lmomco pwm2lmom are.lmom.valid pwm.pp
NULL

last <- \(x) x[[length(x)]]

# x -> cdf -> z
qGamma <- function(x, coef, cdf = cdfgam) {
  pze <- last(coef) # the probability of x <= 0
  npar <- length(coef) - 1
  cdf_res <- cdf(x, coef[1:npar]) # exclude nans
  z <- qnorm(pze + (1 - pze) * cdf_res)
  z[x <= 0] <- -Inf
  z
}

qNormal <- function(x, coef, cdf = cdfnor) {
  cdf(x, coef) |> qnorm()
}

qPearsonIII <- \(x, coef) qGamma(x, coef, cdf = cdfpe3)

qLogLogistic <- \(x, coef) qNormal(x, coef, cdf = cdfglo)


par_nor <- setNames(rep(NA, 2), c("mu", "sigma"))
par_gam <- setNames(rep(NA, 3), c("alpha", "beta", "pzero"))
par_pe3 <- setNames(rep(NA, 4), c("mu", "sigma", "gamma", "pzero"))
par_glo <- setNames(rep(NA, 3), c("xi", "alpha", "kappa"))
# par_gev = setNames(rep(NA, 3), c("xi", "alpha", "kappa"))

.options <- list(
  "Normal"       = list(pel = pelnor, par = parnor, inv = qNormal, param = par_nor),
  "Gamma"        = list(pel = pelgam, par = pargam, inv = qGamma, param = par_gam),
  "PearsonIII"   = list(pel = pelpe3, par = parpe3, inv = qPearsonIII, param = par_pe3),
  "log-Logistic" = list(pel = pelglo, par = parglo, inv = qLogLogistic, param = par_glo)
)

#' cal_spei
#'
#' @param x A monthly time series of certain month, which you want to calculate drought index.
#' @param x_ref the reference time series to fit distribution
#'
#' @param distribution Distribution name, one of `Normal`, `Gamma`, `PearsonIII`, `log-Logistic`.
#' @param fit Fitting method, one of `ub-pwm`, `pp-pwm`, `max-lik`.
#' @param ... ignored
#'
#' @examples
#' cal_spei(wb)
#' cal_spei(wb, fit = "max-lik")
#' cal_spei(wb, distribution = "Normal")
#' cal_spi(wb)
#' @export
cal_spei <- function(x, x_ref = x, distribution = "log-Logistic", fit = "ub-pwm", ...) {
  r <- list(z = x * NA, coef = .options[[distribution]]$param)

  x.mon <- x_ref[!is.na(x_ref)]
  if (distribution %in% c("Gamma", "PearsonIII")) {
    pze <- sum(x.mon <= 0) / length(x.mon)
    x.mon <- x.mon[x.mon > 0]
  }

  # Get functions
  dist <- .options[[distribution]]
  pel_lmom <- dist$pel
  par_lmomco <- dist$par
  inv <- dist$inv

  # Fit distribution parameters
  if (length(x.mon) < 4) {
    return(r)
  }

  x.mon_sd <- sd(x.mon, na.rm = TRUE)
  if (is.na(x.mon_sd) || (x.mon_sd == 0)) {
    return(r)
  }

  # Calculate probability weighted moments based on `lmomco` or `TLMoments`
  fit2 <- fit
  if (!(fit2 %in% c("pp-pwm", "ub-pwm"))) fit2 <- "ub-pwm" # init param for lmom

  pwm <- switch(fit2,
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
  coef <- tryCatch(
    pel_lmom(fortran_vec),
    error = function(e) {
      par_lmomco(lmom)$para
    }
  )

  # Adjust if user chose `log-Logistic` and `max-lik`
  if (distribution == "log-Logistic" && fit == "max-lik") {
    coef <- parglo.maxlik(x.mon, coef)$para
  }
  if (distribution %in% c("Gamma", "PearsonIII")) coef["pzero"] <- pze

  z <- inv(x, coef)
  list(z = z, coef = coef)
}

#' @rdname cal_spei
#' @export
cal_spi <- function(x, x_ref = x, distribution = "Gamma", fit = "ub-pwm", ...) {
  cal_spei(x, x_ref, distribution, fit, ...)
}
