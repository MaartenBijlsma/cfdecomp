
#' Flexible Function Decomposition: decompose any function that returns a vector
#'
#' @param formula.y the \code{\link{formula}} for the multivariable model (see \code{\link{glm}}) for the outcome Y.
#' @param formula.m the \code{\link{formula}} for the multivariable model (see \code{\link{glm}}) for the mediator M.
#' @param mediator the column name of the mediator M.
#' @param group column name of a factor variable containing the group identifier.
#' @param data a data frame containing the variables in the model.
#' @param family.y a description of the error distribution to be used in the model, see \code{\link{family}} for details. For the outcome variable any member of the \code{glm} family can be used.
#' @param family.m a description of the error distribution to be used in the model, see \code{\link{family}} for details. For the mediator, currently \code{gaussian}, \code{binomial} and \code{poisson} are supported.
#' @param bs.size the number of bootstrap iterations to be performed.
#' @param mc.size the number of Monte Carlo iterations to be performed (more = more MC error reduction).
#' @param FUN.y a function to compute the statistics which can be applied to all data subsets, this function should return a vector and should be ran on pred_y (simulated y values in the natural course or counterfactual) and optional additional columns.
#' @param alpha the alpha level used to construct confidence intervals (0.05 = 95 percent confidence interval).
#' @param cluster.sample set to TRUE if data are clustered in the long format (i.e. multiple rows per individual or other cluster).
#' @param cluster.name the name (as a character) of the column containing the cluster identifiers.
#' @param cluster.mrows for the mediator model, only allows 1 observation per mediator so that the mediator model is not weighted by number of observations. e.g. set to TRUE if the mediator is time constant in longitudinal analysis of long format data.
#' @param sample.resid.y sample.resid if the \code{outcome} is Gaussian, should the simulation sample from the residuals of the linear regression model of the outcome to approximate the empirical distribution of the outcome in the simulation (Monte Carlo integration) (if so, set to \code{TRUE}), or should it sample from a Gaussian distribution with the standard deviation of the outcome? If the true distribution of the continuous outcome is not very Gaussian, the former may be preferred.
#' @param sample.resid.m sample.resid if the \code{mediator} is Gaussian, should the simulation sample from the residuals of the linear regression model of the mediator to approximate the empirical distribution of the mediator in the simulation (Monte Carlo integration) (if so, set to \code{TRUE}), or should it sample from a Gaussian distribution with the standard deviation of the mediator? If the true distribution of the continuous mediator is not very Gaussian, the former may be preferred.
#' @param print.iteration print the bootstrap iteration
#' @param ... further arguments passed to or used by methods.
#'
#' @return \code{out_nc_m} returns the mean level of the mediator under the natural course, which is a value that should be close to the empirically observed value of the mediator for each group. \code{out_nc_quantile} provides the \code{alpha/2} and \code{1-alpha/2} bootstrap quantiles for this mean (AKA bootstrap percentile confidence intervals). \code{out_nc_y} provides the output of the function fed into FUN.y for each bootstrap iteration, with \code{out_nc_quantile_y} providing the \code{alpha/2} and \code{1-alpha/2} bootstrap quantiles of that output. Similarly, \code{out_cf_m}, \code{out_cf_quantile_m},\code{out_cf_y}, and \code{out_cf_quantile_y} provide the corresponding values for the counterfactual scenario where the mediators of the groups are equalized. \code{mediation} and \code{mediation_quantile} are not provided for this function, so should be calculated by the user based on the output. \code{mc_conv_info_m} and \code{mc_conv_info_y} provide information that can help determine the number of Monte Carlo and Bootstrap iterations needed to achieve stability. See the \code{Examples} for more information.
#' @export
#'
#' @examples
#' set.seed(100)
#' # the decomposition functions in our package are computationally intensive
#' # to make the example run quick, I perform it on a subsample (n=250) of the data:
#' cfd.example.sample <- cfd.example.data[sample(250),]
#' # define some function (here one that calculates the mean from the data)
#' # such a function already exists, but this is to demonstrate how to do it for one that
#' # will be implemented in cfd.FUN:
#' mean.fun <- function(data,yname) {
#' x <- data
#' return(mean(x[,yname],na.rm=TRUE))
#' }
#' # test if the function works on normal data:
#' mean.fun(cfd.example.sample,yname="med.pois")
#' # then enter it into cfd.FUN and run:
#' mean.results <- cfd.FUN(formula.y='out.gauss ~ SES + med.gauss + med.binom + age',
#'                           formula.m='med.gauss ~ SES + age',
#'                           mediator='med.gauss',
#'                           group='SES',
#'                           data=cfd.example.sample,
#'                           family.y='gaussian',
#'                           family.m='gaussian',
#'                           FUN.y=mean.fun,
#'                           bs.size=15,
#'                           mc.size=5,
#'                           alpha=0.05,
#'                           print.iteration=TRUE,
#'                           yname="pred_y")
#' # more advanced code demonstrating how to do this with a function that calculates
#' # the age-adjusted rate ratio and life expectancy will hopefully soon be available
#' # in a publication.
#' #' @import stats utils
cfd.FUN <- function (formula.y, formula.m, mediator, group,data,
                     family.y = "binomial",
                     family.m = "binomial",
                     bs.size = 250,
                     mc.size = 50,
                     FUN.y=mean,
                     alpha = 0.05,
                     cluster.sample = FALSE,
                     cluster.name = NA,
                     cluster.mrows=FALSE,
                     sample.resid.y = FALSE,
                     sample.resid.m = FALSE,
                     print.iteration = FALSE,...)
{
  call <- match.call()
  if (missing(data))
    data <- environment(formula.y)
  mf <- match.call(expand.dots = FALSE)
  ini_fit.y <- glm(formula = formula.y, data = data, family = family.y)
  ini_fit.y$model_matrix <- model.matrix(as.formula(formula.y),
                                         data = data)
  fam.y <- family(ini_fit.y)
  if(cluster.mrows==FALSE) {
    ini_fit.m <- glm(formula = formula.m, data = data, family = family.m)
    ini_fit.m$model_matrix <- model.matrix(as.formula(formula.m),
                                           data = data)
  } else if(cluster.mrows==TRUE) {
    data.m <- do.call(rbind,by(data, list(data[,cluster.name]),FUN=function(x) head(x,1)))
    ini_fit.m <- glm(formula = formula.m, data = data.m, family = family.m)
    ini_fit.m$model_matrix <- model.matrix(as.formula(formula.m),
                                           data = data.m)
  }
  fam.m <- family(ini_fit.m)

  out_nc_m <- matrix(NA, nrow = bs.size, ncol = nlevels(data[,
                                                             group]), dimnames = list(1:bs.size, levels(data[, group])))
  out_cf_m <- matrix(NA, nrow = bs.size, ncol = nlevels(data[,
                                                             group]), dimnames = list(1:bs.size, levels(data[, group])))
  temp_nc_m <- matrix(NA, nrow = mc.size, ncol = nlevels(data[,
                                                              group]))
  temp_cf_m <- matrix(NA, nrow = mc.size, ncol = nlevels(data[,
                                                              group]))

  # ! let's determine the length of the matrix needed to store results from FUN.y
  test.data <- data

  pred_mean_y <- predict(ini_fit.y, test.data, type = "response")
  n <- length(pred_mean_y)
  if (family.y[1] == "gaussian" & sample.resid.y ==
      FALSE) {
    sd.ref.y <- sd(ini_fit.y$residuals[test.data[, group] ==
                                         levels(test.data[, group])[1]])
  }
  if (family.y[1] == "gaussian" & sample.resid.y ==
      TRUE) {
    resid.ref.y <- ini_fit.y$residuals[test.data[, group] ==
                                         levels(test.data[, group])[1]]
  }
  if (family.y[1] == "poisson") {
    pred_y <- rpois(n, pred_mean_y)
  }      else if (family.y[1] == "binomial") {
    pred_y <- rbinom(n, 1, pred_mean_y)
  }      else if (family.y[1] == "gaussian" & sample.resid.y ==
                  FALSE) {
    pred_y <- rnorm(n, pred_mean_y, sd = sd.ref.y)
  }      else if (family.y[1] == "gaussian" & sample.resid.y ==
                  TRUE) {
    pred_y <- pred_mean_y + sample(resid.ref.y, n,
                                   replace = TRUE)
  }
  test.data$pred_y <- pred_y
  test.vec <- FUN.y(test.data,...)
  rm(test.data)
  rm(pred_mean_y)
  vec.length <- length(test.vec)

  out_nc_y <- array(NA,dim=c(bs.size,vec.length))
  out_cf_y <- array(NA,dim=c(bs.size,vec.length))
  temp_nc_y <- array(NA,dim=c(mc.size,vec.length))
  temp_cf_y <- array(NA,dim=c(mc.size,vec.length))

  for (i in 1:bs.size) {
    if (print.iteration == TRUE) {
      print(i)
    }
    if (cluster.sample == FALSE) {
      resample_ind <- sample(nrow(data), replace = TRUE)
      data_bs <- data[resample_ind, ]
    } else if (cluster.sample == TRUE) {
      data_bs <- cluster.resample(data, cluster.name, size=length(unique(data[,cluster.name])))
    }
    bs_fit.y <- update(ini_fit.y, data = data_bs)
    coef_bs.y <- coef(bs_fit.y)
    if (family.y[1] == "gaussian" & sample.resid.y ==
        FALSE) {
      sd.ref.y <- sd(bs_fit.y$residuals[data_bs[, group] ==
                                          levels(data_bs[, group])[1]])
    }
    if (family.y[1] == "gaussian" & sample.resid.y ==
        TRUE) {
      resid.ref.y <- bs_fit.y$residuals[data_bs[, group] ==
                                          levels(data_bs[, group])[1]]
    }

    if(cluster.mrows==FALSE) {
      bs_fit.m <- update(ini_fit.m, data = data_bs)
      coef_bs.m <- coef(bs_fit.m)
    } else if(cluster.mrows==TRUE) {
      data_bs.m <- do.call(rbind,by(data_bs, list(data_bs[,cluster.name]),FUN=function(x) head(x,1)))
      bs_fit.m <- update(ini_fit.m, data = data_bs.m)
      coef_bs.m <- coef(bs_fit.m)
    }
    if (family.m[1] == "gaussian" & sample.resid.m ==
        FALSE) {
      sd.ref.m <- sd(bs_fit.m$residuals[data_bs[, group] ==
                                          levels(data_bs[, group])[1]])
    }
    if (family.m[1] == "gaussian" & sample.resid.m == TRUE & cluster.mrows == FALSE) {
      resid.ref.m <- bs_fit.m$residuals[data_bs[, group] ==
                                          levels(data_bs[, group])[1]]
    }
    if (family.m[1] == "gaussian" & sample.resid.m == TRUE & cluster.mrows == TRUE) {
      resid.ref.m <- bs_fit.m$residuals[data_bs.m[, group] ==
                                          levels(data_bs.m[, group])[1]]
    }
    data_bs$truegroup <- data_bs[, group]
    for (ii in 1:mc.size) {
      data_mc <- data_bs
      pred_mean_m <- predict(bs_fit.m, data_mc, type = "response")
      n <- length(pred_mean_m)
      if (family.m[1] == "poisson") {
        data_mc[, mediator] <- rpois(n, pred_mean_m)
      }      else if (family.m[1] == "binomial") {
        data_mc[, mediator] <- rbinom(n, 1, pred_mean_m)
      }      else if (family.m[1] == "gaussian" & sample.resid.m ==
                      FALSE) {
        data_mc[, mediator] <- rnorm(n, pred_mean_m,
                                     sd = sd.ref.m)
      }      else if (family.m[1] == "gaussian" & sample.resid.m ==
                      TRUE) {
        data_mc[, mediator] <- pred_mean_m + sample(resid.ref.m,
                                                    n, replace = TRUE)
      }
      temp_nc_m[ii, ] <- tapply(data_mc[, mediator], list(data_mc[,
                                                                  group]), mean, na.rm = TRUE)
      pred_mean_y <- predict(bs_fit.y, data_mc, type = "response")
      data_mc$pred_mean_y <- pred_mean_y
      if (family.y[1] == "poisson") {
        pred_y <- rpois(n, pred_mean_y)
      }      else if (family.y[1] == "binomial") {
        pred_y <- rbinom(n, 1, pred_mean_y)
      }      else if (family.y[1] == "gaussian" & sample.resid.y ==
                      FALSE) {
        pred_y <- rnorm(n, pred_mean_y, sd = sd.ref.y)
      }      else if (family.y[1] == "gaussian" & sample.resid.y ==
                      TRUE) {
        pred_y <- pred_mean_y + sample(resid.ref.y, n,
                                       replace = TRUE)
      }
      # ! append pred_y to data_mc so that the other info in there
      # can also be used by FUN.y
      data_mc$pred_y <- pred_y
      temp_nc_y[ii,] <- FUN.y(data_mc,...)
      data_mc$pred_y <- NA
      data_mc$pred_mean_y <- NA
      data_mc[, group][which(data_mc[, group] %in% levels(data_mc[,
                                                                  group])[-1])] <- levels(data_mc[, group])[1]
      pred_mean_m <- predict(bs_fit.m, data_mc, type = "response")
      if (family.m[1] == "poisson") {
        data_mc[, mediator] <- rpois(n, pred_mean_m)
      }      else if (family.m[1] == "binomial") {
        data_mc[, mediator] <- rbinom(n, 1, pred_mean_m)
      } else if (family.m[1] == "gaussian" & sample.resid.m ==
                 FALSE) {
        data_mc[, mediator] <- rnorm(n, pred_mean_m,
                                     sd = sd.ref.m)
      }       else if (family.m[1] == "gaussian" & sample.resid.m ==
                       TRUE) {
        data_mc[, mediator] <- pred_mean_m + sample(resid.ref.m,
                                                    n, replace = TRUE)
      }
      temp_cf_m[ii, ] <- tapply(data_mc[, mediator], list(data_mc$truegroup),
                                mean, na.rm = TRUE)
      data_mc[, group] <- data_mc$truegroup
      pred_mean_y <- predict(bs_fit.y, data_mc, type = "response")
      data_mc$pred_mean_y <- pred_mean_y
      if (family.y[1] == "poisson") {
        pred_y <- rpois(n, pred_mean_y)
      }      else if (family.y[1] == "binomial") {
        pred_y <- rbinom(n, 1, pred_mean_y)
      }      else if (family.y[1] == "gaussian" & sample.resid.y ==
                      FALSE) {
        pred_y <- rnorm(n, pred_mean_y, sd = sd.ref.y)
      }      else if (family.y[1] == "gaussian" & sample.resid.y ==
                      TRUE) {
        pred_y <- pred_mean_y + sample(resid.ref.y, n,
                                       replace = TRUE)
      }
      # ! append pred_y to data_mc so that the other info in there
      # can also be used by FUN.y
      data_mc$pred_y <- pred_y
      temp_cf_y[ii,] <- FUN.y(data_mc,...)
      data_mc$pred_y <- NA
      data_mc$pred_mean_y <- NA
    }
    if (i == 1) {
      mc_conv_info_m <- apply(temp_nc_m, 2, conv.mean)
      mc_conv_info_y <- apply(temp_nc_y, 2, conv.mean)
    }
    out_nc_m[i, ] <- apply(temp_nc_m, 2, mean, na.rm = TRUE)
    out_cf_m[i, ] <- apply(temp_cf_m, 2, mean, na.rm = TRUE)
    out_nc_y[i, ] <- apply(temp_nc_y, 2, FUN=mean,na.rm=TRUE)
    out_cf_y[i, ] <- apply(temp_cf_y, 2, FUN=mean,na.rm=TRUE)
  }
  return(list(out_nc_m = out_nc_m,
              out_cf_m = out_cf_m,
              out_nc_quantile_m = apply(out_nc_m, 2, quantile, c(alpha/2, 0.5, 1 - alpha/2),na.rm=TRUE),
              out_cf_quantile_m = apply(out_cf_m, 2, quantile, c(alpha/2, 0.5, 1 - alpha/2),na.rm=TRUE),
              out_nc_y = out_nc_y,
              out_cf_y = out_cf_y,
              out_nc_quantile_y = apply(out_nc_y, 2, quantile, c(alpha/2, 0.5, 1 - alpha/2),na.rm=TRUE),
              out_cf_quantile_y = apply(out_cf_y, 2, quantile, c(alpha/2, 0.5, 1 - alpha/2),na.rm=TRUE),
              mc_conv_info_m = mc_conv_info_m,
              mc_conv_info_y = mc_conv_info_y))
}
