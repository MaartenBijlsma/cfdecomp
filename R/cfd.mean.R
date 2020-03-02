#'
#' @title Mean Decomposition: parametric version
#'
#' @description Decompose the mean difference in outcome Y between groups.
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
#' @param alpha the alpha level used to construct confidence intervals (0.05 = 95 percent confidence interval).
#' @param sample.resid if the \code{mediator} is Gaussian, should the simulation sample from the residuals of the linear regression model of the mediator to approximate the empirical distribution of the mediator in the simulation (Monte Carlo integration) (if so, set to \code{TRUE}), or should it sample from a Gaussian distribution with the standard deviation of the mediator? If the true distribution of the continuous mediator is not very Gaussian, the former may be preferred.
#' @param print.iteration print the bootstrap iteration
#'
#' @return \code{out_nc_m} returns the mean level of the mediator under the natural course, which is a value that should be close to the empirically observed value of the mediator for each group. \code{out_nc_quantile} provides the \code{alpha/2} and \code{1-alpha/2} bootstrap quantiles for this mean (AKA bootstrap percentile confidence intervals). \code{out_nc_y} and \code{out_nc_quantile_y} provide the corresponding values, but then for the outcome variable Y. Similarly, \code{out_cf_m}, \code{out_cf_quantile_m},\code{out_cf_y}, and \code{out_cf_quantile_y} provide the corresponding values for the counterfactual scenario where the mediators of the groups are equalized. \code{mediation} returns the proportion mediated by setting the intervened on mediator to be equal in level to the reference group and \code{mediation_quantile} returns the 1-alpha confidence interval. \code{mc_conv_info_m} and \code{mc_conv_info_y} provide information that can help determine the number of Monte Carlo and Bootstrap iterations needed to achieve stability. See the \code{Examples} for more information.
#' @export
#'
#' @examples
#' \dontrun{
#' mean.results.1 <- cfd.mean(formula.y='out.gauss ~ SES + med.gauss + med.binom + age',
#' formula.m='med.gauss ~ SES + age',
#' mediator='med.gauss',
#' group='SES',
#' data=cfd.example.data,
#' family.y='gaussian',
#' family.m='gaussian',
#' bs.size=250,
#' mc.size=50,
#' alpha=0.05)
#' # the differences between SES groups 1 and 2 were first:
#' mean(mean.results.1$out_nc_y[,2] - mean.results.1$out_nc_y[,1])
#' # and after giving the gaussian mediator of SES group 2 the distribution of the one in group 1
#' # the difference becomes:
#' mean(mean.results.1$out_cf_y[,2] - mean.results.1$out_nc_y[,1])
#' # so the % of the outcome Y that is due to differences between the two SES groups
#' # in the gaussian mediator is
#' mean(1-(mean.results.1$out_cf_y[,2] - mean.results.1$out_nc_y[,1]) /
#' (mean.results.1$out_nc_y[,2] - mean.results.1$out_nc_y[,1]))
#' # we can also get this number, and the one from the comparison of the other SES group
#' # with group 1, straight from the object
#' mean.results.1$mediation
#' # and we can get the 1-alpha CI for each:
#' mean.results.1$mediation_quantile
#' }
#' @import stats
#'
#'
cfd.mean <- function(formula.y,formula.m,mediator,group,
                     data,
                     family.y = 'binomial',
                     family.m = 'binomial',
                     bs.size = 1000,
                     mc.size=50,
                     alpha=0.05,
                     sample.resid=FALSE,
                     print.iteration=FALSE) {

  ## envir check
  call <- match.call()
  if (missing(data))
    data <- environment(formula.y)
  mf <- match.call(expand.dots = FALSE)

  ## get initial fit of outcome y
  ini_fit.y <- glm(formula = formula.y,
                   data=data, family=family.y)
  ini_fit.y$model_matrix <- model.matrix(as.formula(formula.y), data = data)
  fam.y <- family(ini_fit.y)

  ## get initial fit of mediator m
  ini_fit.m <- glm(formula = formula.m,
                   data=data, family=family.m)
  ini_fit.m$model_matrix <- model.matrix(as.formula(formula.m), data = data)
  fam.m <- family(ini_fit.m)

  ## loop over bootstrap samples
  out_nc_m <- matrix(NA,nrow=bs.size,ncol=nlevels(data[,group]),dimnames = list(1:bs.size,levels(data[,group])))
  out_cf_m <- matrix(NA,nrow=bs.size,ncol=nlevels(data[,group]),dimnames = list(1:bs.size,levels(data[,group])))
  out_nc_y <- matrix(NA,nrow=bs.size,ncol=nlevels(data[,group]),dimnames = list(1:bs.size,levels(data[,group])))
  out_cf_y <- matrix(NA,nrow=bs.size,ncol=nlevels(data[,group]),dimnames = list(1:bs.size,levels(data[,group])))

  ## temp matrix to save mc results
  temp_nc_m <- matrix(NA, nrow=mc.size, ncol=nlevels(data[,group]))
  temp_cf_m <- matrix(NA, nrow=mc.size, ncol=nlevels(data[,group]))
  temp_nc_y <- matrix(NA, nrow=mc.size, ncol=nlevels(data[,group]))
  temp_cf_y <- matrix(NA, nrow=mc.size, ncol=nlevels(data[,group]))

  for(i in 1:bs.size) {
    if(print.iteration==TRUE){print(i)}

    ## resampling
    resample_ind <- sample(nrow(data), replace=TRUE)
    data_bs <- data[resample_ind,]

    ## refit models for y and m
    bs_fit.y <- update(ini_fit.y, data = data_bs)
    coef_bs.y <- coef(bs_fit.y)

    bs_fit.m <- update(ini_fit.m, data = data_bs)
    coef_bs.m <- coef(bs_fit.m)
    if(family.m[1]=='gaussian' & sample.resid==FALSE) {sd.ref.m <- sd(bs_fit.m$residuals[data_bs[,group]==levels(data_bs[,group])[1]])}
    if(family.m[1]=='gaussian' & sample.resid==TRUE) {resid.ref.m <- bs_fit.m$residuals[data_bs[,group]==levels(data_bs[,group])[1]]}

    # save group identifier in an additional column since
    # we overwrite it to create counterfactuals
    data_bs$truegroup <- data_bs[,group]

    ## start Monte Carlo loop
    ## ! The parametric version does not equalize mediator values over strata
    for(ii in 1:mc.size){

      # save and overwrite data_mc since we create counterfactuals later
      data_mc <- data_bs

      ## natural course simulation ##
      # simulate mediator
      pred_mean_m <- predict(bs_fit.m,data_mc,type='response')
      n <- length(pred_mean_m)
      if(family.m[1]=='poisson') {data_mc[,mediator] <- rpois(n,pred_mean_m)} else if
      (family.m[1]=='binomial') {data_mc[,mediator] <- rbinom(n,1,pred_mean_m)} else if
      (family.m[1]=='gaussian' & sample.resid==FALSE) {data_mc[,mediator] <- rnorm(n,pred_mean_m,sd=sd.ref.m)} else if
      (family.m[1]=='gaussian' & sample.resid==TRUE) {data_mc[,mediator] <- pred_mean_m + sample(resid.ref.m,n,replace=TRUE)}
      temp_nc_m[ii,] <-  tapply(data_mc[,mediator], list(data_mc[,group]),mean,na.rm=T)

      # simulate outome
      # ! Note that this function predicts means only because it compares means
      # ! Hence the entire distribution of values (of the outcome) is not re-created
      # ! If that is desired, see the function for the quantiles instead.
      pred_mc_nc <- predict(bs_fit.y,data_mc,type='response')
      temp_nc_y[ii,] <-  tapply(pred_mc_nc, list(data_mc[,group]),mean,na.rm=T)

      ## counterfactual simulation ##
      # equalize distribution of mediator
      # this is done by assigning everyone the reference group and then (re)predicting
      data_mc[,group][which(data_mc[,group] %in% levels(data_mc[,group])[-1])] <- levels(data_mc[,group])[1]
      pred_mean_m <- predict(bs_fit.m,data_mc,type='response')
      if(family.m[1]=='poisson') {data_mc[,mediator] <- rpois(n,pred_mean_m)} else if
      (family.m[1]=='binomial') {data_mc[,mediator] <- rbinom(n,1,pred_mean_m)} else if
      (family.m[1]=='gaussian' & sample.resid==FALSE) {data_mc[,mediator] <- rnorm(n,pred_mean_m,sd=sd.ref.m)} else if
      (family.m[1]=='gaussian' & sample.resid==TRUE) {data_mc[,mediator] <- pred_mean_m + sample(resid.ref.m,n,replace=TRUE)}
      temp_cf_m[ii,] <-  tapply(data_mc[,mediator], list(data_mc$truegroup),mean,na.rm=T)

      # set group identifier back to true levels
      data_mc[,group] <- data_mc$truegroup
      # simulate outcome
      pred_mc_cf <- predict(bs_fit.y,data_mc,type='response')
      temp_cf_y[ii,] <-  tapply(pred_mc_cf, list(data_mc$truegroup),mean,na.rm=T)

    }

    # for thef first bootstrap iteration, also save mc information
    # to be used for convergence information
    if(i==1) {
      mc_conv_info_m <- apply(temp_nc_m,2,conv.mean)
      mc_conv_info_y <- apply(temp_nc_y,2,conv.mean)
    }

    ## save results for this bootstrap iteration
    out_nc_m[i,] <- apply(temp_nc_m,2,mean,na.rm=T)
    out_cf_m[i,] <- apply(temp_cf_m,2,mean,na.rm=T)
    out_nc_y[i,] <- apply(temp_nc_y,2,mean,na.rm=T)
    out_cf_y[i,] <- apply(temp_cf_y,2,mean,na.rm=T)

  }
  return(list(out_nc_m=out_nc_m,
              out_cf_m=out_cf_m,
              out_nc_quantile_m=apply(out_nc_m,2,quantile,c(alpha/2,0.5,1-alpha/2)),
              out_cf_quantile_m=apply(out_cf_m,2,quantile,c(alpha/2,0.5,1-alpha/2)),

              out_nc_y=out_nc_y,
              out_cf_y=out_cf_y,
              out_nc_quantile_y=apply(out_nc_y,2,quantile,c(alpha/2,0.5,1-alpha/2)),
              out_cf_quantile_y=apply(out_cf_y,2,quantile,c(alpha/2,0.5,1-alpha/2)),

              mediation=apply(1-(out_cf_y - out_nc_y[,1]) / (out_nc_y - out_nc_y[,1]),2,mean)[-1],
              mediation_quantile=apply(1-(out_cf_y - out_nc_y[,1]) / (out_nc_y - out_nc_y[,1]),2,quantile,probs=c(alpha/2,1-alpha/2),na.rm=T)[,-1],

              mc_conv_info_m=mc_conv_info_m,
              mc_conv_info_y=mc_conv_info_y
  )
  )
}
