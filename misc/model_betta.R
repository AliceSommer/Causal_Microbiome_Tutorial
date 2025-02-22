#' Modelling total diversity with betta
#' 
#' This function tests for heterogeneity of total diversity (observed plus
#' unobserved) across multiple sites. It can account or test for fixed effects
#' that may explain diversity. It returns the significance of the covariates in
#' explaining diversity and a hypothesis test for heterogeneity.
#' 
#' 
#' @param chats A vector of estimates of total diversity at different sampling
#' locations. \samp{breakaway} estimates are suggested in the high-diversity
#' case but not enforced. 
#' @param ses The standard errors in \code{chats}, the diversity estimates. This
#' can either be a vector of standard errors (with the arguments \code{chats} and 
#' \code{X}), or the name of the variable in the dataframe \code{data} that contains 
#' the standard errors (with the arguments \code{formula} and \code{data}). 
#' @param X A numeric matrix of covariates. If not supplied, an intercept-only
#' model will be fit. This is optional with the \code{chats} argument. 
#' @param formula A formula object of the form \eqn{y ~ x | group}. Required with 
#' the \code{data} argument. 
#' @param data A dataframe containing the response, response standard errors, covariates, 
#' and grouping variable. Required with the \code{formula} argument. 
#' @param initial_est (Optional) A vector of length 1 + ncol(X) giving the starting values for the likelihood maximisation search. The first element is the starting estimate for sigma^2_u, and the remaining elements are the starting elements for beta. Defaults to NULL, in which case the starting values outlined in the paper are used. 
#' @return \item{table}{ A coefficient table for the model parameters. The
#' columns give the parameter estimates, standard errors, and p-values,
#' respectively. This model is only as effective as your diversity estimation
#' procedure; for this reason please confirm that your estimates are
#' appropriate and that your model is not misspecified. betta_pic may be
#' useful for this purpose.  } \item{cov}{ Estimated covariance matrix of the
#' parameter estimates.  } \item{ssq_u}{ The estimate of the heterogeneity
#' variance.  } \item{homogeneity}{ The test statistic and p-value for the test
#' of homogeneity.  } \item{global}{ The test statistic and p-value for the
#' test of model explanatory power.  } \item{blups}{ The conditional expected
#' values of the diversity estimates (conditional on the random effects). The
#' authors propose that if the practitioner believes that information from one
#' diversity estimator may inform the others, then using the \samp{condfits} as
#' estimators of total diversity rather than \samp{Chats} may reduce variance
#' of diversity estimates by ``sharing strength'' across the samples.  }
#' \item{blupses}{ The estimated standard deviation (standard errors) in the
#' blups.  }
#' \item{loglikelihood}{ The log likelihood of the fitted model.  }
#' \item{aic}{ The Akaike information criterion for the fitted model. }
#' \item{aicc}{ The finite sample correction of the Akaike information criterion for the fitted model.  }
#' \item{r_squared_wls}{  The weighted R^2 statistic, appropriate for heteroskedastic linear models. }
#' @note Ecologists who are interested in the way species richness varies with
#' covariate information often run a regression-type analysis on the observed
#' diversity using their covariate information as predictors. However, in many
#' settings (especially microbial), rare and unobserved taxa play a hugely
#' important role in explaining the subtleties of the ecosystem, however, a
#' regression analysis on the observed diversity level fails to account for
#' these unobserved taxa. By predicting the total level of diversity (for
#' example, via \code{\link{breakaway}}) and estimating the standard error in
#' the estimate, one can take account of these unobserved, but important, taxa.
#' In order to account for the estimated nature of the response, a mixed model
#' approach is taken, whereby the varying levels of confidence in the estimates
#' contributes to a diagonal but heteroscedastic covariance matrix. Given
#' covariates constitute the fixed effects in the mixed model, and significance
#' of the random effect term ``sigsq_u'' reflects heterogeneity in the sample,
#' that is, variability that cannot be explained by only the covariates. The
#' authors believe this to be the first attempt at modelling total diversity in
#' a way that accounts for its estimated nature.
#' @author Amy Willis
#' 
#' @importFrom stats coef dexp dgeom dnbinom dpois fitted lm model.matrix nls optim pchisq pnorm predict quantile rbeta rbinom rnbinom rnorm runif sd var vcov
#' 
#' @seealso \code{\link{breakaway}}; \code{\link{breakaway_nof1}};
#' \code{\link{apples}}
#' @references Willis, A., Bunge, J., and Whitman, T. (2015). Inference for
#' changes in biodiversity. \emph{arXiv preprint.}
#' 
#' Willis, A. and Bunge, J. (2015). Estimating diversity via frequency ratios.
#' \emph{Biometrics.}
#' @keywords diversity
#' @examples
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' df <- data.frame(chats = c(2000, 3000, 4000, 3000), ses = c(100, 200, 150, 180),
#'                  Cont_var = c(100, 150, 100, 50))
#' 
#' # formula notation 
#' betta(formula = chats ~ Cont_var, ses = ses, data = df)
#' 
#' # direct input 
#' betta(c(2000, 3000, 4000, 3000), c(100, 200, 150, 180), cbind(1, c(100, 150, 100, 
#'     50)))
#' 
#' ## handles missing data
#' betta(c(2000, 3000, 4000, 3000), c(100, 200, 150, NA))
#' 
#' ## A test for heterogeneity of apples diversity estimates vs butterfly estimates
#' betta(c(1552, 1500, 884), c(305, 675, 205), cbind(1, c(0, 0, 1)))
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' @export betta
betta <- function(chats = NULL, ses, X = NULL, 
                  initial_est = NULL, formula = NULL, data = NULL) {
  if (!is.null(formula)) {
    if (is.null(data)) {
      stop("Please include a dataframe that corresponds with your formula.")
    }
  } else {
    if (is.null(chats)) {
      stop("Please include 'ses' along with either the argument 'chats' 
           or the arguments 'formula' and 'data'.")
    }
  }
  if (!is.null(formula)) {
    ses <- data[,deparse(substitute(ses))]
    X <- stats::model.matrix(formula, data)
    chats <- stats::model.response(stats::model.frame(formula = formula, data = data))
  }
  if (isTRUE(is.null(X))) { X <- matrix(rep(1,length(chats)),ncol = 1) }
  
  consider <- !(is.na(chats) | is.na(ses) | apply(is.na(X), 1, sum))
  
  chats_effective <- chats[consider]
  ses_effective <- ses[consider]
  X_effective <- as.matrix(X[consider,])
  
  n <- dim(X_effective)[1]
  p <- dim(X_effective)[2]
  
  likelihood <- function(input) {
    ssq_u <- input[1]
    beta <- input[2:length(input)]
    W <- diag(1/(ssq_u + ses_effective^2))
    -0.5*(sum(log(ssq_u + ses_effective^2) + (chats_effective - X_effective %*% beta)^2/(ssq_u + ses_effective^2))  +  log(det(t(X_effective) %*% W %*% X_effective)))
  }
  
  if (any(is.null(initial_est))) {
    initial_est <- c(var(chats_effective), solve(t(X_effective) %*% X_effective) %*% t(X_effective) %*% chats_effective)
  }
  output <- try(optim(initial_est, 
                      likelihood, 
                      hessian = FALSE, 
                      control = list(fnscale = -1), # fnscale => maximises if -ve
                      lower = c(0, rep(-Inf, p)), 
                      method = "L-BFGS-B"), 
                silent = TRUE)
  
  # while loop
  i = 0
  
  # perturb the initialisation incrementally, up to 200 times
  
  while ("try-error" %in% class(output) & i < 200) {
    i <- i + 1
    
    perturb <- rnorm(n = length(initial_est), 
                     mean = c(0,0), 
                     sd = 0.001 * i * abs(initial_est))
    initial_est_perturbed <- pmax(c(0, rep(-Inf, p)), 
                                  initial_est + perturb)
    output <- try(optim(initial_est_perturbed, 
                        likelihood, 
                        hessian = FALSE, 
                        control = list(fnscale = -1), # fnscale => maximises if -ve
                        lower = c(0, rep(-Inf, p)), 
                        method = "L-BFGS-B"), 
                  silent = TRUE)
  } 
  
  if ("try-error" %in% class(output) ) {
    stop(paste("The starting value and 200 perturbations were not", 
               "enough to find a maximum likelihood solution.", 
               "Please try again with a new choice of `initial_est`."))
  }
  
  ssq_u <- output$par[1]
  beta <- output$par[2:length(output$par)]
  
  W <- diag(1/(ssq_u + ses_effective^2))
  vars <- 1/diag(t(X_effective) %*% W %*% X_effective)
  
  global <- t(beta) %*% (t(X_effective) %*% W %*% X_effective) %*% beta ## global test
  
  Q <- sum((chats_effective - X_effective %*% beta)^2/ses_effective^2)
  R <- diag(ses_effective^2)
  G <- diag(ssq_u, n)
  
  mytable <- list()
  mytable$table <- cbind("Estimates"=beta, 
                         "Standard Errors"=sqrt(vars), 
                         "p-values"=round(2*(1-pnorm(abs(beta/sqrt(vars)))), 3))
  mytable$cov <- solve(t(X_effective) %*% W %*% X_effective)
  mytable$ssq_u <- ssq_u
  mytable$homogeneity <- c(Q, 1-pchisq(Q, n-p))
  mytable$global <- c(global, 1-pchisq(global, p-1))
  
  us <-  c(ssq_u*W %*% (chats_effective - X_effective %*% beta))
  blups <- rep(NA, length(chats)) 
  blups[consider] <- c(X_effective %*% beta + us)
  mytable$blups <- blups
  
  if (class(try(getvar(), silent=T)) != "try-error") {
  # get BLUPs SEs
  var_matrix <- matrix(NA, nrow=(n + p), ncol=(n + p))
  var_matrix[1:p, 1:p] <- t(X_effective) %*% solve(R) %*% X_effective
  var_matrix[(p + 1):(n + p), (p + 1):(n + p)] <- solve(R) + MASS::ginv(G)
  var_matrix[1:p, (p + 1):(n + p)] <- t(X_effective) %*% solve(R)
  var_matrix[(p + 1):(n + p), 1:p] <- solve(R) %*% X_effective
  var_matrix_inv <- MASS::ginv(var_matrix)
  blupvars <- rep(NA, length(chats))
  blupvars[consider] <- (cbind(X_effective, diag(1, n)) %*%
                           var_matrix_inv %*%
                           t(cbind(X_effective, diag(1, n)))) %>% diag %>% sqrt %>% c
  mytable$blupses <- blupvars
  }
  
  # For n-dimensional MVN random variable, likelihood is...
  logLhat <- -0.5*(n*log(2*pi) + # -n/2 * log(2 * pi) +
                     sum(log(ssq_u + ses_effective^2) + # -1/2 * log | Sigma| + 
                           (chats_effective - X_effective %*% beta)^2/(ssq_u + ses_effective^2))) # -1/2 (x - mu)^T * Sigma^-1 * (x - mu)
  
  mytable$loglikelihood <- logLhat
  # AIC = 2k - 2 log L(theta = MLE)
  # k = # fitted parameters = 1 + p (one variance term; p regression terms)
  mytable$aic <- -2 * logLhat + 2 * (1 + p)
  
  # AICc = AIC + (2k^2+2k)/(n-k-1)
  mytable$aicc <- mytable$aic + (2*(1 + p)^2 + 2*(1 + p))/(n - (1 + p) - 1)
  
  # R-squared (WLS): Eqn 7 of Willett & Singer, 1988, American Statistician.
  mytable$r_squared_wls <- 1 - sum((chats_effective - X_effective %*% beta)^2) / 
    (sum(chats_effective^2) - n*(mean(chats_effective))^2)
  
  
  return(mytable)
}

#' modelling total diversity with random effects
#' 
#' This function extends betta() to permit random effects modelling.
#' 
#' 
#' @param chats A vector of estimates of total diversity at different sampling
#' locations. Required with the \code{groups} argument, and optionally with the 
#' \code{X} argument.
#' @param ses The standard errors in \code{chats}, the diversity estimates. This
#' can either be a vector of standard errors (with the arguments \code{chats} and 
#' \code{X}), or the name of the variable in the dataframe \code{data} that contains 
#' the standard errors (with the arguments \code{formula} and \code{data}). 
#' @param X A numeric matrix of covariates corresponding to fixed effects. If
#' not supplied, an intercept-only model will be fit. Optional with the \code{chats} 
#' and \code{groups} arguments.
#' @param groups A categorical variable representing the random-effects groups
#' that each of the estimates belong to. Required with the \code{chats} argument and 
#' optionally with the \code{X} argument. 
#' @param formula A formula object of the form \eqn{y ~ x | group}. Required with 
#' the \code{data} argument. 
#' @param data A dataframe containing the response, response standard errors, covariates, 
#' and grouping variable. Required with the \code{formula} argument. 
#' @return \item{table}{ A coefficient table for the model parameters. The
#' columns give the parameter estimates, standard errors, and p-values,
#' respectively. This model is only as effective as your diversity estimation
#' procedure; for this reason please confirm that your estimates are
#' appropriate and that your model is not misspecified. betta_pic may be
#' useful for this purpose.  } \item{cov}{ Estimated covariance matrix of the
#' parameter estimates.  } \item{ssq_u}{ The estimate of the heterogeneity
#' variance.  } \item{ssq_g}{ Estimates of within-group variance. The estimate
#' will be zero for groups with only one observation.  } \item{homogeneity}{
#' The test statistic and p-value for the test of homogeneity.  }
#' \item{global}{ The test statistic and p-value for the test of model
#' explanatory power.  } \item{blups}{ The conditional expected values of the
#' diversity estimates (conditional on the random effects). Estimates of
#' variability for the random effects case are unavailable at this time; please
#' contact the maintainer if needed.  }
#' @author Amy Willis
#' 
#' @importFrom stats coef dexp dgeom dnbinom dpois fitted lm model.matrix nls optim pchisq pnorm predict quantile rbeta rbinom rnbinom rnorm runif sd var vcov
#' @import lme4 
#' 
#' @seealso \code{\link{betta}}; 
#' @references Willis, A., Bunge, J., and Whitman, T. (2015). Inference for
#' changes in biodiversity. \emph{arXiv preprint.}
#' @keywords diversity
#' @examples
#' 
#' df <- data.frame(chats = c(2000, 3000, 4000, 3000), ses = c(100, 200, 150, 180), Cont_var = c(100, 150, 100, 50), groups = c("a", "a", "a", "b"))
#' 
#' # formula notation
#' betta_random(formula = chats ~ Cont_var| groups, ses = ses,  data = df)    
#' 
#' # direct input
#' betta_random(c(2000, 3000, 4000, 3000), c(100, 200, 150, 180), X = cbind(Int = 1, Cont_var = c(100, 150, 100, 50)), groups = c("a", "a", "a", "b"))
#'     
#' ## handles missing data
#' betta_random(c(2000, 3000, 4000, 3000), c(100, 200, 150, NA), groups = c("a", NA, 
#'     "b", "b"))
#' 
#' @export 

betta_random <- function(chats = NULL, ses, X = NULL, groups = NULL, formula = NULL, data = NULL) {
  if (!is.null(formula)) {
    if (is.null(data)) {
      stop("Please include a dataframe that corresponds with your formula.")
    }
  } else {
    if (is.null(chats)) {
      stop("Please include 'ses' along with either the arguments 'chats' and 'groups' 
           or the arguments 'formula' and 'data'.")
    }
    if (is.null(groups)) {
      stop("Please include a vector of group memberships as 'groups'.")
    }
  }
  if (!is.null(formula)) {
    ses <- data[,deparse(substitute(ses))]
    group_var <- lme4:::barnames(lme4::findbars(lme4:::RHSForm(formula)))
    groups <- data[,group_var]
    full_form <- lme4::subbars(formula)
    sm_form <- update(full_form, paste("~.-",group_var))
    X <- stats::model.matrix(sm_form, data)
    chats <- stats::model.response(stats::model.frame(formula = sm_form, data = data))
  }
  if (isTRUE(is.null(X))) { X <- matrix(rep(1,length(chats)),ncol=1) }
  consider <- !(is.na(chats) | is.na(ses) | is.na(groups) | apply(is.na(X),1,sum))
  
  chats_effective <- chats[consider]; ses_effective <- ses[consider]; X_effective <- as.matrix(X[consider,])
  groups_effective <- groups[consider]
  n <- dim(X_effective)[1]; p <- dim(X_effective)[2]; gs <- length(unique(groups_effective))
  
  likelihood <- function(input) {
    ssq_u <- input[1]
    beta <- input[2:(p+1)]
    ssq_group <- input[(p+2):(p+gs+1)]
    group_variance <- ssq_group[groups_effective]
    W <- diag(1/(ssq_u+ses_effective^2+group_variance))
    -0.5*(sum(log(ssq_u+ses_effective^2+group_variance)+(chats_effective - X_effective %*% beta)^2/(ssq_u+ses_effective^2+group_variance)) + log(det(t(X_effective) %*% W %*% X_effective)))
  }
  
  within_groups_start <- c(by(chats_effective, groups_effective, var, simplify = T))
  within_groups_start[which(is.na(within_groups_start))]  <- 0
  
  initial_est <- c(var(chats_effective),
                   solve(t(X_effective) %*% X_effective) %*% t(X_effective) %*% chats_effective,
                   within_groups_start)
  
  output <- optim(initial_est, 
                  likelihood, 
                  hessian=FALSE, 
                  control=list(fnscale=-1),
                  lower=c(0, rep(-Inf, p), rep(0, gs)), 
                  method="L-BFGS-B")
  ssq_u <- output$par[1]
  beta <- output$par[2:(p+1)]
  ssq_group <- output$par[(p+2):(p+gs+1)]
  
  W <- diag(1/(ssq_u+ses_effective^2+ssq_group[groups_effective]))
  vars <- 1/diag(t(X_effective) %*% W %*% X_effective)
  
  global <- t(beta) %*% (t(X_effective) %*% W %*% X_effective) %*% beta ## global test
  
  Q <- sum((chats_effective - X_effective %*% beta)^2/ses_effective^2)
  
  mytable <- list()
  mytable$table <- cbind("Estimates"=beta,"Standard Errors"=sqrt(vars),"p-values"=round(2*(1-pnorm(abs(beta/sqrt(vars)))),3))
  try(rownames(mytable$table) <- colnames(X), silent = T)
  mytable$cov <- solve(t(X_effective) %*% W %*% X_effective)
  mytable$ssq_u <- ssq_u
  mytable$ssq_group <- ssq_group
  mytable$homogeneity <- c(Q,1-pchisq(Q,n-p))
  mytable$global <- c(global,1-pchisq(global,p-1))
  
  us <-  c(ssq_u*W %*% (chats_effective - X_effective %*% beta))
  blups <- rep(NA, length(chats))
  blups[consider] <- c(X_effective %*% beta + us)
  mytable$blups <- blups
  
  return(mytable)
}