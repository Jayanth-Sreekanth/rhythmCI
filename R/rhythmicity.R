#' Rhythmicity analysis using Harmonic Regression and bootstrapping
#'
#' @description
#' `rhythmicity` applies harmonic regression to the variables and calculates confidence intervals
#' using bootstrapping. The intervals are corrected using bonferroni method in the end.
#'
#' @details
#' A matrix with samples across time in rows and with columns of variables are to be provided to perform harmonic regression.
#' After regressing the data with cosine functions, the best fit and the residuals are obtained. The residuals
#' are then sampled with replacement using bootstrapping. In each instance of bootstrapping, the sampled residuals
#' are summed to the best fit yielding a bootstrap estimate of amplitude. The Confidence Interval (CI) for amplitude
#' estimates will be computed using the provided alpha value (also includes Bonferroni correction). Using the CI and fc_cutoff
#' the rhythmic status of each variables are defined as follows: rhythmic (if the lower and upper CI lies above the cutoff),
#' arrhythmic (if the lower and upper CI lies below or is equal to the cutoff) and inconclusive (if the lower and upper CI straddles
#' the cutoff). The obtained mesor (the average value), amplitude, and phase is estimated before bootstrapping.
#'
#' Please note that at least 4 samples are needed per subject. Only the variables with at least 4 samples will be analyzed.
#' The column-wise operation will need more time than matrix-wise operation.
#'
#' @param input Matrix, where columns correspond to variables and rows correspond to samples w.r.t time.
#' @param operation Anyone of the following should be provided: "matrix" or "column". Do not use "matrix" mode when there are missing values
#'   in the matrix, instead use "column" mode.
#' @param seed Is used to reproduce the same results across multiple runs.
#' @param tau The period of oscillation.
#' @param num_boot The number of bootstrap samples to be generated.
#' @param alpha The significance level to estimate the confidence level. Default is 0.05.
#' @param fc_cutoff The amplitude cutoff to define the rhythmic status of variables (genes and other).
#' @return A list for each subject containing the following: a lower (ci.lower) and upper (ci.upper) value of the confidence
#'   interval of the amplitude, the rhythmic status of each variable (rhythmicity), mesor from harmonic regression (hr.mesor),
#'   amplitude from harmonic regression (hr.amplitude), and phase from harmonic regression (hr.phase).
#' @importFrom stringr str_split
#' @importFrom stats lm na.omit fitted residuals coef
#' @importFrom boot boot boot.ci
#' @export
#' @author Jayanth Sreekanth
#' @examples
#'
#' # Load the data that comes along with the package and run rhythmicity.
#' # This data was generated using simphony.
#'
#' data(simdf)
#'
#' suppressWarnings(
#'
#' output_hr <- rhythmicity(input = simdf,
#'                          operation = "matrix",
#'                          seed = 7,
#'                          tau = 24, #value in general for humans
#'                          num_boot = 999,
#'                          alpha = 0.05,
#'                          fc_cutoff = 0.5)
#'
#' )


rhythmicity <- function(input,
                        operation,
                        seed=NA,
                        tau=24,
                        num_boot,
                        alpha=0.05,
                        fc_cutoff)
{
  #creating a list of ids in the df - contains duplicates
  ids <- str_split(rownames(input), pattern = '\\_', simplify = TRUE)[,1]
  #Creating a list of data frames splitted according to ids
  data_log_list <- split.data.frame(input, f=ids)
  #total number of patients
  ttl_sub <- length(data_log_list)

  if (operation == 'matrix') {

    if(any(is.na(input))) {
      stop('Found missing value(s). Cannot peform the analysis in \"matrix\" mode. \n Try \"column\" mode.')
    }

    if (!is.na(seed)) {
      set.seed(seed)
    }

    ot <- lapply(data_log_list, function(x) rhy_func_mat(x, ttl_sub, tau, alpha, fc_cutoff, num_boot))

    return(ot)

  } else if (operation == 'column') {

    if (!is.na(seed)) {
      set.seed(seed)
    }

    ot <- lapply(data_log_list, function(x) rhy_func(x, ttl_sub, tau, alpha, fc_cutoff, num_boot))

    return(ot)

  } else{
    stop(paste("Please specify one of the following operation mode: \"matrix\" or \"column\""))
  }


}


#' Matrix-wise operation - main function
#'
#' @description
#' `rhy_func_mat` performs harmonic regression and returns the CI of amplitude, rhythmic status and mesor.
#' Note that the mesor is estimated before bootstrapping
#'
#' @param df Matrix of value, columns corresponding to variables and rows corresponding to samples w.r.t time
#' @param ttl_sub Total number of subjects
#' @param tau The period of one oscillation
#' @param alpha The significance level to estimate the confidence level. Default is 0.05
#' @param fc_cutoff The amplitude cutoff to define the rhythmic status of variables (genes and other)
#' @param num_boot The number of bootstrap samples to be generated
#' @return A list with the CI of amplitude, rhythmic status and mesor for each variable
#' separated according to subjects
#'
#' @author Jayanth Sreekanth

rhy_func_mat <- function(df, ttl_sub, tau, alpha, fc_cutoff, num_boot) {
  #retrieving time stamps and id from rowname
  tms_df <- as.numeric(str_split(rownames(df), pattern = '\\_', simplify = TRUE)[, 2])
  id <- unique(str_split(rownames(df), pattern = '\\_', simplify = TRUE)[, 1])

  #Excluding genes/columns having less than four time-points
  df_col_exc <- names(which(colSums(!is.na(df)) <= 3))
  #Including genes/columns having more than three time-points
  df_col_inc <- names(which(colSums(!is.na(df)) > 3))

  if (length(df_col_exc) != 0) {
    cat("For ", id, ", the following parameters have fewer than four time points:\n\"",
        paste(df_col_exc, collapse = "\", \""), "\"\n\n", sep = "")
  }

  if (length(df_col_inc) == 0) {
    stop("Can't perform the analysis because none of the parameters have atleast four time points.")
  }

  #Create a df by selecting only the columns having more than three time-points
  df_sltd <- df[, colnames(df) %in% df_col_inc]

  #Harmonic Regression
  input_ts_fit <- lm(df_sltd ~ cos((2*pi*tms_df)/tau) +
                       sin((2*pi*tms_df)/tau), x=TRUE, na.action = na.omit)

  # Checking the determinant of the model matrix
  mdl_det <- zapsmall(input_ts_fit$x)
  if (det(t(mdl_det) %*% mdl_det) == 0) {
    stop("Problem with regression. Basis functions are linearly dependent. \nCheck the input.")
  }

  #Obtaining the fit, residuals, and mesor
  fit_vals <- as.matrix(fitted(input_ts_fit, na.omit = T))
  resi <- residuals(input_ts_fit)
  avg_exp <- coef(input_ts_fit)[1,]
  ampli <- sqrt(coef(input_ts_fit)[2, ]^2 + coef(input_ts_fit)[3, ]^2)
  phase <- atan2(coef(input_ts_fit)[3, ], coef(input_ts_fit)[2, ]) %% (2*pi) #normalize the phase to have it positive

  #Function for bootstrapping
  boot_final <- boot(data = resi, statistic = func_boot_mat, R=num_boot, lm_fitted = fit_vals, tau=tau)

  # Creating a dataframe to store the confidence intervals for patients after bonferroni correction
  boot_ci_bonf <- data.frame(matrix(ncol = 6, nrow = ncol(boot_final$t)))
  colnames(boot_ci_bonf) <- c('ci.lower', 'ci.upper', 'rhythmicity', 'hr.mesor', 'hr.amplitude', 'hr.phase')

  #the loop runs for the number of components present
  for (c in 1:ncol(boot_final$t)) {
    boot_dyna_ci_bonf <- boot.ci(boot_final, conf = (1 - (alpha / (ncol(boot_final$t) * ttl_sub))),
                                 type = 'bca',
                                 index = c)

    rownames(boot_ci_bonf)[c] <- colnames(boot_final$t0)[c]
    boot_ci_bonf[c, 'ci.lower'] <- round(boot_dyna_ci_bonf$bca[4], digits = 3) #assigning values to a cell
    boot_ci_bonf[c, 'ci.upper'] <- round(boot_dyna_ci_bonf$bca[5], digits = 3) #assigning values to a cell
    #calculating the rhythmic feature
    if ((boot_dyna_ci_bonf$bca[4] <= fc_cutoff) & (boot_dyna_ci_bonf$bca[5] <= fc_cutoff)) {
      boot_ci_bonf[c, 'rhythmicity'] <- 'arrhythmic'
    } else if ((boot_dyna_ci_bonf$bca[4] > fc_cutoff) & (boot_dyna_ci_bonf$bca[5] > fc_cutoff)) {
      boot_ci_bonf[c, 'rhythmicity'] <- 'rhythmic'
    } else {
      boot_ci_bonf[c, 'rhythmicity'] <- 'inconclusive'
    }
  }

  boot_ci_bonf <- boot_ci_bonf[order(rownames(boot_ci_bonf)),]
  boot_ci_bonf[match(names(avg_exp), rownames(boot_ci_bonf)), 'hr.mesor'] <- round(avg_exp, digits =3)
  boot_ci_bonf[match(names(ampli), rownames(boot_ci_bonf)), 'hr.amplitude'] <- round(ampli, digits =3)
  boot_ci_bonf[match(names(phase), rownames(boot_ci_bonf)), 'hr.phase'] <- round(phase, digits =3)

  return(boot_ci_bonf)
}

#' Matrix-wise operation - sub function to estimate the amplitude during every bootstrap instance
#'
#' @description
#' `func_boot_mat` adds the sampled residuals with the obtained best fit. Then performs harmonic regression
#' and returns the amplitude for each bootstrap instance.
#'
#'
#' @param df_resi The sampled matrix of residuals in an instance
#' @param n Row number in the matrix
#' @param tau The period of one oscillation
#' @param lm_fitted The fitted value of initial harmonic regression
#' @return An amplitude estimate for the instance
#'
#' @author Jayanth Sreekanth

func_boot_mat <- function(df_resi, n, lm_fitted, tau){
  boot_df <- lm_fitted + df_resi[n, ] #Summing the fitted values with residuals

  ampli_boot_df <- matrix(nrow = 1, ncol = ncol(boot_df)) #creating an empty df with NAs
  colnames(ampli_boot_df) <- colnames(boot_df) #defining proper column names

  # Checking in total if there are more than 3 time-points available
  # because it is tricky while bootstrapping with NA
  if (all(colSums(!is.na(boot_df)) <= 3)) {
    return(ampli_boot_df) #if not then return the empty df for ampli
  } else {

    #Excluding genes/columns having less than four time-points
    boot_col_exc <- names(which(colSums(!is.na(boot_df)) <= 3))
    #Including genes having more than three time-points
    boot_col_inc <- names(which(colSums(!is.na(boot_df)) > 3))
    #Create a df by selecting only the columns having more than three time-points
    boot_df_sltd <- boot_df[, colnames(boot_df) %in% boot_col_inc]
    #retrieving time stamps from rowname
    tms <- as.numeric(str_split(rownames(boot_df_sltd), pattern = '\\_', simplify = TRUE)[, 2]) #determining the time from row names
    #HR
    boot_ts_fit <- lm(boot_df_sltd ~ cos((2*pi*tms)/tau) +
                        sin((2*pi*tms)/tau), x=TRUE, na.action = na.omit)
    #retrieving the coefficients for all the components
    boot_coeff <- t(coef(boot_ts_fit))
    #calculating the amplitude
    boot_ampli <- sqrt(boot_coeff[, 2]^2 + boot_coeff[, 3]^2)

    #Only append the values of the columns that were analysed
    ampli_boot_df[match(names(boot_ampli), colnames(ampli_boot_df))] <- boot_ampli

    return(ampli_boot_df)

  }
}

#' Column-wise operation - main function
#'
#' @description
#' `rhy_func` performs harmonic regression for each variable separately and returns the CI of amplitude,
#' rhythmic status and mesor. Note that the mesor is estimated before bootstrapping
#'
#'
#' @param df Matrix of value, columns corresponding to variables and rows corresponding to samples w.r.t time
#' @param ttl_sub Total number of subjects
#' @param tau The period of one oscillation
#' @param alpha The significance level to estimate the confidence level. Default is 0.05
#' @param fc_cutoff The amplitude cutoff to define the rhythmic status of variables (genes and other)
#' @param num_boot The number of bootstrap samples to be generated
#' @return A list with the CI of amplitude, rhythmic status and mesor for each variable
#' separated according to subjects
#'
#' @author Jayanth Sreekanth

rhy_func <- function(df, ttl_sub, tau, alpha, fc_cutoff, num_boot) {
  #retrieving time stamps from row name
  tms_df <- as.numeric(str_split(rownames(df), pattern = '\\_', simplify = TRUE)[, 2])

  #Excluding genes/columns having less than four time-points
  df_col_exc <- names(which(colSums(!is.na(df)) <= 3))
  #Including genes having more than three time-points
  df_col_inc <- names(which(colSums(!is.na(df)) > 3))

  #Create a df by selecting only the columns having more than three time-points
  df_sltd <- df[, colnames(df) %in% df_col_inc]

  #calling another following function
  boot_ci_bonf  <- sapply(colnames(df_sltd), fit_fun, df_dyn = df_sltd, n_col = ncol(df_sltd),
                          ttl_sub = ttl_sub, tau = tau, alpha = alpha, fc_cutoff = fc_cutoff, num_boot = num_boot)
  #converting a list to df
  boot_ci_bonf <- data.frame(t(boot_ci_bonf))
  #ordering the rows
  boot_ci_bonf <- boot_ci_bonf[order(rownames(boot_ci_bonf)),]

  return(boot_ci_bonf)
}

#' Column-wise operation - sub function
#'
#' @description
#' `fit_fun` performs harmonic regression, calls bootstrapping, calculates the amplitude CI
#'
#' @param nm Column number
#' @param df_dyn The matrix with columns having more than three time-points
#' @param n_col Number of columns in the matrix
#' @param ttl_sub Total number of subjects
#' @param tau The period of one oscillation
#' @param alpha The significance level to estimate the confidence level. Default is 0.05
#' @param fc_cutoff The amplitude cutoff to define the rhythmic status of variables (genes and other)
#' @param num_boot The number of bootstrap samples to be generated
#' @return A list with the CI of amplitude, rhythmic status and mesor for each variable separated according
#' to subjects
#'
#' @author Jayanth Sreekanth

fit_fun <- function(nm, df_dyn, n_col, ttl_sub, tau, alpha, fc_cutoff, num_boot) {
  #retrieving time stamps from row name
  tms_df <- as.numeric(str_split(rownames(df_dyn), pattern = '\\_', simplify = TRUE)[, 2])
  #HR
  input_ts_fit <- lm(df_dyn[,nm] ~ cos((2*pi*tms_df)/tau) +
                       sin((2*pi*tms_df)/tau), x=TRUE, na.action = na.omit)
  #Retrieving the fit, residuals and average expression
  fit_val <- as.matrix(fitted(input_ts_fit, na.omit = T))
  colnames(fit_val) <- nm #assigning the column name
  res <- as.matrix(residuals(input_ts_fit))
  colnames(res) <- nm #assigning the column name
  avg_exp <- round(coef(input_ts_fit)[[1]], digits = 3)
  ampli <- round(sqrt(coef(input_ts_fit)[[2]]^2 + coef(input_ts_fit)[[3]]^2), digits = 3)
  phase <- round((atan2(coef(input_ts_fit)[[3]], coef(input_ts_fit)[[2]]) %% (2*pi)), digits = 3) #normalize the phase to have it positive


  #Bootstrapping
  boot_final <- boot(data = res, statistic = func_boot, R=num_boot, lm_fitted = fit_val, tau=tau)
  #Calculating the bootstrap confidence intervals using Bias corrected and acclerated method
  boot_ci_bonf <- boot.ci(boot_final, conf = (1 - (alpha / (n_col * ttl_sub))), type = 'bca')
  lower_ci <- round(boot_ci_bonf$bca[4], digits = 3) #assigning values to the variable
  upper_ci <- round(boot_ci_bonf$bca[5], digits = 3) #assigning values to the variable
  rhythmicity <- ifelse(((lower_ci <= fc_cutoff) & (upper_ci <= fc_cutoff)), 'arrhythmic', ifelse(
    ((lower_ci > fc_cutoff) & (lower_ci > fc_cutoff)), 'rhythmic', 'inconclusive'))

  return(c('ci.lower' = as.double(lower_ci), 'ci.upper' = as.double(upper_ci),
           'rhythmicity' = rhythmicity,'hr.mesor' = as.double(avg_exp),
           'hr.amplitude' = as.double(ampli), 'hr.phase' = as.double(phase)))
}


#' Column-wise operation to estimate the amplitude during every bootstrap instance
#'
#' @description
#' `func_boot` adds the sampled residuals with the obtained best fit. Then performs harmonic regression
#' and returns the amplitude for each bootstrap instance.
#'
#'
#' @param df_resi The sampled matrix of residuals in an instance
#' @param n Row number in the matrix
#' @param tau The period of one oscillation
#' @param lm_fitted The fitted value of initial harmonic regression
#' @return An amplitude estimate for the instance
#'
#' @author Jayanth Sreekanth

func_boot <- function(df_resi, n, lm_fitted, tau){
  boot_df <- lm_fitted + df_resi[n, ] #Summing the fitted values with residuals
  #retrieving time stamps from row name
  tms <- as.numeric(str_split(rownames(boot_df), pattern = '\\_', simplify = TRUE)[, 2]) #determining the time from row names
  #HR
  boot_ts_fit <- lm(boot_df ~ cos((2*pi*tms)/tau) +
                      sin((2*pi*tms)/tau), x=TRUE)
  #Obtaining the coefficients of the fit
  boot_coeff <- t(coef(boot_ts_fit))
  #Calculating the amplitude using coefficients
  ampli_boot <- sqrt(boot_coeff[, 2]^2 + boot_coeff[, 3]^2)
  names(ampli_boot) <- colnames(boot_df) #properly naming the variable

  return(ampli_boot)
}

