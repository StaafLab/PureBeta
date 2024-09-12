#' Purity Estimation
#'
#' The purity_estimation() function estimates tumour sample purity based on DNA methylation
#' beta values and reference linear models that reflect the correlation between
#' tumour purity and beta values in each analysed CpG. This function can be
#' applied to both single and multiple samples in each run and allows multi-core
#' execution.
#'
#' @param reference_regressions List object containing the parameters of the
#' reference regressions determined using the reference_regression_generator()
#' function (both short and extended versions are valid). The input must
#' at least include the list containing a named vector with the variance of
#' the betas of CpGs used to build the regressions (input$cpg.variance), the slopes,
#' intercepts residual standard error and degrees of freedom of the regression
#' calculated per CpG (input$reg.slopes, input$reg.intercepts, input$reg.RSE
#' and input$df as matrices). If the prediction intervals for each regression
#' are intended to be calculated through bootstrapping instead of assuming t
#' distrbution the extended version of the regressions is required.
#'
#' @param beta_values A matrix with CpGs as rows and analysed samples (or an
#' individual sample) as columns with the uncorrected beta values from the CpGs
#' of the samples whose purities are intended to be estimated. The values must
#' be numeric, the rows must be names with the CpG ID, and the columns with the
#' sample IDs. An example of the required format is available in the
#' example_betas_to_correct matrix.
#'
#' @param alpha Default = 0.7. The alpha value used to determine the width of
#' the prediction intervals for each CpG and sample determined from the
#' reference regressions. The usage of a excessively wide or narrow interval can
#' increase the prediction error.
#'
#' @param slope_threshold Default = 0.2. Minimum slope allowed per regression
#' for them to be taken into consideration in the sample purity estimation. he
#' inclusion of regressions below the threshold can increase the prediction error,
#' and therefore will be considered uninformative an ignored.
#'
#' @param variance_threshold Default = 0.05. CpG beta value variance cutoff to
#' filter reference data. If the CpGs are not variable enough PureBeta can not
#' obtain information for the estimation. These CpGs are non-informative and can
#' increase the execution time and the prediction error.
#'
#' @param proportion_to_interval Default = 0.96. Percentage of the maximum 1 -
#' purity coverage detected to include in the estimated the 1-Purity interval.
#' A lower value will generate wider intervals, and a higher one narrower ranges.
#'
#' @param cores Default = 1. Number of cores to be used to run the function in
#' parallel.
#'
#' @param extended_output Default = FALSE. Set this argument to TRUE if the user
#' wants to obtain the CpGs used (informative for the estimation) for the purity
#' estimation of each analysed sample in the output list.
#'
#' @param assume_t_distribution Default = TRUE. Set this argument to FALSE if the user
#' does not want to use the t statistic to calculate prediction intervals from each of the
#' refernce regression. In this case the interval will be determined through bootstrapping,
#' a non-parametric strategy that may signifcantly increase the execution time.
#'
#' @param Boost_N The number of times the that the values will be bootstrapped to generate
#' the prediction interval should be entered here if assume_t_distribution = FALSE. While an
#' excessively low number will generate an unreliable output, choosing an excesively high value
#' will significantly decrease the function's execution speed.
#'
#' @returns extended_output List containing a data frame with the predicted 1 - Purity values
#' (output$`Estimated_1mPurities`). It contains the identified estimates (in very
#' exceptional cases it could be different to 1) the estimated 1-Purity values
#' and intervals predicted per each sample. If more than one estimates are
#' obtained, one independent line per estimate will be created in the data frame.
#' A list with the CpGs used for the purity prediction of each sample after the
#' filtering steps is also included in the output list (output$`Used_CpGs`) if
#' the extended_output = TRUE option is selected.
#'
#' @export
#'
#' @examples
#'
#' # Using the default parameters
#' purity_estimation(reference_regressions = output_from_reference_regression_generator,
#'                   beta_values = example_betas_to_correct)
#'
#' # Specifying new parameters
#' purity_estimation(reference_regressions = output_from_reference_regression_generator,
#'                   beta_values = example_betas_to_correct,
#'                   alpha = 0.75,
#'                   slope_threshold = 0.25,
#'                   variance_threshold = 0.06,
#'                   proportion_to_interval = 0.93,
#'                   cores = 5,
#'                   extended_output = TRUE)
#'
#' # Using bootstrapping
#' purity_estimation(reference_regressions = output_from_reference_regression_generator,
#'                   beta_values = example_betas_to_correct,
#'                   cores = 5,
#'                   extended_output = TRUE,
#'                   assume_t_distribution = FALSE,
#'                   Boots_N = 1000)
#'
purity_estimation <- function(

  reference_regressions,
  beta_values,
  alpha = 0.7,
  slope_threshold = 0.2,
  variance_threshold = 0.05,
  proportion_to_interval = 0.96,
  cores = 1,
  assume_t_distribution = TRUE,
  Boots_N = NULL,
  extended_output = FALSE
) {

  #
  # INITIALIZING VARIABLES
  #

  list_of_predicted_intervals <- list() #Create a list to append all the predicted purity intervals
  list_of_used_cpgs <- list() #Create a list to append the used CpGs for each prediction
  output_list <- list() #List to append the final output


  #
  # FILTERING REFERENCE REGRESSIONS BASED ON VARIANCE
  #

  #Generate a vector with the CpGs to filter
  cpgs_to_keep <- names(reference_regressions$cpg.variance[reference_regressions$cpg.variance >= variance_threshold])

  #Filtering regression objects
  reference_regressions$reg.slopes <- reference_regressions$reg.slopes[cpgs_to_keep,]
  reference_regressions$reg.intercepts <- reference_regressions$reg.intercepts[cpgs_to_keep,]
  reference_regressions$reg.RSE <- reference_regressions$reg.RSE[cpgs_to_keep,]
  my_df <- beta_values[rownames(beta_values) %in% cpgs_to_keep,]

  #
  # QC OF THE REGRESSIONS
  #

  # QC of reference regerssions to avoid errors due to flexmix. Check if there are NAs in wrong 
  # positions or degrees of freedom that are equal to 0
  check_na_last <- function(vec) {
  
    # Check if all the elements are not NA or NaN
    if(sum(is.na(vec)) == length(vec)) {
      return(FALSE) # Valid if all NAs are at the end
    }
    
    # Find the last non-NA element index
    last_non_na_index <- max(which(!is.na(vec)), na.rm = TRUE)
    
    # Check if there are any NAs before the last non-NA element
    if (any(is.na(vec[1:last_non_na_index]))) {
      return(FALSE) # Invalid if any NAs are found before the last non-NA
    }
    
    return(TRUE) # Valid if all NAs are at the end
  }

  # Check if th
  check_df <- function(vec) {

    # Check if ant of the elements is equal to 0
    if (0 %in% vec) {
      return(FALSE)
    } else {
      return(TRUE)
    }

  }

  # QC remove regressions with potential errors
  qc_slope_NA <- apply(reference_regressions$reg.slopes, check_na_last, MARGIN=1)
  qc_intercept_NA <- apply(reference_regressions$reg.intercepts, check_na_last, MARGIN=1)
  qc_df_NA <- apply(my_df, check_na_last, MARGIN=1)
  qc_df_0 <- apply(my_df, check_df, MARGIN=1)

  # Remove problematic CpGs from the regression list
  cpgs_to_keep <- rownames(reference_regressions$cpg.slope)[qc_slope_NA & qc_intercept_NA & qc_df_NA & qc_df_0]
  print(cpgs_to_keep)


  #Filtering regression objects
  reference_regressions$reg.slopes <- reference_regressions$reg.slopes[cpgs_to_keep,]
  reference_regressions$reg.intercepts <- reference_regressions$reg.intercepts[cpgs_to_keep,]
  reference_regressions$reg.RSE <- reference_regressions$reg.RSE[cpgs_to_keep,]
  my_df <- beta_values[rownames(beta_values) %in% cpgs_to_keep,]

  #
  # CONFIGURING PARALLELIZATION
  #

  # Printing the number of cores to be used
  cat("\nUsing", cores, "core(s)\n")

  # Creating clusters to run the script in  parallel
  cl <- makeCluster(cores)

  #Registering the clusters
  registerDoSNOW(cl)

  # Export all the functions in the package to the defined cores
  parallel::clusterExport(cl = cl,
                  varlist = unclass(lsf.str(envir = asNamespace("PureBeta"), all = TRUE)),
                  envir = as.environment(asNamespace("PureBeta")))


  #
  # RUNNING THE ANALYSIS WITH A PROGRESS BAR
  #

  # Printing command line message
  cat("\nRunning the analysis...\n\n")

if (assume_t_distribution) {

  print("T_dist")

  # Getting the names of the samples to analyse
  samples <- colnames(beta_values)

  # Defining the progress bar
  p_bar <- txtProgressBar(min = 0,
                          max = length(samples),
                          style = 3,
                          width = 80)

  # Creating a function to follow the execution of the script
  progress <- function(n) setTxtProgressBar(p_bar, n)
  opts <- list(progress = progress)

  # Running the sourced functions in parallel for each sample. The execution level will be followed through a progress bar
  out_list <- foreach(s = samples, .packages = "Kendall", .options.snow = opts) %dopar% {

    # Defining an empty matrix with the cpg ids as rownames to add the all the 1-Purity predicted intervals for all
    # the CpGs of a sample
    interval_mat <- matrix(ncol=2, nrow=length(rownames(beta_values)))
    rownames(interval_mat) <- rownames(beta_values)


    # Predicting all the 1-Purity intervals for each CpG of each sample and append them to the empty interval_mat
    for (cpg in rownames(beta_values)) {

      # The following if statement will be used to take into account only CpGs included into the
      # refernce regression dataset
      if (cpg %in% rownames(reference_regressions$reg.slopes)) {

        interval_mat[cpg,] <- predicting_purity(beta=beta_values[cpg, s],
                                                slopes=reference_regressions$reg.slopes[cpg, ],
                                                intercepts=reference_regressions$reg.intercepts[cpg, ],
                                                RSE=reference_regressions$reg.RSE[cpg, ],
                                                degrees_of_freedom=reference_regressions$reg.df[cpg, ],
                                                slope_threshold=slope_threshold,
                                                alpha=alpha,
                                                assume_t_distribution = assume_t_distribution)

      }
    }

    # Calculate the 1-Purity estimate and interval for the sample analysed.
    # The results with be shown in list named with the sample id
    list(name = s,
         value = purity_coverage(
         pred_purity_confidence=interval_mat,
         interval_threshold=100*(1-proportion_to_interval)),
         cpgs = rownames(na.omit(interval_mat))
    )
  }

  # Append the list defined for each sample to the list containing the predicted values for all the samples.
  # The sample id is used to identify each element of the list
  list_of_predicted_intervals <- setNames(lapply(out_list, function(x) x$value), sapply(out_list, function(x) x$name))
  list_of_used_cpgs <- setNames(sapply(out_list, function(x) x$cpgs), sapply(out_list, function(x) x$name))

} else {

    # Getting the names of the samples to analyse
  samples <- colnames(beta_values)

  # Defining the progress bar
  p_bar <- txtProgressBar(min = 0,
                          max = length(samples),
                          style = 3,
                          width = 80)

  # Creating a function to follow the execution of the script
  progress <- function(n) setTxtProgressBar(p_bar, n)
  opts <- list(progress = progress)

  # Running the sourced functions in parallel for each sample. The execution level will be followed through a progress bar
  out_list <- foreach(s = samples,
                      .packages = "Kendall",
                      .options.snow = opts) %dopar% {


    # Defining an empty matrix with the cpg ids as rownames to add the all the 1-Purity predicted intervals for all
    # the CpGs of a sample
    interval_mat <- matrix(ncol=2, nrow=length(rownames(beta_values)))
    rownames(interval_mat) <- rownames(beta_values)

    # Predicting all the 1-Purity intervals for each CpG of each sample and append them to the empty interval_mat
    for (cpg in rownames(beta_values)) {

      # The following if statement will be used to take into account only CpGs included into the
      # refernce regression dataset
      if (cpg %in% rownames(reference_regressions$reg.slopes)) {

        interval_mat[cpg,] <- predicting_purity(beta=beta_values[cpg, s],
                                                 slopes=reference_regressions$reg.slopes[cpg, ],
                                                 intercepts=reference_regressions$reg.intercepts[cpg, ],
                                                 RSE=reference_regressions$reg.RSE[cpg, ],
                                                 degrees_of_freedom=reference_regressions$reg.df[cpg, ],
                                                 slope_threshold=slope_threshold,
                                                 populations=reference_regressions$cpg.populations[cpg, ],
                                                 original_betas = reference_regressions$betas.original[cpg, ],
                                                 original_purities = reference_regressions$purities,
                                                 B = Boots_N,
                                                 alpha = alpha,
                                                 assume_t_distribution = assume_t_distribution
        )

      }
    }


    # Calculate the 1-Purity estimate and interval for the sample analysed.
    # The results with be shown in list named with the sample id
    list(name = s,
        cpgs = rownames(na.omit(interval_mat)),
        intervals = interval_mat,
        value = purity_coverage(pred_purity_confidence=interval_mat,
                                interval_threshold=100*(1-proportion_to_interval))
     )
  }


  # Append the list defined for each sample to the list containing the predicted values for all the samples.
  # The sample id is used to identify each element of the list
  list_of_predicted_intervals <- setNames(lapply(out_list, function(x) x$value), sapply(out_list, function(x) x$name))
  list_of_used_cpgs <- setNames(sapply(out_list, function(x) x$cpgs), sapply(out_list, function(x) x$name))
  list_of_cpg_intervals <- setNames(lapply(out_list, function(x) x$intervals), sapply(out_list, function(x) x$name))

}



  #
  # GENERATING OUTPUT
  #

  # Create a vector with the column names of the output dataframe
  cols <- c("#sample", "num_of_est", "estimate_1-purity", "low_bound", "top_bound")

  # Creating a dataframe with the columns below
  output_df <- data.frame(matrix(nrow=0, ncol=length(cols)))

  # Appending the values of the list of predicted intervals
  for (sample in names(list_of_predicted_intervals)) {

    # A different row will be appended per each detected estimate per sample
    for (num in length(list_of_predicted_intervals[[sample]][["1-Pur_estimates"]])) {

      # Creating vector with the data to append
      row <- c(sample,
               length(list_of_predicted_intervals[[sample]][["1-Pur_estimates"]]), # This will indicate the number of estimates detected
               list_of_predicted_intervals[[sample]][["1-Pur_estimates"]][num],
               list_of_predicted_intervals[[sample]][["interval(s)"]][[num]][1],
               list_of_predicted_intervals[[sample]][["interval(s)"]][[num]][2]
      )

      # Appending row to dataframe
      output_df <- rbind(output_df, row)

    }

  }

  # Stop clusters used in parallelization
  stopCluster(cl)


  #Setting column names
  colnames(output_df) <- cols


  # Adding purity estimates
  output_list$"Estimated_1-Purities" = output_df

  # Adding uised CpGs for purity estimation if the user requires so.
  if (extended_output) {
    output_list$"Used_CpGs" = list_of_used_cpgs
  }



  cat("\n\n=================\n")
  cat ("PROCESS FINISHED")
  cat("\n=================\n\n")

  # Returning output list
  return(output_list)

}
