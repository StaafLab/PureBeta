#' Reference regression generator
#'
#' This function generates the reference regressions' parameters and variance of
#' the reference CpGs required in subsequent steps of the PureBeta
#' workflow based on the Staaf & Aine methylation beta value correction strategy.
#' This function allows multi-core execution.
#'
#' @param beta_values A matrix with CpGs as rows and analysed samples as columns
#' with the uncorrected beta values from the CpGs of the samples that are
#' intended to be used to build the reference regressions. The values must be
#' numeric, the rows must be names with the CpG ID, and the columns with the
#' sample IDs. An example of the required format is available in the
#' example_betas_reference matrix.
#'
#' @param tumour_purities Named vector containing the purity values of
#' of the samples whose DNA methylation beta values are intended to be used to
#' build the reference regressions. The vector must be named with the sample ID,
#' which must match with the sample IDs from the matrix containing the beta
#' values. An example of the required format is available in the
#' example_purities_reference vector
#'
#' @param set_seed  Default = FALSE. A seed for the FlexMix package to detect
#' the different CpG methylation patterns can be used by setting this argument
#' to TRUE.
#'
#' @param seed_num Default = 2000. The seed to be used when set_seed = TRUE can
#' be specified here.
#'
#' @param cores  Default = 1. Number of cores to be used to run the function in
#' parallel.
#'
#' @param regression_details Default = FALSE. The user can set this argument to
#' TRUE to obtain an extended version of the output reference list containing
#' the original beta values, the corrected tumour and microenvironment betas
#' after applying the Staaf & Aine beta correction for whole cohorts in the
#' reference data, and the methylation pattern or population assigned to each
#' sample used to buld the reference regressions per CpG.
#'
#' @returns A list with the parameters of the computed reference regressions. The
#' list contains the variance of the betas of CpGs used to build the regressions
#' (output$cpg.variance), the slopes, intercepts residual standard error and
#' degrees of freedom of the regression calculated per CpG (output$reg.slopes,
#' output$reg.intercepts, output$reg.RSE and output$df). If the extended output
#' list is required it will also include the original uncorrected beta values
#' (output$betas.original), the corrected tumour beta values (output$betas.tumour),
#' the corrected microenvironment beta values (output$betas.microenvironment) and
#' the methylation pattern or population assigned to each sample used to buld
#' the reference regressions per CpG (output$cpg.populations).
#'
#' @export
#'
#' @examples
#'
#' # Using the default parameters
#' reference_regression_generator(beta_values = example_betas_reference,
#'                              tumour_purities = example_purities_reference)
#'
#' # Specifying new parameters
#' reference_regression_generator(beta_values = example_betas_reference,
#'                                tumour_purities = example_purities_reference,
#'                                set_seed = TRUE,
#'                                seed_num = 1,
#'                                cores = 5,
#'                                regression_details = TRUE)
#'
reference_regression_generator <- function(

  beta_values, # Matrix of beta values of CpGs of certain samples
  tumour_purities, # Purity values of the samples to be analysed
  set_seed = FALSE, # Boolean to specify if a seed should be set
  seed_num = 2000, # Seed to be used if set_seed == TRUE
  cores = 1, # Cores to be used to run the function
  regression_details = FALSE # Include the corrected betas based on the Staaf & Aine method and the populations identified for each sample in each CpG in the output

) {

  # ==============
  # ERROR HANDLING
  # ==============

  # Checking if the the samples of betas and purities match

  if (
    !identical(sort(colnames(beta_values)), sort(names(tumour_purities)))
  ) {
    stop("The samples of the beta matrix and purity vector do not match. Review the input data.")
  }


  # ===========================
  # CONFIGURING PARALLELIZATION
  # ===========================

  cat("\nUsing", cores,"core(s)\n\n")

  #Creating the cluster to run the process in parallel
  cl <- makeCluster(cores)
  registerDoParallel(cl)

  #Making sure that all packages have access to the flexmix package. Using invisible()
  #to avoid printing anything to the terminal.
  invisible(clusterEvalQ(cl, {library("flexmix")}))

  # Export all the functions in the package to the defined cores
  parallel::clusterExport(cl = cl,
                  varlist = unclass(lsf.str(envir = asNamespace("PureBeta"), all = TRUE)),
                  envir = as.environment(asNamespace("PureBeta")))


  # ================================
  # DETERMINING VARIANCE OF EACH CPG
  # ================================

  #Create a vector with the variance of each cpg (row)
  cpg_variance <- apply(beta_values, 1, var)
  names(cpg_variance) <- rownames(beta_values)


  # =======================
  # CALCULATING REGRESSIONS
  # =======================

  cat("\nCalculating reference regressions\n\n")

  # Setting seed if necessary
  if (set_seed == TRUE) {

    #Adding seed to each row of the beta value dataframe
    betaRun <- cbind(seed=seed_num:seed_num+nrow(beta_values), beta_values)

  } else {

    betaRun <- beta_values

  }

  #Storing sample names
  betaNames <- colnames(beta_values)

  # Initializing progress bar and specifying options
  pbo <- pboptions(type = "txt", char="=", txt.width=80)

  #Running the analysis in parallel with a progress bar
  #Using supress warnings to avoid the function to raise a warning when an
  #identified methylation population generates a regression with perfect fit.
  #This does not happen due to errors
  res <- suppressWarnings(pbapply(cl = cl, #Clusters to run the process
                 MARGIN = 1, #Apply the function to the rows
                 FUN = adjustBeta, #Function to correct betas
                 purity=tumour_purities, #Purity values
                 snames=betaNames, #Sample names
                 seed=set_seed, #Specify if the seed has been added to the data or not
                 betaRun)) #Beta values+the added seed

  # Ensuring the cluster is stopped properly
  stopCluster(cl)

  # =================
  # GENERATING OUTPUT
  # =================

  cat("\nGenerating output...\n\n")


  # Generating output list. Short and extended versions
  if (regression_details == TRUE) {

    # Creating a list to add the results.
    result_list <- list(
      betas.original = do.call("rbind",lapply(res,function(x) x$y.orig)), #Original beta values
      betas.tumor = do.call("rbind",lapply(res,function(x) x$y.tum)), #Corrected tumor beta values
      betas.microenvironment = do.call("rbind",lapply(res,function(x) x$y.norm)), #Corrected microenvironment beta values
      cpg.populations =  do.call("rbind",lapply(res,function(x) x$groups)), #Methylation patterns (populations) of each CpG
      cpg.variance = cpg_variance, # Variance value of the beta values of each CpG
      reg.slopes = do.call("rbind",lapply(res,function(x) x$res.slopes)), #Slopes of the populations
      reg.intercepts = do.call("rbind",lapply(res,function(x) x$res.int)), #Intercepts of the populations
      reg.RSE = do.call("rbind",lapply(res,function(x) x$res.rse)), #Residual standard error
      reg.df = do.call("rbind",lapply(res,function(x) x$res.df)) #Degrees of freedom of the reversed regressions
    )

  } else {

    # Creating a list to add the results.
    result_list <- list(
      cpg.variance = cpg_variance, # Variance value of the beta values of each CpG
      reg.slopes = do.call("rbind",lapply(res,function(x) x$res.slopes)), #Slopes of the populations
      reg.intercepts = do.call("rbind",lapply(res,function(x) x$res.int)), #Intercepts of the populations
      reg.RSE = do.call("rbind",lapply(res,function(x) x$res.rse)), #Residual standard error
      reg.df = do.call("rbind",lapply(res,function(x) x$res.df)) #Degrees of freedom of the reversed regressions
    )

  }


  cat("\n\n=================\n")
  cat ("PROCESS FINISHED")
  cat("\n=================\n\n")


  # Returning output
  return(result_list)

}
