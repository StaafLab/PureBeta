#' beta_correction_for_cohorts
#'
#' This functions applies the original Staaf & Aine DNA methylation beta value
#' correction strategy to adjust the effect of sample purities in cohorts with
#' known sample purities and CpG beta values (Staaf & Aine, Plos One, 2022).
#' This function allows mult-core execution.
#'
#' @param beta_values A matrix with CpGs as rows and tumour samples as columns
#' with the uncorrected beta values from the CpGs of the samples that are
#' intended to be corrected. The values must be numeric, the rows must be named
#' with the CpG ID, and the columns with the sample IDs. An example of the
#' required format is available in the example_betas_reference matrix.
#'
#' @param tumour_purities Named vector containing the sample purity values of
#' of the samples whose DNA methylation beta values are intended to be corrected.
#' The vector must be named with the sample ID, which must match with the sample
#' IDs from the matrix containing the beta values. An example of the required
#' format is available in the example_purities_reference vector.
#'
#' @param set_seed Default = FALSE. A seed for the FlexMix package to detect
#' the different CpG methylation patterns can be used by setting this argument
#' to TRUE.
#'
#' @param seed_num Default = 2000. The seed to be used when set_seed = TRUE can
#' be specified here.
#'
#' @param cores Default = 1. Number of cores to be used to run the function in
#' parallel.
#'
#' @return List containing the original uncorrected beta values
#' (output$betas.original), the corrected tumour beta values (output$betas.tumour)
#' and the corrected microenvironment beta values (output$betas.microenvironment).
#'
#' @export
#'
#' @examples
#'
#' # Using the default parameters
#' beta_correction_for_cohorts (beta_values = example_betas_reference,
#'                              tumour_purities = example_purities_reference)
#'
#' # Specifying new parameters
#' beta_correction_for_cohorts (beta_values = example_betas_reference,
#'                              tumour_purities = example_purities_reference,
#'                              set_seed = TRUE,
#'                              seed_num = 1,
#'                              cores = 5)
beta_correction_for_cohorts <- function(

  beta_values, # Matrix of beta values of CoGs of certain samples
  tumour_purities, # Purity values of the samples to be analysed
  set_seed = FALSE, # Boolean to specify if a seed should be set
  seed_num = 2000, # Seed to be used if set_seed == TRUE
  cores = 1 # Cores to be used to run the function

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

  # Ensuring the cluster is stopped properly
  on.exit(stopCluster(cl), add = TRUE)


  # ======================
  # CORRECTING BETA VALUES
  # ======================

  cat("\nCorrecting betas based on a cohort\n\n")

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
  #identified methylation population generates a regerssion with perfect fit.
  #This does not happen due to errors
  res <- suppressWarnings(pbapply(cl = cl, #ClusterS to run the process
                 MARGIN = 1, #Apply the function to the rows
                 FUN = adjustBeta, #Function to correct betas
                 purity=tumour_purities, #Purity values
                 snames=betaNames, #Sample names
                 seed=set_seed, #Specify if the seed has been added to the data or not
                 X = betaRun)) #Beta values+the added seed


  # =================
  # GENERATING OUTPUT
  # =================

  cat("\nGenerating output...\n\n")

  # Generating output list. Short and extended versions

  result_list <- list(
      betas.original = do.call("rbind",lapply(res,function(x) x$y.orig)), #Original beta values
      betas.tumour = do.call("rbind",lapply(res,function(x) x$y.tum)), #Corrected tumor beta values
      betas.microenvironment = do.call("rbind",lapply(res,function(x) x$y.norm)) #Corrected microenvironment beta values
    )


  cat("\n\n=================\n")
  cat ("PROCESS FINISHED")
  cat("\n=================\n\n")

  # Returning output
  return(result_list)

}
