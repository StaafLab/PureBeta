#' Reference based beta correction
#'
#' This function adjusts CpG beta values for tumor cells and inferred
#' normal cells using reference regressions and estimated purities. This can
#' be carried out refitting the regressions to include the new data points
#' (betas + estimated purities) or using the original reference regressions. Unlike
#' beta_correction_for_cohorts(), this function does not require the usage
#' of a full cohort of samples, as it is single sample and single CpG applicable.
#' This function allows multi-core execution.
#'
#' @param betas_to_correct A matrix with CpGs as rows and analysed samples (or an
#' individual sample) as columns with the uncorrected beta values from the CpGs
#' of the samples that are intended to be corrected. The values must be numeric,
#' the rows must be named with the CpG ID, and the columns with the sample IDs.
#' An example of the required format is available in the example_betas_to_correct
#' matrix.
#'
#' @param purities_samples_to_correct The output of the purity_estimation
#' function in the original format should be entered here. If the user intends
#' to use any alternative format a dataframe with sample IDs in the first column
#' and sample purity could be eneterd here after setting the purities_purebeta_format
#' argument to FALSE (IMPORTANT! The user MUST enter sample purity values, not 1-Purity
#' values).
#'
#' @param purities_purebeta_format Default = TRUE. If the user wanted to use a
#' different input format for sample purity (see purities_samples_to_correct
#' for the format especifications) this argument should be set to FALSE.
#'
#' @param only_certain_CpGs Default = FALSE. If the beta correction has to be
#' applied only to certain CpGs and not to all the ones included in the matrix
#' provided as the betas_to_correct argument this argument should be set to
#' TRUE.
#'
#' @param CpGs_to_correct_vector. A vector of the CpG IDs to be corrected when
#' only_certain_CpGs = TRUE has been specified.
#'
#' @param refitting Default = FALSE. This argument should be set to TRUE if the
#' user wants to refit the reference regression to also take into account the
#' betas and PureBeta estimated sample purity values as additional data points.
#' This could be advisable when having a non-fully representative reference
#' dataset and a significant number of samples whose purity has been estimated
#' using PureBeta. Else, the reference regressions will be directly used for the
#' beta correction.
#'
#' @param reference_regressions The output of the reference regression generator
#' should be entered here if the refitting argument has NOT been set to TRUE.
#' Else, this argument should be ignored (both short and extended versions are
#' valid). The input list must at least include the list contains the a named
#' vector with the variance of the betas of CpGs used to build the regressions
#' (input$cpg.variance), the slopes, intercepts residual standard error and
#' degrees of freedom of the regression calculated per CpG (input$reg.slopes,
#' input$reg.intercepts, input$reg.RSE and input$df as matrices).
#'
#' @param reference_betas A matrix with CpGs as rows and analysed samples (or an
#' individual sample) as columns with the uncorrected beta values from the CpGs
#' of the samples that are intended to be used to as reference data should be
#' entered here if the refitting argument has been set to TRUE. The values must
#' be numeric, the rows must be names with the CpG ID, and the columns with the
#' sample IDs. An example of the required format is available in the
#' example_betas_reference matrix.
#'
#' @param reference_purities Named vector containing the sample purity values of
#' of the samples whose DNA methylation beta values are intended to be used as
#' reference data should been entered here if the refitting argument has been set
#' to TRUE. The vector must be named with the sample ID, which must match with
#' the sample IDs from the matrix containing the beta values. An example of the
#' required format is available in the example_purities_reference vector.
#'
#' @param use_seed Default = FALSE. A seed for the FlexMix package to detect
#' the different CpG methylation patterns can be used by setting this argument
#' to TRUE. This argument will only be used if the refitting argument has been
#' set to TRUE.
#'
#' @param seed_num Default = 2000. The seed to be used when set_seed = TRUE can
#' be specified here.
#'
#' @param cores Default = 1. Number of cores to be used to run the function in
#' parallel.
#'
#' @returns List with the corrected betas for the tumour (output$`Corrected_tumor`)
#' and microenvironment (output$`Corrected_microenvironment`) when refitting =
#' FALSE. If refitting = TRUE has been selected the corrected betas for the
#' tumour and microenvironment (output$`Corrected_betas`) will be available in
#' addition to the parameters of the refitted new regressions
#' (output$`Regression_parameters`).
#'
#' @export
#'
#' @examples
#'
#' # Using the non-refitting approach for all the CpGs
#' reference_based_beta_correction(betas_to_correct = example_betas_to_correct,
#'                                 purities_samples_to_correct = purity_estimation_output,
#'                                 only_certain_CpGs = FALSE,
#'                                 refitting = FALSE,
#'                                 reference_regressions = reference_regression_generator_output,
#'                                 cores = 5)
#'
#'# Using the non-refitting approach for certain CpGs
#' reference_based_beta_correction(betas_to_correct = example_betas_to_correct,
#'                                 purities_samples_to_correct = purity_estimation_output,
#'                                 only_certain_CpGs = TRUE,
#'                                 CpGs_to_correct_vector = c("cg09248054", "cg08231710"),
#'                                 refitting = FALSE,
#'                                 refernce_regressions = reference_regression_generator_output,
#'                                 cores = 5)
#'
#'
#' # Using the refitting approach
#' reference_based_beta_correction(betas_to_correct = example_betas_to_correct,
#'                                 purities_samples_to_correct = purity_estimation_output,
#'                                 only_certain_CpGs = FALSE,
#'                                 refitting = TRUE,
#'                                 reference_betas = example_betas_reference,
#'                                 reference_purities = example_purities_reference,
#'                                 set_seed = TRUE,
#'                                 seed_num = 1,
#'                                 cores = 5)
#'
#'
reference_based_beta_correction <- function(

  betas_to_correct,
  purities_samples_to_correct,
  purities_purebeta_format = TRUE,
  only_certain_CpGs = FALSE,
  CpGs_to_correct_vector,
  refitting,
  reference_regressions,
  reference_betas,
  reference_purities,
  use_seed = TRUE,
  seed_num = 2000,
  cores = 1

) {


  # ===========================
  # CONFIGURING PARALLELIZATION
  # ===========================

  cat("\nUsing", cores,"core(s)\n\n")

  #Creating the cluster to run the process in parallel
  cl <- makeCluster(cores)


  # =============================================
  # CORRECT BETAS REFITTING REFERENCE REGRESSIONS
  # =============================================

  if (refitting == TRUE) {


    registerDoParallel(cl)

    #Making sure that all cores have access to the flexmix package and the
    #Using invisible() to avoid printing anything to the terminal
    invisible(clusterEvalQ(cl, {library("flexmix")}))

    # Export all the functions in the package to the defined cores
    parallel::clusterExport(cl = cl,
                  varlist = unclass(lsf.str(envir = asNamespace("PureBeta"), all = TRUE)),
                  envir = as.environment(asNamespace("PureBeta")))



    # PROCESSING PREDICTED PURITIES

    cat("\nPreprocessing the data...\n\n")

    if (purities_purebeta_format == TRUE) {

      # Getting predicted 1-purities from purity_estimation output
      predicted_1mPurities <- purities_samples_to_correct$`Estimated_1-Purities`

      # Removing samples with more than one estimates (if any)
      if (nrow(predicted_1mPurities[which(predicted_1mPurities[,2]!=1),])!=0) {

        #Calculate the number of samples to remove
        samples_to_remove <- nrow(predicted_1mPurities[which(predicted_1mPurities[,2]!=1),]) / 2

        #Print warining message
        cat("\n", samples_to_remove, "samples have more than one predicted purity. Samples removed from the beta correction.\n")

        #Filtering samples with more than one purity values
        predicted_1mPurities <- predicted_1mPurities[which(predicted_1mPurities[,2]==1),]
      }

      # Transforming the predicted_purities dataframe into a vector
      predicted_purities_vec <- 1 - as.numeric(predicted_1mPurities[,3]) # Using purity, no 1-Purity
      names(predicted_purities_vec) <- predicted_1mPurities[,1]

    } else {

      # Reassigning variable names
      predicted_Purities <- purities_samples_to_correct

      # Transforming the predicted_purities dataframe into a vector
      predicted_purities_vec <- as.numeric(predicted_Purities[,2]) # Using Purity
      names(predicted_purities_vec) <- predicted_Purities[,1]

    }


    # FILTERING CPGS

    cat("\nChecking CpGs to be corrected...\n\n")

    # Use only the specified CpGs if that option has been selected
    if (only_certain_CpGs) {

      #Keeping only CpGs of interest
      betas_to_correct <- betas_to_correct[CpGs_to_correct_vector,]

    }


    # Checking if the CpGs are included in the reference data
    if (sum(!(rownames(betas_to_correct) %in% rownames(reference_betas))) != 0) {

      # Printing warning message
      cat("\n",  sum(!(rownames(betas_to_correct) %in% rownames(reference_betas))), "CpG(s) is/are not included into the reference cohort, so it/they can not be corrected.\n\n")

      # Filtering not included CpGs
      betas_to_correct <- betas_to_correct[rownames(betas_to_correct) %in% rownames(reference_betas),]
    }

    # Remove CpGs from the cohort dataset that are not included into the data to correct to speed up the process.
    reference_betas <- reference_betas[names(reference_betas) %in% names(betas_to_correct)]

    #Sorting the cohort betas dataframe based on the rownames of betas_to_correct
    reference_betas <- reference_betas[rownames(betas_to_correct)]


    # MERGING REFERENCE AND PREDICTED VALUES TO REFIT THE REFERENCE REGRESSIONS

    # Creating a single purity vector
    purities <- c(reference_purities, predicted_purities_vec)

    # Creating a single betas dataframe
    betas <- cbind(reference_betas, betas_to_correct)

    #Removing sample purities not included into the beta dataset. It generates errors
    purities <- purities[colnames(betas)]


    # RUNNING BETA CORRECTION

    cat("\nCorrecting betas refitting the reference regressions...\n\n")

    #Acding seed if necessary
    if (use_seed) {

      #Adding seed to each row of the beta value dataframe
      betaRun <- cbind(seed=seed_num:seed_num+nrow(betas),betas)

    } else {

      betaRun <- betas

    }

    #Storing sample names
    betaNames <- colnames(betas)

    print(betaRun)
    print(purities)

    # Initializing progress bar and specifying options
    pbo <- pboptions(type = "txt", char="=", txt.width=80)

    #Running the analysis in parallel with a progress bar
    #Using supress warnings to avoid the function to raise a warning when an
    #identified methylation population generates a regression with perfect fit.
    #This does not happen due to errors
    res <- suppressWarnings(pbapply(cl = cl, #Clusters to run the process
                   MARGIN = 1, #Apply the function to the rows
                   FUN = adjustBeta, #Function to correct betas
                   purity=purities, #Purity values
                   snames=betaNames, #Sample names
                   seed=use_seed, #Specify if the seed has been added to the data or not
                   betaRun #Beta values+the added seed
    ))

    # Stop clusters used in parallelization
    stopCluster(cl)

    # GENERATING RESULT LIST

    cat("\n\nGenerating output...\n\n")


    # Creating a list to add the results
    result_list <- list(
      betas.original = do.call("rbind",lapply(res,function(x) x$y.orig[names(predicted_purities_vec)])), #Original beta values
      betas.tumor = do.call("rbind",lapply(res,function(x) x$y.tum[names(predicted_purities_vec)])), #Corrected tumor beta values
      betas.microenvironment = do.call("rbind",lapply(res,function(x) x$y.norm[names(predicted_purities_vec)])) #Corrected microenvironment beta values
    )

    # Creating a list to add the parameters of the correction regressions
    reg_list <- list(
      cpg.populations =  do.call("rbind",lapply(res,function(x) x$groups)), #Methylation patterns (populations) of each CpG
      reg.slopes = do.call("rbind",lapply(res,function(x) x$res.slopes)), #Slopes of the populations
      reg.intercepts = do.call("rbind",lapply(res,function(x) x$res.int)), #Intercepts of the populations
      reg.RSE = do.call("rbind",lapply(res,function(x) x$res.rse)), #Residual standard error
      reg.df = do.call("rbind",lapply(res,function(x) x$res.df)) #Degrees of freedom of the regressions
    )

    # Filtering results. Keeping only CpGs that were intended to be corrected
    lapply(names(result_list), function(n) {
      result_list[[n]][,names(predicted_purities_vec)]
    })



    cat("\n=================\n")
    cat ("PROCESS FINISHED")
    cat("\n=================\n")

    # Return output
    return(
      list(
        "Corrected_betas" = result_list,
        "Regression_parameters" = reg_list
      )
    )

  # =====================================================
  # CORRECT BETAS WITHOUT REFITTING REFERENCE REGRESSIONS
  # =====================================================


  } else {

    # Registering clusters
    registerDoSNOW(cl)

    # Exporting all the functions of the package to the clusters
    parallel::clusterExport(cl = cl,
                  varlist = unclass(lsf.str(envir = asNamespace("PureBeta"), all = TRUE)),
                  envir = as.environment(asNamespace("PureBeta")))


    # PROCESSING PREDICTED PURITIES

    cat("\nPreprocessing the data...\n\n")

    if (purities_purebeta_format == TRUE) {

      # Getting predicted 1-purities from purity_estimation output
      predicted_1mPurities <- purities_samples_to_correct$`Estimated_1-Purities`

      # Removing samples with more than one estimates (if any)
      if (nrow(predicted_1mPurities[which(predicted_1mPurities[,2]!=1),])!=0) {

        #Calculate the number of samples to remove
        samples_to_remove <- nrow(predicted_1mPurities[which(predicted_1mPurities[,2]!=1),]) / 2

        #Print warining message
        cat("\n", samples_to_remove, "samples have more than one predicted purity. Samples removed from the beta correction.\n")

        #Filtering samples with more than one purity values
        predicted_1mPurities <- predicted_1mPurities[which(predicted_1mPurities[,2]==1),]
      }

      # Transforming the predicted_purities dataframe into a vector
      predicted_purities_vec <- 1 - as.numeric(predicted_1mPurities[,3]) # Using purity, no 1-Purity
      names(predicted_purities_vec) <- predicted_1mPurities[,1]

    } else {

      # Reassigning variable names
      predicted_Purities <- purities_samples_to_correct

      # Transforming the predicted_purities dataframe into a vector
      predicted_purities_vec <- as.numeric(predicted_Purities[,2]) # Using Purity
      names(predicted_purities_vec) <- predicted_Purities[,1]

    }

  # PROCESSING REGRESSION PARAMETERS

  # Reformatting reference regression parameters
  my_slopes <- reference_regressions$reg.slopes
  my_intercepts <- reference_regressions$reg.intercepts
  my_RSE <- reference_regressions$reg.RSE
  my_df <- reference_regressions$reg.df

  # Filtering reference regressions to account for possible errors caused by FlexMix

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

  # Check if the degrees of freedom are 0
  check_df <- function(vec) {

    # Check if ant of the elements is equal to 0
    if (0 %in% vec) {
      return(FALSE)
    } else {
      return(TRUE)
    }

  }

  # QC remove regressions with potential errors
  qc_slope_NA <- apply(my_slopes, check_na_last, MARGIN=1)
  qc_intercept_NA <- apply(my_intercepts, check_na_last, MARGIN=1)
  qc_df_NA <- apply(my_df, check_na_last, MARGIN=1)
  qc_df_0 <- apply(my_df, check_df, MARGIN=1)

  # Remove problematic CpGs from the regression list
  cpgs_to_keep <- rownames(reference_regressions$reg.slopes)[qc_slope_NA & qc_intercept_NA & qc_df_NA & qc_df_0]


  #Filtering regression objects and Beta values to correct
  my_slopes <- my_slopes[cpgs_to_keep,]
  my_intercepts <- my_intercepts[cpgs_to_keep,]
  my_RSE <- my_RSE[cpgs_to_keep,]
 
  # FILTERING CPGS

  cat("\nChecking CpGs...\n\n")

  # Use only the specified CpGs if that option has been selected
  if (only_certain_CpGs) {

    #Keeping only CpGs of interest
    betas_to_correct <- betas_to_correct[CpGs_to_correct_vector,]

  }


  # Checking if the CpGs are included in the reference regressions
  if (sum(!(rownames(betas_to_correct) %in% rownames(my_slopes))) != 0) {

    # Printing warning message
    cat("\n",  sum(!(rownames(betas_to_correct) %in% rownames(my_slopes))), "CpG(s) is/are not included into the refernce cohort, so it/they can not be corrected.\n\n")

    # Filtering not included CpGs
    betas_to_correct <- betas_to_correct[rownames(betas_to_correct) %in% rownames(my_slopes),]
  }

  # Remove CpGs from the regressions that are not included into the data to correct to speed up the process.
  my_slopes <- my_slopes[rownames(my_slopes) %in% rownames(betas_to_correct),]
  my_intercepts <- my_intercepts[rownames(my_intercepts) %in% rownames(betas_to_correct),]


  # CORRECTING BETAS BASED ON REFERENCE REGRESSIONS

  cat("\nCorrecting betas without refitting reference regressions...\n\n")


  # Configure progress bar
  p_bar <- txtProgressBar(min=0,
                          max=nrow(betas_to_correct),
                          style=3,
                          width=80)


  # Creating a function to follow the execution of the script
  progress <- function(n) setTxtProgressBar(p_bar, n)
  opts <- list(progress = progress)

  # Correcting betas through a parallelized for loop
  output <- foreach(cpg = rownames(betas_to_correct), .packages = "Kendall", .options.snow = opts) %dopar% {

    # ASSIGN REGRESSION TO THE DIFFERENT SAMPLES FOR EACH CPG

    identified_regressions <- identify_regression(
      vec_betas = as.numeric(betas_to_correct[cpg,]),
      vec_estimated_1mPurity = 1 - predicted_purities_vec,
      vec_slopes = as.numeric(my_slopes[cpg,]),
      vec_intercepts = as.numeric(my_intercepts[cpg,])
    )


    # CORRECTING BETAS BASED ON THE IDENTIFIED REGRESSIONS

    return(
      list(

        "Tumour" =  correcting_betas(
          slopes_vec = as.numeric(identified_regressions$Slope),
          intercepts_vec = as.numeric(identified_regressions$Intercept),
          distances_vec = as.numeric(identified_regressions$Distance),
          to_correct = "Tumor"
        ),
        "Microenvironment" = correcting_betas(
          slopes_vec = as.numeric(identified_regressions$Slope),
          intercepts_vec = as.numeric(identified_regressions$Intercept),
          distances_vec = as.numeric(identified_regressions$Distance),
          to_correct = "Microenvironment"
        )
      )
    )

    }

    # Stop clusters used in parallelization
    stopCluster(cl)

    # PROCESSING OUTPUT

    cat("\n\nGenerating output...\n\n")

    #Converting output list into dataframe of corrected tumor betas
    corrected_tumor <- as.data.frame(do.call(rbind, lapply(output, function(item) item[[1]])))
    colnames(corrected_tumor) <- colnames(betas_to_correct)
    rownames(corrected_tumor) <- rownames(betas_to_correct)

    #Converting output list into dataframe of corrected tumor betas
    corrected_microenvironment <- as.data.frame(do.call(rbind, lapply(output, function(item) item[[2]])))
    colnames(corrected_microenvironment) <- colnames(betas_to_correct)
    rownames(corrected_microenvironment) <- rownames(betas_to_correct)


    cat("\n\n=================\n")
    cat ("PROCESS FINISHED")
    cat("\n=================\n\n")

    return(
      list(
        "Corrected_tumour" = corrected_tumor,
        "Corrected_microenvironment"  = corrected_microenvironment
      )
    )
  }
}
