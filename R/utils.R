# Internal functions used in the PureBeta package

#
# Function to correct beta values based on a cohort of samples with known
# purities. This function is also used to determine reference regressions
# for each CpG.
#

adjustBeta <- function(

    methylation,
    purity,
    snames,
    nmax = 3,
    nrep = 3,
    seed = TRUE

  ) {

    #If a seed is provided as the first element of the methylation vector (seed=TRUE)
    #set seed for the cluster determination and remove it from the methylation vector
    if(seed) {
      set.seed(as.integer(methylation[1]))
      methylation<-methylation[-1]
    }

    #Defining variables for clustering and regression
    x <- as.numeric(purity)
    x2 <- 1-as.numeric(purity)
    y <- as.numeric(methylation)

    #Calculate global correlation between beta and purity
    gl.corr <- suppressWarnings(cor(x,y))
    gl.corr[is.na(gl.corr)] <- 0
    gl.corr <- round(gl.corr,3)

    #Add small gaussian noise to x (avoid errors when large number of zero samples)
    y2<-y+rnorm(length(y),mean=0,sd=.005)

    #Modeling the methylation patterns (CpG populations)
    model <- stepFlexmix(y2 ~ x,k = 1:nmax, nrep = nrep,verbose = FALSE)
    model <- getModel(model, "BIC")
    cl <- clusters(model)
    
    #Make sure clusters are numbered 1 to 3, odd cases exist where one pop has zero members from flexMix
    #can rename because flexmix object not used after this
    cl <- as.integer(factor(cl))



    ## CALCULATING REGRESSIONS

    #Determining the b_vs_pur regressions for the clusters identified
    b_vs_pur <- lapply(1:nmax,function(z) {
      if(z %in% cl) {
        lm(y[cl==z]~x[cl==z]) #Beta VS 1-Purity regression

      } else { NA } #If less than 3 populations were detected NA is added

    })

    #Determining the b_vs_1mp regressions for the clusters identified
    b_vs_1mp <- lapply(1:nmax,function(z) {
      if(z %in% cl) {
        lm(y[cl==z]~x2[cl==z]) #Beta VS 1-Purity regression
      } else { NA } #If less than 3 populations were detected NA is added
    })


    ## DETERMINING THE FUNCTION'S OUTPUT FROM THE CALCULATED REGRESSIONS

    output_list <- list() #Creating a list to store the data

    #Adding the population to which each CpG belongs to the list
    output_list$groups <- cl

    #Adding the number of population to which each CpG belongs to the list
    output_list$n.groups <- length(levels(factor(cl)))

    #Adding the correlation between betas and purity to the list
    output_list$glob.corr <- gl.corr

    #Adding corrected microenvironment betas to the list
    output_list$y.norm <- sapply(
      X = 1:nmax,
      FUN = function(z) {
        if (z %in% cl) {
          n_vals <- coefficients(b_vs_pur[[z]])[1]+residuals(b_vs_pur[[z]])
          names(n_vals)<-snames[cl==z]
          n_vals

        } else {
          NULL
        }
      })
    output_list$y.norm <- unlist(output_list$y.norm)[snames]

    #Adding corrected tumor betas to the list
    output_list$y.tum <- sapply(
      X = 1:nmax,
      FUN = function(z) {
        if (z %in% cl) {
          n_vals <- coefficients(b_vs_1mp[[z]])[1]+residuals(b_vs_1mp[[z]])
          names(n_vals)<-snames[cl==z]
          n_vals

        } else {
          NULL
        }
      })
    output_list$y.tum <- unlist(output_list$y.tum)[snames]


    ##in very rare instances flexmix calls 3 populations but one groups has zero members. -> has parameters for non-existant pop in output?!
    #output_list$res.int <- round(as.numeric(unlist(lapply(slot(model,"components"),function(z) slot(z[[1]],"parameters")$coef[1]))),3)
    #output_list$res.slopes <- round(as.numeric(unlist(lapply(slot(model,"components"),function(z) slot(z[[1]],"parameters")$coef[2]))),3)

    #Getting line parameters from b_vs_1mp: INTERCEPTS
    output_list$res.int <- sapply(
      X = 1:nmax,
      FUN = function(z) {
        if (z %in% cl) {
          int <- coefficients(b_vs_1mp[[z]])[1]
          round(as.numeric(int),3)

        } else {
          NA
        }
      },
      simplify = TRUE)

    #Getting line parameters from b_vs_1mp: SLOPES
    output_list$res.slopes <- sapply(
      X = 1:nmax,
      FUN = function(z) {
        if (z %in% cl) {
          int <- coefficients(b_vs_1mp[[z]])[2]
          round(as.numeric(int),3)

        } else {
          NA
        }
      },
      simplify = TRUE)


    #Getting line parameters from b_vs_1mp: RESIDUAL STANDARD ERROR
    output_list$res.rse <- sapply(
      X = 1:nmax,
      FUN = function(z) {
        if (z %in% cl) {
          RSE <- summary(b_vs_1mp[[z]])$sigma
          round(as.numeric(RSE),6)

        } else {
          NA
        }
      },
      simplify = TRUE)

    #Getting line parameters from b_vs_1mp: RESIDUAL STANDARD ERROR
    output_list$res.df <- sapply(
      X = 1:nmax,
      FUN = function(z) {
        if (z %in% cl) {
          df.residual(b_vs_1mp[[z]])

        } else {
          NA
        }
      },
      simplify = TRUE)

    #Cap the corrected betas to 0 and 1
    output_list$y.tum[output_list$y.tum > 1] <- 1
    output_list$y.tum[output_list$y.tum < 0] <- 0
    output_list$y.norm[output_list$y.norm > 1] <- 1
    output_list$y.norm[output_list$y.norm < 0] <- 0

    #Round the betas to three decimal places
    output_list$y.tum <- round(output_list$y.tum,3)
    output_list$y.norm <- round(output_list$y.norm,3)

    #Adding the rounded original betas to the list
    output_list$y.orig <- round(methylation,3)

    return(output_list)

  }


#
# Function to determine a purity interval for each CpG's beta value of each
# sample
#

predicting_purity <- function(

    beta,
    slopes,
    intercepts,
    RSE,
    degrees_of_freedom,
    slope_threshold,
    alpha,
    populations = NULL,
    original_betas = NULL,
    original_purities = NULL,
    B = NULL,
    assume_t_distribution

  ) {

    if (assume_t_distribution) {

      # Identifying the regressions whose residual standard error is over the threshold
      # and whose slope is too close to 0 to be informative
      to_be_ignored <- (slopes > -slope_threshold & slopes < slope_threshold) | is.na(slopes)

      # Removing the values that do not meet the conditions previously analysed.
      slopes <- slopes[!to_be_ignored]
      intercepts <- intercepts[!to_be_ignored]
      RSE <- RSE[!to_be_ignored]
      degrees_of_freedom <- degrees_of_freedom[!to_be_ignored]

      # Deleting unnecessary variables to improve memory usage
      rm(to_be_ignored, pops_to_keep)

      # Check that all the regressions are not uninformative and execute the following code
      # to estimate the 1-Purity value. If they weren't NA would be assigned to the 1-Purity value,
      # so that CpG would not be used to estimate the final 1-Purity of the sample

      if (length(slopes)!=0) {

        #Creating a vector to store the regressions in which the beta value analysed is included using booleans
        matches <- ((beta - intercepts) / slopes >= 0) & ((beta - intercepts) / slopes <= 1)

        # If the length of the matches vector when matches==TRUE is one, a unique CpG methylation pattern (regression)
        # can be assigned to each beta value. If this value is different to one no significant association to a pattern can
        # be made
        if (sum(matches)==1) {

          #Store the identified population (indexes of matches vector) into a variable
          identified_pop <- which(matches)

          #Predicting the 1-purity value from the inverse regression
          predicted_1_minus_p <- (beta - intercepts[identified_pop]) /slopes[identified_pop]

          # Get the prediction interval of the beta identified
          beta_int <- c(
            beta - qt((1-alpha)/2, degrees_of_freedom[identified_pop]) * RSE[identified_pop],
            beta + qt((1-alpha)/2, degrees_of_freedom[identified_pop]) * RSE[identified_pop]
          )


          #In some cases if the degrees of freedom are too low the result of qt(), the quantile of the Student's
          #t-distribution can be undefined (NaN). Those cases are very uncommon,  but they can generate bugs in the
          #code. Therefore, if that happened, the one_minus_purity interval would be set to (NA, NA) in order to avoid
          #possible bugs

          if (anyNA(beta_int)) {

            one_minus_purity <- c(NA,NA)

          } else {

            #Getting the prediction interval of the 1 - purity value. Extrapolate from beta interval.
            one_minus_purity <- (beta_int - intercepts[identified_pop]) / slopes[identified_pop]
            one_minus_purity <- sort(one_minus_purity)

            one_minus_purity[one_minus_purity < 0] <- 0
            one_minus_purity[one_minus_purity > 1] <- 1

          }

        } else {
          #If the CpG can not be assigned to a single population (methylation pattern) the 1-purity value will not be assigned.
          one_minus_purity <- c(NA,NA)
        }

        #If all the regressions are uninformative NA value will be assigned to one_minus_purity.
      } else {
        one_minus_purity <- c(NA,NA)
      }

      #Setting the interval's endpoints to 0 or 1 if one of the limits lower or higher than those values.
      #The NAs are ignored
      if (!(is.na(one_minus_purity[1] & is.na(one_minus_purity[2])))) {
        one_minus_purity[1] <- max(0, one_minus_purity[1])
        one_minus_purity[2] <- min(1, one_minus_purity[2])
      }

      #Returning the 1-purity vector
      return(one_minus_purity)

    } else {

      #Identifying the regressions whose slope is too close to 0 to be informative
      to_be_ignored <- (slopes > -slope_threshold & slopes < slope_threshold) | is.na(slopes)

      # Removing the values that do not meet the conditions previously analysed.
      slopes <- slopes[!to_be_ignored]
      intercepts <- intercepts[!to_be_ignored]
      RSE <- RSE[!to_be_ignored]
      degrees_of_freedom <- degrees_of_freedom[!to_be_ignored]

      # Filtering original betas
      pops_to_keep <- which(!to_be_ignored)
      original_betas <- original_betas[populations %in% pops_to_keep]
      original_purities <- original_purities[populations %in% pops_to_keep]

      # Filtering populations
      populations <- match(populations, sort(pops_to_keep))[!is.na(match(populations, sort(pops_to_keep)))]

      # Deleting unnecessary variables to improve memory usage
      rm(to_be_ignored, pops_to_keep)

      # Check that all the regressions are not uninformative and execute the following code
      # to estimate the 1-Purity value. If they weren't NA would be assigned to the 1-Purity value,
      # so that CpG would not be used to estimate the final 1-Purity of the sample

      if (length(slopes)!=0) {

        #Creating a vector to store the regressions in which the beta value analysed is included using booleans
        # The matches are bounded to the 0-1 interval.
        matches <- ((beta - intercepts) / slopes >= 0) & ((beta - intercepts) / slopes <= 1)

        # If the length of the matches vector when matches==TRUE is one, a unique CpG methylation pattern (regression)
        # can be assigned to each beta value. If this value is different to one no significant association to a pattern can
        # be made
        if (sum(matches)==1) {

          #Store the identified population (indexes of matches vector) into a variable
          identified_pop <- which(matches)

          # Remove unnecesary variables
          rm(matches)

          ### INTERVAL PREDICTION THROUGH BOOTSTRAPPING

          # Getting values to perform booustrapping
          my_values <- data.frame("purities"= 1 - unname(original_purities[populations == identified_pop]), 
                                  "betas" = unname(original_betas[populations == identified_pop]))

          # Generate purity value to predict interval
          to_predict <- data.frame("purities"=(beta - intercepts[identified_pop]) /slopes[identified_pop])

          # Generate empty vector to store predicted values
          pred <- rep(NaN, B)

          # Calculate values to estimate prediction interval.
          # Take into account the error from estimatying a values and random variation of individual values sampled to determine prediction interval
          for (i in 1:B) {
            new_data <- my_values[sample(nrow(my_values), nrow(my_values), replace = TRUE),]
            fit.l <- lm(betas ~ purities, data = my_values)
            eps <- sample(residuals(fit.l), size = 1)
            pred[i] <- predict(fit.l, to_predict, type = "response") + eps
          }

          rm(my_values, to_predict)
          
          # Determining interval endpoints based on the chosen alpha
          bound_1 <- (1 - alpha)/2
          bound_2 <- 1 - (1 - alpha)/2

          # Calculating botstrapping based prediction interval
          beta_int <- c(
            unname(quantile(pred, probs = bound_1)),
            unname(quantile(pred, probs = bound_2))
          )


          #In some cases if the degrees of freedom are too low the result of qt(), the quantile of the Student's
          #t-distribution can be undefined (NaN). Those cases are very uncommon,  but they can generate bugs in the
          #code. Therefore, if that happened, the one_minus_purity interval would be set to (NA, NA) in order to avoid
          #possible bugs

          if (anyNA(beta_int)) {

            one_minus_purity <- c(NA,NA)

          } else {

            #Getting the prediction interval of the 1 - purity value. Extrapolate from beta interval.
            one_minus_purity <- (beta_int - intercepts[identified_pop]) / slopes[identified_pop]
            one_minus_purity <- sort(one_minus_purity)

            one_minus_purity[one_minus_purity < 0] <- 0
            one_minus_purity[one_minus_purity > 1] <- 1

          }

        } else {
          #If the CpG can not be assigned to a single population (methylation pattern) the 1-purity value will not be assigned.
          one_minus_purity <- c(NA,NA)
        }

        #If all the regressions are uninformative NA value will be assigned to one_minus_purity.
      } else {
        one_minus_purity <- c(NA,NA)
      }

      #Setting the interval's endpoints to 0 or 1 if one of the limits lower or higher than those values.
      #The NAs are ignored
      if (!(is.na(one_minus_purity[1] & is.na(one_minus_purity[2])))) {
        one_minus_purity[1] <- max(0, one_minus_purity[1])
        one_minus_purity[2] <- min(1, one_minus_purity[2])
      }

      #Returning the 1-purity vector
      return(one_minus_purity)

    }


  }


#
# Function to determine the coverage of the 1-Purity intervals determined using
# reference regressions
#

purity_coverage <- function(

    pred_purity_confidence,
    interval_threshold,
    min_endpoint = 0,
    max_endpoint = 1

  ) {

  #Creating a list to store the output
  output_list = list()

  #Removing NA values from the predicted purity dataframe and formatting the values
  #to three decimal positions
  pred_purity_confidence <- format(round(na.omit(pred_purity_confidence),3),nsmall=3)


  # =============================================
  # DETERMINE THE COVERAGE OF THE PURITY INTERVAL
  # =============================================

  #Creating a vector to store the sections with the maximum coverage
  coverage_per_section <- c()

  # Divide the (min_endpoint,max_endpoint) interval in sections of length 0.001.
  #Generating a list with the each section of the interval as keys (the keys
  #are written as characters)
  sections <- seq(min_endpoint,max_endpoint,by=0.001)
  coverage_per_section <- setNames(rep(0,length(sections)), as.character(sections))


  #Check if the sections of the 0-1 purity interval are covered in each interval and determine the coverage (how
  #many times is each section includedin the intervals) per section

  for (cpg in rownames(pred_purity_confidence)) {

    pos <- as.character(seq(pred_purity_confidence[cpg,1], pred_purity_confidence[cpg,2], by=0.001))
    included_section <- names(coverage_per_section) %in% pos
    coverage_per_section[included_section] <- coverage_per_section[included_section] + 1

  }


  #Correcting the overrepresentation of purity values between 0.8 and 1. Fitting linear regression and using the residuals
  coverage_per_section <- setNames(residuals(lm(unname(coverage_per_section)~as.numeric(names(coverage_per_section)))),names(coverage_per_section))

  #Smoothening the plot using spline
  smooth <- smooth.spline(x=as.numeric(names(coverage_per_section)),
                          y=unname(coverage_per_section),
                          n=30)


  #Predict values per each section using the smoothed function
  smoothed_coverage_values <- predict(smooth, newdata=list(x=sections))$y
  names(smoothed_coverage_values) <- as.character(sections)

  #Getting the maximum corrected coverage value
  max_ccov <- max(smoothed_coverage_values)



  # CORRECTING THE CORRECTION

  # In some cases, when there is only one peak between 0.75 and 1 because
  # the noise peak and the actual peak are mixed, the correction method may
  # generate the predicted 1-Purity to be 0. In order to deal with that, the
  # intercept of the regression used to correct the data will be set to the first corrected coverage value when
  # the 1-P is predicted to be 0 and the 1-P recalculated based on that.

  if (sections[which(smoothed_coverage_values == max_ccov)] == 0) {

    # Determining the intercept. The first coverage value of the 0-1 1-Purity range.
    inter <- unname(coverage_per_section["0"])

    #Correcting the overrepresentation of purity values between 0.8 and 1. Fitting linear regression and using the residuals and fixed intercept
    # The intercept fixation is done this way to avoid R from interpreting inter as another variable of the linear model it is
    # building, so inter is subtracted to the coverage, and then the intercept is set to 0.
    corrected_coverage <- setNames(residuals(lm(I(unname(coverage_per_section) - inter) ~ as.numeric(names(coverage_per_section)) + 0)),names(coverage_per_section))

    #Smoothening the plot using spline
    smooth <- smooth.spline(x=as.numeric(names(corrected_coverage)),
                            y=unname(corrected_coverage),
                            n=30)


    #Predict values per each section using the smoothed function
    smoothed_coverage_values <- predict(smooth, newdata=list(x=sections))$y
    names(smoothed_coverage_values) <- as.character(sections)

    #Getting the maximum corrected coverage value
    max_ccov <- max(smoothed_coverage_values)

  }

  # ======================================================
  # DETERMINE THE MAXIMUM COVERAGE ESTIMATES AND INTERVALS
  # ======================================================

  ## GETTING THE ESTIMATES

  #Appending the max value(s), the 1-Purity estimate(s), to the output_list
  output_list[["1-Pur_estimates"]] <- sections[which(smoothed_coverage_values == max_ccov)]

  ## GETTING THE INTERVALS

  # Get the intervals with the maximum coverage. The minimum coverage threshold is the maximum value minus the interval_
  #threshold percentage selected (default value is 10%)
  selected_values <- sections[which(smoothed_coverage_values >= max_ccov*(1-interval_threshold/100))]

  #Defining variables to store parameters of the intervals
  start_val <- selected_values[1] # Start point of the interval
  ref_val <- NULL # Value to use as the refernce for the ref value of the loop
  end_val <- NULL # End point of the interval

  interval_list <- list() # A list to store the intervals detected

  # Iterate through the values over the coverage threshold (except the first element, as
  # it has already been assigned to start_val)
  for (val in selected_values[-1]) {

    #If ref_val is not defined assign it to the start_value
    if (is.null(ref_val)) {
      ref_val <- start_val
    }

    # Compare the current value with ref_val to check if they are contiguous (the differnece between them is 0.001)
    if (format(val - ref_val,2)==0.001) {
      #If they were contiguous the current value would be assigned to ref_val
      ref_val <- val
    } else {
      #If they are not contiguous the reference value will be selected as the endpoint of the interval
      end_val <- ref_val
      #The interval is added to the interval list
      interval_list <- c(interval_list, list(c(start_val, end_val)))


      #Reestarting the loop after identifying an interval
      start_val <- val # Assign the current value to start_val
      end_val <- NULL # Assign NULL to end_val
      ref_val <- NULL # Assign NULL to ref_val

    }
  }

  #Add the last selected interval to the interval list at the end of the for loop
  if (!is.null(ref_val)) {
    end_val <- ref_val #Assign the ref_val to end_val if ref_val is defined
  } else {end_val <- start_val}

  #Add the interval to the interval list
  interval_list <- c(interval_list, list(c(start_val, end_val)))


  # ==============================================
  # REMOVING THE INTERVALS WITHOUT A MAXIMUM VALUE
  # ==============================================

  # Checking if any of the maxs identified are inside the detected intervals. If they are not, the interval will be removed
  interval_list <- Filter(function(interval) {
    length(intersect(seq(interval[1], interval[2], by = 0.001), format(output_list[["1-Pur_estimates"]], 3))) != 0
  }, interval_list)

  #Adding the intervals detected to the output list
  output_list[["interval(s)"]] <- interval_list


  #The maximum coverage interval list will be returned
  return(output_list)

  }


#
# Function to identify reference regression (methylation pattern or population)
# which each CpG of each sample belongs to
#

identify_regression <- function(

    vec_betas,
    vec_estimated_1mPurity,
    vec_slopes,
    vec_intercepts

  ) {

    # Checking if the arguments (vectors and simple numeric arguments) are numeric
    if (!is.numeric(vec_betas) | !is.numeric(vec_slopes) | !is.numeric(vec_intercepts)) {
      stop("Betas, slopes and intercepts must be numeric to identify the regression.")
    }

    # Initializing dataframe to store the distances of each sample's CpG to each population
    distances_df <- data.frame(
      Beta = vec_betas,
      Purity = as.numeric(vec_estimated_1mPurity),
      Distance_1 = rep(NA, length(vec_betas)),
      Distance_2 = rep(NA, length(vec_betas)),
      Distance_3 = rep(NA, length(vec_betas))
    )

    # Filling the matrix using apply
    distances_matrix <- t(apply(distances_df, 1, function(row) {
    
        # If beta is NA, return NA for all distances
        if (is.na(row["Beta"])) {
          return(c(Distance_1 = NA, Distance_2 = NA, Distance_3 = NA))
        }
    
        # Compute distances only if slopes and intercepts are available
        distance_1 <- if (!is.na(vec_slopes[1]) & !is.na(vec_intercepts[1])) {
          row["Beta"] - (vec_slopes[1] * row["Purity"] + vec_intercepts[1])
        } else { NA }
    
        distance_2 <- if (!is.na(vec_slopes[2]) & !is.na(vec_intercepts[2])) {
          row["Beta"] - (vec_slopes[2] * row["Purity"] + vec_intercepts[2])
        } else { NA }
    
        distance_3 <- if (!is.na(vec_slopes[3]) & !is.na(vec_intercepts[3])) {
          row["Beta"] - (vec_slopes[3] * row["Purity"] + vec_intercepts[3])
        } else { NA }
    
        return(c(Distance_1 = distance_1, Distance_2 = distance_2, Distance_3 = distance_3))
    }))

    #Determining the population (vector index) with the lowest absolute Euclidean distance. If the distances are equal the first
    #population with be chosen by default
    pop_identified <- apply(distances_matrix, 1, function(row) {
        if (any(is.na(row))) {
            return(NA)
        } else {
            return(which.min(abs(row)))
      }
    })

    #Generating and returning output dataframe with the parameters of the identified regressions
    output_df <- data.frame(
      Sample = names(vec_estimated_1mPurity), # Adding sample names
      Slope = ifelse(is.na(pop_identified), NA, vec_slopes[pop_identified]),
      Intercept = ifelse(is.na(pop_identified), NA, vec_intercepts[pop_identified]),
      Distance = ifelse(is.na(pop_identified), NA, distances_matrix[cbind(1:nrow(distances_df), pop_identified)])
    )

    return(output_df)

  }


#
# Function to correct beta values based on precomputed reference regressions
#

correcting_betas <- function(

    slopes_vec,
    intercepts_vec,
    distances_vec,
    to_correct

  ) {


#    if (!is.numeric(slopes_vec) | !is.numeric(intercepts_vec)) {
#      stop("Slope and intercept must be numeric to correct betas.")
#    }

    if (to_correct=="Tumor") {

       # Only perforn beta correction if all the distances to regressions are not NA
       if (all(is.na(distances_vec))) {
        return(NA)
       } else {
          #The tumor beta value will be obtained using the intercept and the calculated distance.
          tum_betas <- intercepts_vec + distances_vec
    
          #The maximum possible value will always be kept below or equal to 1 and minimum to 0
          tum_betas <- sapply(tum_betas, function(x) if(x > 1) {1} else {x})
          tum_betas <- sapply(tum_betas, function(x) if(x < 0) {0} else {x})
    
          return(tum_betas)
       }

    } else if (to_correct=="Microenvironment") {

       # Only perforn beta correction if all the distances to regressions are not NA
       if (all(is.na(distances_vec))) {
        return(NA)
       } else {
           
          #The microenvironment beta value will be obtained using the intercept and slope when 1-P=1 and the calculated distance.
          #The minimum possible value will always be kept below or equal to 1
          env_betas <- intercepts_vec + slopes_vec + distances_vec
    
          #The maximum possible value will allways be kept below or equal to 1 and minimum to 0
          env_betas <- sapply(env_betas, function(x) if(x > 1) {1} else {x})
          env_betas <- sapply(env_betas, function(x) if(x < 0) {0} else {x})
           
       }

      return(env_betas)

    }
  }
