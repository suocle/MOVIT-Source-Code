#' Obtain Y-Y plot R^2 from IV and OLS regressions
#'
#' @param  X A length n vector or a size (n x 1) matrix / data frame of exposure
#' @param  Y A length n vector or a size (n x 1) matrix / data frame of outcome
#' @param  Z A length n vector (in case of only one instrument) or a size (n x k) matrix / data frame of instrument(s), k is the number of instruments
#' @param  cov A length n vector (in case of only one covariate) or a size (n x q) matrix / data frame of covariate(s), q is the number of covariates
#' @return The Y-Y plot R^2 between the residuals of the IV and OLS regressions
#' @export
yy.r2 <- function(X, Y, Z, cov = NULL) {
  # Check if inputs are vectors or matrices/data frames and convert to data frame if necessary
  if (is.vector(X)) X <- data.frame(X = X)
  if (is.vector(Y)) Y <- data.frame(Y = Y)
  if (is.vector(Z)) Z <- data.frame(Z = Z)
  if (!is.null(cov) && is.vector(cov)) cov <- data.frame(cov = cov)
  
  # Check if X, Y, Z, and cov have the same number of rows
  n <- nrow(X)
  if (n != nrow(Y) || n != nrow(Z) || (!is.null(cov) && n != nrow(cov))) {
    stop("All inputs must have the same number of rows")
  }
  
  # Rename columns of Z and cov if necessary
  if (ncol(Z) > 1) {
    colnames(Z) <- paste0("Z", 1:ncol(Z))
  } else {
    colnames(Z) <- "Z"
  }
  
  if (!is.null(cov) && ncol(cov) > 1) {
    colnames(cov) <- paste0("cov", 1:ncol(cov))
  } else if (!is.null(cov)) {
    colnames(cov) <- "cov"
  }
  
  # Combine all data into one data frame
  df <- cbind(Y, X, Z, cov)
  
  # Construct formulas for IV and OLS regressions
  IV_formula <- as.formula(paste("Y ~", paste(c("X", colnames(Z), colnames(cov)), collapse = " + "), "- 1"))
  OLS_formula <- as.formula(paste("Y ~", paste(c("X", colnames(cov)), collapse = " + ")))
  
  fit.IV <- ivreg(IV_formula, data = df)
  fit.OLS <- lm(OLS_formula, data = df)
  
  # Calculate and return squared correlation of residuals
  r2 <- cor(summary(fit.IV)$residuals, summary(fit.OLS)$residuals)^2
  return(r2)
}

#' Compute the Y-Y plot R^2 cutoff for given parameters
#'
#' @param alpha The association between the instrument and the exposure
#' @param sigmaX2 Noise level for exposure equation (default is 1)
#' @param sigmaY2 Noise level for outcome equation (default is 1)
#' @param sigmaXY Covariance between the exposure and outcome noise (default is 0.1)
#' @param n Sample size
#' @param n.rep Number of replications used to calculate the cutoff (default is 1000)
#' @param probs Desired quantile for Y-Y plot R^2 (default is 0.05), can be a vector
#' @param sigmaZ2 Variance of the instrument (default is 1)
#' @return The Y-Y plot R^2 cutoff value corresponding to the specified quantile
#' @export
r2.cutoff <- function(alpha, sigmaX2 = 1, sigmaY2 = 1, sigmaXY = 0.1, 
                      n, n.rep = 1000, probs = 0.05, sigmaZ2 = 1) {
  
  # Initialize vector to store R^2 values from simulations
  r2.vec <- rep(NA, n.rep)
  
  # Set the true coefficient for the effect of X on Y
  beta <- 1
  
  # Define the covariance matrix for the multivariate normal distribution
  cov <- matrix(0, 3, 3)
  cov[1,1] <- sigmaZ2   # Variance of Z
  cov[2,2] <- sigmaX2   # Variance of epsX
  cov[3,3] <- sigmaY2   # Variance of epsY
  cov[2,3] <- cov[3,2] <- sigmaXY    # Covariance between epsX and epsY
  
  # Loop over the number of repetitions to perform the simulation
  for (i in 1:n.rep) {
    # Set a random seed for reproducibility
    set.seed(i)
    
    # Generate trivariate normal random variables for Z, epsX, and epsY
    Z_epsX_epsY <- rmvn(n, rep(0, 3), cov)
    
    # Extract Z, epsX, and epsY from the generated data
    Z <- Z_epsX_epsY[,1,drop=FALSE]
    epsX <- Z_epsX_epsY[,2,drop=FALSE]
    epsY <- Z_epsX_epsY[,3,drop=FALSE]
    
    # Generate the exposure variable X and the outcome variable Y
    X <- alpha * Z + epsX
    Y <- beta * X + epsY
    
    # Create a data frame with the generated data
    dat <- data.frame(X, Y, Z)
    
    # Fit the IV and OLS models
    fit.IV <- ivreg(Y ~ X - 1 | Z - 1, data = dat)
    fit.OLS <- lm(Y ~ X + Z - 1, data = dat)
    
    # Calculate the squared correlation of the residuals and store in r2.vec
    r2.vec[i] <- cor(summary(fit.IV)$residuals, summary(fit.OLS)$residuals)^2
  }
  
  # Compute and return the specified quantile of the R^2 values
  return(unname(quantile(r2.vec, probs = probs)))
}

#' Generate synthetic data for instrumental variable analysis
#'
#' @param alpha The association between the instrument (Z) and the exposure (X)
#' @param beta1 The causal effect of exposure (X) on the outcome (Y)
#' @param beta2 The degree of violation for the exclusion restriction (ER) assumption (if not 0)
#' @param sigmaX2 Noise level for the exposure equation (default is 1)
#' @param sigmaY2 Noise level for the outcome equation (default is 1)
#' @param sigmaXY Covariance between the noise in the exposure and outcome equations (default is 0.1)
#' @param n Sample size
#' @param n.rep Number of replications (not used in this function, but included for consistency with similar functions)
#' @param sigmaZ2 Variance of the instrument (default is 1)
#' @param sigmaZY Covariance between the instrument and the noise in the outcome equation (default is 0.1)
#' @param seed Seed for random number generation (default is 06212024)
#' @return A data frame containing the generated synthetic data with columns X, Y, and Z
#' @export
gen.data <- function(alpha, beta1, beta2, sigmaX2 = 1, sigmaY2 = 1, 
                     sigmaXY = 0.1, n, n.rep = 1000, sigmaZ2 = 1, 
                     sigmaZY = 0.1, seed = 06212024) {
  
  # Define the covariance matrix for the multivariate normal distribution
  cov <- matrix(0, 3, 3)
  cov[1,1] <- sigmaZ2   # Variance of Z
  cov[2,2] <- sigmaX2   # Variance of epsX
  cov[3,3] <- sigmaY2   # Variance of epsY
  cov[1,3] <- cov[3,1] <- sigmaZY    # Covariance between Z and epsY
  cov[2,3] <- cov[3,2] <- sigmaXY    # Covariance between epsX and epsY
  
  # Set a random seed for reproducibility
  set.seed(seed)
  
  # Generate trivariate normal random variables for Z, epsX, and epsY
  Z_epsX_epsY <- rmvn(n, rep(0, 3), cov)
  
  # Extract Z, epsX, and epsY from the generated data
  Z <- Z_epsX_epsY[,1, drop = FALSE]
  epsX <- Z_epsX_epsY[,2, drop = FALSE]
  epsY <- Z_epsX_epsY[,3, drop = FALSE]
  
  # Generate the exposure variable X
  X <- alpha * Z + epsX
  
  # Generate the outcome variable Y
  Y <- beta1 * X + beta2 * Z + epsY
  
  # Create a data frame with the generated data
  dat <- data.frame(X, Y, Z)
  
  # Return the generated data
  return(dat)
}
