############Simulations

# Simulate wrapper
simCAT <- function(n = 10,
                   store = FALSE,
                   location = getwd(),
                   items = createTestBank("GPCM"),
                   test = initTest(items),
                   fixed_theta = NULL,
                   theta0 = list(mean = rep(0, items$Q),
                                 covar = diag(items$Q)),
                   prior = theta0)
  
{
  cat("\n", items$model, test$selection, test$estimator, fixed_theta)
  
  if (store) { 
    dir <- paste0(location,"/simulations/",paste(items$model, test$selection, test$estimator, sep="-"),"/",make.names(Sys.time()))
    cat(" - Storing in ", location, "/simulations/", sep='')
    dir.create(dir, TRUE, TRUE)
  }
  
  cat("\n")
  results <- list()
  thetas <- estimates <- variances <- matrix(NA, n, test$items$Q)
  
  start.time <- Sys.time()
  
  for (i in 1:n) {
    # draw a theta, if not specified (draws from prior, which means the prior is always correct - slightly cheating) (not used in reported simulations)
    if (is.null(fixed_theta)) theta <- rmvnorm(1, theta0$mean, theta0$covar)
    else theta <- fixed_theta
    # run the cat
    results[[i]] <- person <- ShadowCAT(initPerson(items, theta = theta, prior = prior$covar), test, FALSE)
    thetas[i,] <- theta
    estimates[i,] <- person$estimate
    variances[i,] <- diag( attr(person$estimate, 'variance') )
    
    if (store) save(person, file = paste0(dir, "/", i, '.CAT'))
    
    # some feedback, time
    if (i > 2) {
      time.spent <- difftime(Sys.time(), start.time, units = 'secs')
      time.left <- (time.spent / i) * (n - i)
      time.spent <- format(.POSIXct(time.spent, tz = "GMT"), "%H:%M:%S")
      time.string <- format(.POSIXct(time.left, tz = "GMT"), "%H:%M:%S")
    } else {
      time.spent <- 'unknown'
      time.string <- 'unknown'
    }
    cat("\r", i, " of ", n, ", deviation: ", paste(round(person$estimate - theta, 2), collapse = ' '), ", time spent: ", format(time.spent, format="%T"), ", time left: ", format(time.string, format = "%T"), "                 ",sep='')
    
  }
  
  deviations <- estimates - thetas
  bias <- colMeans(deviations)
  rmse <- sqrt(colMeans(deviations**2))
  mean_var <- colMeans(variances)
  SE_var <- apply(variances, 2, sd)
  
  summary <- data.frame(bias = bias, RMSE = rmse, mean_var = mean_var, SE_var = SE_var)
  cat('\n')
  print(round(summary,2))
  
  if (store) save(thetas, estimates, deviations, bias, rmse, mean_var, SE_var, test, results, file = paste0(dir, "/!RESULTS!"))
  if (store) save(summary, file = paste0(dir, "/!SUMMARY!"))
  return(invisible(results))
}

require(ShadowCAT)
models = c("GPCM", "GRM", "SM")
priors = list(matrix(c(1, .4,
                            .4, 1), 2, 2),
                   matrix(c(1, .8,
                            .8, 1), 2, 2))

theta_grid = list(NULL, c(-3,-3), c(-2,-2), c(-1,-1), c(0,0), c(1,1), c(2,2), c(3,3))
n = 100
K = 200
Q = 2
between = FALSE
for (model in models) {
    items <- createTestBank(model, K, Q, 7, between)
    test_segall <- initTest(items, estimator = "MAP", selection = "MI", objective = "PD")
    test_shadow <- initTest(items, estimator = "MAP", selection = "Shadow", objective = "PEKL")
    for (j in 1:length(theta_grid)) {
      for (prior in priors){
      simCAT(n, items = items, test = test_segall, fixed_theta = theta_grid[[j]], theta0 = list(mean = c(0,0), covar = prior), store = TRUE)
      #simCAT(n, items = items, test = test_shadow, fixed_theta = theta_grid[[j]], theta0 = list(mean = c(0,0), prior = prior), store = TRUE)
    }
  }
}

# fetch all results
# result_files <- list.files(path = paste0(getwd(),"/simulations"), pattern = '!RESULTS', all.files = FALSE,
#                            full.names = TRUE, recursive = TRUE,
#                            ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
# 
# result_files