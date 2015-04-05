############Simulations

# Simulate wrapper
simCAT <- function(n = 10,
                   store = FALSE,
                   location = getwd(),
                   items = createTestBank("GPCM"),
                   test = initTest(items),
                   theta0 = list(mean = rep(0,items$Q),
                                 covar = diag(items$Q)),
                   prior = theta0)
                   
{
  if (store) { 
    dir <- paste0(location,"/simulations/",paste(items$model, test$selection, test$estimator, sep="-"),"/",make.names(Sys.time()))
    cat("Storing in", dir, "...\n\n")
    dir.create(dir, TRUE, TRUE)
  }
  
  results <- list()
  thetas <- estimates <- matrix(NA, n, test$items$Q)
  
  start.time <- Sys.time()
  
  for (i in 1:n) {
    # draw a theta
    theta = rmvnorm(1, theta0$mean, theta0$covar)

    # run the cat
    results[[i]] <- person <- ShadowCAT(initPerson(items, theta = theta, prior = prior$covar), test)
    thetas[i,] <- theta
    estimates[i,] <- person$estimate
    
    if (store) save(person, file = paste0(dir, "/", i, '.CAT'))
    
    # some feedback, time
    if (i > 2) {
      time.spent <- difftime(Sys.time(), start.time, units = 'secs')
      time.left <- (time.spent / i) * (n - i)
      time.string <- format(.POSIXct(time.left, tz = "GMT"), "%H:%M:%S")
    } else {
      time.string <- 'unknown'
    }
    cat("\r", i, "of", n, ", deviation: ", round(person$estimate - theta, 2), ", estimated time left: ", format(time.string, format = "%T"), "                 ")
    
  }
  
  deviations <- estimates - thetas
  
  cat("\n\nBias;\n")
  print(bias <- colMeans(deviations))
  
  cat("\nMSE;\n")
  print(rmse <- colMeans(deviations**2))
  
  if (store) save(thetas, estimates, deviations, bias, rmse, test, results, file = paste0(dir, "/!RESULTS!"))
  return(invisible(results))
}

n <- 1000
require(ShadowCAT)
items <- createTestBank("3PLM", K = 100, Q = 2, between = FALSE)
test <- initTest(items, estimator = "ML", selection = "MI", objective = "PD")
simCAT(n, items = items, test = test, store = TRUE)