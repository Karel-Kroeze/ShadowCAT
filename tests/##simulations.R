############Simulations

# Simulate wrapper
simCAT <- function(n = 10,
                   store = FALSE,
                   location = paste0(getwd(), '/simulations'),
                   items = createTestBank("GPCM"),
                   test = initTest(items),
                   fixed_theta = NULL,
                   theta0 = list(mean = rep(0, items$Q),
                                 covar = diag(items$Q)),
                   prior = theta0)
  
{
  cat("\n", items$model, test$selection, test$estimator, fixed_theta)
  if (store) { 
    dir <- paste0(location,"/", paste(items$model, test$selection, test$estimator, sep="-"), "/", make.names(Sys.time()))
    cat(" [ Storing in", location, "]")
    dir.create(dir, TRUE, TRUE)
  }
  cat("\n")
  
  # set some timing stuff on the redis host
  redisSet('start', Sys.time())
  redisSet('jobname', paste(test$selection, 'prior', theta0$covar[1,2], ifelse(is.null(fixed_theta), 'NORM', fixed_theta[1]), sep="-"))
  redisDelete('done')
  redisSet('n', n)
  
  # doRedis to the rescue!
  results <- foreach(icount(n), .inorder = FALSE, .packages = c('mvtnorm','MultiGHQuad','ShadowCAT'), .verbose = FALSE) %dopar% {
    # draw a theta, if not specified (draws from prior, which means the prior is always correct - slightly cheating) (not used in reported simulations)
    if (is.null(fixed_theta)) theta <- rmvnorm(1, theta0$mean, theta0$covar)
    else theta <- fixed_theta
    
    # run the cat
    person <- ShadowCAT(initPerson(items, theta = theta, prior = prior$covar), test, FALSE)
    
    # increment redis counter
    redisIncr("done")
    
    person
  }
  
  # set up empty containers
  deviations <- estimates <- thetas <- variances <- matrix(NA, n, items$Q)
  
  # fill from results
  for (i in seq_along(results)) {
    estimates[i,] <- results[[i]]$estimate
    thetas[i,] <- results[[i]]$theta
    variances[i,] <- diag(attr(results[[i]]$estimate, 'variance'))
  }
  
  # compute the rest
  deviations <- estimates - thetas
  bias <- colMeans(deviations)
  rmse <- sqrt(colMeans(deviations**2))
  mean_var <- colMeans(variances)
  SE_var <- apply(variances, 2, sd)
  
  #
  summary <- data.frame(bias = bias, RMSE = rmse, mean_var = mean_var, SE_var = SE_var)
  cat('\n')
  print(round(summary,2))

  #
  if (store) save(results, test, fixed_theta, theta0, prior, file = paste0(dir, "/!RESULTS!"))
  if (store) save(summary, file = paste0(dir, "/!SUMMARY!"))
  return(invisible(results))
}




require(ShadowCAT)
models = c("GPCM" 
           ,"GRM"
           ,"SM"
           )
priors = list(matrix(c(1, .4,
                       .4, 1), 2, 2)
              ,matrix(c(1, .8,
                       .8, 1), 2, 2)
              )

theta_grid = list(
                  c(-3,-3)
                  , c(-2,-2)
                  , c(-1,-1)
                  , c(0,0)
                  , c(1,1)
                  , c(2,2)
                  , c(3,3)
                  )

n = 1
n_norm = 1
K = 200
Q = 2

require(rredis)
require(doRedis)
require(foreach)

redis.host = '128.199.63.229'
startLocalWorkers(n = 1, queue = "jobs", redis.host)
registerDoRedis("jobs", redis.host)

for (model in models) {
  items <- createTestBank(model, K, Q, 7)
  test_segall <- initTest(items, estimator = "MAP", selection = "MI", objective = "PD")
  test_shadow <- initTest(items, estimator = "MAP", selection = "Shadow", objective = "PEKL")
  
  for (prior in priors){
    simCAT(n_norm, items = items, test = test_segall, theta0 = list(mean = c(0,0), covar = prior), store = TRUE, location = paste0(getwd(), '/simulations/unconstrained'))
    simCAT(n_norm, items = items, test = test_shadow, theta0 = list(mean = c(0,0), covar = prior), store = TRUE, location = paste0(getwd(), '/simulations/unconstrained'))
    for (j in 1:length(theta_grid)) {
      simCAT(n, items = items, test = test_segall, fixed_theta = theta_grid[[j]], theta0 = list(mean = c(0,0), covar = prior), store = TRUE, location = paste0(getwd(), '/simulations/unconstrained'))
      simCAT(n, items = items, test = test_shadow, fixed_theta = theta_grid[[j]], theta0 = list(mean = c(0,0), covar = prior), store = TRUE, location = paste0(getwd(), '/simulations/unconstrained'))
    }
  }
}

removeQueue("jobs")
