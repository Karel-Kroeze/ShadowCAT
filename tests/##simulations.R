############Simulations

# Simulate wrapper
simCAT <- function(n = 10,
                   store = FALSE,
                   location = paste0(getwd(), '/simulations'),
                   test = initTest(createTestBank("GPCM")),
                   fixed_theta = NULL,
                   theta0 = list(mean = rep(0, items$Q),
                                 covar = diag(items$Q)),
                   prior = theta0,
                   exposure.control = FALSE,
                   constrained = FALSE)
  
{
  cat("\n", test$items$model, test$selection, test$estimator, fixed_theta)
  if (store) { 
    dir <- paste0(location,"/", paste(test$items$model, test$selection, test$estimator, sep="-"), "/", make.names(Sys.time()))
    cat(" [ Storing in", location, "]")
    dir.create(dir, TRUE, TRUE)
  }
  cat("\n")
  
  # set some timing stuff on the redis host
  redisSet('start', Sys.time())
  redisDelete('done')
  redisSet('n', n)
  
  # give this job a unique name
  jobname <- paste(test$items$model,
                   test$selection,
                   theta0$covar[1,2],
                   ifelse(is.null(fixed_theta), 'NORM', fixed_theta[1]),
                   ifelse(constrained, "constrained", "unconstrained"),
                   ifelse(exposure.control, "exposure controlled", "not exposure controlled"), sep="_")
  redisSet('jobname', jobname)
  
  if(exposure.control) {
    # start/reset exposure administrative variables on the redis server.
    redisMulti()
    raw0 <- charToRaw('0')
    redisSet('exposure:total', raw0)
    redisSet('exposure:feasible', raw0)
    redisHMSet('exposure:exposure', setNames(as.list(rep(raw0,test$items$K)), as.character(1:test$items$K)))
    redisHMSet('exposure:eligible', setNames(as.list(rep(raw0,test$items$K)), as.character(1:test$items$K)))
    redisExec()
  }
  
  # doRedis to the rescue!
  results <- foreach(icount(n), .inorder = FALSE, .packages = c('mvtnorm','MultiGHQuad','ShadowCAT'), .verbose = TRUE) %do% {
    # draw a theta, if not specified (draws from prior, which means the prior is always correct - slightly cheating) (not used in reported simulations)
    if (is.null(fixed_theta)) {
      theta <- rmvnorm(1, theta0$mean, theta0$covar)
    } else { 
      theta <- fixed_theta
    }
    
    # exposure control constraints
    if(exposure.control){ 
      total <- as.numeric(redisGet("exposure:total"))
      feasible <- as.numeric(redisGet("exposure:feasible"))
      eligible <- as.numeric(unlist(redisHGetAll("exposure:eligible")))
      exposure <- as.numeric(unlist(redisHGetAll("exposure:exposure")))
      test <- createExposureConstraint(test, exposure, eligible, feasible, total, .25)
    }
    
    # run the cat
    person <- ShadowCAT(initPerson(test$items, theta = theta, prior = prior$covar), test, FALSE)
    
    # exposure control administration
    if(exposure.control){
      redisMulti()
      redisIncr('exposure:total')
      if(test$exposure$feasible) redisIncr('exposure:feasible')
      for (i in (1:test$items$K)[test$exposure$eligible]) redisHIncrBy('exposure:eligible', as.character(i), 1)
      for (i in person$administered) redisHIncrBy('exposure:exposure', as.character(i), 1)
      redisExec()
    }
    
    # increment redis counter
    redisIncr('done')
    redisSet('update', Sys.time())
    
    # return person object
    person
  }
  
  # set up empty containers
  deviations <- estimates <- thetas <- variances <- matrix(NA, n, test$items$Q)
  
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
  
  # create a summary data-frame
  summary <- data.frame(bias = bias, RMSE = rmse, mean_var = mean_var, SE_var = SE_var)
  cat('\n')
  print(round(summary,2))
  
  # store it on the redis server, as well as on disk.
  redisHSet('summary', jobname, summary)
  # redisHSet('results', jobname, results)
  if (store) save(results, test, fixed_theta, theta0, prior, file = paste0(dir, "/!RESULTS.Rdata"))
  if (store) save(summary, file = paste0(dir, "/!SUMMARY.Rdata"))
  return(invisible(results))
}




require(ShadowCAT)
models = c(
    "3PLM",
     "GPCM", 
     "GRM",
  "SM")

priors = list(
    matrix(c(1, .4,
             .4, 1), 2, 2),
  matrix(c(1, .8,
           .8, 1), 2, 2)
)

theta_grid = list(
    c(-3,-3), 
    c(-2,-2),
    c(-1,-1),
    c(0,0),
    c(1,1),
    c(2,2),
  c(3,3)
)

n = 250
n_norm = 1000
K = 200
Q = 2
N = 30 # test length

# stop rule
stoprule <- list(type = 'length', n = N)

# constraints
# chisq distributed time for all items.
time <- pmax( rnorm(K, 300, 90), 60)

# three sets of overlapping content specs;
# ..|......
# ....|....
# ......|..
contentA <- factor(rep(c('alpha', 'beta'), times = c(floor(K/3), K - floor(K/3))))   
contentB <- factor(rep(c('rincewind', 'twoflower'), times = c(floor(K/2), K - floor(K/2))))   
contentC <- factor(rep(c('charly', 'wonka'), times = c(K - floor(K/3), floor(K/3))))   

# a bunch of enemies
enemies <- setNames(lapply(1:10, function (x) {
  x <- rep(0, K)
  x[sample(K, 2)] <- 1
  x
}), as.character(1:10))

characteristics <- data.frame(time, contentA, contentB, contentC, enemies = enemies)
constraints <- list(
  list('time', '><', c(mean(time)*.5*N, mean(time)*1.2*N)),
  list('contentA/alpha',     '><', floor(c(N/3-3, N/3+3))),
  list('contentA/beta',      '><', floor(c(N/3*2-3, N/3*2+3))),
  list('contentB/rincewind', '><', floor(c(N/2-3, N/2+3))),
  list('contentB/twoflower', '><', floor(c(N/2-3, N/2+3))),
  list('contentC/charly',    '><', floor(c(N/3*2-3, N/3*2+3))),
  list('contentC/wonka',     '><', floor(c(N/3-3, N/3+3))),
  list('enemies.1', '<', 1),
  list('enemies.2', '<', 1),
  list('enemies.3', '<', 1),
  list('enemies.4', '<', 1),
  list('enemies.5', '<', 1),
  list('enemies.6', '<', 1),
  list('enemies.7', '<', 1),
  list('enemies.8', '<', 1),
  list('enemies.9', '<', 1),
  list('enemies.10', '<', 1)
  )

require(rredis)
require(doRedis)
require(foreach)

redis.host = '128.199.63.229'
startLocalWorkers(n = 1, queue = "simulations", redis.host)
registerDoRedis("test", redis.host)

for (model in models) {
  items <- createTestBank(model, K, Q, 5)
  test_segall <- initTest(items, estimator = "MAP", selection = "MI", objective = "PD", stop = stoprule)
  test_shadow <- initTest(items, estimator = "MAP", selection = "Shadow", objective = "PEKL", stop = stoprule)
  test_shadow_content <- initTest(items, estimator = "MAP", selection = "Shadow", objective = "PEKL",
                                  stop = stoprule, constraints = list(constraints = constraints, characteristics = characteristics))  
  
  for (prior in priors){
    # normal distributed
    simCAT(n_norm, test = test_segall, theta0 = list(mean = c(0,0), covar = prior),
           store = TRUE, location = paste0(getwd(), '/simulations/segall'))
    simCAT(n_norm, test = test_shadow, theta0 = list(mean = c(0,0), covar = prior),
           store = TRUE, location = paste0(getwd(), '/simulations/shadow-unconstrained'))
    simCAT(n_norm, test = test_shadow_content, theta0 = list(mean = c(0,0), covar = prior), constrained = TRUE,
           store = TRUE, location = paste0(getwd(), '/simulations/shadow-constrained'))
    simCAT(n_norm, test = test_shadow_content, theta0 = list(mean = c(0,0), covar = prior), constrained = TRUE,
           exposure.control = TRUE, store = TRUE, location = paste0(getwd(), '/simulations/shadow-constrained-exposure'))
    
    
    for (j in 1:length(theta_grid)) {
      # theta grid.
      simCAT(n, test = test_segall, fixed_theta = theta_grid[[j]], theta0 = list(mean = c(0,0), covar = prior),
             store = TRUE, location = paste0(getwd(), '/simulations/segall'))
      simCAT(n, test = test_shadow, fixed_theta = theta_grid[[j]], theta0 = list(mean = c(0,0), covar = prior),
             store = TRUE, location = paste0(getwd(), '/simulations/shadow-unconstrained'))
      simCAT(n, test = test_shadow_content, fixed_theta = theta_grid[[j]], theta0 = list(mean = c(0,0), covar = prior), constrained = TRUE,
             store = TRUE, location = paste0(getwd(), '/simulations/shadow-constrained'))
      simCAT(n, test = test_shadow_content, fixed_theta = theta_grid[[j]], theta0 = list(mean = c(0,0), covar = prior), constrained = TRUE,
             exposure.control = TRUE, store = TRUE, location = paste0(getwd(), '/simulations/shadow-constrained-exposure'))
    }
  }
}

removeQueue("jobs")
