#' Comprehensive testing wrapper.
#' 
#' Test a CAT run under various conditions.
#' 
#' @param models
#' @param estimators
#' @param dims
#' @param size
#' @param length
testCAT <- function(models = c('3PLM','GRM','GPCM','SM'),
                    estimators = c('ML','MAP','EAP'),
                    objectives = c('PEKL','TPFI','DPFI','TFI','DFI'),
                    selectors = c('Shadow','MI'),
                    dims = c(1, 2, 3), between = c(TRUE, FALSE),
                    covars = NULL,
                    K = 100,
                    start = list(type = 'random', n = 5),
                    stop = list(type = 'length', n = 30)) {
  
  total <- length(models) * length(dims) * length(between) * length(estimators) * length(objectives) * length(selectors)
  i <- 0
  cat("Trying", total, "different permutations of options.\n")
  
  # progress bar
  #pb <- txtProgressBar(0, total, style = 3)
  
  for (model in models) {
    for (Q in dims) {
      for (within in between) {
        
        # assumes Q is sequential from 1:X    
        items <- createTestBank(model, K, Q, M = 4, between = within)
        
        for (estimator in estimators) {
          for (objective in objectives) {
            for (selection in selectors) {
              # output
              i <- i + 1
              if (within) { bet <- 'between' } else { bet <- 'within' }
              cat("\r", i, model, Q, 'dimensions', bet, estimator, objective, selection, "                        ")
              
    
              person <- initPerson(items, rnorm(Q), diag(Q))
              test <- initTest(items,start,stop,estimator,objective,selection)
              test <- createConstraints(test)
              
              # go on when fail
              tryCatch({
                person <- ShadowCAT(person, test)
              }, error = function(e) {
                cat("\n\n", model, Q, 'dimensions', bet, estimator, objective, selection)
                cat("\n FAILED: ", e$message)
                cat("\n IN: ", as.character(e$call), "items\n\n\n")
              }
              , warning = function(w) {
                cat("\n\n", model, Q, 'dimensions', bet, estimator, objective, selection)
                cat("\n WARNING: ", w$message)
                cat("\n IN: ", as.character(w$call), "items\n\n\n")
              }
              )
            }
          }
        }
      }
    }
  }  
}

#' Test Estimators
#' 
#' Test the various estimators
#' 
#' @param models
#' @param estimators
#' @param dims
#' @param between
#' @param theta
#' @param prior
#' @param N
#' @param K
#' @param detail
#' @return NULL
#' @export
testEstimators <- function(models = c("3PLM","GRM","GPCM","SM"),
                           estimators = c("ML", "MAP", "EAP"),
                           dims = 1:3,
                           between.dim = c(TRUE, FALSE),
                           theta = lapply(dims, function(x) rep(0,x)),
                           prior = lapply(dims, function(x) diag(x)), 
                           N = 100, 
                           K = 20,
                           detail = FALSE, silent = FALSE) {
  
  
  require(mvtnorm)
  
  results <- list()
  for (q in seq_along(dims)){
    # set dimension, draw thetas
    Q <- dims[q]
    thetas <- rmvnorm(N, theta[[q]], prior[[q]])
    
    results[[as.character(Q)]] <- list()
    
    for (between in between.dim){
      results[[as.character(Q)]][[as.character(between)]] <- list()
      for (model in models) {
        results[[as.character(Q)]][[as.character(between)]][[model]] <- list()
        
        # create itembank
        items <- createTestBank(model, K, Q, 4, between)
        
        for (estimator in estimators) {     
          
          # create test, create persons, answer questions
          test <- initTest(items, estimator = estimator)
          persons <- apply(thetas, 1, function(x) {
            answer(initPerson(items, x, prior[[q]]), test, 1:items$K)
          })
          
          # estimate latent traits
          persons <- lapply(persons, function(x) {
            estimate(x, test)
          })
          
          # store results
          results[[as.character(Q)]][[as.character(between)]][[model]][[estimator]] <- persons
          
          # summary statistics
          deviation <- sapply(persons, function(x) {
            x$estimate - x$theta
          })
          
          # cast into matrix
          deviation <- matrix(deviation, ncol = Q)
          
          # bias
          bias <- colMeans(deviation)
          
          # RMSE
          RMSE <- sqrt(colMeans(deviation**2))
          
          # output
          if (! silent){
            between.text <- "within"
            if (between) between.text <- "between"
            cat(Q,'dimensional',between.text,model,K,'items, for',N,'respondents,',estimator,'estimates.\n')
            cat('bias\n')
            print(bias)
            cat('RMSE\n')
            print(RMSE)
            cat('max deviation\n')
            print(max(abs(deviation)))
            cat('\n')
          }
        }
      }
    }
  }
  return(invisible(results))
}

#' Create a test Itembank
#' 
#' Quick and simple itembanks for testing purposes.
#' @param model
#' @param K
#' @param Q
#' @param M
#' @return ShadowCAT.itembank
#' @export
createTestBank <- function(model, K = 50, Q = 1, M = 4, between = FALSE){
  # 3PLM is dichotomous by definition
  if (model == "3PLM") M <- 1 
  
  # make sure the number of items is divisible by the number of dimensions
  if (between) K <- ceiling(K/Q) * Q 
  
  # set up alpha, very rough uniform from .3 to 1.5
  alpha <- matrix(runif(K * Q, .3, 1.5), K, Q)
  
  # if between, force items to load on one dimension each.
  if (between){
    set = K / Q
    for (i in 1:Q){
      alpha[((i-1)*set+1):(i*set), (1:Q)[-i]] <- 0    
    }
  }
  
  # spread polytomous items cats -2 to +2.
  spread <- seq(-2,2,length.out=M)
  
  # base loading for items
  beta <- matrix(rnorm(K), K, 1)
  
  # apply spread for polytomous, betas are strictly monotously increasing because the spread is.
  # apply transposes the matrix...
  if (M > 1) beta <- t(apply(beta, 1, function(x) x + spread))
  
  # reparameterize GPCM
  if (model == "GPCM") {
    # make betas
    eta <- beta
    for (i in 1:M) {
      # rolling sum, apply over all items
      beta[,i] <- apply(eta[,1:i, drop=FALSE], 1, sum)
    }
  }
  
  # create Itembank object
  items <- initItembank(model, alpha, beta, silent = TRUE)  
  
  # return
  return(invisible(items))
}


#' Test Information Measures
#' 
#' @param models
#' @param dims
#' @param methods
#' @param K
#' @export
testFI <- function(models = c("3PLM", "GPCM", "SM", "GRM"), dims = 1:3, methods = c("FI", "PFI"), K = 1
                   ,show.items = FALSE){
  for (Q in dims) {
    for (model in models) {
      # create item(s) and person
      items <- createTestBank(model, K, Q)
      person <- initPerson(items, prior = diag(Q))
      
      # item feedback
      if (show.items) {
        cat("item\n")
        print(items$pars)
      }
      
      for (method in methods) {
        
        # init test
        test <- initTest(items, objective = method)
        
        # get information
        info <- FI(test, person)
        
        # sum over items
        total <- apply(info, c(1,2), sum)
        
        # apply prior
        if (method == "PFI") total <- total + solve(person$prior)
        
        # feedback
        cat("\n+++++",method,"for",Q,"dimensional",model,"(theta = 0 vector).\n")
        cat("info\n")
        print(total)
        tryCatch({ # suppress singularity errors
          covar <- solve(total)
          cat("(posterior) expected variance\n")
          print(diag(covar))
          cat("det (size of (posterior) credibility region)\n")
          print(det(covar))
          cat("\n")
        }, error = function(e) {
          cat(e$message)
          if (grepl("singular", e$message)) cat("\nThis is likely because too few items were used.")
        })
      }
    }  
  }
}

#' Test probability functions
#' 
#' @param models
#' @param alpha
#' @param beta
#' @param Q
#' @importFrom plot3D persp3D
#' @export
testProb <- function(models=c("3PLM","GPCM","SM","GRM"), Q=1, alpha=1, beta=0){
  
  # theta range
  theta <- seq(-4,4,length.out = 100)
  
  # set plot pars.
  cols = 1 + (length(models) > 1)
  rows = ceiling(length(models) / cols)
  par(mfrow = c(rows, cols))
  
  if (Q == 1){ # one dimension
    
    for (model in models) {
      # item
      BETA <- matrix(beta, nrow = 1)
      
      # 3PLM is dichotomous -> first beta only
      if (model == "3PLM") BETA <- BETA[1,1]
      
      # reparameterize GPCM
      if (model == "GPCM") {
        # make betas
        eta <- BETA
        for (i in 1:length(eta)) {
          # rolling sum, apply over all items
          BETA[,i] <- apply(eta[,1:i, drop=FALSE], 1, sum)
        }
      }
      
      item <- initItembank(model, alpha, BETA, silent=T)
      test <- initTest(item)
      
      # get probs
      p <- matrix(0,100,item$M+1)
      for (i in seq_along(theta)){
        p[i,] <- prob(test,theta=theta[i])$P
      }
      
      # plot it!
      matplot(theta,p,type='l',main=model)
    }
  } else {
    Q <- 2 # higher dimensions are fairly pointless for plotting
    theta <- seq(-4, 4, length.out = 25) # more respectable size
    theta.grid <- as.matrix(expand.grid(theta, theta))
    
    # make alpha 2 dim
    if (length(alpha) != 2) alpha  <- matrix(1, nrow = 1, ncol = 2)
    
    # do the work!    
    for (model in models) {
      # item
      BETA <- matrix(beta, 1)
      
      # 3PLM is dichotomous -> first beta only
      if (model == "3PLM") BETA <- BETA[1,1]
      
      # reparameterize GPCM
      if (model == "GPCM") {
        # make betas
        eta <- BETA
        for (i in 1:length(eta)) {
          # rolling sum, apply over all items
          BETA[,i] <- apply(eta[,1:i, drop=FALSE], 1, sum)
        }
      }
      
      item <- initItembank(model, alpha, BETA, silent=T)
      test <- initTest(item)
      
      # get probs
      p <- matrix(0,nrow(theta.grid),item$M+1)
      for (i in 1:nrow(theta.grid)){
        p[i,] <- prob(test,theta=theta.grid[i,])$P
      }
      
      # plot it!
      persp3D(theta, theta, matrix(p[,1], nrow = length(theta)),theta = 45, phi = 45,xlab='theta 1',ylab='theta 2',zlab='P',facets=T,col='#FF3434',border='black',colkey=F,contour=F,zlim=c(0,1),main=model)
      # other cats.
      for (i in (1:length(BETA))+1){
        persp3D(theta, theta, matrix(p[,i], nrow = length(theta)),theta = 45, phi = 45,xlab='theta 1',ylab='theta 2',zlab='P',facets=T,col='#FF9034',border='black',colkey=F,contour=F,zlim=c(0,1),add=T)
      }
    }
  }
  par(mfrow = c(1,1))
}

#' Test analytical derivatives
#' 
#' Compare analytical derivatives with their numerical counterparts
#' @param models
#' @param dims
#' @param eval.points
#' @param K
#' @param estimator
#' @importFrom numDeriv jacobian hessian
#' @export
testDeriv <- function(models = c("3PLM", "GPCM", "GRM", "SM"), dims = 1:3, eval.points = c(-2,0,2), K = 10, estimator = "ML") {
  for (Q in dims){
    for (model in models) {
      # create items, test and person. Obtain responses pattern.
      items <- createTestBank(model, K, Q)
      test <- initTest(items, estimator = estimator)
      person <- answer(initPerson(items), test, 1:K)
      
      # for each evaluation point
      for (eval.point in eval.points) {
        eval.point <- rep(eval.point, Q)
        PROB <- prob(test, person, eval.point, TRUE)
        
        # output
        cat("\n=================================================\n",
            "derivatives for",Q,"dim",model,"at theta (vector)",eval.point,"\n")
        
        analytical <- list(jacobian = PROB$d1, hessian = PROB$d2)
        numerical <- list(jacobian = jacobian(ShadowCAT:::LL,eval.point,"complex",test = test, person = person),
                          hessian = hessian(ShadowCAT:::LL,eval.point,"complex",test = test, person = person))
        
        if(isTRUE(all.equal(analytical, numerical,tolerance = 1e-3))) {
          cat("match\n")
        } else {
          print(all.equal(analytical, numerical,tolerance = 1e-3))
        }
      } 
    }
  }
}

#' Plot derivatives
#' 
#' @param models
#' @param d
#' @param estiamtor
#' @param theta
#' @param K
#' @param M
#' @export
plotDeriv <- function(models = c("3PLM", "GRM", "SM", "GPCM"),
                      d=1, estimator = 'ML', prior = NULL,
                      theta = seq(-6,6, length.out = 100),
                      K = 10, M = 4) {
  # set plot pars.
  cols = 1 + (length(models) > 1)
  rows = ceiling(length(models) / cols)
  par(mfrow = c(rows, cols))
  
  for (model in models){
    # init items, test and person (with responses)
    items <- createTestBank(model,K=K,Q=1,M=M)
    test <- initTest(items, estimator = estimator)
    person <- answer(initPerson(items, 0, prior = prior), test, 1:K)
    
    # get values from numerical and analytical, as well as LL values.
    y <- matrix(0,length(theta), 3)
    for (i in seq_along(theta)){
      PROB <- prob(test, person, theta = theta[i], deriv = TRUE)
      if (d == 1){
        numeric <- jacobian(LL, theta[i], test = test, person = person)
        y[i,] <- c(PROB$LL, PROB$d1, numeric)
      } else {
        numeric <- hessian(LL, theta[i], test = test, person = person)
        y[i,] <- c(PROB$LL, PROB$d2, numeric)
      }
    }
    
    # plot it!
    matplot(theta, y, type='l', col=c('black','red','black'), lty=c(1,2,2), main = model, xlab = paste0('Analytical = red (', paste(person$responses,collapse=','), ")", collapse=""))
    abline(h = 0, col='grey')
    
  }
  par(mfrow = c(1,1))
}

testKL <- function(model = "3PLM", Q = 1) {
  # set up objects
  items <- createTestBank(model, Q = Q)
  person <- initPerson(items)
  test <- initTest(items, estimator = "EAP")
  
  # get some answers
  person <- answer(person, test, sample(items$K, 10))
  
  # special case GPCM
  if (model == "GPCM"){
    data <- GPCM
    items <- initItembank(model = "GPCM", alpha = data[[1]], beta = data[[2]])
    person <- initPerson(items, prior = data[[3]])
    test <- initTest(items, estimator = "EAP")
    person <- answer(person, test, c(1,2,20,21,190,191))
  }
  
  # get fisher information
  FI_all <- FI(test, person)
  so_far <- apply(FI_all[,,person$administered, drop = FALSE], c(1,2), sum)
  FI <- apply(FI_all[,,person$available, drop = FALSE], 3, function(x) det(so_far + x))
  
  # get KL information
  KL <- PEKL(test, person)
  
  # print simplified results
  cat("\n", cor(FI, KL))
  return(invisible(list(FI = FI, KL = KL)))
}
