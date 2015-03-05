#' Compare ShadowCAT estimates to mirt estimates
#' @param models
#' @param estimators
#' @param dims
#' @param K
#' @param N
#' @importFrom mvtnorm rmvnorm
#' @export
compareEstimators <- function(models = c("3PLM", "GRM", "GPCM"), estimators = c("ML", "MAP", "EAP"), dims = c(1,2,3), K = 20, N = 1000) {
  # loop over models
  for (model in models){
    # loop over dims
    for (Q in dims){
      # set M to 1 for 3PLM, 4 otherwise
      M <- switch(model, "3PLM" = 1, 4)
      
      # create some data
      items <- createTestBank(model, K = K, Q = Q)
      test <- initTest(items)
      
      # create parameter data.frame
      df <- with(items$pars, data.frame(a = alpha, d = -beta))
      colnames(df) <- c(paste0('a',1:Q), paste0('d',1:M))
      
      # 3PL is special...
      if (model == "3PLM") {
        colnames(df) <- c(paste0('a',1:Q), 'd')
        df$d = df$d * rowSums(items$pars$alpha)
      }
      
      # load into mirt
      mod <- generate.mirt_object(df, itemtype = switch(model, "3PLM" = "4PL", "GRM" = "graded", "GPCM" = "gpcm"))
      
      # get a new response pattern
      person <- answer(person = initPerson(items, rmvnorm(1,rep(0,Q),diag(Q))),test,1:items$K)
      
      for (estimator in estimators){      
        # get estimates
        mirt.est <- as.numeric(fscores(mod, method = estimator, response.pattern = person$responses, scores.only = TRUE)[,paste0('F',1:Q)])
        test$estimator <- estimator
        shadow.est <- estimate(person, test)$estimate
        
        # output
        cat(model,Q,estimator,"\n")
        print(mirt.est)
        print(shadow.est)
      }      
    }
  }
}