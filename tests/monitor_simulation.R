# This rredis package keeps giving erros, will fix this later
if (FALSE) {
require(rredis)

# connect
host <- '128.199.63.229'
redisConnect(host)

# interval
interval = 5

# for plot
x <- numeric(0)

while(TRUE) {
  job <- redisGet('jobname')
  cat("\n\nCurrent job:", job, "\n")
  start <- max(as.numeric(redisGet('done')), 0) # we may drop in halfway thru
  n <- as.numeric(redisGet('n'))
  start_time <- redisGet('start')
  last <- start
  
  Sys.sleep(interval)
    
  while(redisGet('jobname') == job) {
    done <- as.numeric(redisGet('done'))
    
    # plot
    x <- c(x, max(done-last, 0))
    last <- done
    xlab <- paste("Last ", round( length(x) * interval / 60, 1), " minutes.")
    plot(x*60/interval, type='l', col ='steelblue', ylim = c(0,max(x*60/interval)), ylab = 'rate p/m', xlab = xlab, xaxt='n')
    
    if (length(x) > 100) x <- x[(length(x)-100):length(x)]
    
    # cat
    cur_time <- Sys.time()
    time.spent <- difftime(cur_time, start_time, tz = "GMT", units = "sec")
    time.per.unit <- time.spent / done
    time.left <- time.per.unit * (n - done)
    
    #cat(time.per.unit)
    time.spent <- format(.POSIXct(time.spent, tz = "GMT"), "%H:%M:%S")
    time.left <- format(.POSIXct(time.left, tz = "GMT"), "%H:%M:%S")
    
    cat("\rSimulation ", done, " of ", n, ".", " Time taken: ",
        format(time.spent, format = "%T"), ", time left: ",
        format(time.left, format = "%T"), ", rate: ",
        round(60 / as.numeric(time.per.unit), 0), " p/m.                  ",
        sep='')
    Sys.sleep(interval)
  }
  Sys.sleep(min(5 * interval, 10))
}

}
