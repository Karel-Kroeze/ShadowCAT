# package
require(ShadowCAT)

# Load in the example data.
data(GPCM)

# create itembank
items <- initItembank("GPCM", alpha = GPCM[[1]], beta = GPCM[[2]])
person <- initPerson(items = items, theta = c(.5,.5,.5), prior = GPCM[[3]])
test <- initTest(items, estimator = "MAP", objective = "PFI", start = list(type = 'fixed', indeces = c(1,2,15,16,182,183)))

person <- answer(person, test, c(1,2,15,16,182,183))
person$responses <- rep(c(1,0), 3)
estimate(person, test)$estimate

x <- FI(test, person)
y <- apply(x[,,person$administered],c(1,2),sum)
det(y)
z <- apply(x[,,person$available],3,function(x) det(x + y))
z
which(z == max(z))


test$estimator <- "EAP"
estimate(person, test)$estimate
adaptIntegrate(ShadowCAT:::LL, test$lowerBound, test$upperBound, fDim = 1, test = test, person = person)

for (est in c("EAP", "MAP", "ML")) {
  test$estimator <- est
  print(estimate(person, test)$estimate)
}

