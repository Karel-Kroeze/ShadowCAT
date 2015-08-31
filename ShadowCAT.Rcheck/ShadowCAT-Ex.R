pkgname <- "ShadowCAT"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('ShadowCAT')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("answer")
### * answer

flush(stderr()); flush(stdout())

### Name: answer
### Title: Simulate response(s)
### Aliases: answer

### ** Examples

items <- createTestBank("GPCM")
test <- initTest(items)
person <- initPerson(items)

# simulates responses to all questions, and returns a vector response pattern.
answer(person, test)

# simulates responses to the specified question indeces, and returns an updated person object.
answer(person, test, sample(test$items$K, 5))



cleanEx()
nameEx("createConstraints")
### * createConstraints

flush(stderr()); flush(stdout())

### Name: createConstraints
### Title: Creates a constraints object.
### Aliases: createConstraints

### ** Examples

# set up a simple itembank and test.
items <- createTestBank("GPCM")
test <- initTest(items, selection = "Shadow" , objective = "PEKL")

# set up some dummy characteristics.
content <- sample(c('algebra','physics','calculus'), items$K, TRUE)
time <- rnorm(items$K)
exclusive <- rep(0, items$K)
exclusive[sample(items$K, 4)] <- 1

# bind them in a data.fame
characteristics <- data.frame(content, time, exclusive)

# set up the constraints
constraints <- list(
  list(name = 'content/algebra',
       op = '><',
       target = c(5,10)),
  list(name = 'content/physics',
       op = '><',
       target = c(2,5)),
  list(name = 'time',
       op = '<',
       target = 20),
  list(name = 'exclusive',
       op = '<',
       target = 2))

# update the test object
test$constraints <- createConstraints(test, characteristics, constraints)

# or do it all at once;
test2 <- initTest(items, constraints = list(characteristics = characteristics, constraints = constraints))

# results are identical (initTest uses createConstraints internally);
all.equal(test$constraints, test2$constraints)



cleanEx()
nameEx("estimate")
### * estimate

flush(stderr()); flush(stdout())

### Name: estimate
### Title: Latent trait estimation
### Aliases: estimate

### ** Examples

# create a basic test + person
items <- createTestBank("GPCM")
test <- initTest(items, estimator = "ML")
person <- initPerson(items)

# answer a few items
person <- answer(person, test, sample(items$K, 10))

# obtain estimates
ML <- estimate(person, test)$estimate
test$estimator <- "MAP"
MAP <- estimate(person, test)$estimate
test$estimator <- "EAP"
EAP <- estimate(person, test)$estimate
ML; MAP; EAP

# access variance
attr(ML, "variance")

# Note that EAP takes considerably more time when dimensionality is higher...
items5d <- createTestBank("GPCM", Q=5)
test5dEAP <- initTest(items, estimator = "EAP")
test5dMAP <- initTest(items, estimator = "MAP")
person5d <- answer(initPerson(items), test, sample(items5d$K, 10))

system.time(estimate(person5d, test5dEAP))
system.time(estimate(person5d, test5dMAP))



### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
