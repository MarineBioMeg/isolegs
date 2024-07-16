####################
# Helper Functions #

# Unscaling Function
unscale <- function(z, center = attr(z, "scaled:center"), scale = attr(z, "scaled:scale")) {
  if(!is.null(scale))  z <- sweep(z, 2, scale, `*`)
  if(!is.null(center)) z <- sweep(z, 2, center, `+`)
  structure(z,
            "scaled:center"   = NULL,
            "scaled:scale"    = NULL,
            "unscaled:center" = center,
            "unscaled:scale"  = scale
  )
}


# Expected Fish Immigration #
Immigration <- function(pars, density) {
  # Expected number of immigrants for a given starting density
  a <- pars[1]
  b <- pars[2]
  lambda <- pars[3]
  out <- lambda + a * density * exp(-b * density)
  return(out)
}
# 
# LikelihoodSSE <- function(predicted, observed) {
#   # Sum of squared errors
#   # We'll use just the predicted immigrants
#   
#   
#   out <- sum((predicted - observed)^2)
#   return(out)
# }

LikelihoodSSE <- function(pars, density, immigrants) {
  # Expected number of immigrants for a given starting density
  a <- pars[1]
  b <- pars[2]
  lambda <- pars[3]
  predicted <- lambda + a * density * exp(-b * density)
  observed <- immigrants
  out <- sum((predicted - observed)^2)
  return(out)
}

WeightedLikelihoodSSE <- function(pars, densityChinook, densitySteelhead, immigrants) {
  # Expected number of immigrants for a given starting density, divided into con- and inter-specific
  a <- pars[1]
  b <- pars[2]
  lambda <- pars[3]
  eta0 <- pars[4]
  eta1 <- pars[5]
  eta2 <- pars[6]
  dc <- eta0 + eta1*densityChinook + eta2*densitySteelhead
  predicted <- lambda + a * dc * exp(-b * dc)
  observed <- immigrants
  out <- sum((predicted - observed)^2)
  return(out)
}

LikelihoodPois <- function(pars, density, immigrants) {
  # Expected number of immigrants for a given starting density
  a <- pars[1]
  b <- pars[2]
  lambda <- pars[3]
  predicted <- lambda + a * density * exp(-b * density)
  observed <- immigrants
  out <- -sum(dpois(observed, predicted, log=TRUE))
  return(out)
}

WeightedLikelihoodPois <- function(pars, densityChinook, densitySteelhead, immigrants) {
  # Expected number of immigrants for a given starting density, divided into con- and inter-specific
  a <- pars[1]
  b <- pars[2]
  lambda <- pars[3]
  eta0 <- pars[4]
  eta1 <- pars[5]
  eta2 <- pars[6]
  dc <- eta0 + eta1*densityChinook + eta2*densitySteelhead
  predicted <- lambda + a * dc * exp(-b * dc)
  observed <- immigrants
  out <- -sum(dpois(observed, predicted, log=TRUE))
  return(out)
}

############################
# Initialize, read data ####
############################
###################
# Fish Immigation Model

setwd("~/Dropbox/PostDoc/Fish/")
setwd("C:/Users/Spencer/Dropbox/PostDoc/Fish/")
dat <- read.csv("EarlyMarkRecapClean.csv")
dat$density <- dat$num.mark / dat$pool.area

pools <- unique(dat$pool)
dat$dc <- rep(NA, 138) # combined density of both species
dat$do <- rep(NA, 138) # density of opposing species

restored <- dat[dat$habitat=="restored",]
unrestored <- dat[dat$habitat=="unrestored",]

for (i in 1:68) {
  restored$dc[i] <- sum(restored[(restored$date==restored$date[i]&restored$pool==restored$pool[i]),]$num.mark) / restored$pool.area[i]
  restored$do[i] <-sum(restored[(restored$date==restored$date[i]&restored$pool==restored$pool[i]&restored$species!=restored$species[i]),]$num.mark) / restored$pool.area[i]
}

for (i in 1:70) {
  unrestored$dc[i] <- sum(unrestored[(unrestored$date==unrestored$date[i]&unrestored$pool==unrestored$pool[i]),]$num.mark) / unrestored$pool.area[i]
  unrestored$do[i] <- sum(unrestored[(unrestored$date==unrestored$date[i]&unrestored$pool==unrestored$pool[i]&unrestored$species!=unrestored$species[i]),]$num.mark) / unrestored$pool.area[i]
  
}
############################
# Bayesian Version (in JAGS)
############################
require(BEST)
require(rjags)
require(scatterplot3d)

params <- c("alpha0", "alpha1", "alpha2", "prediction")
params2 <- c("beta0", "beta1", "beta2", "eta0", "eta1", "eta2", "prediction")


jags.data <- list(N.cells = length(x), density = x[ ,1], immigrants = y)
jm <- jags.model("model.txt", data = jags.data, n.chains = 5, n.adapt = 1000)
update(jm, n.iter = 10000)
jm.sample <- jags.samples(jm, variable.names = params, n.iter = 1000, thin = 1)
plot(as.mcmc.list(jm.sample$beta0), main = "Beta_0")
plot(as.mcmc.list(jm.sample$beta1), main = "Beta_1")
plot(as.mcmc.list(jm.sample$beta2), main = "Beta_2")

predictions <- summary(as.mcmc.list(jm.sample$prediction))
#prds <- data.frame(sc50 = x[,1], predictions$quantiles)
prds <- data.frame(sc50 = x, predictions$quantiles)
prds <- prds[order(prds[, 1]), ]

plot(x, y, cex = 1, col = "darkgrey", pch = 19, ylab = "# of Immigrants", 
     xlab = "Chinook Initial Density", bty = "l")
lines(prds[, 1], prds[, 2], lwd = 2)
lines(prds[, 1], prds[, 4], lwd = 2, col = "red")
lines(prds[, 1], prds[, 6], lwd = 2)

legend("topleft", legend = c("95% P.I.", "lambda_i"), col = c("black", "red"), 
       lwd = c(2, 2))

###################
# Here we will recreate the below in full Bayesian framework
#png("chinookBaseModel.png", width = 600, height = 400)

# Toggle comments for immigration/emigration, and change dc for combined density
setEPS()
#postscript("chinookEmigration.eps")
postscript("combinedDensity.eps")
par(mfrow = c(2,2))

# Chinook x Chinook, unrestored
#x <- scale(unrestored$density[unrestored$species=="chinook"]); y <- unrestored$immig[unrestored$species=="chinook"]
#x <- scale(unrestored$density[unrestored$species=="chinook"]); y <- unrestored$emig[unrestored$species=="chinook"]
x <- scale(unrestored$dc[unrestored$species=="chinook"]); y <- unrestored$immig[unrestored$species=="chinook"]

jags.data <- list(N.cells = length(x), density = x[ ,1], immigrants = y)
jm <- jags.model("model.txt", data = jags.data, n.chains = 3, n.adapt = 10000)
update(jm, n.iter = 10000)
jm.sample <- jags.samples(jm, variable.names = params, n.iter = 10000, thin = 10)
predictions <- summary(as.mcmc.list(jm.sample$prediction))
#prds <- data.frame(sc50 = x[,1], predictions$quantiles)
prds <- data.frame(sc50 = x, predictions$quantiles)
prds <- prds[order(prds[, 1]), ]

plot(x, y, cex = 1, col = "darkgrey", pch = 19, ylab = "# of Chinook immigrants", 
     xlab = "Scaled Initial Combined Density", bty = "l", main = "Unrestored Pools")
lines(prds[, 1], prds[, 2], lwd = 2)
lines(prds[, 1], prds[, 4], lwd = 2, col = "red")
lines(prds[, 1], prds[, 6], lwd = 2)

#dev.off()
#png(filename = "gelmanCUW.png", width=600, height = 400)
# par(mfrow = c(3,2))
# plot(as.mcmc.list(jm.sample$beta0), main = "Beta_0", auto.layout = FALSE)
# plot(as.mcmc.list(jm.sample$beta1), main = "Beta_1", auto.layout = FALSE)
# plot(as.mcmc.list(jm.sample$beta2), main = "Beta_2", auto.layout = FALSE)
# par(mfrow=c(1,3))
# gelman.plot(as.mcmc.list(jm.sample$beta0), auto.layout = FALSE)
# gelman.plot(as.mcmc.list(jm.sample$beta1), auto.layout = FALSE)
# gelman.plot(as.mcmc.list(jm.sample$beta2), auto.layout = FALSE)
# 
# gelman.diag(as.mcmc.list(jm.sample$beta0))
# gelman.diag(as.mcmc.list(jm.sample$beta1))
# gelman.diag(as.mcmc.list(jm.sample$beta2))

#legend("topleft", legend = c("95% P.I.", "lambda_i"), col = c("black", "red"), 
#       lwd = c(2, 2))

# Chinook x Chinook, restored
#x <- scale(restored$density[restored$species=="chinook"]); y <- restored$immig[restored$species=="chinook"]
#x <- scale(restored$density[restored$species=="chinook"]); y <- restored$emig[restored$species=="chinook"]
x <- scale(restored$dc[restored$species=="chinook"]); y <- restored$immig[restored$species=="chinook"]
jags.data <- list(N.cells = length(x), density = x[ ,1], immigrants = y)
jm <- jags.model("model.txt", data = jags.data, n.chains = 3, n.adapt = 1000)
update(jm, n.iter = 1000)
jm.sample <- jags.samples(jm, variable.names = params, n.iter = 1000, thin = 1)
predictions <- summary(as.mcmc.list(jm.sample$prediction))
#prds <- data.frame(sc50 = x[,1], predictions$quantiles)
prds <- data.frame(sc50 = x, predictions$quantiles)
prds <- prds[order(prds[, 1]), ]

plot(x, y, cex = 1, col = "darkgrey", pch = 19, ylab = "# of Chinook Immigrants", 
     xlab = "Scaled Initial Combined Density", bty = "l", main = "Restored Pools")
lines(prds[, 1], prds[, 2], lwd = 2)
lines(prds[, 1], prds[, 4], lwd = 2, col = "red")
lines(prds[, 1], prds[, 6], lwd = 2)



x <- scale(unrestored$dc[unrestored$species=="steelhead"]); y <- unrestored$immig[unrestored$species=="steelhead"]

jags.data <- list(N.cells = length(x), density = x[ ,1], immigrants = y)
jm <- jags.model("model.txt", data = jags.data, n.chains = 3, n.adapt = 1000)
update(jm, n.iter = 100000)
jm.sample <- jags.samples(jm, variable.names = params, n.iter = 100000, thin = 10)
predictions <- summary(as.mcmc.list(jm.sample$prediction))
#prds <- data.frame(sc50 = x[,1], predictions$quantiles)
prds <- data.frame(sc50 = x, predictions$quantiles)
prds <- prds[order(prds[, 1]), ]

plot(x, y, cex = 1, col = "darkgrey", pch = 19, ylab = "# of Steelhead immigrants", 
     xlab = "Scaled Initial Combined Density", bty = "l", main = "Unrestored Pools")
lines(prds[, 1], prds[, 2], lwd = 2)
lines(prds[, 1], prds[, 4], lwd = 2, col = "red")
lines(prds[, 1], prds[, 6], lwd = 2)

x <- scale(restored$dc[restored$species=="steelhead"]); y <- restored$immig[restored$species=="steelhead"]
jags.data <- list(N.cells = length(x), density = x[ ,1], immigrants = y)
jm <- jags.model("model.txt", data = jags.data, n.chains = 3, n.adapt = 1000)
update(jm, n.iter = 1000)
jm.sample <- jags.samples(jm, variable.names = params, n.iter = 1000, thin = 1)
predictions <- summary(as.mcmc.list(jm.sample$prediction))
#prds <- data.frame(sc50 = x[,1], predictions$quantiles)
prds <- data.frame(sc50 = x, predictions$quantiles)
prds <- prds[order(prds[, 1]), ]

plot(x, y, cex = 1, col = "darkgrey", pch = 19, ylab = "# of Steelhead Immigrants", 
     xlab = "Scaled Initial Combined Density", bty = "l", main = "Restored Pools")
lines(prds[, 1], prds[, 2], lwd = 2)
lines(prds[, 1], prds[, 4], lwd = 2, col = "red")
lines(prds[, 1], prds[, 6], lwd = 2)

dev.off()

# par(mfrow=c(1,3))
# gelman.plot(as.mcmc.list(jm.sample$alpha0), auto.layout = FALSE)
# gelman.plot(as.mcmc.list(jm.sample$alpha1), auto.layout = FALSE)
# gelman.plot(as.mcmc.list(jm.sample$alpha2), auto.layout = FALSE)
# 
# gelman.diag(as.mcmc.list(jm.sample$alpha0))
# gelman.diag(as.mcmc.list(jm.sample$alpha1))
# gelman.diag(as.mcmc.list(jm.sample$alpha2))
#legend("topleft", legend = c("95% P.I.", "lambda_i"), col = c("black", "red"), 
#       lwd = c(2, 2))

# Chinook x steelhead, unrestored
#x <- scale(unrestored$do[unrestored$species=="chinook"]); y <- unrestored$immig[unrestored$species=="chinook"]
#x <- scale(unrestored$do[unrestored$species=="chinook"]); y <- unrestored$emig[unrestored$species=="chinook"]

jags.data <- list(N.cells = length(x), density = x[ ,1], immigrants = y)
jm <- jags.model("model.txt", data = jags.data, n.chains = 3, n.adapt = 1000)
update(jm, n.iter = 1000)
jm.sample <- jags.samples(jm, variable.names = params, n.iter = 1000, thin = 1)
predictions <- summary(as.mcmc.list(jm.sample$prediction))
#prds <- data.frame(sc50 = x[,1], predictions$quantiles)
prds <- data.frame(sc50 = x, predictions$quantiles)
prds <- prds[order(prds[, 1]), ]

plot(x, y, cex = 1, col = "darkgrey", pch = 19, ylab = "# of Chinook Emigrants", 
     xlab = "Scaled Initial Steelhead Density", bty = "l", main = "Unrestored Pools")
lines(prds[, 1], prds[, 2], lwd = 2)
lines(prds[, 1], prds[, 4], lwd = 2, col = "red")
lines(prds[, 1], prds[, 6], lwd = 2)
# par(mfrow=c(1,3))
# gelman.plot(as.mcmc.list(jm.sample$alpha0), auto.layout = FALSE)
# gelman.plot(as.mcmc.list(jm.sample$alpha1), auto.layout = FALSE)
# gelman.plot(as.mcmc.list(jm.sample$alpha2), auto.layout = FALSE)
# 
# gelman.diag(as.mcmc.list(jm.sample$alpha0))
# gelman.diag(as.mcmc.list(jm.sample$alpha1))
# gelman.diag(as.mcmc.list(jm.sample$alpha2))



#legend("topleft", legend = c("95% P.I.", "lambda_i"), col = c("black", "red"), 
#       lwd = c(2, 2))

# Chinook x steelhead, restored
#x <- scale(restored$do[restored$species=="chinook"]); y <- restored$immig[restored$species=="chinook"]
x <- scale(restored$do[restored$species=="chinook"]); y <- restored$emig[restored$species=="chinook"]
jags.data <- list(N.cells = length(x), density = x[ ,1], immigrants = y)
jm <- jags.model("model.txt", data = jags.data, n.chains = 3, n.adapt = 1000)
update(jm, n.iter = 1000)
jm.sample <- jags.samples(jm, variable.names = params, n.iter = 1000, thin = 1)
predictions <- summary(as.mcmc.list(jm.sample$prediction))
#prds <- data.frame(sc50 = x[,1], predictions$quantiles)
prds <- data.frame(sc50 = x, predictions$quantiles)
prds <- prds[order(prds[, 1]), ]

plot(x, y, cex = 1, col = "darkgrey", pch = 19, ylab = "# of Chinook Emigrants", 
     xlab = "Scaled Initial Steelhead Density", bty = "l", main = "Restored Pools")
lines(prds[, 1], prds[, 2], lwd = 2)
lines(prds[, 1], prds[, 4], lwd = 2, col = "red")
lines(prds[, 1], prds[, 6], lwd = 2)
dev.off()
# par(mfrow=c(1,3))
# gelman.plot(as.mcmc.list(jm.sample$alpha0), auto.layout = FALSE)
# gelman.plot(as.mcmc.list(jm.sample$alpha1), auto.layout = FALSE)
# gelman.plot(as.mcmc.list(jm.sample$alpha2), auto.layout = FALSE)
# 
# gelman.diag(as.mcmc.list(jm.sample$alpha0))
# gelman.diag(as.mcmc.list(jm.sample$alpha1))
# gelman.diag(as.mcmc.list(jm.sample$alpha2))
# par(mfrow=c(1,2))

###################



# Chinook x steelhead, unrestored
x <- scale(unrestored$do[unrestored$species=="chinook"]); y <- unrestored$immig[unrestored$species=="chinook"]
jags.data <- list(N.cells = length(x), density = x[ ,1], immigrants = y)
jm <- jags.model("model.txt", data = jags.data, n.chains = 3, n.adapt = 1000)
update(jm, n.iter = 1000)
jm.sample <- jags.samples(jm, variable.names = params, n.iter = 1000, thin = 1)
predictions <- summary(as.mcmc.list(jm.sample$prediction))
#prds <- data.frame(sc50 = x[,1], predictions$quantiles)
prds <- data.frame(sc50 = x, predictions$quantiles)
prds <- prds[order(prds[, 1]), ]

plot(x, y, cex = 1, col = "darkgrey", pch = 19, ylab = "# ofChinook Immigrants", 
     xlab = "Scaled Initial Heterospecific Density", bty = "l", main = "Unrestored Pools")
lines(prds[, 1], prds[, 2], lwd = 2)
lines(prds[, 1], prds[, 4], lwd = 2, col = "red")
lines(prds[, 1], prds[, 6], lwd = 2)
# par(mfrow=c(1,3))
# gelman.plot(as.mcmc.list(jm.sample$alpha0), auto.layout = FALSE)
# gelman.plot(as.mcmc.list(jm.sample$alpha1), auto.layout = FALSE)
# gelman.plot(as.mcmc.list(jm.sample$alpha2), auto.layout = FALSE)
# 
# gelman.diag(as.mcmc.list(jm.sample$alpha0))
# gelman.diag(as.mcmc.list(jm.sample$alpha1))
# gelman.diag(as.mcmc.list(jm.sample$alpha2))
# par(mfrow=c(1,1))
#legend("topleft", legend = c("95% P.I.", "lambda_i"), col = c("black", "red"), 
#       lwd = c(2, 2))

# Chinook x steelhead, restored
x <- scale(restored$do[restored$species=="chinook"]); y <- restored$immig[restored$species=="chinook"]
jags.data <- list(N.cells = length(x), density = x[ ,1], immigrants = y)
jm <- jags.model("model.txt", data = jags.data, n.chains = 3, n.adapt = 1000)
update(jm, n.iter = 1000)
jm.sample <- jags.samples(jm, variable.names = params, n.iter = 1000, thin = 1)
predictions <- summary(as.mcmc.list(jm.sample$prediction))
#prds <- data.frame(sc50 = x[,1], predictions$quantiles)
prds <- data.frame(sc50 = x, predictions$quantiles)
prds <- prds[order(prds[, 1]), ]

plot(x, y, cex = 1, col = "darkgrey", pch = 19, ylab = "# of Chinook Immigrants", 
     xlab = "Scaled Initial Heterospecific Density", bty = "l", main = "Restored Pools")
lines(prds[, 1], prds[, 2], lwd = 2)
lines(prds[, 1], prds[, 4], lwd = 2, col = "red")
lines(prds[, 1], prds[, 6], lwd = 2)
par(mfrow=c(1,3))
gelman.plot(as.mcmc.list(jm.sample$alpha0), auto.layout = FALSE)
gelman.plot(as.mcmc.list(jm.sample$alpha1), auto.layout = FALSE)
gelman.plot(as.mcmc.list(jm.sample$alpha2), auto.layout = FALSE)

gelman.diag(as.mcmc.list(jm.sample$alpha0))
gelman.diag(as.mcmc.list(jm.sample$alpha1))
gelman.diag(as.mcmc.list(jm.sample$alpha2))


#dev.off()
#legend("topleft", legend = c("95% P.I.", "lambda_i"), col = c("black", "red"), 
#       lwd = c(2, 2))

####
# Same for Steelhead
png("SteelheadBaseModel.png", width = 600, height = 400)

setEPS()
postscript("steelheadEmigration.eps")
par(mfrow = c(2,2))

# Steelhead x Steelhead, unrestored
x <- scale(unrestored$density[unrestored$species=="steelhead"]); y <- unrestored$emig[unrestored$species=="steelhead"]
jags.data <- list(N.cells = length(x), density = x[ ,1], immigrants = y)
jm <- jags.model("model.txt", data = jags.data, n.chains = 3, n.adapt = 1000)
update(jm, n.iter = 1000)
jm.sample <- jags.samples(jm, variable.names = params, n.iter = 1000, thin = 1)
predictions <- summary(as.mcmc.list(jm.sample$prediction))
#prds <- data.frame(sc50 = x[,1], predictions$quantiles)
prds <- data.frame(sc50 = x, predictions$quantiles)
prds <- prds[order(prds[, 1]), ]

plot(x, y, cex = 1, col = "darkgrey", pch = 19, ylab = "# of Steelhead Emigrants", 
     xlab = "Scaled Initial Steelhead Density", bty = "l", main = "Unrestored Pools")
lines(prds[, 1], prds[, 2], lwd = 2)
lines(prds[, 1], prds[, 4], lwd = 2, col = "red")
lines(prds[, 1], prds[, 6], lwd = 2)

#legend("topleft", legend = c("95% P.I.", "lambda_i"), col = c("black", "red"), 
#       lwd = c(2, 2))

# Steelhead x Steelhead, restored
x <- scale(restored$density[restored$species=="steelhead"]); y <- restored$emig[restored$species=="steelhead"]
jags.data <- list(N.cells = length(x), density = x[ ,1], immigrants = y)
jm <- jags.model("model.txt", data = jags.data, n.chains = 3, n.adapt = 1000)
update(jm, n.iter = 1000)
jm.sample <- jags.samples(jm, variable.names = params, n.iter = 1000, thin = 1)
predictions <- summary(as.mcmc.list(jm.sample$prediction))
#prds <- data.frame(sc50 = x[,1], predictions$quantiles)
prds <- data.frame(sc50 = x, predictions$quantiles)
prds <- prds[order(prds[, 1]), ]

plot(x, y, cex = 1, col = "darkgrey", pch = 19, ylab = "# of Steelhead Emigrants", 
     xlab = "Scaled Initial Steelhead Density", bty = "l", main = "Restored Pools")
lines(prds[, 1], prds[, 2], lwd = 2)
lines(prds[, 1], prds[, 4], lwd = 2, col = "red")
lines(prds[, 1], prds[, 6], lwd = 2)

#legend("topleft", legend = c("95% P.I.", "lambda_i"), col = c("black", "red"), 
#       lwd = c(2, 2))

# Steelhead x chinook, unrestored
x <- scale(unrestored$do[unrestored$species=="steelhead"]); y <- unrestored$emig[unrestored$species=="steelhead"]
jags.data <- list(N.cells = length(x), density = x[ ,1], immigrants = y)
jm <- jags.model("model.txt", data = jags.data, n.chains = 3, n.adapt = 1000)
update(jm, n.iter = 1000)
jm.sample <- jags.samples(jm, variable.names = params, n.iter = 1000, thin = 1)
predictions <- summary(as.mcmc.list(jm.sample$prediction))
#prds <- data.frame(sc50 = x[,1], predictions$quantiles)
prds <- data.frame(sc50 = x, predictions$quantiles)
prds <- prds[order(prds[, 1]), ]

plot(x, y, cex = 1, col = "darkgrey", pch = 19, ylab = "# of Steelhead Emigrants", 
     xlab = "Scaled Initial Chinook Density", bty = "l", main = "Unrestored Pools")
lines(prds[, 1], prds[, 2], lwd = 2)
lines(prds[, 1], prds[, 4], lwd = 2, col = "red")
lines(prds[, 1], prds[, 6], lwd = 2)

#legend("topleft", legend = c("95% P.I.", "lambda_i"), col = c("black", "red"), 
#       lwd = c(2, 2))

# Steelhead x chinook, restored
x <- scale(restored$do[restored$species=="steelhead"]); y <- restored$emig[restored$species=="steelhead"]
jags.data <- list(N.cells = length(x), density = x[ ,1], immigrants = y)
jm <- jags.model("model.txt", data = jags.data, n.chains = 3, n.adapt = 10000)
update(jm, n.iter = 10000)
jm.sample <- jags.samples(jm, variable.names = params, n.iter = 10000, thin = 10)
predictions <- summary(as.mcmc.list(jm.sample$prediction))
#prds <- data.frame(sc50 = x[,1], predictions$quantiles)
prds <- data.frame(sc50 = x, predictions$quantiles)
prds <- prds[order(prds[, 1]), ]

plot(x, y, cex = 1, col = "darkgrey", pch = 19, ylab = "# of Steelhead Emigrants", 
     xlab = "Scaled Initial Chinook Density", bty = "l", main = "Restored Pools")
lines(prds[, 1], prds[, 2], lwd = 2)
lines(prds[, 1], prds[, 4], lwd = 2, col = "red")
lines(prds[, 1], prds[, 6], lwd = 2)
dev.off()
#legend("topleft", legend = c("95% P.I.", "lambda_i"), col = c("black", "red"), 
#       lwd = c(2, 2))

par(mfrow=c(1,2))
# Chinook x steelhead, unrestored
x <- scale(unrestored$do[unrestored$species=="steelhead"]); y <- unrestored$immig[unrestored$species=="steelhead"]
jags.data <- list(N.cells = length(x), density = x[ ,1], immigrants = y)
jm <- jags.model("model.txt", data = jags.data, n.chains = 3, n.adapt = 10000)
update(jm, n.iter = 10000)
jm.sample <- jags.samples(jm, variable.names = params, n.iter = 10000, thin = 10)
predictions <- summary(as.mcmc.list(jm.sample$prediction))
#prds <- data.frame(sc50 = x[,1], predictions$quantiles)
prds <- data.frame(sc50 = x, predictions$quantiles)
prds <- prds[order(prds[, 1]), ]

plot(x, y, cex = 1, col = "darkgrey", pch = 19, ylab = "# Steelhead Immigrants", 
     xlab = "Scaled Initial Heterospecific Density", bty = "l", main = "Unrestored Pools", ylim = c(0,30))
lines(prds[, 1], prds[, 2], lwd = 2)
lines(prds[, 1], prds[, 4], lwd = 2, col = "red")
lines(prds[, 1], prds[, 6], lwd = 2)
par(mfrow=c(1,3))
gelman.plot(as.mcmc.list(jm.sample$alpha0), auto.layout = FALSE)
gelman.plot(as.mcmc.list(jm.sample$alpha1), auto.layout = FALSE)
gelman.plot(as.mcmc.list(jm.sample$alpha2), auto.layout = FALSE)

gelman.diag(as.mcmc.list(jm.sample$alpha0))
gelman.diag(as.mcmc.list(jm.sample$alpha1))
gelman.diag(as.mcmc.list(jm.sample$alpha2))
par(mfrow=c(1,1))
#legend("topleft", legend = c("95% P.I.", "lambda_i"), col = c("black", "red"), 
#       lwd = c(2, 2))

# Steelhead x Chinook, restored
x <- scale(restored$do[restored$species=="steelhead"]); y <- restored$immig[restored$species=="steelhead"]
jags.data <- list(N.cells = length(x), density = x[ ,1], immigrants = y)
jm <- jags.model("model.txt", data = jags.data, n.chains = 3, n.adapt = 10000)
update(jm, n.iter = 10000)
jm.sample <- jags.samples(jm, variable.names = params, n.iter = 10000, thin = 10)
predictions <- summary(as.mcmc.list(jm.sample$prediction))
#prds <- data.frame(sc50 = x[,1], predictions$quantiles)
prds <- data.frame(sc50 = x, predictions$quantiles)
prds <- prds[order(prds[, 1]), ]

plot(x, y, cex = 1, col = "darkgrey", pch = 19, ylab = "# of Steelhead Immigrants", 
     xlab = "Scaled Initial Heterospecific Density", bty = "l", main = "Restored Pools", ylim = c(0,30))
lines(prds[, 1], prds[, 2], lwd = 2)
lines(prds[, 1], prds[, 4], lwd = 2, col = "red")
lines(prds[, 1], prds[, 6], lwd = 2)
par(mfrow=c(1,3))
gelman.plot(as.mcmc.list(jm.sample$beta0), auto.layout = FALSE)
gelman.plot(as.mcmc.list(jm.sample$beta1), auto.layout = FALSE)
gelman.plot(as.mcmc.list(jm.sample$beta2), auto.layout = FALSE)

gelman.diag(as.mcmc.list(jm.sample$beta0))
gelman.diag(as.mcmc.list(jm.sample$beta1))
gelman.diag(as.mcmc.list(jm.sample$beta2))



# Now to figure out the weighted density shit

params2 <- c("beta0", "beta1", "beta2", "eta0", "eta1", "eta2", "prediction")

# Chinook weighted, Unrestored
x <- scale(unrestored$density[unrestored$species=="chinook"]); y <- unrestored$immig[unrestored$species=="chinook"]
#x <- unrestored$density[unrestored$species=="chinook"]; y <- unrestored$immig[unrestored$species=="chinook"]
xc <- scale(unrestored$do[unrestored$species=="chinook"])

jags.data <- list(N.cells = length(x), denC = x[ , 1], denS = xc[ , 1], immigrants = y)
jm <- jags.model("modelCombined.txt", data = jags.data, n.chains = 5, n.adapt = 300000)
update(jm, n.iter = 1500000)
jm.sample0 <- jags.samples(jm, variable.names = params2, n.iter = 100000, thin = 100)
png(filename = "betaCUW.png", width = 600, height = 400)
par(mfrow = c(3,2))
plot(as.mcmc.list(jm.sample0$beta0), main = "Beta_0", auto.layout = FALSE)
plot(as.mcmc.list(jm.sample0$beta1), main = "Beta_1", auto.layout = FALSE)
plot(as.mcmc.list(jm.sample0$beta2), main = "Beta_2", auto.layout = FALSE)
dev.off()
#par(mfrow=c(2,2))
png(filename = "etaCUW.png", width = 600, height = 400)
par(mfrow = c(3,2))
plot(as.mcmc.list(jm.sample0$eta0), main = "Eta_0", auto.layout = FALSE)
plot(as.mcmc.list(jm.sample0$eta1), main = "Eta_1", auto.layout = FALSE)
plot(as.mcmc.list(jm.sample0$eta2), main = "Eta_2", auto.layout = FALSE)
dev.off()
png(filename = "gelmanCUW.png", width=600, height = 400)
par(mfrow=c(2,3))
gelman.plot(as.mcmc.list(jm.sample0$beta0), auto.layout = FALSE)
gelman.plot(as.mcmc.list(jm.sample0$beta1), auto.layout = FALSE)
gelman.plot(as.mcmc.list(jm.sample0$beta2), auto.layout = FALSE)
gelman.plot(as.mcmc.list(jm.sample0$eta0), auto.layout = FALSE)
gelman.plot(as.mcmc.list(jm.sample0$eta1), auto.layout = FALSE)
gelman.plot(as.mcmc.list(jm.sample0$eta2), auto.layout = FALSE)
dev.off()
gelman.diag(as.mcmc.list(jm.sample0$beta0))
gelman.diag(as.mcmc.list(jm.sample0$beta1))
gelman.diag(as.mcmc.list(jm.sample0$beta2))
gelman.diag(as.mcmc.list(jm.sample0$eta0))
gelman.diag(as.mcmc.list(jm.sample0$eta1))
gelman.diag(as.mcmc.list(jm.sample0$eta2))

# Attempt to plot this mess
scatterplot3d(x, xc, y, pch = 19, highlight.3d = TRUE, type = "h", xlab = "Chinook Density", ylab = "Steelhead Density", zlab = "Chinook Immigrants", main = "Combined Density Model")


plot3d(x, xc, y, pch = 19, xlab = "Chinook Density", ylab = "Steelhead Density", zlab = "Chinook Immigrants", main = "Combined Density Model", col = "red", size = 3)

# Chinook weighted, Restored
x <- scale(restored$density[restored$species=="chinook"]); y <- restored$immig[restored$species=="chinook"]
#x <- unrestored$density[unrestored$species=="chinook"]; y <- unrestored$immig[unrestored$species=="chinook"]
xc <- scale(restored$do[restored$species=="chinook"])

jags.data <- list(N.cells = length(x), denC = x[ , 1], denS = xc[ , 1], immigrants = y)
jm <- jags.model("modelCombined.txt", data = jags.data, n.chains = 5, n.adapt = 300000)
update(jm, n.iter = 1500000)
jm.sample1 <- jags.samples(jm, variable.names = params2, n.iter = 100000, thin = 100)
png(filename = "betaCRW.png", width = 600, height = 400)
par(mfrow = c(3,2))
plot(as.mcmc.list(jm.sample1$beta0), main = "Beta_0", auto.layout = FALSE)
plot(as.mcmc.list(jm.sample1$beta1), main = "Beta_1", auto.layout = FALSE)
plot(as.mcmc.list(jm.sample1$beta2), main = "Beta_2", auto.layout = FALSE)
#par(mfrow=c(2,2))
dev.off()
png(filename = "etaCRW.png", width = 600, height = 400)
par(mfrow=c(3,2))
plot(as.mcmc.list(jm.sample1$eta0), main = "Eta_0", auto.layout = FALSE)
plot(as.mcmc.list(jm.sample1$eta1), main = "Eta_1", auto.layout = FALSE)
plot(as.mcmc.list(jm.sample1$eta2), main = "Eta_2", auto.layout = FALSE)
dev.off()
png(filename = "gelmanCRW.png", width = 600, height = 400)
par(mfrow=c(2,3))
gelman.plot(as.mcmc.list(jm.sample1$beta0), auto.layout = FALSE)
gelman.plot(as.mcmc.list(jm.sample1$beta1), auto.layout = FALSE)
gelman.plot(as.mcmc.list(jm.sample1$beta2), auto.layout = FALSE)
gelman.plot(as.mcmc.list(jm.sample1$eta0), auto.layout = FALSE)
gelman.plot(as.mcmc.list(jm.sample1$eta1), auto.layout = FALSE)
gelman.plot(as.mcmc.list(jm.sample1$eta2), auto.layout = FALSE)
dev.off()

gelman.diag(as.mcmc.list(jm.sample1$beta0))
gelman.diag(as.mcmc.list(jm.sample1$beta1))
gelman.diag(as.mcmc.list(jm.sample1$beta2))
gelman.diag(as.mcmc.list(jm.sample1$eta0))
gelman.diag(as.mcmc.list(jm.sample1$eta1))
gelman.diag(as.mcmc.list(jm.sample1$eta2))


# steelhead weighted, Unrestored
x <- scale(unrestored$density[unrestored$species=="steelhead"]); y <- unrestored$immig[unrestored$species=="steelhead"]
#x <- unrestored$density[unrestored$species=="steelhead"]; y <- unrestored$immig[unrestored$species=="steelhead"]
xc <- scale(unrestored$do[unrestored$species=="steelhead"])

jags.data <- list(N.cells = length(x), denC = x[ , 1], denS = xc[ , 1], immigrants = y)
jm <- jags.model("modelCombined.txt", data = jags.data, n.chains = 5, n.adapt = 500000)
update(jm, n.iter = 200000)
jm.sample2 <- jags.samples(jm, variable.names = params2, n.iter = 100000, thin = 100)
png(filename = "betaSUW.png", width = 600, height = 400)
par(mfrow = c(3,2))
plot(as.mcmc.list(jm.sample2$beta0), main = "Beta_0", auto.layout = FALSE)
plot(as.mcmc.list(jm.sample2$beta1), main = "Beta_1", auto.layout = FALSE)
plot(as.mcmc.list(jm.sample2$beta2), main = "Beta_2", auto.layout = FALSE)
dev.off()
png(filename = "etaSUW.png", width = 600, height = 400)
par(mfrow=c(3,2))
plot(as.mcmc.list(jm.sample2$eta0), main = "Eta_0", auto.layout = FALSE)
plot(as.mcmc.list(jm.sample2$eta1), main = "Eta_1", auto.layout = FALSE)
plot(as.mcmc.list(jm.sample2$eta2), main = "Eta_2", auto.layout = FALSE)
dev.off()
png(filename = "gelmanSUW.png", width = 600, height = 400)
par(mfrow=c(2,3))
gelman.plot(as.mcmc.list(jm.sample2$beta0), auto.layout = FALSE)
gelman.plot(as.mcmc.list(jm.sample2$beta1), auto.layout = FALSE)
gelman.plot(as.mcmc.list(jm.sample2$beta2), auto.layout = FALSE)
gelman.plot(as.mcmc.list(jm.sample2$eta0), auto.layout = FALSE)
gelman.plot(as.mcmc.list(jm.sample2$eta1), auto.layout = FALSE)
gelman.plot(as.mcmc.list(jm.sample2$eta2), auto.layout = FALSE)
dev.off()

gelman.diag(as.mcmc.list(jm.sample2$beta0))
gelman.diag(as.mcmc.list(jm.sample2$beta1))
gelman.diag(as.mcmc.list(jm.sample2$beta2))
gelman.diag(as.mcmc.list(jm.sample2$eta0))
gelman.diag(as.mcmc.list(jm.sample2$eta1))
gelman.diag(as.mcmc.list(jm.sample2$eta2))

# steelhead weighted, Restored
x <- scale(restored$density[restored$species=="steelhead"]); y <- restored$immig[restored$species=="steelhead"]
#x <- unrestored$density[unrestored$species=="steelhead"]; y <- unrestored$immig[unrestored$species=="steelhead"]
xc <- scale(restored$do[restored$species=="steelhead"])

jags.data <- list(N.cells = length(x), denC = x[ , 1], denS = xc[ , 1], immigrants = y)
jm <- jags.model("modelCombined.txt", data = jags.data, n.chains = 5, n.adapt = 500000)
update(jm, n.iter = 200000)
jm.sample3 <- jags.samples(jm, variable.names = params2, n.iter = 100000, thin = 100)
png(filename = "betaSRW.png", width = 600, height = 400)
par(mfrow = c(3,2))
plot(as.mcmc.list(jm.sample3$beta0), main = "Beta_0", auto.layout = FALSE)
plot(as.mcmc.list(jm.sample3$beta1), main = "Beta_1", auto.layout = FALSE)
plot(as.mcmc.list(jm.sample3$beta2), main = "Beta_2", auto.layout = FALSE)
dev.off()
png(filename = "etaSRW.png", width = 600, height = 400)
par(mfrow=c(3,2))
plot(as.mcmc.list(jm.sample3$eta0), main = "Eta_0", auto.layout = FALSE)
plot(as.mcmc.list(jm.sample3$eta1), main = "Eta_1", auto.layout = FALSE)
plot(as.mcmc.list(jm.sample3$eta2), main = "Eta_2", auto.layout = FALSE)
dev.off()
png(filename = "gelmanSRW.png", width = 600, height = 400)
par(mfrow=c(2,3))
gelman.plot(as.mcmc.list(jm.sample3$beta0), auto.layout = FALSE)
gelman.plot(as.mcmc.list(jm.sample3$beta1), auto.layout = FALSE)
gelman.plot(as.mcmc.list(jm.sample3$beta2), auto.layout = FALSE)
gelman.plot(as.mcmc.list(jm.sample3$eta0), auto.layout = FALSE)
gelman.plot(as.mcmc.list(jm.sample3$eta1), auto.layout = FALSE)
gelman.plot(as.mcmc.list(jm.sample3$eta2), auto.layout = FALSE)
dev.off()


gelman.diag(as.mcmc.list(jm.sample3$beta0))
gelman.diag(as.mcmc.list(jm.sample3$beta1))
gelman.diag(as.mcmc.list(jm.sample3$beta2))
gelman.diag(as.mcmc.list(jm.sample3$eta0))
gelman.diag(as.mcmc.list(jm.sample3$eta1))
gelman.diag(as.mcmc.list(jm.sample3$eta2))



# store "competition" parameters
fit.pars <- matrix(NA, nrow = 4, ncol = 6)

summary(restored$density[restored$species=="chinook"])
summary(restored$density[restored$species=="steelhead"])

summary(unrestored$density[unrestored$species=="chinook"])
summary(unrestored$density[unrestored$species=="steelhead"])

a<-summary(unrestored$immig[unrestored$species=="steelhead"])
b<-summary(restored$immig[restored$species=="steelhead"])
t.test(a,b, alternative = "t")

a<-unrestored$immig[unrestored$species=="steelhead"]
b<-restored$immig[restored$species=="steelhead"]

a<-unrestored$immig[unrestored$species=="chinook"]
b<-restored$immig[restored$species=="chinook"]

summary(restored$immig)
summary(unrestored$immig)
t_test_chinook <- BESTmcmc(a,b)
t.test(a,b)

setEPS()
postscript("best_chinook.eps")
plotAll(t_test_chinook)
dev.off()

setEPS()
postscript("best_steelhead.eps")
plotAll(t_test)
dev.off()
summary(t_test_chinook)

par(mfrow=c(2,2))
x <- scale(unrestored$density[unrestored$species=="chinook"]); y <- unrestored$immig[unrestored$species=="chinook"]
#x <- unrestored$density[unrestored$species=="chinook"]; y <- unrestored$immig[unrestored$species=="chinook"]
xc <- unrestored$do[unrestored$species=="chinook"]
plot(x, y, xlab = "Chinook Initial Density", ylab = "Chinook Immigration", main = "Unrestored Pools", bty = "l", pch = 19)

#unRestChinook <- unrestored$density[unrestored$species=="chinook"]
#LikelihoodSSE(Immigration(c(1,0.5,2),unRestChinook), unrestored$immig[unrestored$species=="chinook"])

pars0 <- c(1,0.5,2)
pars1 <- c(1,0.5,2,0,0.5,0.5)
fit <- optim(pars1, fn=WeightedLikelihoodSSE, densityChinook = x, densitySteelhead = xc, immigrants =y)
fit.pars[1, ] <- fit$par
#fit <- optim(pars0, fn=LikelihoodSSE, density = x, immigrants = y)
fit <- optim(pars0, fn=LikelihoodPois, density = x, immigrants = y)

lines(seq(min(x),max(x),by=0.1), Immigration(fit$par, seq(min(x),max(x),by=0.1)))


x <- restored$density[restored$species=="chinook"]; y <- restored$immig[restored$species=="chinook"]
xc <- restored$do[restored$species=="chinook"]
plot(x, y, xlab = "Chinook Initial Density", ylab = "Chinook Immigration", main = "Restored Pools", bty = "l", pch = 19)

fit <- optim(pars1, fn=WeightedLikelihoodSSE, densityChinook = x, densitySteelhead = xc, immigrants =y)
fit.pars[2, ] <- fit$par

#fit <- optim(pars0, fn=LikelihoodSSE, density = x, immigrants = y)
fit <- optim(pars0, fn=LikelihoodPois, density = x, immigrants = y)
lines(seq(min(x),max(x),by=0.1), Immigration(fit$par, seq(min(x),max(x),by=0.1)))


x <- unrestored$dc[unrestored$species=="chinook"]; y <- unrestored$immig[unrestored$species=="chinook"]
plot(x, y, xlab = "Total Initial Density", ylab = "Chinook Immigration", main = "Unrestored Pools", bty = "l", pch = 19)
fit <- optim(pars0, fn=LikelihoodSSE, density = x, immigrants = y)
lines(seq(min(x),max(x),by=0.1), Immigration(fit$par, seq(min(x),max(x),by=0.1)))

x <- restored$dc[restored$species=="chinook"]; y <- restored$immig[restored$species=="chinook"]
plot(x, y, xlab = "Total Initial Density", ylab = "Chinook Immigration", main = "Restored Pools", bty = "l", pch = 19)
fit <- optim(pars0, fn=LikelihoodPois, density = x, immigrants = y)
lines(seq(min(x),max(x),by=0.1), Immigration(fit$par, seq(min(x),max(x),by=0.1)))


# Now Steelhead
x <- unrestored$density[unrestored$species=="steelhead"]; y <- unrestored$immig[unrestored$species=="steelhead"]
xc <- unrestored$do[unrestored$species=="steelhead"]

plot(x, y, xlab = "steelhead Initial Density", ylab = "Steelhead Immigration", main = "Unrestored Pools", bty = "l", pch = 19)

fit <- optim(pars1, fn=WeightedLikelihoodSSE, densityChinook = xc, densitySteelhead = x, immigrants =y)
fit.pars[3, ] <- fit$par

fit <- optim(pars0, fn=LikelihoodPois, density = x, immigrants = y)
lines(seq(min(x),max(x),by=0.1), Immigration(fit$par, seq(min(x),max(x),by=0.1)))

x <- restored$density[restored$species=="steelhead"]; y <- restored$immig[restored$species=="steelhead"]
xc <- restored$do[restored$species=="steelhead"]
plot(x, y, xlab = "steelhead Initial Density", ylab = "Steelhead Immigration", main = "Restored Pools", bty = "l", pch = 19)

fit <- optim(pars1, fn=WeightedLikelihoodSSE, densityChinook = xc, densitySteelhead = x, immigrants =y)
fit.pars[4, ] <- fit$par

fit <- optim(pars0, fn=LikelihoodPois, density = x, immigrants = y)
lines(seq(min(x),max(x),by=0.1), Immigration(fit$par, seq(min(x),max(x),by=0.1)))

x <- unrestored$dc[unrestored$species=="steelhead"]; y <- unrestored$immig[unrestored$species=="steelhead"]
plot(x, y, xlab = "Total Initial Density", ylab = "Steelhead Immigration", main = "Unrestored Pools", bty = "l", pch = 19)
fit <- optim(pars0, fn=LikelihoodPois, density = x, immigrants = y)
lines(seq(min(x),max(x),by=0.1), Immigration(fit$par, seq(min(x),max(x),by=0.1)))

x <- restored$dc[restored$species=="steelhead"]; y <- restored$immig[restored$species=="steelhead"]
plot(x, y, xlab = "Total Initial Density", ylab = "Steelhead Immigration", main = "Restored Pools", bty = "l", pch = 19)
fit <- optim(pars0, fn=LikelihoodPois, density = x, immigrants = y)
lines(seq(min(x),max(x),by=0.1), Immigration(fit$par, seq(min(x),max(x),by=0.1)))

# Against the opposite species only
# x <- unrestored$density[unrestored$species=="steelhead"]; y <- unrestored$immig[unrestored$species=="chinook"]
# plot(x, y, xlab = "steelhead Initial Density", ylab = "Chinook Immigration", main = "Unrestored Pools", bty = "l", pch = 19)
# fit <- optim(pars0, fn=LikelihoodSSE, density = x, immigrants = y)
# lines(seq(min(x),max(x),by=0.1), Immigration(fit$par, seq(min(x),max(x),by=0.1)))
# 
# x <- restored$density[restored$species=="steelhead"]; y <- restored$immig[restored$species=="chinook"]
# plot(x, y, xlab = "steelhead Initial Density", ylab = "Chinook Immigration", main = "Restored Pools", bty = "l", pch = 19)
# fit <- optim(pars0, fn=LikelihoodSSE, density = x, immigrants = y)
# lines(seq(min(x),max(x),by=0.1), Immigration(fit$par, seq(min(x),max(x),by=0.1)))
# 
# x <- unrestored$density[unrestored$species=="chinook"]; y <- unrestored$immig[unrestored$species=="steelhead"]
# plot(x, y, xlab = "Chinook Initial Density", ylab = "Steelhead Immigration", main = "Unrestored Pools", bty = "l", pch = 19)
# fit <- optim(pars0, fn=LikelihoodSSE, density = x, immigrants = y)
# lines(seq(min(x),max(x),by=0.1), Immigration(fit$par, seq(min(x),max(x),by=0.1)))
# 
# 
# x <- restored$density[restored$species=="chinook"]; y <- restored$immig[restored$species=="steelhead"]
# plot(x, y, xlab = "Chinook Initial Density", ylab = "Steelhead Immigration", main = "Restored Pools", bty = "l", pch = 19)
# fit <- optim(pars0, fn=LikelihoodSSE, density = x, immigrants = y)
# lines(seq(min(x),max(x),by=0.1), Immigration(fit$par, seq(min(x),max(x),by=0.1)))
