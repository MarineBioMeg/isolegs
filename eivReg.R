# Steelhead vs Chinook Density

# libraries
require(deming) #implements TLS and Deming regression

setwd("~/Dropbox/PostDoc/Fish/") # or wherever you've stored your data
isodar <- read.csv("lowerIsodarAll.csv")

# subset the data
x1 <- isodar[isodar$species=="steelhead",]$d.untreat # steelhead density in untreated reaches
x2 <- isodar[isodar$species=="chinook",]$d.untreat # chinook density in untreated reaches
y1 <- isodar[isodar$species=="steelhead",]$d.treat # steelhead, treated
y2 <- isodar[isodar$species=="chinook",]$d.treat # chinook, treated

# recreate the manuscript figure 3
plot(x1, y1, xlim = c(0,1.5), ylim = c(0,2), xlab = "unrestored density", ylab = "restored density")
points(x2, y2, pch=19)
legend(x=1, y=1, legend = c("Chinook", "Steelhead"), lty=c(1,2), pch = c(19, 1), bty="n")

# deming() implements total least squares regression, which accounts for 
# observation errors in both dependent and independent variables as long as
# errors in both variables are equally distributed. Otherwise it behaves 
# similarly to lm()

dreg1 <- deming(y1 ~ 0+x1) # original manuscript figure 3 imposed intercept of 0, so I follow suit
dreg2 <- deming(y2 ~ 0+x2) # as above
abline(dreg2)
abline(dreg1, lty=2)
