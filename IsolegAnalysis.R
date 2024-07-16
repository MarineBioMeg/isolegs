####Selectivity Analysis - Isolegs
##M. Malone
##July 16, 2024

#Set working directory

#load libaries
#packages
library(cowplot)
library(ggplot2)
library(lme4)
library(rockchalk)
library(equatiomatic)
library(cowplot)
library(ggplot2)

#load 2014 and 2015 data
bigtreat<-read.csv("20142015_HSI_isolegs.csv")

#Chinook
# x, y and z coordinates
x <- bigtreat$total_C
z <- bigtreat$P.C.deep
y <- bigtreat$total_S

fit <- glm(z ~ x + y, family = "quasibinomial")
summary(fit)

extract_eq(fit, wrap = TRUE, use_coefs = TRUE)
##This is the equation for deriving Chinook isoleg

#Chinook 3d fig
plotPlane(fit, "x", "y", pch=16, col=rgb(0,0,1,.1), drawArrows=TRUE, alength=0, 
          x1lab = "Chinook Density", ylab = "Chinook Selectivity", x2lab = "Steelhead Density",  
          acol="grey", alty=1,alwd=1, theta=25, phi=0, main = "Chinook Selectivity Deep", ticktype = "detailed")


#Steelhead
# x, y and z coordinates
x <- bigtreat$total_C
z <- bigtreat$P.S.deep
y <- bigtreat$total_S

fit <- glm(z ~ x + y, family = "quasibinomial")
summary(fit)

extract_eq(fit, wrap = TRUE, use_coefs = TRUE)
##This is the equation for deriving steelhead isoleg


#steelhead 3d fig
plotPlane(fit, "x", "y", pch=16, col=rgb(0,0,1,.1), drawArrows=TRUE, alength=0, 
          x1lab = "Chinook Density", ylab = "Steelhead Selectivity", x2lab = "Steelhead Density",  
          acol="grey", alty=1,alwd=1, theta=25, phi=0, main = "Steelhead Selectivity Deep", ticktype = "detailed")

####2d graphs
###Treated Steelhead Only Preference Deep and Shallow
steeldeep<-read.csv("steelpreftreat.csv") #selectivity of s for deep habitats
###Treated Chinook Preference Deep and Shallow
chindeep<-read.csv("chinookpreftreat.csv") #selectivity of c for deep habitats

chinpgraph<-ggplot(data = chindeep, aes(x = total_C, y = P.C.deep), size = 3) + 
  geom_point() +
  #geom_smooth(method = "lm", se = TRUE) + #use this if you want linear model
  geom_smooth(method = "glm",method.args = list(family = "quasibinomial"), se = TRUE) + #use this if you want binomial
  # xlim(0,6) +
  #ylim(0,1) +
  #geom_hline(yintercept= 0.5, col = "red") +
  xlab(" ") +
  ylab("Chinook Selectivity") +
  #  ggtitle("Treated")+
  theme_classic(base_size = 15)
chinpgraph

chinpgraph2<-ggplot(data = chindeep, aes(x = total_S, y = P.C.deep), size = 3) + 
  geom_point() +
  #geom_smooth(method = "lm", se = TRUE) +
  geom_smooth(method = "glm",method.args = list(family = "quasibinomial"), se = TRUE) + #use this if you want binomial
  # xlim(0,6) +
  #ylim(0,1) +
  #geom_hline(yintercept= 0.5, col = "red") +
  xlab(" ") +
  ylab(" ") +
  #  ggtitle("Treated")+
  theme_classic(base_size = 15)
chinpgraph2

steelpgraph<-ggplot(data = steeldeep, aes(x = total_S, y = P.S.deep), size = 3) + 
  geom_point() +
  #geom_smooth(method = "lm", se = TRUE) +
  geom_smooth(method = "glm",method.args = list(family = "quasibinomial"), se = TRUE) + #use this if you want binomial
  # xlim(0,6) +
  #ylim(0,1) +
  # geom_hline(yintercept= 0.5, col = "red") +
  xlab("Steelhead Density") +
  ylab(" ") +
  #  ggtitle("Treated")+
  theme_classic(base_size = 15)
steelpgraph

steelpgraph2<-ggplot(data = steeldeep, aes(x = total_C, y = P.S.deep), size = 3) + 
  geom_point() +
  #geom_smooth(method = "lm", se = TRUE) +
  geom_smooth(method = "glm",method.args = list(family = "quasibinomial"), se = TRUE) + #use this if you want binomial
  # xlim(0,6) +
  #ylim(0,1) +
  # geom_hline(yintercept= 0.5, col = "red") +
  xlab("Chinook Density") +
  ylab("Steelhead Selectivity") +
  #  ggtitle("Treated")+
  theme_classic(base_size = 15)
steelpgraph2

png("log_preference_graphs_95CI_Sept2023.png", width = 8.0, height = 6.5, units = 'in', res = 600)
plot_grid(chinpgraph, chinpgraph2,  steelpgraph2, steelpgraph, labels = c('A', 'B', 'C', 'D'), label_size = 12 )
dev.off()
