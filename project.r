library(RColorBrewer)
library(rgdal)
library(INLA)
library(SpatialEpi)
library(ggplot2)

#load the map
data(pennLC)
pmap = pennLC$spatial.polygon

#load the data
CHD = read.csv("C:/Users/niuzc/Desktop/CHD_18_Penn_2020.csv", sep = ",", header = TRUE)
CHD$fips = paste0("0",as.character(CHD$fips))
#compute the expected values of CHD
CHD$CHD_exp = sum(CHD$CHD_obs)/sum(CHD$population)*CHD$population

#check the order of data and the map
all.equal(sapply(slot(pmap, "polygons"), function(x) {
  slot(x, "ID")
}), CHD$county_name)


#explore the variables
summary(CHD$CHD_rate)
sd(CHD$CHD_rate)
summary(CHD$CHD_SMR)
sd(CHD$CHD_SMR)
summary(CHD$high_chol_rate)
sd(CHD$high_chol_rate)
summary(CHD$smoke_rate)
sd(CHD$smoke_rate)
summary(CHD$low_edu_rate)
sd(CHD$low_edu_rate)
summary(CHD$family_income)
sd(CHD$family_income)
summary(CHD$PM25)
sd(CHD$PM25)

#map the observed data
ColorProb <- brewer.pal(7, "RdYlGn")[7:1]
quantile(CHD$CHD_obs,probs=seq(0,1,1/7))
plot(pmap, col = ColorProb[as.numeric(
  cut(CHD$CHD_obs,breaks = quantile(CHD$CHD_obs,probs=seq(0,1,1/7)),
      include.lowest=TRUE))])
title("Observed", cex = 0.75)
legend(x = "topright", fill = ColorProb[7:1],
       legend = c(">22208", "10634-22208", "7708-10634", "4535-7708",
                  "3063-4535", "2445-3063", "<2445"), cex = 0.65,
       inset = 0.03, title = "Observed")

#map the expected data
plot(pmap, col = ColorProb[as.numeric(
  cut(CHD$CHD_exp,breaks = quantile(CHD$CHD_obs,probs=seq(0,1,1/7)),
      include.lowest=TRUE))])
title("Expected", cex = 0.75)
legend(x = "topright", fill = ColorProb[7:1],
       legend = c(">22208", "10634-22208", "7708-10634", "4535-7708",
                  "3063-4535", "2445-3063", "<2445"), cex = 0.65,
       inset = 0.03, title = "Expected")

#map the rate
ColorSMR <- brewer.pal(7, "BrBG")[7:1]
quantile(CHD$CHD_rate,probs=seq(0,1,1/7))
plot(pmap, col = ColorSMR[as.numeric(
  cut(CHD$CHD_rate,breaks = quantile(CHD$CHD_rate,probs=seq(0,1,1/7)),
      include.lowest=TRUE))])
title("Rate", cex = 0.75)
legend(x = "topright", fill = ColorSMR[7:1],
       legend = c(">9.1", "8.8-9.1", "8.4-8.8", "8-8.4",
                  "7.6-8", "7.04-7.6", "<7.04"), cex = 0.65,
       inset = 0.03, title = "Rate")

#compute the SMR
CHD$CHD_SMR = CHD$CHD_obs/CHD$CHD_exp
#map the SMR
ColorSMR <- brewer.pal(7, "BrBG")[7:1]
plot(pmap, col = ColorSMR[as.numeric(
  cut(CHD$CHD_SMR, c(-0.1, 1/1.5, 1/1.25, 1/1.1,1.1, 1.25, 1.5, 100)))])
title("SMR", cex = 0.75)
legend(x = "topright", fill = ColorSMR[7:1],
       legend = c(">1.5", "1.25-1.50", "1.10-1.25","0.91-1.10", "0.80-0.91", "0.66-0.80","<0.66"),
       cex = 0.65, inset = 0.03, title = "SMR")
#plot the SMRs verus the estimated standard errors
plot(sqrt(CHD$CHD_SMR/CHD$CHD_exp)~CHD$CHD_SMR,
     xlab="SMR",ylab="Estimated standard error")
abline(coef = c(0,1),col="grey")

#compute smoothed SMR with Poisson-Lognormal model
smooth_data = data.frame(O=CHD$CHD_obs, E=CHD$CHD_exp,id=1:length(CHD$CHD_obs))
smooth_result = inla(O ~ 1 + f(id, model="iid"),data=smooth_data,
                     family="poisson", E=E,
                     control.predictor = list(compute = TRUE),
                     control.compute= list(return.marginals=TRUE))
result = data.frame("Median"=c(smooth_result$summary.fixed$`0.5quant`,
                               1/sqrt(smooth_result$summary.hyper$`0.5quant`)),
                    "CI_lower"=c(smooth_result$summary.fixed$`0.025quant`,
                                 1/sqrt(smooth_result$summary.hyper$`0.975quant`)),
                    "CI_upper"=c(smooth_result$summary.fixed$`0.975quant`,
                                 1/sqrt(smooth_result$summary.hyper$`0.025quant`)))
rownames(result) = c("beta","sigma")
knitr::kable(round(result,2), "simple")
#map the medians of the relative risk
plot(pmap, col = ColorSMR[as.numeric(
  cut(smooth_result$summary.fitted.values$`0.5quant`,
      c(-0.1, 1/1.5, 1/1.25, 1/1.1,1.1, 1.25, 1.5, 100)))])
title("Smoothed Rate", cex = 0.75)
legend(x = "topright", fill = ColorSMR[7:1],
       legend = c(">1.5", "1.25-1.50", "1.10-1.25",
                  "0.91-1.10", "0.80-0.91", "0.66-0.80","<0.66"),
       cex = 0.65, inset = 0.03, title = "Smoothed Rate")
#plot the RR versus SMR
plot(smooth_result$summary.fitted.values$`0.5quant`~CHD$CHD_SMR,
     xlab="SMR",ylab="Relative Risk",
     main="Relative Risk vs SMR")
abline(coef = c(0,1),col="red")
#plot the standard deviations of RR and SMR
plot(smooth_result$summary.fitted.values$sd~sqrt(CHD$CHD_SMR/CHD$CHD_exp),
     xlab="SMR",ylab="Relative Risk",
     main="SD of Relative Risk vs SMR")
abline(coef = c(0,1),col="red")

#compute smoothed SMR with Poisson-Lognormal-Spatial model
formula = O ~ 1 +f(id, model="bym2", graph="C:/Users/niuzc/Desktop/penn.graph", scale.model=T, constr=T,
                   hyper=list(phi=list(prior="pc", param=c(0.5, 0.5), initial=1),
                              prec=list(prior="pc.prec", param=c(0.3,0.01), initial=5)))
smooth_result2 <- inla(formula,data=smooth_data,
                       family="poisson", E=E,
                       control.predictor = list(compute = TRUE),
                       control.compute=list(return.marginals=TRUE))
result2 = data.frame("Median"=c(smooth_result2$summary.fixed$`0.5quant`,
                               1/(smooth_result2$summary.hyper$`0.5quant`[1]),
                               smooth_result2$summary.hyper$`0.5quant`[2]),
                    "CI_lower"=c(smooth_result2$summary.fixed$`0.025quant`,
                                 1/(smooth_result2$summary.hyper$`0.975quant`[1]),
                                 smooth_result2$summary.hyper$`0.025quant`[2]),
                    "CI_upper"=c(smooth_result2$summary.fixed$`0.975quant`,
                                 1/(smooth_result2$summary.hyper$`0.025quant`[1]),
                                 smooth_result2$summary.hyper$`0.975quant`[2]))
rownames(result2) = c("beta","total_var","prop_spatial")
knitr::kable(round(result2,2), "simple")
#map the medians of the relative risk
plot(pmap, col = ColorSMR[as.numeric(
  cut(smooth_result2$summary.fitted.values$`0.5quant`,
      c(-0.1, 1/1.5, 1/1.25, 1/1.1,1.1, 1.25, 1.5, 100)))])
title("Spatial Smoothed Rate", cex = 0.75)
legend(x = "topright", fill = ColorSMR[7:1],
       legend = c(">1.5", "1.25-1.50", "1.10-1.25",
                  "0.91-1.10", "0.80-0.91", "0.66-0.80","<0.66"),
       cex = 0.65, inset = 0.03, title = "Spatial Smoothed Rate")
#plot the spatial smoothed rate vs non-spatial smoothed rate
plot(smooth_result2$summary.fitted.values$`0.5quant`~smooth_result$summary.fitted.values$`0.5quant`,
     xlab="nonspatial Poisson-Lognormal",
     ylab="spatial Poisson-Lognormal",
     main="")
abline(coef = c(0,1),col="red")
#plot the spatial smoothed rate versus SMR
plot(smooth_result2$summary.fitted.values$`0.5quant`~CHD$CHD_SMR,
     xlab="SMR",ylab="spatial smoothed rate",
     main="spatial smoothed rate vs SMR")
abline(coef = c(0,1),col="red")

#(map) and plot other variables vs SMR of CHD
#high cholesterol rate ** X1
ggplot(CHD, aes(x = high_chol_rate, y = CHD_SMR)) + geom_point() + labs(y = "SMR")
quantile(CHD$high_chol_rate,probs=seq(0,1,1/7))
plot(pmap, col = ColorProb[as.numeric(
  cut(CHD$high_chol_rate,breaks = quantile(CHD$high_chol_rate,probs=seq(0,1,1/7)),
      include.lowest=TRUE))])
title("high_chol_rate (%)", cex = 0.75)
legend(x = "topright", fill = ColorProb[7:1],
       legend = c(">36.5", "35.5-36.5", "34.8-35.5", "33.9-34.8",
                  "33.1-33.9", "31.8-33.1", "<31.8"), cex = 0.65,
       inset = 0.03, title = "Observed")
#smoke rate ** X2
ggplot(CHD, aes(x = smoke_rate, y = CHD_SMR)) + geom_point() + labs(y = "SMR")
quantile(CHD$smoke_rate,probs=seq(0,1,1/7))
plot(pmap, col = ColorProb[as.numeric(
  cut(CHD$smoke_rate,breaks = quantile(CHD$smoke_rate,probs=seq(0,1,1/7)),
      include.lowest=TRUE))])
title("smoke_rate (%)", cex = 0.75)
legend(x = "topright", fill = ColorProb[7:1],
       legend = c(">20.7", "20.1-20.7", "19.4-20.1", "18.8-19.4",
                  "17.7-18.8", "16.5-17.7", "<16.5"), cex = 0.65,
       inset = 0.03, title = "Observed")
#low education rate X3
ggplot(CHD, aes(x = low_edu_rate, y = CHD_SMR)) + geom_point() + labs(y = "SMR")
#family_income ** X4
ggplot(CHD, aes(x = family_income, y = CHD_SMR)) + geom_point() + labs(y = "SMR")
#PM2.5 level * X5
ggplot(CHD, aes(x = PM25, y = CHD_SMR)) + geom_point() + labs(y = "SMR")

#Poisson-Lognormal non-spatial model a covariate
#preprare the data
regression_data = data.frame(id=1:length(CHD$CHD_obs),obs=CHD$CHD_obs,exp=CHD$CHD_exp,rate=CHD$CHD_rate,
                             chol=CDH$high_chol_rate,smoke=CDH$smoke_rate,edu=CDH$low_edu_rate,
                             inc=CDH$family_income,pm=CDH$PM25)
#~high_chol_rate(X1)
fit1X1 = inla(obs ~ 1+ chol + f(id, model="iid"),data=regression_data,
                     family="poisson", E=exp,
                     control.predictor = list(compute = TRUE),
                     control.compute= list(return.marginals=TRUE))
exp(fit1X1$summary.fixed[,3:5])
sigma1=1/sqrt(fit1X1$summary.hyperpar[3:5])
c(as.numeric(exp(-1.96*sigma1[2])),as.numeric(exp(1.96*sigma1[2])))
#~smoke_rate(X2)
fit1X2 = inla(obs ~ 1+ smoke + f(id, model="iid"),data=regression_data,
              family="poisson", E=exp,
              control.predictor = list(compute = TRUE),
              control.compute= list(return.marginals=TRUE))
exp(fit1X2$summary.fixed[,3:5])
sigma2=1/sqrt(fit1X2$summary.hyperpar[3:5])
c(as.numeric(exp(-1.96*sigma2[2])),as.numeric(exp(1.96*sigma2[2])))
#low_edu_rate(X3)
fit1X3 = inla(obs ~ 1+ edu + f(id, model="iid"),data=regression_data,
              family="poisson", E=exp,
              control.predictor = list(compute = TRUE),
              control.compute= list(return.marginals=TRUE))
exp(fit1X3$summary.fixed[,3:5])
sigma3=1/sqrt(fit1X3$summary.hyperpar[3:5])
c(as.numeric(exp(-1.96*sigma3[2])),as.numeric(exp(1.96*sigma3[2])))
#~family_income(X4)
fit1X4 = inla(obs ~ 1+ inc + f(id, model="iid"),data=regression_data,
              family="poisson", E=exp,
              control.predictor = list(compute = TRUE),
              control.compute= list(return.marginals=TRUE))
exp(fit1X4$summary.fixed[,3:5])
sigma4=1/sqrt(fit1X4$summary.hyperpar[3:5])
c(as.numeric(exp(-1.96*sigma4[2])),as.numeric(exp(1.96*sigma4[2])))
#~PM2.5(X5)
fit1X5 = inla(obs ~ 1+ pm + f(id, model="iid"),data=regression_data,
              family="poisson", E=exp,
              control.predictor = list(compute = TRUE),
              control.compute= list(return.marginals=TRUE))
exp(fit1X5$summary.fixed[,3:5])
sigma5=1/sqrt(fit1X5$summary.hyperpar[3:5])
c(as.numeric(exp(-1.96*sigma5[2])),as.numeric(exp(1.96*sigma5[2])))

#Poisson-Lognormal spatial model with a covariate
#~high_chol_rate(X1)
formula1 <- obs ~ 1 + chol + f(id, model = "bym2", graph = "C:/Users/niuzc/Desktop/penn.graph",
                              scale.model = T, constr = T, rankdef = 1, 
                              hyper = list(phi = list(prior = "pc",param = c(0.5, 0.5), initial = 1),
                                           prec = list(prior = "pc.prec",param = c(0.5/0.8753691, 0.01), initial = 5)))
fit2X1 <- inla(formula1,data=regression_data,
                     family="poisson", E=exp,
                     control.predictor = list(compute = TRUE),
                     control.compute=list(return.marginals=TRUE))
exp(fit2X1$summary.fixed[,3:5])
fit2X1$summary.hyper[,3:5]
#~smoke_rate(X2)
formula2 <- obs ~ 1 + smoke + f(id, model = "bym2", graph = "C:/Users/niuzc/Desktop/penn.graph",
                               scale.model = T, constr = T, rankdef = 1, 
                               hyper = list(phi = list(prior = "pc",param = c(0.5, 0.5), initial = 1),
                                            prec = list(prior = "pc.prec",param = c(0.5/0.8753392, 0.01), initial = 5)))
fit2X2 <- inla(formula2,data=regression_data,
               family="poisson", E=exp,
               control.predictor = list(compute = TRUE),
               control.compute=list(return.marginals=TRUE))
exp(fit2X2$summary.fixed[,3:5])
fit2X2$summary.hyper[,3:5]
#~low_edu_rate(X3)
formula3 <- obs ~ 1 + edu + f(id, model = "bym2", graph = "C:/Users/niuzc/Desktop/penn.graph",
                                scale.model = T, constr = T, rankdef = 1, 
                                hyper = list(phi = list(prior = "pc",param = c(0.5, 0.5), initial = 1),
                                             prec = list(prior = "pc.prec",param = c(0.5/0.7871205, 0.01), initial = 5)))
fit2X3 <- inla(formula3,data=regression_data,
               family="poisson", E=exp,
               control.predictor = list(compute = TRUE),
               control.compute=list(return.marginals=TRUE))
exp(fit2X3$summary.fixed[,3:5])
fit2X3$summary.hyper[,3:5]
#~family_income(X4)
formula4 <- obs ~ 1 + inc + f(id, model = "bym2", graph = "C:/Users/niuzc/Desktop/penn.graph",
                              scale.model = T, constr = T, rankdef = 1, 
                              hyper = list(phi = list(prior = "pc",param = c(0.5, 0.5), initial = 1),
                                           prec = list(prior = "pc.prec",param = c(0.5/0.850831, 0.01), initial = 5)))
fit2X4 <- inla(formula4,data=regression_data,
               family="poisson", E=exp,
               control.predictor = list(compute = TRUE),
               control.compute=list(return.marginals=TRUE))
exp(fit2X4$summary.fixed[,3:5])
fit2X4$summary.hyper[,3:5]
#~PM2.5(X5)
formula5 <- obs ~ 1 + pm + f(id, model = "bym2", graph = "C:/Users/niuzc/Desktop/penn.graph",
                              scale.model = T, constr = T, rankdef = 1, 
                              hyper = list(phi = list(prior = "pc",param = c(0.5, 0.5), initial = 1),
                                           prec = list(prior = "pc.prec",param = c(0.5/0.8067265, 0.01), initial = 5)))
fit2X5 <- inla(formula5,data=regression_data,
               family="poisson", E=exp,
               control.predictor = list(compute = TRUE),
               control.compute=list(return.marginals=TRUE))
exp(fit2X5$summary.fixed[,3:5])
fit2X5$summary.hyper[,3:5]
