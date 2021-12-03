                                                              ##HOMEASSIGNMENT 1,2,3 (Zoltan)## 

#######_______Initial preparations - installment and loading of packages and custom functions from course literature (exercise 1-14) for specific operations within the assignments__#

# package installment 
install.packages("tidyverse")
install.packages("psych")
install.packages("gridExtra")
install.packages("lm.beta")
install.packages("car")
install.packages("lmtest")
install.packages("sandwich")
install.packages("boot")
install.packages("lmboot")
install.packages("cAIC4")
install.packages("r2glmm")
install.packages("lme4")
install.packages("lmerTest")
install.packages("MuMIn")
install.packages("optimx")
install.packages("gsheet")

# package loading
library(influence.ME) 
library(lattice)
library(magrittr)
library(gsheet)
library(gridExtra)
library(lmtest)
library(sandwich)
library(cAIC4)
library(optimx)
library(car)
library(lm.beta)
library(tidyverse)
library(MuMIn)
library(lmerTest)
library(boot)
library(lme4)
library(r2glmm)
library(lmboot)
library(psych)

# custom functions loaded

bs_to_boot <- function(model, data, indices) {
  d <- data[indices, ] 
  fit <- lm(formula(model), data = d)
  return(coef(fit))
}

adjR2_to_boot <- function(model, data, indices) {
  d <- data[indices, ] 
  fit <- lm(formula(model), data = d)
  return(summary(fit)$adj.r.squared)
}

confint.boot <- function(model, data = NULL, R = 1000) {
  if (is.null(data)) {
    data = eval(parse(text = as.character(model$call[3])))
  }
  
  boot.ci_output_table = as.data.frame(matrix(NA, nrow = length(coef(model)),
                                              ncol = 2))
  row.names(boot.ci_output_table) = names(coef(model))
  names(boot.ci_output_table) = c("boot 2.5 %", "boot 97.5 %")
  results.boot = results <- boot(data = data, statistic = bs_to_boot,
                                 R = 1000, model = model)
  for (i in 1:length(coef(model))) {
    boot.ci_output_table[i, ] = unlist(unlist(boot.ci(results.boot,
                                                      type = "bca", index = i))[c("bca4", "bca5")])
  }
  return(boot.ci_output_table)
}

wild.boot.confint <- function(model, data = NULL, B = 1000) {
  if (is.null(data)) {
    data = eval(parse(text = as.character(model$call[3])))
  }
  wild_boot_estimates = wild.boot(formula(model), data = data,
                                  B = B)
  result = t(apply(wild_boot_estimates[[1]], 2, function(x) quantile(x,
                                                                     probs = c(0.025, 0.975))))
  return(result)
}

error_plotter <- function(mod, col = "black", x_var = NULL){
  mod_vars = as.character(mod$call[2])
  data = as.data.frame(eval(parse(text = as.character(mod$call[3]))))
  y = substr(mod_vars, 1, as.numeric(gregexpr(pattern ='~',mod_vars))-2)
  x = substr(mod_vars, as.numeric(gregexpr(pattern ='~',mod_vars))+2, nchar(mod_vars))
  data$pred = predict(mod)
  if(x == "1" & is.null(x_var)){x = "response_ID"
  data$response_ID = 1:nrow(data)} else if(x == "1"){x = x_var}
  plot(data[,y] ~ data[,x], ylab = y, xlab = x)
  abline(mod)
  for(i in 1:nrow(data)){
    clip(min(data[,x]), max(data[,x]), min(data[i,c(y,"pred")]), max(data[i,c(y,"pred")]))
    abline(v = data[i,x], lty = 2, col = col)
  }
}

coef_table = function(model) {
  require(lm.beta)
  mod_sum = summary(model)
  mod_sum_p_values = as.character(round(mod_sum$coefficients[,
                                                             4], 3))
  mod_sum_p_values[mod_sum_p_values != "0" & mod_sum_p_values !=
                     "1"] = substr(mod_sum_p_values[mod_sum_p_values != "0" &
                                                      mod_sum_p_values != "1"], 2, nchar(mod_sum_p_values[mod_sum_p_values !=
                                                                                                            "0" & mod_sum_p_values != "1"]))
  mod_sum_p_values[mod_sum_p_values == "0"] = "<.001"
  mod_sum_table = cbind(as.data.frame(round(cbind(coef(model),
                                                  confint(model), c(0, lm.beta(model)$standardized.coefficients[c(2:length(model$coefficients))])),
                                            2)), mod_sum_p_values)
  names(mod_sum_table) = c("b", "95%CI lb", "95%CI ub", "Std.Beta",
                           "p-value")
  mod_sum_table["(Intercept)", "Std.Beta"] = "0"
  return(mod_sum_table)
}

stdCoef.merMod <- function(object) {
  sdy <- sd(getME(object,"y"))
  sdx <- apply(getME(object,"X"), 2, sd)
  sc <- fixef(object)*sdx/sdy
  se.fixef <- coef(summary(object))[,"Std. Error"]
  se <- se.fixef*sdx/sdy
  return(data.frame(stdcoef=sc, stdse=se))
}

# STEP 1...........................Retrieving dataset1 



Data <- read.csv("https://tinyurl.com/yxm5rd89")



# STEP 2.............................Inspecting data and correcting errors



describe(Data) 
summary(Data)

# data exclusion
#______ a max value of 55 within pain data is identified --> the scale interval is from 1-10 (error)
#_______ a min value of 4.2 within STAI_trait data is identified --> the scale interval is from 20-80 (error)
#_______ ID88,ID34 is identified and excluded from data-set

Data1 <- Data[-c(88,34), ] 
view(Data1) 

# visual inspection of data-variables individually, using histograms and box-plots to detect possible outliers - normality of data-variables also investigated with Shapiro-Wilks 

# pain
Data1%>%ggplot() +aes(x = pain) +geom_histogram(binwidth=0.5, fill="#69b3a2", color="#e9ecef", alpha=0.9)
boxplot(Data1$pain) 
shapiro.test(Data1$pain) #________assumption of normality is violated based on W = 0.95673, p-value = 7.98e-05
#_________histogram and boxplot does not show any outliers. The histogram contrary to the normality test, makes me think the data is normally distributed --> no further action will be made for now.

# age
Data1 %>%ggplot() +aes(x = age) +geom_histogram(binwidth=0.5, fill="#69b3a2", color="#e9ecef", alpha=0.9)
boxplot(Data1$age) 
shapiro.test(Data1$age) #____________assumption of normality met, W = 0.99031, p-value = 0.3541
#___________histogram and boxplot do not show any outliers

# pain_cat
Data1%>%ggplot() +aes(x = pain_cat) +geom_histogram(binwidth=0.5, fill="#69b3a2", color="#e9ecef", alpha=0.9)
boxplot(Data1$pain_cat) 
shapiro.test(Data1$pain_cat) #_______________assumption of normality met, W = 0.98913, p-value = 0.2623
#________________histogram and boxplot show 2 outliers relative to the 1.5 interquartile range (ID8,ID80)

# IQ
Data1%>%ggplot() +aes(x = IQ) +geom_histogram(binwidth=3, fill="#69b3a2", color="#e9ecef", alpha=0.9) 
boxplot(Data1$IQ) 
shapiro.test(Data1$IQ) #___________assumption of normality met, W = 0.98457, p-value = 0.07611
#_________histogram and boxplot show 7 outliers relative to the 1.5 interquartile range: (ID113,ID89,ID24,ID72,ID74,ID124,ID53)

# STAI_trait
Data1%>%ggplot() +aes(x = STAI_trait) +geom_histogram(binwidth=0.5, fill="#69b3a2", color="#e9ecef", alpha=0.9)
boxplot(Data1$STAI_trait) 
shapiro.test(Data1$STAI_trait) #____________assumption of normality met, W = 0.98658, p-value = 0.1326
#____________histogram and boxplot show 1 outlier relative to the 1.5 interquartile range: (ID65)

# mindfulness
Data1 %>%ggplot() +aes(x = mindfulness) +geom_histogram(binwidth=0.2, fill="#69b3a2", color="#e9ecef", alpha=0.9) 
boxplot(Data1$mindfulness) 
shapiro.test(Data1$mindfulness) #___________assumption of normality met, W = 0.98337, p-value = 0.05464
#___________histogram and boxplot show 2 outliers relative to the 1.5 interquartile range: (ID80,ID91)

# cortisol_serum
Data1 %>%ggplot() +aes(x = cortisol_serum) +geom_histogram(binwidth=0.2, fill="#69b3a2", color="#e9ecef", alpha=0.9) 
boxplot(Data1$cortisol_serum) 
shapiro.test(Data1 $cortisol_serum) #__________________assumption of normality met, W = 0.98964, p-value = 0.2995
#______________histogram and boxplot do not show any outliers

# cortisol_saliva
Data1%>%ggplot() +aes(x = cortisol_saliva) +geom_histogram(binwidth=0.2, fill="#69b3a2", color="#e9ecef", alpha=0.9)
boxplot(Data1$cortisol_saliva) 
shapiro.test(Data1$cortisol_saliva) #____________assumption of normality met, W = 0.99193, p-value = 0.5163
#____________histogram and boxplot do not show any outliers

# sex data is grouped by sex and visually inspected in box-plot relative to pain - not numerical data 
Data1 %>% group_by(sex) %>% count()
# female =85, male =73
Data1 %>% ggplot() + aes(x = sex, y = pain) + geom_boxplot() 
# males and females largely report the same amount of pain, however females report a smaller minimum and maximum than males. 


#___________Due to the fact that IQ is not included in the regression models in assignment 1, I will not take further action to exclude the data for now --> the data also met the assumption of normally --> no data exclusion  
#___________Participant (ID80) accounts for 2 outliers. He is high in mindfulness and low in The Pain Catastrophizing Scale measures, which seems like a plausible relationship of character traits--> no data exclusion 




# STEP 3.................................Building models




# Regression model1 (age+sex as predictors of pain). 
REG1 <- lm(pain ~ sex + age, data = Data1)
summary(REG1) #Multiple R-squared:  0.08535,	Adjusted R-squared:  0.07355, F-statistic: 7.232,  p-value: 0.0009935
confint(REG1)
#               2.5 %      97.5 %
#(Intercept)  6.3481166 10.13590824
#sexmale     -0.1669602  0.75910695
#age         -0.1343980 -0.04026618


# Regressionmodel2 with age, sex, STAI_trait, pain_catastrophizing, mindfulness, and cortisol measures as predictors of pain. 
REG2 = lm(pain ~ age + sex +  STAI_trait + pain_cat + mindfulness + cortisol_serum + cortisol_saliva, data = Data1)
summary(REG2)  #Multiple R-squared:  0.5395,	Adjusted R-squared:  0.518,  F-statistic:  25.1,  p-value: < 2.2e-16
confint(REG2)
#                    2.5 %      97.5 %
#  (Intercept)   -2.98333660  2.96807841
#age             -0.06146302  0.02079715
#sexmale         -0.20301300  0.50411954
#STAI_trait      -0.07228581  0.02379824
#pain_cat         0.08767065  0.18443507
#mindfulness     -0.45439578 -0.03990979
#cortisol_serum  -0.22391900  0.57090815
#cortisol_saliva  0.07213227  0.90617760


# Visualization of the relationship of pain and the predictive variables individually and identifying cases with high leverage by creating scatter plots    

Scat1 = Data1 %>% ggplot() + aes(x = age, y = pain) + geom_point(shape = 21, fill = "red2", size = 2) + geom_smooth(method = "lm")
Scat1

Scat2 = Data1 %>% ggplot() + aes(x = sex, y = pain) + geom_point(shape = 21, fill = "red2", size = 3) + geom_smooth(method = "lm")
Scat2

Scat3 = Data1 %>% ggplot() + aes(x = STAI_trait, y = pain) + geom_point(shape = 21, fill = "red2", size = 2) + geom_smooth(method = "lm")
Scat3

Scat4 = Data1 %>% ggplot() + aes(x = pain_cat, y = pain) + geom_point(shape = 21, fill = "red2", size = 2) + geom_smooth(method = "lm")
Scat4

Scat5 = Data1 %>% ggplot() + aes(x = mindfulness, y = pain) + geom_point(shape = 21, fill = "red2", size = 2) + geom_smooth(method = "lm")
Scat5

Scat6 = Data1 %>% ggplot() + aes(x = cortisol_serum, y = pain) + geom_point(shape = 21, fill = "red2", size = 2) + geom_smooth(method = "lm")
Scat6

Scat7 = Data1 %>% ggplot() + aes(x = cortisol_saliva, y = pain) + geom_point(shape = 21, fill = "red2", size = 2) + geom_smooth(method = "lm")
Scat7

# Sccatter plots can visually indicate whether there is a positive or negative relationship with pain the predictive variable

grid.arrange(Scat1, Scat2, nrow = 2) # predictive variables in Regression model 1 (REG1)

grid.arrange(Scat1, Scat2, Scat3, Scat4, Scat5, Scat6, Scat7,  nrow = 3) # predictive variables in Regression model 2 (REG2) 




# STEP 4..................... ......Diagnostics 





# Regression model 1 (REG1)

#Identifying extreme cases with Cook's distance and high levarage (using residual - leverage plot) 

REG1 %>% plot(which = 4) # Cook's distance plot --> none of the values identified is likely to have problematic effect to the model - 
REG1 %>% plot(which = 5) # 3 cases of high residual error and high leverage - row 8,23,47 is identified -  (Residuals - Levevage plot)
#Cook's distance > 1 or values > 4/160= 0,25 is problematic --> the 3 cases of high residual error and high leverage is  evaluated from Cook's distance plot
Data1 %>%slice(c(8, 23, 47)) # data seems fine --> no exclusion

# Testing for normality of residuals with QQ plot
REG1 %>% plot(which = 2)
describe(residuals(REG1))  #skew and kurtosis --> Between +-1-,,,all good
hist(residuals(REG2), breaks=20, main="Residuals for regression model 2", density=10, angle=50, ylim=c(0,30), xlim=c(-5,5), ylab="Frequency", xlab="Residuals", las=1, border="gray20", col="gray70", labels=TRUE)

# Linearity 
REG1 %>%residualPlots() #Looks good

# Homoscedasticty
REG1 %>%plot(which = 3) 
REG1 %>%ncvTest() #Looks good 
REG1 %>% bptest() #Looks good 

# No multicollinearity
REG1 %>% vif() #looks good 



# Regression model 2 (REG2)



#Identifying extreme cases with Cook's distance and high levarage (using residual - leverage plot) 

REG2 %>% plot(which = 4) # Cook's distance plot --> none of the values identified is likely to have problematic effect to the model - 
REG2 %>% plot(which = 5) # 3 cases of high residual error and high leverage - row 46, 73, 85 is identified -  (Residuals - Levevage plot)
#Cook's distance > 1 or values > 4/160= 0,25 is problematic --> the 3 cases of high residual error and high leverage is  evaluated from Cook's distance plot
Data1 %>%slice(c(47, 74, 86)) # --> data seem good - no exclusion  

#testing for normality of residuals with QQ plot
REG2 %>% plot(which = 2)
describe(residuals(REG2))  #skew and kurtosis --> Between +-1-,,,all good
hist(residuals(REG2), breaks=20, main="Residuals for regression model 2", density=10, angle=50, ylim=c(0,30), xlim=c(-5,5), ylab="Frequency", xlab="Residuals", las=1, border="gray20", col="gray70", labels=TRUE)

#Linearity 
REG2 %>%residualPlots() #Looks good

#Homoscedasticty
REG2 %>%plot(which = 3) 
REG2 %>%ncvTest() #Looks good 
REG2 %>% bptest() #Looks good 

#No multicollinearity
REG2 %>%vif() # --> Cortisol Saliva and Cortisol_Serum have variance inflation factors # cortisol_serum = 4,79 - cortisol_saliva = 5,07
# cortisol_serum is thought to be a more reliable biological indicator of stress compared to cortisol_saliva --> Exclusion from model.



#............Rebuilding REG2 into --> REG_V2


REG2_V2 = lm(pain ~ age + sex +  STAI_trait + pain_cat + mindfulness + cortisol_serum, data = Data1)
summary(REG2_V2) # Multiple R-squared:  0.523,	Adjusted R-squared:  0.5041, F-statistic: 27.59, p-value: < 2.2e-16
confint(REG2_V2)

#                     2.5 %      97.5 %
# (Intercept)    -1.25909508  4.20679714
# age            -0.07834225 -0.00246924
# sexmale        -0.20302123  0.51418310
# STAI_trait     -0.05649792  0.03755381
# pain_cat        0.06715857  0.15550824
# mindfulness    -0.48608132 -0.06895363
# cortisol_serum  0.35243253  0.78286713

# Diagnosis of regression model 2 version 2

REG2_V2 %>% plot(which = 4) # Cook's distance plot --> none of the values identified is likely to have problematic effect to the model - 
REG2_V2 %>% plot(which = 5) # 3 cases of high residual error and high leverage - row 46, 73, 85 is identified -  (Residuals - Levevage plot)
#Cook's distance > 1 or values > 4/160= 0,25 is problematic --> the 3 cases of high residual error and high leverage is  evaluated from Cook's distance plot
Data1 %>%slice(c(47, 65, 86))

#testing for normality of residuals with QQ plot
REG2_V2 %>% plot(which = 2)
describe(residuals(REG2_V2))  #skew and kurtosis --> Between +-1-,,,all good
hist(residuals(REG2_V2), breaks=20, main="Residuals for regression model 2", density=10, angle=50, ylim=c(0,30), xlim=c(-5,5), ylab="Frequency", xlab="Residuals", las=1, border="gray20", col="gray70", labels=TRUE)

#Linearity 
REG2_V2 %>%residualPlots() #Looks good

#Homoscedasticty
REG2_V2 %>%plot(which = 3) 
REG2_V2 %>%ncvTest() #Looks good 
REG2_V2 %>% bptest() #Looks good 

#No multicollinearity
REG2_V2 %>%vif() #Looks good



# STEP 5................. Model evaluation 




# Regression model1 (age+sex as predictors of pain). 
REG1 <- lm(pain ~ sex + age, data = Data1)
summary(REG1)$adj.r.squared #Multiple R-squared:  0.08535,	Adjusted R-squared:  0.07355, F-statistic: 7.232,  p-value: 0.0009935
AIC(REG1)


# Regressionmodel2 with age, sex, STAI_trait, pain_catastrophizing, mindfulness, and cortisol measures as predictors of pain. 
REG2_V2 = lm(pain ~ age + sex +  STAI_trait + pain_cat + mindfulness + cortisol_serum, data = Data1)
summary(REG2_V2)$adj.r.squared
AIC(REG2_V2)


ANOVA_REG1_REG2_V2 = anova(REG1, REG2_V2)
ANOVA_REG1_REG2_V2 # F = 29.58

AIC(REG1) #574.1267
AIC(REG2_V2) #475.7032 

#REG2_V2 is the best model based on AIC and F values --> Significant difference between REG1 and REG2_V2  (98,42 points.! Adding additional predictors to the model improve the predictive abilities of the model

confint(REG1)
lm.beta(REG1)

Confint(REG2_V2)
lm.beta(REG2_V2)

SUMMARY_COF_table1 = coef_table(REG1)
SUMMARY_COF_table2 = coef_table(REG2_V2)

SUMMARY_COF_table1
#b 95%CI lb 95%CI ub Std.Beta p-value
#(Intercept)  8.24     6.35    10.14        0   <.001
#sexmale      0.30    -0.17     0.76      0.1    .208
#age         -0.09    -0.13    -0.04    -0.28   <.001

SUMMARY_COF_table2
#b 95%CI lb 95%CI ub Std.Beta p-value
#(Intercept)     1.47    -1.26     4.21        0    .288
#age            -0.04    -0.08     0.00    -0.13    .037
#sexmale         0.16    -0.20     0.51     0.05    .393
#STAI_trait     -0.01    -0.06     0.04    -0.03    .691
#pain_cat        0.11     0.07     0.16     0.38   <.001
#mindfulness    -0.28    -0.49    -0.07    -0.18    .009
#cortisol_serum  0.57     0.35     0.78     0.34   <.001


#Equation - REG2_V2: 𝑌 = 𝑏0 + 𝑏1 ∗ X1 + 𝑏2 ∗ X2 +...+ bn * Xn,
# y(pain)= 1.47 + (-0.04)*age + (0.15)*sex + (-0.009)*STAI_trait + (0.11)*pain_cat + (-0.28)*mindfulness + (0.57)*cortisol_serum 

