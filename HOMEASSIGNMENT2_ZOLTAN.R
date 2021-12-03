                                                                    ##HOMEASSIGNMENT 1 (Zoltan)

#######_______Initial preparations - loading of packages and custom functions from course literature (exercise 1-14) for specific operations within the assignments__#

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

# STEP 1...........................Retrieving data  

Data = read_csv("https://tinyurl.com/yxm5rd89")




# STEP 2.............................Inspecting data and correcting errors

view(data)

describe(Data) 
summary(Data)

Data1 <- Data[-c(88,34), ] 
view(Data1) 

# STEP3.......................Building research model 

REG3 <- lm( pain ~ age + sex + STAI_trait + pain_cat + cortisol_serum +  mindfulness + weight + IQ + household_income, data = Data1)
REG3
summary(REG3) 
confint(REG3)



# STEP 4........................MODEL Diagnostics 




REG3 %>% plot(which = 4) # Cook's distance plot --> none of the values identified is likely to have problematic effect to the model - 
REG3 %>% plot(which = 5) # 3 cases of high residual error and high leverage - row 46,84,85 is identified -  (Residuals - Levevage plot)
#Cook's distance > 1 or values > 4/160= 0,25 is problematic --> the 3 cases of high residual error and high leverage is  evaluated from Cook's distance plot
Data1 %>%slice(c(46, 84, 85)) # data seems fine --> no exclusion

# Testing for normality of residuals with QQ plot and historigram 
REG3 %>% plot(which = 2)
describe(residuals(REG3))  #skew and kurtosis --> Between +-1-,,,all good
hist(residuals(REG3), breaks=20, main="Residuals for regression model 3", density=10, angle=50, ylim=c(0,30), xlim=c(-5,5), ylab="Frequency", xlab="Residuals", las=1, border="gray20", col="gray70", labels=TRUE)

# Linearity 
REG3%>% residualPlots() #Looks good 

#3.Homoscedasticity 
REG3%>% plot(which=3) #Looks good 
REG3 %>%ncvTest() #Looks good 
REG3w %>%bptest() #Looks good 

#4.No multicolinealrity 
REG3 %>% vif() #Looks good 
#------------------------
REG3 %>% summary()
AIC(REG3) #484,01




# STEP 5.......................Backwards Regression




BACK_REG3 = step(REG3, direction= "backward") 
summary(BACK_REG3)

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)     1.27627    1.33113   0.959   0.3392    
#age            -0.04116    0.01819  -2.263   0.0250 *  
#  pain_cat        0.11359    0.02038   5.573 1.10e-07 ***
#  cortisol_serum  0.53383    0.09982   5.348 3.18e-07 ***
#  mindfulness    -0.26852    0.10478  -2.563   0.0114 *  
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 1.069 on 153 degrees of freedom
#Multiple R-squared:  0.5199,	Adjusted R-squared:  0.5073 
#F-statistic: 41.41 on 4 and 153 DF,  p-value: < 2.2e-16

# Rebuilding backwards model, without variables indicated

BACK_REG3_V2= lm(pain~ age + pain_cat + cortisol_serum + mindfulness, data= Data1)
BACK_REG3_V2
summary(BACK_REG3_V2)
confint(BACK_REG3_V2)


# Diagnosis 



BACK_REG3_V2 %>% plot(which = 4) # Cook's distance plot --> none of the values identified is likely to have problematic effect to the model - 
BACK_REG3_V2 %>% plot(which = 5) # 3 cases of high residual error and high leverage - row 46,102,115 is identified -  (Residuals - Levevage plot)
#Cook's distance > 1 or values > 4/160= 0,25 is problematic --> the 3 cases of high residual error and high leverage is  evaluated from Cook's distance plot
Data1 %>%slice(c(46, 102, 115)) # data seems fine --> no exclusion

# Testing for normality of residuals with QQ plot and historigram 
BACK_REG3_V2 %>% plot(which = 2)
describe(residuals(BACK_REG3_V2))  #skew and kurtosis --> Between +-1-,,,all good
hist(residuals(BACK_REG3_V2), breaks=20, main="Residuals for regression model 3", density=10, angle=50, ylim=c(0,30), xlim=c(-5,5), ylab="Frequency", xlab="Residuals", las=1, border="gray20", col="gray70", labels=TRUE)

# Linearity 
BACK_REG3_V2%>% residualPlots() #Looks good 

#3.Homoscedasticity 
BACK_REG3_V2 %>%plot(which=3) #Looks good 
BACK_REG3_V2 %>%ncvTest() #Looks good 
BACK_REG3_V2 %>%bptest() #Looks good 

#4.No multicolinealrity 
BACK_REG3_V2 %>% vif() #Looks good 
#----------------------------------
BACK_REG3_V2 %>% summary()
AIC(BACK_REG3_V2) #534.73




#..................STEP 6 comparision of models.............

# Retrieve dataset 2 
Data2 = read.csv("https://tinyurl.com/87v6emky")


REG2_V2 = lm(pain ~ age + sex +  STAI_trait + pain_cat + mindfulness + cortisol_serum, data = Data2)
REG2_V2
#Equation - REG2_V2: 𝑌 = 𝑏0 + 𝑏1 ∗ X1 + 𝑏2 ∗ X2 +...+ bn * Xn,
# y(pain)= 1.51 + (-0.05)*age + (0.44)*sex + (-0.02)*STAI_trait + (0.11)*pain_cat + (-0.15)*mindfulness + (0.64)*cortisol_serum 

BACK_REG3_V2= lm(pain~ age + pain_cat + cortisol_serum + mindfulness, data= Data2)
BACK_REG3_V2
#Equation - REG2_V2: 𝑌 = 𝑏0 + 𝑏1 ∗ X1 + 𝑏2 ∗ X2 +...+ bn * Xn,
# y(pain)= 1.42 + (-0.05)*age + (0.11)*pain_cat + (0.57) (-0.14)*mindfulness 

summary(REG2_V2)
summary(BACK_REG3_V2)

summary(REG2_V2)$adj.r.squared # =0.451
summary(BACK_REG3_V2)$adj.r.squared # =0.438

AIC(REG2_V2) #532.79
AIC(BACK_REG3_V2) #534.74
AIC(REG2_V2)-AIC(BACK_REG3_V2) # 532.79 - 534.73 = -1.94
532.79 - 534.73

anova(REG2_V2, BACK_REG3_V2)
#Model 1: pain ~ age + sex + STAI_trait + pain_cat + mindfulness + cortisol_serum
#Model 2: pain ~ age + pain_cat + cortisol_serum + mindfulness
#Res.Df    RSS Df Sum of Sq      F  Pr(>F)  
#1    153 236.81                              
#2    155 245.77 -2   -8.9648 2.8961 0.05827 .
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# Predictive evaluation 

Theory = predict(REG2_V2, newdata = Data2) 
Data2 <- cbind(Data2, Theory) #Adding predicted values as a new column to data_sample_2

Backward = predict(BACK_REG3_V2, newdata = Data2) #Running predictions on new dataset using the backward model
Data2 <- cbind(Data2, Backward) #Adding predicted values as a new column to data_sample_2

View(Data2) #Inspecting the data and checking to see if predicted values were added to the dataset

#Visual inspection of prediction accuracy 
BACK_Pre_scatter = Data2 %>% 
  ggplot() +
  aes(x = Backward, y = pain) + geom_point() 
BACK_Pre_scatter

Theory_Pre_scatter = Data2 %>%
  ggplot() +
  aes(x = Theory, y = pain) + geom_point() 
Theory_Pre_scatter

RSS_theory = sum((Data2$pain - Data2$Theory)^2) 
RSS_backward = sum((Data2$pain - Data2$Backward)^2) 

RSS_theory # = 236.8
RSS_backward # =245.77

summary(REG2_V2) # Adjusted R-squared:  0.4513, F-statistic:  22.8 on 6 and 153 DF,  p-value: < 2.2e-16
summary(BACK_REG3_V2) #Adjusted R-squared:  0.4379, F-statistic: 31.96 on 4 and 155 DF,  p-value: < 2.2e-16



Table_REG2_V2 = coef_table(REG2_V2)
summary(REG2_V2)

Table_BACK_REG3_V2 = coef_table(BACK_REG3_V2)
summary(BACK_REG3_V2)



