#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x##x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x
# R code example for shrimp catch reconstruction
# Example of GLMs fit for Sao Paulo state catches
# The construction of coefficient of variation, means, and amplitude indices is omitted in this file
# Author: Matheus Lourenço
#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x

# packages
# install.packages("lmtest") # just lmtest package is required
library(lmtest)
rm(list = ls()) #clean environment

# Change the directory to the folder that contains the script (R_code_Reconstruction_example.R) and the spreadsheet (Data_SP.csv)
setwd("C:/Matheus/Universidade/Doutorado") # set working directory to source file location

## Reading the dataset
data_SP<-read.csv("data_SP.csv",sep=",", dec=".")
#------------------------------------------------------------------------------------------

#------------
# pink
#------------

# Base years for reconstruction # Square root transformation
data_rosa_SP<- sqrt(data_SP[c(33:62),c(1:4,5,8:19,2)]) #Remove the influence of extreme values
row.names(data_rosa_SP)<-NULL

#Model with all variables: Inverse Gaussian identity
fitall<-glm(c.rosa~ tsm_med_chuv_SP+
              c.total+tsm_med_SP+
              tsm_med_seco_SP+
              prec_med_SP+
              prec_med_seco_SP+
              prec_med_chuv_SP+
              prec_amp_seco_SP+
              tsm_coef_SP+
              prec_coef_SP+
              tsm_amp_seco_SP+
              prec_amp_chuv_SP+
              tsm_amp_chuv_SP+
              sinano+ano2,
            family = inverse.gaussian(link="identity"), 
            data = data_rosa_SP)

fitini<-glm(c.rosa ~1,
            family = inverse.gaussian(link = "identity"), 
            data = data_rosa_SP)

step(fitini,direction = "forward",scope = formula(fitall))

mod_rosa_gauss_iden<-glm(formula = c.rosa ~ 
                           c.total + sinano + 
                           tsm_amp_seco_SP, 
                         family = inverse.gaussian(link = "identity"), 
                         data = data_rosa_SP)
AIC(mod_rosa_gauss_iden) #159.18


#Model with all variables: Inverse Gaussian log
fitall<-glm(c.rosa~ tsm_med_chuv_SP+
              c.total+tsm_med_SP+
              tsm_med_seco_SP+
              prec_med_SP+
              prec_med_seco_SP+
              prec_med_chuv_SP+
              prec_amp_seco_SP+
              tsm_coef_SP+
              prec_coef_SP+
              tsm_amp_seco_SP+
              prec_amp_chuv_SP+
              tsm_amp_chuv_SP+
              sinano+ano2,
            family = inverse.gaussian(link="log"), 
            data = data_rosa_SP)

fitini<-glm(c.rosa ~1,
            family = inverse.gaussian(link = "log"), 
            data = data_rosa_SP)

step(fitini,direction = "forward",scope = formula(fitall))

mod_rosa_gauss_log<-glm(formula = c.rosa ~ c.total +
                          sinano + ano2,
                        family = inverse.gaussian(link = "log"), 
                        data = data_rosa_SP)

AIC(mod_rosa_gauss_log) #167.26


#Model with all variables: Gamma identity
fitall<-glm(c.rosa~ tsm_med_chuv_SP+
              c.total+tsm_med_SP+
              tsm_med_seco_SP+
              prec_med_SP+
              prec_med_seco_SP+
              prec_med_chuv_SP+
              prec_amp_seco_SP+
              tsm_coef_SP+
              prec_coef_SP+
              tsm_amp_seco_SP+
              prec_amp_chuv_SP+
              tsm_amp_chuv_SP+
              sinano+ano2,
            family = Gamma(link = "identity"), 
            data = data_rosa_SP)

fitini<-glm(c.rosa ~1,
            family = Gamma(link = "identity"), 
            data = data_rosa_SP)

step(fitini,direction = "forward",scope = formula(fitall))


mod_rosa_gamma_iden<-glm(formula = c.rosa ~ c.total + 
                           sinano, 
                         family = Gamma(link = "identity"), 
                         data = data_rosa_SP)
AIC(mod_rosa_gamma_iden) #165.23


#Model with all variables: Gamma log
fitall<-glm(c.rosa~ tsm_med_chuv_SP+
              c.total+tsm_med_SP+
              tsm_med_seco_SP+
              prec_med_SP+
              prec_med_seco_SP+
              prec_med_chuv_SP+
              prec_amp_seco_SP+
              tsm_coef_SP+
              prec_coef_SP+
              tsm_amp_seco_SP+
              prec_amp_chuv_SP+
              tsm_amp_chuv_SP+
              sinano+ano2,
            family = Gamma(link = "log"), 
            data = data_rosa_SP)

fitini<-glm(c.rosa ~1,
            family = Gamma(link = "log"), 
            data = data_rosa_SP)

step(fitini,direction = "forward",scope = formula(fitall))

mod_rosa_gamma_log<- glm(formula = c.rosa ~ 
                           c.total + sinano + 
                           prec_med_chuv_SP, 
                         family = Gamma(link = "log"), 
                         data = data_rosa_SP)

AIC(mod_rosa_gamma_log) #169.90


#Best models
par(mfrow=c(2,2))
plot(mod_rosa_gauss_iden)

#Residual homoscedasticity test
require("lmtest")
bptest(mod_rosa_gauss_iden)

#Normality of residuals
shapiro.test(resid(mod_rosa_gauss_iden))

par(mfrow=c(2,2))
plot(mod_rosa_gamma_iden)

#Homoscedasticity test of residuals
require("lmtest")
bptest(mod_rosa_gamma_iden)

#Normality of residuals
shapiro.test(resid(mod_rosa_gamma_iden))

#Best model (pink) #removing leverage points
best_mod_rosa<-glm(formula = c.rosa ~ 
                     c.total + sinano + 
                     tsm_amp_seco_SP, 
                   family = inverse.gaussian(link = "identity"), 
                   data = data_rosa_SP[c(-13),])
par(mfrow=c(2,2))
plot(best_mod_rosa) #AIC 148.2

#Residual homoscedasticity test
require("lmtest")
bptest(best_mod_rosa) 

#Residual normality
shapiro.test(resid(best_mod_rosa))
#-----------------------------------------------------------------------------------

#------------
#White
#------------

#Base years for reconstruction # Square root transformation
data_branco_SP<-sqrt(data_SP[c(35:45,51:62),c(1:4,7,8:19,2)]) 
row.names(data_branco_SP)<-NULL


#Model with all variables: Inverse Gaussian identity
fitall<-glm(c.branco~ tsm_med_chuv_SP+
              c.total+tsm_med_SP+
              tsm_med_seco_SP+
              prec_med_SP+
              prec_med_seco_SP+
              prec_med_chuv_SP+
              prec_amp_seco_SP+
              tsm_coef_SP+
              prec_coef_SP+
              tsm_amp_seco_SP+
              prec_amp_chuv_SP+
              tsm_amp_chuv_SP+
              sinano+ano2,
            family = inverse.gaussian(link="identity"), 
            data = data_branco_SP)

fitini<-glm(c.branco ~1,
            family = inverse.gaussian(link = "identity"), 
            data = data_branco_SP)

step(fitini,direction = "forward",scope = formula(fitall))

mod_branco_gauss_iden<-glm(formula = c.branco ~ c.total + 
                             prec_med_chuv_SP, 
                           family = inverse.gaussian(link = "identity"), 
                           data = data_branco_SP)

AIC(mod_branco_gauss_iden) #108.95


#Model with all variables: Inverse Gaussian log
fitall<-glm(c.branco~ tsm_med_chuv_SP+
              c.total+tsm_med_SP+
              tsm_med_seco_SP+
              prec_med_SP+
              prec_med_seco_SP+
              prec_med_chuv_SP+
              prec_amp_seco_SP+
              tsm_coef_SP+
              prec_coef_SP+
              tsm_amp_seco_SP+
              prec_amp_chuv_SP+
              tsm_amp_chuv_SP+
              sinano+ano2,
            family = inverse.gaussian(link="log"), 
            data = data_branco_SP)

fitini<-glm(c.branco ~1,
            family = inverse.gaussian(link = "log"), 
            data = data_branco_SP)

step(fitini,direction = "forward",scope = formula(fitall))

mod_branco_gauss_log<-glm(formula = c.branco ~ 
                            c.total + prec_med_chuv_SP +
                            prec_amp_chuv_SP + 
                            sinano, 
                          family = inverse.gaussian(link = "log"),
                          data = data_branco_SP)
AIC(mod_branco_gauss_log) #110.59


#Model with all variables: Gamma identity
fitall<-glm(c.branco~ tsm_med_chuv_SP+
              c.total+tsm_med_SP+
              tsm_med_seco_SP+
              prec_med_SP+
              prec_med_seco_SP+
              prec_med_chuv_SP+
              prec_amp_seco_SP+
              tsm_coef_SP+
              prec_coef_SP+
              tsm_amp_seco_SP+
              prec_amp_chuv_SP+
              tsm_amp_chuv_SP+
              sinano+ano2,
            family = Gamma(link = "identity"), 
            data = data_branco_SP)

fitini<-glm(c.branco ~1,
            family = Gamma(link = "identity"), 
            data = data_branco_SP)

step(fitini,direction = "forward",scope = formula(fitall))

mod_branco_gamma_iden<-glm(formula = c.branco ~ c.total +
                             prec_med_chuv_SP,
                           family = Gamma(link = "identity"), 
                           data = data_branco_SP)
AIC(mod_branco_gamma_iden) #109.32


#Model with all variables: Gamma log
fitall<-glm(c.branco~ tsm_med_chuv_SP+
              c.total+tsm_med_SP+
              tsm_med_seco_SP+
              prec_med_SP+
              prec_med_seco_SP+
              prec_med_chuv_SP+
              prec_amp_seco_SP+
              tsm_coef_SP+
              prec_coef_SP+
              tsm_amp_seco_SP+
              prec_amp_chuv_SP+
              tsm_amp_chuv_SP+
              sinano+ano2,
            family = Gamma(link = "log"), 
            data = data_branco_SP)

fitini<-glm(c.branco ~1,
            family = Gamma(link = "log"), 
            data = data_branco_SP)

step(fitini,direction = "forward",scope = formula(fitall))

mod_branco_gamma_log<- glm(formula = c.branco ~ c.total +
                             ano2 + tsm_med_SP, 
                           family = Gamma(link = "log"), 
                           data = data_branco_SP)

AIC(mod_branco_gamma_log) #110.35


#Best models
par(mfrow=c(2,2))
plot(mod_branco_gauss_iden)

#Homoscedasticity test of residuals
require("lmtest")
bptest(mod_branco_gauss_iden)

#Normality of residuals
shapiro.test(resid(mod_branco_gauss_iden))

par(mfrow=c(2,2))
plot(mod_branco_gamma_iden)

#Homoscedasticity test of residuals
require("lmtest")
bptest(mod_branco_gamma_iden)

#Normality of residuals
shapiro.test(resid(mod_branco_gamma_iden))

#Best model (white) - Removing leverage points
best_mod_branco<- glm(formula = c.branco ~ c.total +
                        prec_med_chuv_SP,
                      family = Gamma(link = "identity"), 
                      data = data_branco_SP[c(-4),])
par(mfrow=c(2,2))
plot(best_mod_branco) #AIC 101.8

#Residual homoscedasticity test
require("lmtest")
bptest(best_mod_branco)

#Residual normality
shapiro.test(resid(best_mod_branco))



#------------
# Seabob
#------------

#Base years for reconstruction #square root transformation
data_7barbas_SP<-sqrt(data_SP[c(33:62),c(1:4,6,8:19,2)]) 
row.names(data_7barbas_SP)<-NULL


#Model with all variables using Inverse Gaussian distribution with identity link
fitall<-glm(c.7barbas~ tsm_med_chuv_SP+
              c.total+tsm_med_SP+
              tsm_med_seco_SP+
              prec_med_SP+
              prec_med_seco_SP+
              prec_med_chuv_SP+
              prec_amp_seco_SP+
              tsm_coef_SP+
              prec_coef_SP+
              tsm_amp_seco_SP+
              prec_amp_chuv_SP+
              tsm_amp_chuv_SP+
              sinano+ano2,
            family = inverse.gaussian(link="identity"), 
            data = data_7barbas_SP)

fitini<-glm(c.7barbas ~1,
            family = inverse.gaussian(link = "identity"), 
            data = data_7barbas_SP)

step(fitini,direction = "forward",scope = formula(fitall))

mod_7barbas_gauss_iden<-glm(formula = c.7barbas ~ c.total + 
                              sinano + prec_med_chuv_SP + 
                              tsm_amp_seco_SP +
                              tsm_med_chuv_SP, 
                            family = inverse.gaussian(link = "identity"), 
                            data = data_7barbas_SP)

AIC(mod_7barbas_gauss_iden) #142.36


#Model with all variables using Inverse Gaussian distribution with log link
fitall<-glm(c.7barbas~ tsm_med_chuv_SP+
              c.total+tsm_med_SP+
              tsm_med_seco_SP+
              prec_med_SP+
              prec_med_seco_SP+
              prec_med_chuv_SP+
              prec_amp_seco_SP+
              tsm_coef_SP+
              prec_coef_SP+
              tsm_amp_seco_SP+
              prec_amp_chuv_SP+
              tsm_amp_chuv_SP+
              sinano+ano2,
            family = inverse.gaussian(link="log"), 
            data = data_7barbas_SP)

fitini<-glm(c.7barbas ~1,
            family = inverse.gaussian(link = "log"), 
            data = data_7barbas_SP)

step(fitini,direction = "forward",scope = formula(fitall))

mod_7barbas_gauss_log<-glm(formula = c.7barbas ~ c.total + 
                             prec_amp_chuv_SP + ano2 + 
                             sinano +
                             tsm_amp_chuv_SP,
                           family = inverse.gaussian(link = "log"), 
                           data = data_7barbas_SP)

AIC(mod_7barbas_gauss_log) #147.56



#Model with all variables using Gamma distribution with identity link
fitall<-glm(c.7barbas~ tsm_med_chuv_SP+
              c.total+tsm_med_SP+
              tsm_med_seco_SP+
              prec_med_SP+
              prec_med_seco_SP+
              prec_med_chuv_SP+
              prec_amp_seco_SP+
              tsm_coef_SP+
              prec_coef_SP+
              tsm_amp_seco_SP+
              prec_amp_chuv_SP+
              tsm_amp_chuv_SP+
              sinano+ano2,
            family = Gamma(link= "identity"), 
            data = data_7barbas_SP)

fitini<-glm(c.7barbas ~1,
            family = Gamma(link= "identity"), 
            data = data_7barbas_SP)

step(fitini,direction = "forward",scope = formula(fitall))

mod_7barbas_gamma_iden<- glm(formula = c.7barbas ~ c.total + 
                               sinano + prec_med_chuv_SP, 
                             family = Gamma(link = "identity"),
                             data = data_7barbas_SP)
AIC(mod_7barbas_gamma_iden) #141.91


#Model with all variables using Gamma distribution with log link
fitall<-glm(c.7barbas~ tsm_med_chuv_SP+
              c.total+tsm_med_SP+
              tsm_med_seco_SP+
              prec_med_SP+
              prec_med_seco_SP+
              prec_med_chuv_SP+
              prec_amp_seco_SP+
              tsm_coef_SP+
              prec_coef_SP+
              tsm_amp_seco_SP+
              prec_amp_chuv_SP+
              tsm_amp_chuv_SP+
              sinano+ano2,
            family = Gamma(link= "log"), 
            data = data_7barbas_SP)

fitini<-glm(c.7barbas ~1,
            family = Gamma(link= "log"), 
            data = data_7barbas_SP)

step(fitini,direction = "forward",scope = formula(fitall))

mod_7barbas_gamma_log<- glm(formula = c.7barbas ~ c.total +
                              ano2 + sinano + prec_amp_chuv_SP, 
                            family = Gamma(link = "log"), 
                            data = data_7barbas_SP)

AIC(mod_7barbas_gamma_log) #143.95


#Best models
par(mfrow=c(2,2))
plot(mod_7barbas_gauss_iden)

#Residual homoscedasticity test
require("lmtest")
bptest(mod_7barbas_gauss_iden)

#Residual normality assessment
shapiro.test(resid(mod_7barbas_gauss_iden))

par(mfrow=c(2,2))
plot(mod_7barbas_gamma_iden)

#Residual homoscedasticity test
require("lmtest")
bptest(mod_7barbas_gamma_iden)

#Normality of residuals
shapiro.test(resid(mod_7barbas_gamma_iden))

#Best model seabob: removing leverage points
best_mod_7barbas<-glm(formula = c.7barbas ~ c.total + 
                        sinano + prec_med_chuv_SP, 
                      family = Gamma(link = "identity"),
                      data = data_7barbas_SP[c(-13),])
par(mfrow=c(2,2))
plot(best_mod_7barbas) #AIC 119.9

#Homoscedasticity test of residuals
require("lmtest")
bptest(best_mod_7barbas)

#Normality of residuals
shapiro.test(resid(best_mod_7barbas))
#------------------------------------------------------------------------------------

#----------------------------
#Predictions with the models
#----------------------------
pred_frame <-sqrt(data_SP) # data.frame prediction


#Predictions with the pink model
pred_modrosa<-data.frame(predict(best_mod_rosa, newdata=pred_frame,interval="prediction", 
                                 level=0.95, type = "response",se.fit = TRUE)) # IC 95%
pred_modrosa$up<-pred_modrosa$fit+qnorm(0.975)*pred_modrosa$se.fit
pred_modrosa$lw<-pred_modrosa$fit-qnorm(0.975)*pred_modrosa$se.fit

pred_modrosa<- pred_modrosa^2 #Converting back to the original scale

plot(pred_frame$ano,pred_modrosa$fit,type="l",ylim = c(0,max(pred_modrosa$up)))
par(new=T)

plot(pred_frame$ano,pred_modrosa$up,type="l",lty=2,ylim = c(0,max(pred_modrosa$up)))
par(new=T)

plot(pred_frame$ano,pred_modrosa$lw,type="l",lty=2,ylim = c(0,max(pred_modrosa$up)))

par(new=T)

plot(data_SP$ano,data_SP$c.rosa,type="l",lwd=2,lty=1,ylim = c(0,max(pred_modrosa$up)))
dev.off()
#------------------------------------------------------------------------------------

#Predictions with the white model

pred_modbranco<-data.frame(predict(best_mod_branco, newdata=pred_frame,interval="prediction", 
                                   level=0.95,type = "response", se.fit = TRUE)) # IC 95%
pred_modbranco$up<-pred_modbranco$fit+qnorm(0.975)*pred_modbranco$se.fit
pred_modbranco$lw<-pred_modbranco$fit-qnorm(0.975)*pred_modbranco$se.fit
pred_modbranco<- pred_modbranco^2 #Converting back to the original scale


plot(data_SP$ano,pred_modbranco$fit,type="l",ylim = c(0,max(pred_modbranco$fit)))
par(new=T)

plot(data_SP$ano,pred_modbranco$up,type="l",lty=2,ylim = c(0,max(pred_modbranco$fit)))
par(new=T)

plot(data_SP$ano,pred_modbranco$lw,type="l",lty=2,ylim = c(0,max(pred_modbranco$fit)))

par(new=T)

plot(data_SP$ano,data_SP$c.branco,type="l",lwd=2,lty=1,ylim = c(0,max(pred_modbranco$fit)))
dev.off()
#-----------------------------------------------------------------------------------
#Predictions with the seabob model

pred_mod7barbas<-data.frame(predict(best_mod_7barbas, newdata=pred_frame,interval="prediction", 
                                    level=0.95,type = "response", se.fit = TRUE)) # IC 95%
pred_mod7barbas$up<-pred_mod7barbas$fit+qnorm(0.975)*pred_mod7barbas$se.fit
pred_mod7barbas$lw<-pred_mod7barbas$fit-qnorm(0.975)*pred_mod7barbas$se.fit
pred_mod7barbas<- pred_mod7barbas^2 #Converting back to the original scale


plot(data_SP$ano,pred_mod7barbas$fit,type="l",ylim = c(0,max(pred_mod7barbas$up)))
par(new=T)

plot(data_SP$ano,pred_mod7barbas$up,type="l",lty=2,ylim = c(0,max(pred_mod7barbas$up)))
par(new=T)

plot(data_SP$ano,pred_mod7barbas$lw,type="l",lty=2,ylim = c(0,max(pred_mod7barbas$up)))

par(new=T)

plot(data_SP$ano,data_SP$c.7barbas,type="l",lwd=2,lty=1,ylim = c(0,max(pred_mod7barbas$up)))
dev.off()
#-----------------------------------------------------------------------------------

summary(best_mod_rosa) #line 13 
data_rosa_SP<-(data_rosa_SP^2)

#Plotting pink Reconstruction
{jpeg("Serie_recons_rosa_SP.jpg",width=15,height=10,units="cm",res=250,quality=120)
  par(mar=c(4,4,0.5,0), cex.lab=1.2,cex.axis=1.2,mfrow=c(1,1))
  plot(data_SP$ano, rep(0,length(data_SP$ano)),type="n", ylim=c(0,2800), 
       ylab="Catch (t)",xlab="Year",xaxt="n", bty="n",xlim = c(1946,2011)) 
  # polygon
  polygon(c(data_SP$ano[c(1:66)],rev(data_SP$ano[c(1:66)])), 
          c(pred_modrosa$lw[c(1:66)],rev(pred_modrosa$up[c(1:66)])),
          col="tomato1", border=NA)
  #mean line
  lines(data_SP$ano[1:66],abs(pred_modrosa$fit[1:66]),lwd=2)
  
  
  library(plotrix)
  draw.circle(data_rosa_SP$ano[c(8)],
              data_rosa_SP$c.rosa[c(8)],radius = 1,lwd=2)
  library(plotrix)
  draw.circle(data_rosa_SP$ano[c(13)],
              data_rosa_SP$c.rosa[c(13)],radius = 1,lwd=2)
  
  #Reported Capture Used
  lines(data_SP$ano[33:62],data_SP$c.rosa[33:62],lwd=2,lty=1,col="blue")
  
  # legend(x=1974.6,y=0.98*(max(pred_modrosa$up)),fill = c("tomato1"),lty=1,merge = TRUE,legend = c("Reconstru?da"),
  #        bty="n",cex=1.1, seg.len = 0.65,lwd=2)
  # 
  # legend(x=1975,y=0.9*(max(pred_modrosa$up)),col = c("dimgray"),lty=1,merge = TRUE,
  #        legend = c("Reportada incorporada"),
  #        bty="n",cex=1.1, seg.len = 1.2,lwd=2)
  # 
  # legend(x=1975,y=0.90*(max(pred_modrosa$up)),col = c("black"),lty=2,merge = TRUE,
  #        legend = c("Reportada n?o incorporada"),
  #        bty="n",cex=1.1, seg.len = 1.2,lwd=2)
  # 
  axis(1,seq(1946,2011,4))
  dev.off()
}
#------------------------------------------------------------------------------------------------------

summary(best_mod_branco) #line 4
data_branco_SP<-(data_branco_SP^2) 

#Plotting White Reconstruction
{jpeg("Serie_recons_branco_SP.jpg",width=15,height=10,units="cm",res=250,quality=120)
  
  par(mar=c(4,4,0,0), cex.lab=1.2,cex.axis=1.2,mfrow=c(1,1))
  plot(data_SP$ano, rep(0,length(data_SP$ano)),type="n", ylim=c(0,max(pred_modbranco$up)), 
       ylab="Catch (t)",xlab="year",xaxt="n", bty="n",xlim = c(1946,2011)) 
  # polygon
  polygon(c(data_SP$ano[1:66],rev(data_SP$ano[1:66])),
          c(pred_modbranco$lw[1:66],rev(pred_modbranco$up[1:66])),
          col="cyan2", border=NA)
  #mean line
  lines(data_SP$ano[1:66],(abs(pred_modbranco$fit[1:66])),lwd=2)
  
  #Incorporated Reported Catch
  #points(1981,4,pch=20,col="black")
  lines(data_SP$ano[c(35:49)],data_SP$c.branco[c(35:49)],lwd=2,lty=1,col="blue")
  lines(data_SP$ano[c(51:62)],data_SP$c.branco[c(51:62)],lwd=2,lty=1,col="blue")
  
  library(plotrix)
  draw.circle(data_branco_SP$ano[c(4)],
              data_branco_SP$c.branco[c(4)],radius = 1,lwd=2)
  
  draw.circle(data_SP$ano[c(47)],
              data_SP$c.branco[c(47)],radius = 2.5,lwd=2)
  
  # legend(x=1974.6,y=1.1*(max(pred_modbranco$up)),fill = c("cyan2"),lty=1,merge = TRUE,legend = c("Reconstru?da"),
  #        bty="n",cex=1.1, seg.len = 0.65,lwd=2)
  # 
  # legend(x=1975,y=1*(max(pred_modbranco$up)),col = c("dimgray"),lty=1,merge = TRUE,
  #        legend = c("Reportada incorporada"),
  #        bty="n",cex=1.1, seg.len = 1.2,lwd=2)
  # 
  # legend(x=1975,y=0.9*(max(pred_modbranco$up)),col = c("dimgray"),lty=2,merge = TRUE,
  #        legend = c("Reportada n?o incorporada"),
  #        bty="n",cex=1.1, seg.len = 1.2,lwd=2)
  
  axis(1,seq(1946,2011,4))
  dev.off()
}
#----------------------------------------------------------------------------------------------------

summary(best_mod_7barbas) #line 13
data_7barbas_SP<-(data_7barbas_SP^2) 

#Plotting seabob Reconstruction
{jpeg("Serie_recons_7barbas_SP.jpg",width=15,height=10,units="cm",res=250,quality=120)
  
  par(mar=c(4,4,0.8,0), cex.lab=1.2,cex.axis=1.2,mfrow=c(1,1))
  plot(data_SP$ano, rep(0,length(data_SP$ano)),type="n", ylim=c(0,max(pred_mod7barbas$up)), 
       ylab="Catch (t)",xlab="Year",xaxt="n", bty="n",xlim = c(1946,2011)) 
  # polygon
  polygon(c(data_SP$ano[1:66],rev(data_SP$ano[1:66])), 
          c(pred_mod7barbas$lw[1:66],rev(pred_mod7barbas$up[1:66])),
          col="springgreen2", border=NA)
  # mean line
  lines(data_SP$ano[1:66],abs((pred_mod7barbas$fit[1:66])),lwd=2)
  
  #Reported catch incorporated
  lines(data_SP$ano[c(33:62)],data_SP$c.7barbas[c(33:62)],lwd=2,lty=1,col="blue")
  
  
  draw.circle(data_7barbas_SP$ano[c(13)],
              data_7barbas_SP$c.7barbas[c(13)],radius = 1,lwd=2)
  
  # legend(x=1970.6,y=1.1*(max(pred_mod7barbas$up)),fill = c("springgreen2"),lty=1,merge = TRUE,legend = c("Reconstru?da"),
  #        bty="n",cex=1.1, seg.len = 0.65,lwd=2)
  # 
  # legend(x=1971,y=1*(max(pred_mod7barbas$up)),col = c("dimgray"),lty=1,merge = TRUE,
  #        legend = c("Reportada incorporada"),
  #        bty="n",cex=1.1, seg.len = 1.2,lwd=2)
  # 
  # legend(x=1945,y=0.88*(max(pred_mod7barbas$up)),col = c("black"),lty=2,merge = TRUE,
  #        legend = c("Reportada n?o incorporada"),
  #        bty="n",cex=1.1, seg.len = 1.2,lwd=2)
  axis(1,seq(1946,2011,4))
  dev.off()
}

#------------------------------------------------------------------------------------
#Adding reconstructed pink series to data
data_SP$c.rosa.recon <- (pred_modrosa$fit)
#------------------------------------------------------------------------------------

#Adding reconstructed white series to data
data_SP$c.branco.recon<-(pred_modbranco$fit)
#------------------------------------------------------------------------------------

#Adding reconstructed seabob series to data
data_SP$c.7barbas.recon <-(pred_mod7barbas$fit)
#------------------------------------------------------------------------------------

#Adding reconstructed series to data and an output data.frame
# outfile  <- "data_SP.csv"
# write.table(data_SP, file = outfile, append = TRUE,dec=".",sep = ",",
#             row.names = FALSE, col.names =!file.exists(outfile)) 
#------------------------------------------------------------------------------------

#End of reconstructions
#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x#x

