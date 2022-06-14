require(mgcv)
require(DHARMa)
require(stringr)
require(Metrics)
require(CombMSC)
require(rgdal)
require(corrplot)

setwd("E:/KIT_Forschung/141_Callie_Tibet/GAMs")
shp <- readOGR(".","NDVI_explanatory_variables_56000_masked")

# get attribute table of shapefile with response and predictor variables
data2 <- shp@data

# check NAs
table(is.na(data2))
# remove all rows containing NAs
data3 <- data2[complete.cases(data2),]

corrplot(cor(data3[,-c(1:9)]))

# set too high height variance to 50
#var <- data2[,6]
#boxplot(var)
#var[var>100] <- 100
#data2[,6] <- var

# check correlations


plot(data3$temp, data3$prec)

colnames(data3)
head(data3)

## prepare dataset (drop all variables that are neither reference or predictor variable)

data4 <- data3[,c(3:22, 24)]

# prepare dataframe containing only considered variables
names <- colnames(data4)[8:21]

##############################
####### GAMs
##############################


#################################
## calculate best subset model for dNBR with 5 variables
#################################

# automatically select variables and smoothing parameters


# get all possible combinations of five variables (available are 14 explanatory variables in total)
# permafrost had to be dropped because it consists of categoric variables (those have to be
# treated differently in the model - could be addapted but it is quite complicated)

combs <- subsets(14,5,seq(1,14,1))



# create empty lists to store results
aic_res_p30 <- list()
aic_res_p50 <- list()

set.seed(23)

nrow(combs)

# to speed up the process I ran this loop in four parts, that is, I opened Rstudio four times and
# ran in the first instance the models 1:500, in the second instance 501:1000, in the third 1001:1500
# and in the fourth 1501:2002.

i=398

for (i in 1:500){
  
  # take currently selected variable combination
  selvars <- combs[i,]
  # take corresponding variable names
  col_sel <- names[selvars]
  
  # calculate GAM model for current variable combination and using perc30 as reference
  
  # prepare formular for model call
  form_p30 <- paste0("perc30 ~ s(", col_sel[1], ") + s(", col_sel[2], ") + s(", col_sel[3], ") + s(", col_sel[4],") + s(", col_sel[5],")")
  
  # run model with dNBR as response
  b3 <- mgcv::gam(as.formula(form_p30),
                  data = data4, method="REML", select=TRUE, family=binomial)
  
  # save AIC of current model
  aic_res_p30[[i]] <- AIC(b3)

  # calculate GAM model for current variable combination and using perc50 as reference
  
  # prepare formular for model call
  form_p50 <- paste0("perc50 ~ s(", col_sel[1], ") + s(", col_sel[2], ") + s(", col_sel[3], ") + s(", col_sel[4],") + s(", col_sel[5],")")
  
  # run model with dNBR as response
  b5 <- mgcv::gam(as.formula(form_p50),
                  data = data4, method="REML", select=TRUE, family=binomial)
  
  # save AIC of current model
  aic_res_p50[[i]] <- AIC(b5)

  
  
print(i)
}


# unlist the results
aics_p30 <- do.call(rbind, aic_res_p30)
aics_p50 <- do.call(rbind, aic_res_p50)

# combine the obtained AIC values (indicating model performance) with the selected predictors
res_fin_2_preds <- data.frame(combs[1:500,], aics_p30, aics_p50)
res_fin_2_preds_p1 <- res_fin_2_preds 
save(res_fin_2_preds_p1, file = "best_subset_GAM_p_1_4.RData")


# now load the results from the other three R instances
load("best_subset_GAM_p_1_4.RData")
load("best_subset_GAM_p_2_4.RData")
load("best_subset_GAM_p_3_4.RData")
load("best_subset_GAM_p_4_4.RData")

# make sure that names are all the same
names(res_fin_2_preds_p1) <- c("v1", "v2", "v3", "v4", "v5", "aic_p30", "aic_p50")
names(res_fin_2_preds_p2) <- c("v1", "v2", "v3", "v4", "v5", "aic_p30", "aic_p50")
names(res_fin_2_preds_p3) <- c("v1", "v2", "v3", "v4", "v5", "aic_p30", "aic_p50")
names(res_fin_2_preds_p4) <- c("v1", "v2", "v3", "v4", "v5", "aic_p30", "aic_p50")

# combine results of all four R instances (this is the same as running the complete loop
# in one instances - but this would take longer)
fin_res_all <- rbind(res_fin_2_preds_p1, res_fin_2_preds_p2, res_fin_2_preds_p3, res_fin_2_preds_p4)

# save final results
save(fin_res_all, file = "GAM_best_subset_NDVI_5_vars_results.RData")

# check whether it worked
head(fin_res_all)

# now sort the models according to AIC to identify the best combinations
# models with lowest AIC indicate best performance
sorted_fin_res_all <- fin_res_all[order(fin_res_all$aic_p30),] 
head(sorted_fin_res_all)

# Variables v1-v5 indicate the variables that lead to the model with
# highest performance accuracy in this case 1, 2, 4, 5, 7

# check which predictors lead to best models
names[c(1,2,4,11,14)]

# now run the models with best results and check how they look
b4 <- mgcv::gam(perc30 ~
                  
                s(temp)+
                s(temp_10)+
                s(prec)+
                s(hf)+
                s(pika),  
                
                data = data4, method="REML", select=TRUE, family=binomial)

summary(b4)
plot(data4$perc30, b4$fitted.values, pch=16, col=rgb(0, 0, 1, 0.25), ylab="fitted RdNBR", xlab="response RdNBR", main="all predictor variables")
abline(0,1)
