
library(glmnet)
library(randomForest)
library(gbm)
library(ROCR)
library(ggfortify)

#After iterations of prin component analysis
base.pca <- prcomp(final_base_table[,c( 15:23, 26:35)], center = TRUE,scale. = TRUE)
summary(base.pca)

autoplot(base.pca, data = final_base_table[,c(  15:23, 26:35)],
         loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 3)

#Check for collinearity
data.frame(colnames(final_base_table))
cor_cols = c(15:35)

res <- cor(final_base_table[,cor_cols])

final_base_table_mod_prep1 =
  final_base_table %>% 
  ungroup() %>%
  dplyr::select(-patient_id,
                -stroke_dt,
                -index_dt,
                -had_death,
                -had_stroke,
                -death,
                -min_nyha_dt
                )

final_base_table_mod_prep1$ino_cnt = as.numeric(final_base_table_mod_prep1$ino_cnt)
final_base_table_mod_prep1$bmi = as.numeric(final_base_table_mod_prep1$bmi)

final_base_table_mod_prep1$region = as.character(final_base_table_mod_prep1$region)
final_base_table_mod_prep1$region[which(final_base_table_mod_prep1$region == "Other/Unknown")] = "Unknown"
final_base_table_mod_prep1$region = as.factor(final_base_table_mod_prep1$region)


final_base_table_mod_prep1 %>% dplyr::select(race) %>% distinct()

data.frame(colnames(final_base_table_mod_prep1))

cols <- c(3,seq(6,28))
final_base_table_mod_prep1[,cols] <- lapply(final_base_table_mod_prep1[,cols], factor)

sapply(final_base_table_mod_prep1, class)
#################################################
#Begin Double Cross Validation
#
#Testing Random Forests and Elastic Net regression
#
#Model selection to choose best 1. # of predictor vars
#     and 2. alpha and lamda values 
#
#################################################

lambdalist = 0:200/2000
alphalist = c(0, 0.05, 0.1, 0.15,.2,.4)
alpha_count = length(alphalist)
#Default mtry for classification is sqrt(predictors)
mtrylist = c(2,3,4,5,6,7,8,9)
m_count = length(mtrylist)

x.matrix = model.matrix(had_strk_or_dth~., data=final_base_table_mod_prep1)[,-1]
y = as.numeric(as.matrix(final_base_table_mod_prep1[,18]))
n = dim(x.matrix)[1]

# define the cross-validation splits 
set.seed(15)
ncv = 10
groups = c(rep(1:ncv,floor(n/ncv)), 1:(n%%ncv))
cvgroups = sample(groups,n)

allpredictedCV = matrix(rep(0,n*2), ncol=2)

for(j in 1:ncv){

 
  print(j)

  #Choose validation set
  group_assesm = (cvgroups == j)
 
  #Assign training data
  trainx = x.matrix[!group_assesm,]
  trainy = as.factor(y[!group_assesm])
  
  #Assign testing data
  testx = x.matrix[group_assesm,]
  testy = as.factor(y[group_assesm])
  
  #Assign groups within training data
  trainx.n = dim(trainx)[1] 
  
  if ((trainx.n%%ncv) == 0) {
    groups.select = rep(1:ncv,floor(trainx.n/ncv))
  } else {
    #account for different-sized input matrices
    groups.select = c(rep(1:ncv,floor(trainx.n / ncv)),(1:(trainx.n%%ncv)))
  }
  
  cvgroups.selection = sample(groups.select, trainx.n)
  
  #END- Assign groups within training data
  
  alllambdabest = rep(NA,alpha_count)
  allcvbest.net = rep(NA,alpha_count)
  rfor_predictions = matrix(rep(0,trainx.n*m_count), ncol=m_count)
  rfor.cv = rep(0,m_count)
 
  #################################################
  #Begin model selection
  #################################################
  
 for (a in 1:alpha_count) {
    
    #First, elastic net choosing best alpha and lamda
    cvfit.net = cv.glmnet(trainx, trainy, lambda=lambdalist, alpha = alphalist[a], 
                          nfolds=ncv, foldid=cvgroups.selection, family = "binomial")
    
    #Cool best lambda plot
    plot(cvfit.net$lambda, cvfit.net$cvm) 
    abline(v=cvfit.net$lambda[order(cvfit.net$cvm)[1]], col = "red")
    
    allcvbest.net[a] = cvfit.net$cvm[order(cvfit.net$cvm)[1]]
    alllambdabest[a] = cvfit.net$lambda[order(cvfit.net$cvm)[1]]
 }
  
  #Second, random forests, choosing best number of predictor vars
  for(m in 1:m_count){
    
    for(h in 1:ncv){
      
      rfgroup = (cvgroups.selection == h)
    
      rforest = randomForest(trainy[!rfgroup] ~ ., data = trainx[!rfgroup,],
                             mtry = mtrylist[m], importance = T, ntree = 500)
      
      rfor_predictions[rfgroup,m] = predict(rforest, newdata = trainx[rfgroup,],type="prob")[,2]
    }
  }

  #For Enet model assessment 
  
  #Choose best alpha and best lambda for elastic net
  whichmodel = order(allcvbest.net)[1]
  bestalpha = alphalist[whichmodel]
  bestlambda = alllambdabest[whichmodel]
  
  #Run the best model on all training data
  bestmodel.enet = glmnet(trainx, trainy, alpha = bestalpha, lambda=bestlambda, family = "binomial")
  
  #For Random Forests model assessment
  
  #Get CV for each mtry value in random forests
  for (rf in 1:m_count) {
    
    y_i <- as.numeric(as.character(trainy))
    u_i <- rfor_predictions[,rf]
    
    deviance.contribs <- (y_i * log(u_i)) + ((1-y_i) * log(1 - u_i))
    deviance.contribs <- (y_i * log(u_i)) + ((1-y_i) * log(1 - u_i))
    deviance <- -2 * sum(deviance.contribs)
    rfor.cv[rf] = deviance/length(rfor_predictions[,rf])
  }
  
  #Find mtry with lowest cv
  whichmodel2 = order(rfor.cv)[1]
  bestmtry = mtrylist[whichmodel2]
  
  #Run the best model on all training data
  bestmodel.rf = randomForest(trainy ~ ., data = trainx,
                              mtry = bestmtry, importance = T,  ntree = 500) 
  
  #################################################
  #End model selection
  #################################################
  
  #Store predictions on test set from best enent model
  allpredictedCV[group_assesm,1] = predict(bestmodel.enet, newx = testx, s = bestlambda, type= "response")

  #Store predictions on test set from random forest model
  allpredictedCV[group_assesm,2] = predict(bestmodel.rf, newdata = testx, type= "prob")[,2]
}

rfor.cv
allcvbest.net
#Assess elastic net model
y_i <- y
u_i <- allpredictedCV[,1]

deviance.contribs <- (y_i * log(u_i)) + ((1-y_i) * log(1 - u_i))
deviance.contribs <- (y_i * log(u_i)) + ((1-y_i) * log(1 - u_i))
deviance.contribs[which(is.na(deviance.contribs))] = 0
deviance <- -2 * sum(deviance.contribs)
CV.enet = deviance/length(y_i)

CV.enet

#Assess random forest model
y_i <- as.numeric(y)
u_i <- allpredictedCV[,2]

deviance.contribs <- (y_i * log(u_i)) + ((1-y_i) * log(1 - u_i))
deviance.contribs <- (y_i * log(u_i)) + ((1-y_i) * log(1 - u_i))
deviance.contribs[which(is.na(deviance.contribs))] = 0
deviance <- -2 * sum(deviance.contribs)
CV.rfor = deviance/length(y_i)

CV.rfor

#################################################
#Elastic net is winner based on deviance, but it is close enough
#that we will proceed with both models
#
#We will use ten-fold cross validation to assess both models
#
#################################################


#################################################
#Testing Elastic Net
#################################################

#Default mtry for classification is sqrt(predictors)
lambdalist = 0:200/2000

# define the cross-validation splits 
set.seed(10)
ncv = 10
groups = c(rep(1:ncv,floor(n/ncv)), 1:(n%%ncv))
cvgroups = sample(groups,n)

factor_y = as.factor(y)

predicted_enet = matrix(rep(0,n), ncol=1)

#Run CV for all values of mtry
for(j in 1:ncv){
  
  print(j)
  
  #Choose validation set
  group_assesm = (cvgroups == j)
  
  #Assign training data
  trainx = x.matrix[!group_assesm,]
  trainy = as.factor(y[!group_assesm])
  
  #Assign testing data
  testx = x.matrix[group_assesm,]
  testy = as.factor(y[group_assesm])
  
  #Assign groups within training data
  trainx.n = dim(trainx)[1] 
  
  #First, elastic net choosing best alpha and lamda
  bestmodel.enet = glmnet(trainx, trainy, alpha = bestalpha, lambda=bestlambda, family = "binomial")
  
  predicted_enet[group_assesm] = predict(bestmodel.enet, newx = testx, s = bestlambda, type= "response")
}

pred <- prediction(predicted_enet ,factor_y)

roc.perf = performance(pred, measure = "tpr", x.measure = "fpr")
auc.perf = performance(pred, measure = "auc")
auc.perf@y.values
plot(roc.perf)



#################################################
#Testing Random Forest
#################################################

mtrylist = c(5)
m_count = length(mtrylist)

# define the cross-validation splits 
set.seed(10)
ncv = 10

groups = c(rep(1:ncv,floor(n/ncv)), 1:(n%%ncv))
cvgroups = sample(groups,n)

factor_y = as.factor(y)

rfor_predictions = matrix(rep(0,n*m_count), ncol=m_count)
rfor.cv = rep(0,1)

#Run CV for a single mtry value of 5
for(b in 1:m_count){
  
  for(c in 1:ncv){
    
    rfgroup = (cvgroups == c)
    
    finalforest = randomForest(factor_y[!rfgroup] ~ ., data = x.matrix[!rfgroup,],
                               mtry = mtrylist[b], importance = T, ntree = 1000)
    
    rfor_predictions[rfgroup,b] = predict(finalforest, newdata = x.matrix[rfgroup,], type= "prob")[,2]
  }
}

pred <- prediction(rfor_predictions[,1],factor_y)

roc.perf = performance(pred, measure = "tpr", x.measure = "fpr")
auc.perf = performance(pred, measure = "auc")
auc.perf@y.values
plot(roc.perf)

pred_df= as.data.frame(rfor_predictions[,1])
names(pred_df)[1] = 'pred_val'
conf_pred = pred_df %>% mutate(class_prod = if_else(pred_val>= .27,1,0)) %>% dplyr::select(class_prod)

table(predicted = conf_pred$class_prod, actual = y)

#Calculate optimal cutoff value
opt.cut = function(perf, pred){
  cut.ind = mapply(FUN=function(x, y, p){
    d = (x - 0)^2 + (y-1)^2
    
    ind = which(d == min(d))
    
    c(sensitivity = y[[ind]], specificity = 1-x[[ind]], 
      cutoff = p[[ind]])},
    
    perf@x.values, perf@y.values, pred@cutoffs)
}
print(opt.cut(roc.perf, pred))


#Last, but not least, fit the final model
bestmodel_final = randomForest(factor_y ~ ., data = x.matrix,
                               mtry = 5, importance = T, ntree = 5000, cutoff = c(.5,.5)) 

#Plot variable importance
varImpPlot(bestmodel_final, main = "Random Forests- Variable Importance")
