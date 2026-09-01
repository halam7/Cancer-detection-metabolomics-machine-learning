install.packages("keras")
install.packages("PRROC")
install.packages("readxl")
install.packages("glmnet")
install.packages("caret")
install.packages("rempsyc")
install.packages("corrr")
install.packages('writexl')
install.packages("reshape2")
install.packages("ggplot2")
library(reshape2)
library(ggplot2)
library(writexl)
library(gtsummary)
library(rempsyc)
library(readxl)
library(keras)
library(caret)
library(PRROC)
library(effectsize)
library(flextable)
library(broom)
library(dplyr)
library(keras)
library(caret)
library(tensorflow)
library(reticulate)
library(tidyr)
library(car)



###################### Confusion matrix function #####################

confusion_matrix <- function (true.values, predicted.values)
{
  true <- true.values
  predicted <- predicted.values
  
  TP <- sum(true == 1 & predicted == 1)
  TN <- sum(true == 0 & predicted == 0)
  FP <- sum(true == 0 & predicted == 1)
  FN <- sum(true == 1 & predicted == 0)
  
  accuracy <- (TP+TN) / (TP + TN + FP + FN)
  sensiticity <- TP / (TP + FN)
  specificity <- TN / (TN + FP)
  f1_score <- 2*TP / (2*TP + FP + FN)
  precision <- TP/ (TP + FP)
  
  roc_curve <- PRROC::roc.curve(scores.class0 = predicted[true==1], 
                                scores.class1 = predicted[true==0], curve = T)
  #plot(roc_curve)
  roc_auc <- roc_curve$auc
  
  n_metrics <- 6
  output <- matrix(ncol=n_metrics, nrow=1)
  output[,1] <- accuracy
  output[,2] <- precision 
  output[,3] <- sensiticity
  output[,4] <- f1_score
  output[,5] <- specificity
  output[,6] <- roc_auc
  
  colnames(output) <- c("Accuracy","Precision", "Sensitivity","F1_score","Specificity", "ROC_AUC")
  
  #return(output)
  print(output)
  
}

#######################################################################################################

setwd("~/Desktop/BSTT527/Homework/")
dat <- read_excel("Metab_Final.xlsx")
str(dat)
dim(dat)
table(dat$type) ## Is the dataset balanced? 313N; 208P


## Turn the label into numeric 0 for no cancer; 1 for cancer
dat$type <- as.numeric (as.factor(dat$type)) - 1
table(dat$type)
names(dat)
dat<- dat [,-1]
dim(dat)
names(dat)

## removed number at end of column name and added M for metabolite b/c starting with a number
## was causing errors

names(dat) <- sub(x = names(dat),"_.*","")
colnames(dat) <- paste0('M_', colnames(dat))

## the spaces,-,',and , in column names are causing issues esp in random forest
names(dat) <- gsub(" ","_", colnames(dat))
names(dat) <- gsub("-","_", colnames(dat))
names(dat) <- gsub(",","_", colnames(dat))
names(dat) <- gsub("'","_", colnames(dat))
names(dat) <- gsub(":","_", colnames(dat))

names(dat)
dim(dat)
colnames(dat)[which(colnames(dat) == 'M_type')] <- 'Y'
dat$Y<-as.numeric(dat$Y)
write_xlsx(dat, "Metabolomic_Names") ### names are changed start with M, Y numeric

## Any NAs in the dataset? NO

sum(is.na(dat))


############################# Exploratory data analysis  ##################################
## all are numeric
dat <- read_excel("Metabolomic_Names")
dat$Y <-as.factor (dat$Y)
table (dat$Y)
str(dat)
rho <- cor(dat[, -1])

library(corrr)
x<- dat[, -1] %>%
  correlate() %>%
  shave() %>%
  stretch() %>% 
  arrange(r) %>%
  fashion()

x$r <- as.numeric (x$r)
highx<- (subset(x, abs(r) > 0.69))
highx <- highx[order(highx$r, decreasing = TRUE),]
rownames(highx) <- 1:nrow(highx)
colnames(highx) <- c("Molecule 1", "Molecule 2", "r")
write_xlsx(highx, "correlation_high.xlsx")
highx<-as.data.frame(highx)
dim(highx)

# Create a histogram
hist(highx$r, 
     breaks = 100,
     main = "Histogram of Correlation", 
     xlab = "Values", 
     col = "blue", 
     border = "black")

install.packages("ggplot2")
library(ggplot2)

ggplot(highx$r)
# Create a histogram using ggplot2
ggplot(x, aes(x = values)) +
  geom_histogram(binwidth = 0.1, fill = "blue", color = "black") +
  labs(title = "Histogram of r", x = "Values", y = "Frequency")

############################# Exploratory data analysis t-test  ##################################
setwd("~/Desktop/BSTT527/Homework/")
dat <- read_excel("Metabolomic_Names")
dat$Y <-as.factor (dat$Y)

t.test.results <- nice_t_test(
  data = dat,
  response = names(dat)[2:148],
  group = "Y",
  warning = FALSE,
)
as.data.frame(t.test.results)
sorted <- t.test.results[order(abs(t.test.results$t), decreasing = TRUE),]
class(sorted)

bonf <- sorted %>% filter(sorted$p < (0.05/14196))

my_table <- nice_table(sorted,
                       title =  "Metabolites",
)
my_table


means <- function(data) {
  data %>%
    group_by(Y) %>%
    summarise(across(everything(), list(
      mean = ~ mean(.x)
    ))) %>%
    unnest(cols = everything())
}

meansdat <-means(dat)
meansdat <- as.data.frame(t(meansdat))
meansdat$diffmean <- as.numeric(meansdat$V2) - as.numeric (meansdat$V1)
meansdat$status = ifelse(meansdat$diffmean > 0, 'higher in cancer', 
                         ifelse(meansdat$diffmean < 0, 'lower in cancer', 'zero'))



############################ Outliers

library(car)
setwd("~/Desktop/BSTT527/Homework/")
dat <- read_excel("Metabolomic_Names")
dat$Y <-as.factor (dat$Y)

datten <- dat %>% dplyr::select (Y, M_Kynurenate, M_Neopterin, M_SAM, M_NMN , M_2_PG_3_PG,
                                 M_Isocitrate, M_Lactate, M_GlcNAc6p, M_Succinate, M_Uridine)

datten <- setNames(datten, c("C", "Kynurenate", "Neopterin", "SAM", "NMN" , "2_PG_3_PG",
                                 "Isocitrate", "Lactate", "GlcNAc6p", "Succinate", "Uridine"))
dattens <- datten %>% mutate_at(c(2,3,4,5,6,7,8,9,10,11), scale)

?Boxplot ## can check rownumer to see if outliers are coming from same samples
Boxplot(dat$M_GlcNAc6p ~ Y,  data = dat)


# Create grouped boxplots
dattens$C<- as.factor (dattens$C)
str(dattens)
dattens.m <- melt(dattens, id.var = "C")
ggplot(data = dattens.m, aes(x=variable, y=value)) + geom_boxplot(aes(fill=C))

############################# Exploratory data normality  ###############################
dat_log <- dat %>% mutate(across(-Y, log))
str(dat_log)

test_normality <- function(data) {
  data %>%
    group_by(Y) %>%
    summarise(across(everything(), list(
      shapiro_p_value = ~ shapiro.test(.x)$p.value
    ))) %>%
    unnest(cols = everything())
}

normality <-test_normality(dat)
normality <- as.data.frame(t(normality))
colnames(normality) 

normal <- normality %>% filter(V1 > 0.05) %>% filter(V2 > 0.05) 

normalitylog <-test_normality(dat_log)
normalitylog <- as.data.frame(t(normalitylog))
colnames(normalitylog) 

normallog <- normalitylog %>% filter(V1 > 0.05) %>% filter(V2 > 0.05) 

### only 3 out of 147 normally distributed; 13 after log transform

test_skewness <- function(data) {
  data %>%
    group_by(Y) %>%
    summarise(across(everything(), list(
      skewness = ~ skewness(.x)
    ))) %>%
    unnest(cols = everything())
}

skewness <-test_skewness(dat)
skewness <- as.data.frame(t(skewness))
colnames(skewness) 

unskewed <- skewness %>% filter(V1 < 0.5 & V1 >-1/2) %>% filter(V2 < 0.5 & V2 >-1/2)

skewnesslog <-test_skewness(dat_log)
skewnesslog <- as.data.frame(t(skewnesslog))
colnames(skewnesslog) 

unskewedlog <- skewnesslog %>% filter(V1 < 0.5 & V1 >-1/2) %>% filter(V2 < 0.5 & V2 >-1/2)
dim(unskewedlog)
### only 12 of 147 not skewed; 56 of 147 after log


######################### Split data: 2/3 train, 1/3 test ###############################

set.seed(111)
index <- sample(1:nrow(dat), 2/3*nrow(dat) )
dat_train <- dat [index, ]
dat_test <- dat [-index, ]

####################### Multiple logistic Regression #######################

setwd("~/Desktop/BSTT527/Homework/")
dat <- read_excel("Metabolomic_Names")
str(dat)
set.seed(111)
index <- sample(1:nrow(dat), 2/3*nrow(dat) )
dat_train <- dat [index, ]
dat_test <- dat [-index, ]
dat_train_log <- dat_train %>% mutate(across(-Y, log)) ### needed to converge
dat_test_log <- dat_test %>% mutate(across(-Y, log))

log.fit <- glm(Y~., data = dat_train_log, family = binomial)
summary (log.fit)

ld.vars <- attributes(alias(log.fit)$Complete)$dimnames[[1]]


VIF <- vif(log.fit)
VIF[order(abs(VIF), decreasing = TRUE)]

problog <- predict(log.fit, dat_test_log, type = c("response"))
predlog.y <- as.numeric (problog > 0.5)
metric.log <- confusion_matrix(dat_test$Y, predlog.y)

######## with ten only
dat10_train <- datten [index, ]
dat10_test <- datten [-index, ]
importance(rf.fit)
str(dat10_test)

log.fit <- glm(Y~., data = dat10_train, family = binomial)
summary (log.fit)

ld.vars <- attributes(alias(log.fit)$Complete)$dimnames[[1]]

vif(log.fit)

problog10 <- predict(log.fit, dat10_test, type = c("response"))
metric.log10 <- confusion_matrix(dat10_test$Y, predlog10.y)


####################### LDA #################################
setwd("~/Desktop/BSTT527/Homework/")
dat <- read_excel("Metabolomic_Names")
str(dat)
set.seed(111)
index <- sample(1:nrow(dat), 2/3*nrow(dat) )
dat_train <- dat [index, ]
dat_test <- dat [-index, ]
dat_train_log <- dat_train %>% mutate(across(-Y, log)) ### needed to converge
dat_test_log <- dat_test %>% mutate(across(-Y, log))

library(MASS)
lda.fit.log <- lda(Y~., data = dat_train_log)
names(lda.fit.log)
# Scatter plot using the first two discriminant components 
par(mfrow=c(1,1))
plot(lda.fit) 

lda.pred <- predict(lda.fit.log, dat_test_log)
lda.pred$class	# predicted outcome classes
lda.pred$post	  # posterior probabilities
metric.lda.log <- confusion_matrix(dat_test_log$Y, lda.pred$class)

####################### QDA #############################

qda.fit <- qda(Y~., data = dat_train) 
qda.pred <- predict(qda.fit, dat_test)
qda.pred$class
metric.qda <- confusion_matrix(dat_test$Y, qda.pred$class)
## Error in qda.default(x, grouping, ...) : some group is too small for 'qda'

######## ten metabolite model
qda10.fit <- qda(Y~., data = dat10_train)
qda10.pred <- predict(qda10.fit, dat10_test)
qda10.pred$class
metric.qda10 <- confusion_matrix(dat10_test$Y, qda.pred$class)

######### log
qda.fit.log <- qda(Y~., data = dat_train_log) 
qda.pred.log <- predict(qda.fit.log, dat_test_log)
qda.pred.log$class
metric.qda.log <- confusion_matrix(dat_test_log$Y, qda.pred$class)
## Error in qda.default(x, grouping, ...) : some group is too small for 'qda'


####################### NAIVE-BAYES #############################
library(e1071)

nb.fit <- naiveBayes(Y~., data = dat_train)
nb.pred <- predict(nb.fit, dat_test)
metric.nb <- confusion_matrix(dat_test$Y, nb.pred)

nb10.fit <- naiveBayes(Y~., data = dat10_train)
nb10.pred <- predict(nb10.fit, dat10_test)
metric.nb10 <- confusion_matrix(dat10_test$Y, nb10.pred)

nb.fit.log <- naiveBayes(Y~., data = dat_train_log)
nb.pred.log <- predict(nb.fit.log, dat_test_log)
metric.nb.log <- confusion_matrix(dat_test_log$Y, nb.pred.log)


####################### KNN #####################################

library(class)
library(e1071)
xtrain <- (model.matrix(Y~., dat_train)[, -1])
dim(xtrain)

ytrain <- as.factor(dat_train$Y)
length(ytrain)

# testing data 
xtest <- (model.matrix(Y~., dat_test)[, -1])

ytest <- as.factor(dat_test$Y)
length(ytest)

# train class is the response variable in the data set
train_class<-ytrain
length(train_class)
test_class<-ytest
length(test_class)

knn.cross <- tune.knn(x = xtrain, y = train_class, k = 1:30,tunecontrol=tune.control(sampling = "cross"), cross=10)
summary(knn.cross)

K<-3
knnfit<-knn(xtrain, xtest, train_class, k = K, 
            prob=TRUE) 

metric.knn <- confusion_matrix(dat_test$Y, knnfit)

####################### SVM #####################################
## linear ten metabolite model
library(e1071)
dat10_train <- datten [index, ]
dat10_test <- datten [-index, ]
str(dat10_train)
str(dat10_test)
dat10_train$Y <- as.factor(dat10_train$Y)
dat10_test$Y <- as.factor(dat10_test$Y)

tune.out <- tune(svm , Y ~ ., data = dat10_train  , kernel = "linear",
                 ranges = list(cost = c(0.01, 0.1, 1, 5)))
summary(tune.out) # not much difference choose 0.01

svm.fit <- svm(Y ~ ., data = dat10_train , kernel = "linear",
               cost = 0.1, scale = FALSE)

svm.pred10 <- predict(svm.fit , dat10_test)

metric.svmlinear <- confusion_matrix(dat10_test$Y, svm.pred10)

## polynomial
tune.out <- tune(svm , Y ~ ., data = dat10_train,
                kernel = "polynomial",
                ranges = list(cost = c(0.01, 0.1 , 1, 10, 100)))

summary(tune.out)

svm.fit.poly <- svm(Y ~ ., data = dat10_train , kernel = "polynomial",
               cost = 10, scale = FALSE)

svm.poly.pred10 <- predict(svm.fit.poly, dat10_test)

metric.svm.poly <- confusion_matrix(dat10_test$Y, svm.poly.pred10)

#####################################log transformed SVM 
setwd("~/Desktop/BSTT527/Homework/")
dat <- read_excel("Metabolomic_Names")
str(dat)
set.seed(111)
index <- sample(1:nrow(dat), 2/3*nrow(dat) )
dat_train <- dat [index, ]
dat_test <- dat [-index, ]
dat_train_log <- dat_train %>% mutate(across(-Y, log)) ### needed to converge
dat_test_log <- dat_test %>% mutate(across(-Y, log))
dat_train_log$Y <- as.factor(dat_train_log$Y)
dat_test_log$Y <- as.factor(ddat_test_log$Y)

tune.out.log <- tune(svm , Y ~ ., data = dat_train_log  , kernel = "linear",
                 ranges = list(cost = c(0.01, 0.1, 1, 5)))
summary(tune.out.log) # not much difference choose 0.01
svm.fit.log <- svm(Y ~ ., data = dat_train_log , kernel = "linear",
               cost = 0.01, scale = FALSE)

svm.pred.log <- predict(svm.fit.log , dat_test_log)
metric.svmlinear <- confusion_matrix(dat_test_log$Y, svm.pred.log)

########### poly

tune.out.logp <- tune(svm , Y ~ ., data = dat_train_log  , kernel = "polynomial",
                     ranges = list(cost = c(0.01, 0.1, 1, 5)))
summary(tune.out.logp) # not much difference choose 0.01
svm.fit.logp <- svm(Y ~ ., data = dat_train_log , kernel = "polynomial",
                   cost = 5, scale = FALSE)

svm.pred.logp <- predict(svm.fit.logp , dat_test_log)
metric.svmpoly <- confusion_matrix(dat_test_log$Y, svm.pred.logp)

############ radial

set.seed (1)
tune.out <- tune(svm , Y ~ ., data = dat_train_log,
                 kernel = "radial",
                 ranges = list(
                   cost = c(0.1 , 1, 10, 100, 1000) ,
                   gamma = c(0.5, 1, 2, 3, 4)
                 )
)
summary(tune.out) ## gamma = 0.5, cost = 0.1

svmfit.rad.log <- svm(Y ~ ., data = dat_train_log, kernel = "radial",
              gamma = 0.5, cost = 0.1)

svm.pred.rad.log <- predict(svmfit.rad.log , dat_test_log)
metric.svm.rad <- confusion_matrix(dat_test_log$Y, svm.pred.rad.log)


####################### TREES ###################################
## Random Forest
setwd("~/Desktop/BSTT527/Homework/")
dat <- read_excel("Metabolomic_Names")
str(dat)
set.seed(111)
index <- sample(1:nrow(dat), 2/3*nrow(dat) )
dat_train <- dat [index, ]
dat_test <- dat [-index, ]
library(randomForest)

set.seed(1)

rf.fit=randomForest(Y~.,
                    data=dat_train,
                    ntree=1000,     
                    mtry=12, 
                    importance=TRUE)


rf.fit
rf.pred = predict(rf.fit, dat_test, type="class")
rf.pred <- as.numeric (rf.pred > 0.5)

imp <- as.data.frame (importance(rf.fit))
imp_sorted_desc <- imp[order(-imp$IncNodePurity), ]
top_features <- imp_sorted_desc[1:20, ]

varImpPlot (rf.fit) 


metric.rf <- confusion_matrix(dat_test$Y, rf.pred)

#################################### Random Forest 10 #####################
rf.fit10=randomForest(Y~.,
                    data=dat10_train,
                    ntree=1000,     
                    mtry=3, 
                    importance=TRUE)


rf.fit10
rf.pred10 = predict(rf.fit10, dat10_test, type="class")

importance(rf.fit10) 
varImpPlot(rf.fit10)  

metric.rf <- confusion_matrix(dat10_test$Y, rf.pred10)

######################################### Bagging
bag.fit=randomForest(Y~.,
                       data=dat_train,
                       ntree=1000,      # number of trees/bootstrap models
                       mtry=147,        # for bagging mtry = total number of features = p  
                       importance=TRUE)

bag.pred = predict(bag.fit, dat_test, type="class")
bag.pred <- as.numeric (bag.pred > 0.5)
importance(bag.fit) 
varImpPlot(bag.fit)  

bagimp <- as.data.frame (importance(bag.fit))
bagimp_sorted_desc <- bagimp[order(-bagimp$IncNodePurity), ]

metric.rf <- confusion_matrix(dat_test$Y, bag.pred)

############################## Boosting
library(gbm)

# Important Notes:
# Although gbm has the ability to recognize Y variable type, 
# gbm prefers 0,1 numerical valued Y variable
str(dat_train)

boost.fit=gbm(Y~.,
                data=dat_train,
                distribution="bernoulli", # distribution of Y
                n.trees=1000,             # number of iterations
                shrinkage=0.01,           # learning rate or step-size reduction ranges (0.001,0.1), small is slower 
                cv.folds=0,               # >1 performs CV and returns CV error 
                interaction.depth=3       # Interaction or tree depth 
)

boost.fit
summary.gbm(boost.fit)


#Predictions using predict() 
boost.pred=as.numeric (predict(boost.fit, dat_test, type='response') > 0.5)
metric.boost <- confusion_matrix(dat_test$Y, boost.pred)

############################### BOOST log -transformed ##############################
setwd("~/Desktop/BSTT527/Homework/")
dat <- read_excel("Metabolomic_Names")
str(dat)
set.seed(111)
index <- sample(1:nrow(dat), 2/3*nrow(dat) )
dat_train <- dat [index, ]
dat_test <- dat [-index, ]
dat_train_log <- dat_train %>% mutate(across(-Y, log)) ### needed to converge
dat_test_log <- dat_test %>% mutate(across(-Y, log))

## Boosting
library(gbm)

# Important Notes:
# Although gbm has the ability to recognize Y variable type, 
# gbm prefers 0,1 numerical valued Y variable
str(dat_train)

boost.fitlt=gbm(Y~.,
              data=dat_train_log,
              distribution="bernoulli", # distribution of Y
              n.trees=1000,             # number of iterations
              shrinkage=0.01,           # learning rate or step-size reduction ranges (0.001,0.1), small is slower 
              cv.folds=0,               # >1 performs CV and returns CV error 
              interaction.depth=3       # Interaction or tree depth 
)

boost.fitlt
summary.gbm(boost.fitlt)


#Predictions using predict() 
boost.predlt=as.numeric (predict(boost.fitlt, dat_test_log, type='response') > 0.5)
metric.boostlt <- confusion_matrix(dat_test_log$Y, boost.predlt)

################################## Boosting 10 meta
library(gbm)
dat10_train <- datten [index, ]
dat10_test <- datten [-index, ]
str(dat10_train)
str(dat10_test)
dat10_train$Y <- as.numeric(dat10_train$Y)
dat10_test$Y <- as.numeric(dat10_test$Y)
# Important Notes:
# Although gbm has the ability to recognize Y variable type, 
# gbm prefers 0,1 numerical valued Y variable

boost.fit10=gbm(Y~.,
                data=dat10_train,
                distribution="bernoulli", # distribution of Y
                n.trees=1000,             # number of iterations
                shrinkage=0.01,           # learning rate or step-size reduction ranges (0.001,0.1), small is slower 
                cv.folds=0,               # >1 performs CV and returns CV error 
                interaction.depth=3       # Interaction or tree depth 
)

boost.fit10
summary.gbm(boost.fit10)


#Predictions using predict() 
boost.pred10=as.numeric (predict(boost.fit10, dat10_test, type='response') > 0.5)
metric.boost10 <- confusion_matrix(dat10_test$Y, boost.pred10)

######################### LASSO ###############################################################
setwd("~/Desktop/BSTT527/Homework/")
dat <- read_excel("Metabolomic_Names")
str(dat)
set.seed(111)
index <- sample(1:nrow(dat), 2/3*nrow(dat) )
dat_train <- dat [index, ]
dat_test <- dat [-index, ]

xtrain <- (model.matrix(Y~., dat_train)[, -1])
ytrain <- as.factor(dat_train$Y)

xtest <- (model.matrix(Y~., dat_test)[, -1])
ytest <- as.factor(dat_test$Y)

library(glmnet)
set.seed(1)
lassocv.fit <- cv.glmnet(xtrain, ytrain, alpha = 1, nfolds = 10, family="binomial") # By default, 10-fold CV.
plot(lassocv.fit)
print(lassocv.fit)

best.lambda <- lassocv.fit$lambda.min # lambda that results in the smallest cv error.
problas <- predict(lassocv.fit, s = "lambda.min", newx = xtest)
sorted <-as.matrix (coef(lassocv.fit, s = "lambda.min"))
sorted[order(abs(sorted[,1]), decreasing = TRUE),]
predlas.y <- as.numeric (problas > 0.5)

metric.lasso <- confusion_matrix(ytest, predlas.y)

##### find lambda so only 10 nonzero remain and make a dataset with only those predictors
sorted10 <-as.matrix (coef(lassocv.fit, s = 0.0939))
sorted10[order(abs(sorted10[,1]), decreasing = TRUE),]
lassoten <- c("M_Kynurenate", "M_Neopterin", "M_SAM", "M_NMN" , "M_2_PG_3_PG",
              "M_Isocitrate", "M_Lactate", "M_GlcNAc6p", "M_Succinate",  "M_Uridine")

datten <- dat %>% dplyr::select (Y, M_Kynurenate, M_Neopterin, M_SAM, M_NMN , M_2_PG_3_PG,
                                 M_Isocitrate, M_Lactate, M_GlcNAc6p, M_Succinate, M_Uridine)
str(datten)

lambda_seq <- lassocv.fit$lambda

# Find the lambda that gives 10 non-zero coefficients
for (lambda in lambda_seq) {
  coef_count <- sum(coef(lassocv.fit, s = lambda) != 0)
  if (coef_count == 11) { # Including the intercept
    lambda_10 <- lambda
    break
  }
}

# Print the lambda value
coefficients <- coef(lassocv.fit, s = lambda_10)
coef_df <- as.data.frame(as.matrix(coefficients))
coef_df <- coef_df[-1, , drop = FALSE] # Remove intercept

# Sort by absolute value in decreasing order
coef_df$coefficients <- as.numeric(coef_df[, 1])
coef_df$abs_value <- abs(coef_df$coefficients)
sorted_coef_df <- coef_df[order(-coef_df$abs_value), ]

# Print the top 15 predictors
top_10_predictors <- head(sorted_coef_df, 10)
print(top_10_predictors)

##ridge 
ridgecv.fit <- cv.glmnet(xtrain, ytrain, alpha = 0, nfolds = 10, family="binomial") # By default, 10-fold CV.
plot(ridgecv.fit)
best.lambda <- ridgecv.fit$lambda.min # lambda that results in the smallest cv error.
sort (abs((coef(ridgecv.fit, s = "lambda.min"))))

probrid <- predict(ridgecv.fit, s=best.lambda, newx = xtest)
predrid.y <- as.numeric (probrid > 0.5)

metric.lasso <- confusion_matrix(ytest, predrid.y)


#################################### ANN ####################################################
setwd("~/Desktop/BSTT527/Homework/")
dat <- read_excel("Metabolomic_Names")

library(nnet) # can only fit one hidden layer NN - need to specify the number of units
# also can use as a tuning parameter
# also have lambda -- similar to ridge regression -- penalty term
library(caret)

x <- dat[, -1]
y <- dat$Y
n <- nrow(dat)

set.seed (111)
set.seed(111)
index <- sample(1:nrow(dat), 2/3*nrow(dat) )
dat_train <- dat [index, ]
dat_test <- dat [-index, ]

nnetFit <- nnet(x[index ,], y[index],
                size = 15,  # with 15 hidden units
                decay = 0.01,
                linout = TRUE,
                trace = FALSE,
                preProc = c("center", "scale"),
                #Expand the number of iterations to find parameter estimates.
                maxit = 500,
                #and the number of parameters used by the model
                MaxNWts = 15 * (ncol(x) + 1) + 15 + 1)

pred.y <- predict(nnetFit, x[-index ,])
pred.y <- as.numeric (pred.y > 0.5)
metric.ann <- confusion_matrix(y[-index], pred.y)

## weight decay and the number of hidden units are tuning parameters
nnetGrid <- expand.grid(.decay = c(0, 0.01, .1),
                        .size = c(1:10), .bag = FALSE)
set.seed(1)

ctrl <- trainControl(method = "cv", number = 10)
nnetTune <- train(x[index ,], y[index],
                  method = "avNNet",
                  tuneGrid = nnetGrid,
                  trControl = ctrl,
                  preProc = c("center", "scale"),
                  linout = TRUE,
                  trace = FALSE,
                  MaxNWts = 10 * (ncol(x) + 1) + 10 + 1,
                  maxit = 500)
nnetTune
plot(nnetTune)

nnetPred <- predict(nnetTune, x[-index ,])

