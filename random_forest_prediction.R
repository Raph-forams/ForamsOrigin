# The following analysis was conducted in R version 4.0.2 using Rstudio version 1.2.1335
# Code written by Christiane Hassenrück
### Approach 1: #####
# using glm coefficients of logistic regression against distance to shore
# to classify benthic and planktonic forams

# We acknowledge the following tutorials, which provided guidance for this analysis:
# https://rpubs.com/markloessi/498787
# https://www.guru99.com/r-random-forest-tutorial.html

# set working directory
#setwd("YOUR_PATH")

# loading required R packages
require(tidyverse) # version 1.3.0
require(caret) # version 6.0-90
require(caTools) # version 1.18.0
require(randomForest) # version 4.6-14

# reading OTU occurences
# Contains occurrences of each OTUs at every station coded in presence/absence, together with the distance of each TARA station from the coast. 
# It is based on the data available in the Supplementary Material 2.  
# First Column is the distance from the coast. 
OTU <- read.table("OTUs_forams_occurences.txt", sep = " ", h = T) 

# Subset distance from coast
ENV <- OTU %>% select(Distance_from_coast) 

# Remove distance from coast from the occurece file
OTU <- OTU %>% select(-Distance_from_coast,-Station) 

# distribution of sampling sites (distance to coast)
quantile(ENV$Distance_from_coast, seq(0, 1, 0.05))

# filtering OTUs that don't occur in at least 5 sample
OTU.filt <- OTU[, apply(OTU, 2, function(x) sum(x > 0) >= 5)]
dim(OTU.filt)

# Performing the GLM
glm.coef <- map_dfr(1:ncol(OTU.filt), function(X) {
  mod <- glm(OTU.filt[, X] ~ sqrt(Distance_from_coast), data = ENV, family = binomial(link='logit'))
  mod.summary <- summary(mod)$coefficients
  mod.summary %>% 
    as_tibble() %>% 
    select(1, 4) %>% 
    flatten_dbl() %>% 
    set_names(c("intercept", "slope", "p_intercept", "p_slope"))
}) %>% 
  mutate(OTU = colnames(OTU.filt)) %>% 
  column_to_rownames(var = "OTU")

# Reading OTU metadata
# These data are included as part of the Supplementary Material 2.
# It classifies each OTU as either planktonic or benthic, or to be determined.
otu.meta <- read.table("OTUs_classification.txt", h = T, sep = "\t", stringsAsFactors = F, row.names = 1)
otu.meta <- otu.meta[rownames(glm.coef), ]
all.equal(rownames(otu.meta), rownames(glm.coef))
glm.coef <- data.frame(otu.meta, glm.coef)
glm.coef$Ecology <- as.factor(glm.coef$Ecology)
table(glm.coef$Ecology)

# data exploration
par(mfrow = c(2, 2))
for(i in 3:6) {
  boxplot(glm.coef[, i] ~ Ecology, main = colnames(glm.coef)[i], data = glm.coef)
}

# For the random forest approach, we will focus only on the model coefficients, not the p values
plot(
  glm.coef$intercept,
  glm.coef$slope,
  col = as.numeric(glm.coef$Ecology),
  pch = 16
)
legend(
  "bottomleft",
  legend = levels(glm.coef$Ecology),
  pch = 16,
  pt.cex = 2,
  col = 1:3
)

# We use the OTUs already classified as benthic or planktonic to train and validate the model

# split classified dataset into training and test data set to evaluate model accuracy
split <- sample.split(glm.coeff.class, SplitRatio = 0.7)
glm.coeff.train <- glm.coeff.class[split == TRUE, ]
glm.coeff.test <- glm.coeff.class[split == FALSE, ]
table(glm.coeff.train$Ecology)

# train random forest model
glm.coef.rf <- randomForest(
  x = glm.coeff.train[, 3:4],
  y = glm.coeff.train$Ecology,
  ntree = 1001,
  mtry = 2
)

# prediction using testing dataset
glm.coeff.test$group.predicted <- predict(glm.coef.rf, newdata = glm.coeff.test[, 3:4])

# confusion matrix for testing dataset
cm.test <- confusionMatrix(
  data = factor(glm.coeff.test$group.predicted, levels = c("Benthic", "Planktonic")),
  reference = factor(glm.coeff.test$Ecology, levels = c("Benthic", "Planktonic")),
  mode = "everything"
)
cm.test

# quick and dirty visualization
plot(
  glm.coeff.train$intercept,
  glm.coeff.train$slope,
  col = as.numeric(glm.coeff.train$Ecology),
  pch = 16
)
plot(
  glm.coeff.test$intercept,
  glm.coeff.test$slope,
  col = as.numeric(glm.coeff.test$Ecology),
  bg = as.numeric(glm.coeff.test$group.predicted),
  pch = 21,
  lwd = 2
)
# black is benthic

# apply rf model to unknown foram OTUs
glm.coeff.unknown <- droplevels(glm.coef[glm.coef$Ecology == "To be determined", ])
glm.coeff.unknown$group.predicted <- predict(glm.coef.rf, newdata = glm.coeff.unknown[, 3:4])
plot(
  glm.coeff.unknown$intercept,
  glm.coeff.unknown$slope,
  col = as.numeric(glm.coeff.unknown$group.predicted),
  pch = 16
)

# write output to file
write.table(glm.coeff.unknown, "results_glm-rf.txt")


### Approach 2: #####
# using OTU presence/absence directly to classify benthic and planktonic forams
# based on the occurrence depending on distance to shore
# using random forest implemented in mlr3

# loading required packages
require(mlr3) # version 0.8.0
require(mlr3learners) # version 0.4.2
require(mlr3tuning) # version 0.4.0
require(ranger) # version 0.12.1

# initial consistency checks and data preparation for mlr3
all.equal(rownames(otu.meta), colnames(OTU.filt))
otu.classif <- data.frame(Ecology = otu.meta$Ecology, t(OTU.filt))
otu.classif$Ecology[otu.classif$Ecology == "To be determined"] <- NA
otu.classif$Ecology <- as.factor(otu.classif$Ecology)
otu.classif.known <- otu.classif[!is.na(otu.classif$Ecology), ]

# define task
task = TaskClassif$new(
  "forams",
  backend = otu.classif.known, 
  target = "Ecology"
)

# define learner
learner <- lrn("classif.ranger", num.trees = 501)
learner$predict_type = "prob"

# inner resampling strategy
resampling <- rsmp("cv", folds = 5)

# assessment
measure = msr("classif.acc")

# hyperparameter range
param_set = paradox::ParamSet$new(
  params = list(
    paradox::ParamInt$new("mtry", lower = 2, upper = 100)
  )
)

# configure autotuner
terminator <- trm("none")
tuner = tnr("grid_search", resolution = 100)
at = AutoTuner$new(
  learner = learner, 
  resampling = resampling,
  measure = measure,
  search_space = param_set, 
  terminator = terminator, 
  tuner = tuner
)

# implement tuning
resampling_outer = rsmp("holdout", iters = 10)
rr = resample(task = task, learner = at, resampling = resampling_outer, store_models = TRUE)

# inspect
rr$aggregate(msr("classif.acc"))
rr$score(msr("classif.acc"))
rr$score(msr("classif.ce"))
rr$prediction()$confusion
rr_out_summary <- map_dfr(1:length(rr$learners), function(X) rr$learners[[X]]$tuning_result)

# first look at output of parameter tuning
plot(0, 0, xlim = c(2, 100), ylim = c(0.85, 1), type = "n")
map(1:length(rr$learners), function(X) {
  lines(
    smooth.spline(
      rr$learners[[X]]$archive$data()$mtry,
      rr$learners[[X]]$archive$data()$classif.acc,
      spar = 0.5
    ),
    col = X,
    lwd = 2
  )
})
abline(v = 12)
# result: 501 trees with mtry 12

# update learner
learner2 <- lrn("classif.ranger", num.trees = 501, mtry = 12)
learner2$predict_type = "prob"

# split train and test data
train_set = sample(task$nrow, 0.7 * task$nrow)
test_set = setdiff(seq_len(task$nrow), train_set)

# training
learner2$train(task, row_ids = train_set)

# validation
test_prediction <- learner2$predict(task, row_ids = test_set)
print(test_prediction)
test_prediction$confusion
measures <- msrs(c("classif.acc", "classif.ce", "classif.fbeta"))
test_prediction$score(measures)

# set classification task for unknown lineages
otu.classif.unknown <- otu.classif[is.na(otu.classif$Ecology), ]

# define task
task.unkown = TaskClassif$new(
  "forams_unkown",
  backend = otu.classif.unknown, 
  target = "Ecology"
)
predict.unknown <- learner2$predict(task.unkown)
predict.unknown.df <- data.frame(
  row.names = rownames(otu.classif.unknown),
  as.data.table(predict.unknown)
)

# write output
write.table(predict.unknown.df, "otu_pa_rf_prediction.txt", sep = "\t", quote = F)

