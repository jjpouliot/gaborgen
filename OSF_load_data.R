# To Load the data in the (data frame and list forms), model summaries and loo; run the following line in R
load("data_summaries_loo.RData")

# labels do not match labels in the manuscript, for example see the comments here next to the summaries of the relevant parameters
model008_fit_summary # Model 1 in manuscript, just adaptation 
model001_fit_summary # Model 2, cue by block 
model003_fit_summary # Model 3, learning model
model002_fit_summary # Model S1, learning model 1 learning rate per participant
model004_fit_summary # Model S2, learning model beta priors on learning rates
model005_fit_summary # Model S3, arousal rating prediction
model006_fit_summary # Model S4, centered-arousal rating prediction
model007_fit_summary # Models S5, ratings arousal + expectacy + intercation (arousal, expectancy)

# The ELPD loo cross-validation information can be found in these files
model008_fit_loo
model001_fit_loo
model003_fit_loo

model002_fit_loo
model004_fit_loo
model005_fit_loo
model006_fit_loo
model007_fit_loo