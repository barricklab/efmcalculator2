library(mgcv)
hi <- read.csv('test.csv')
gam_model <- gam(Mutation_Rate ~ s(RBP_Length, bs = "cr", k = 3) + s(TBD_length, bs = "cr", k = 3), data = hi, method = "REML")

new_data <- data.frame(RBP_Length = 10, TBD_length = 151)

predictions <- predict(gam_model, newdata = new_data, type = "response")

summary(gam_model)

print(predictions)