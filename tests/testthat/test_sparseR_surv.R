context("sparseR survival specification")

data(iris)

# Add another unbalanced factor
iris$Group <- factor(sample(c('A', 'B'), nrow(iris), replace = TRUE))

# Add a nzv variable
iris$NotUseful <- 2

# Add a binary variable
iris$BV <- rbinom(nrow(iris), 1, prob = .5)

# Add an unbalanced binary variable
iris$UBV <- rbinom(nrow(iris), 1, prob = .02)

# Create time-to-event data
y<-iris$time <- survival::Surv(rnorm(nrow(iris), 5)*iris$Sepal.Length, rbinom(nrow(iris), 1, prob = .5))

test_that("Different vals of k and poly work (matrix)", {

  X2 <- bake(sparseR_prep(Sepal.Width ~ ., data = iris), iris)
  X3 <- bake(sparseR_prep(Sepal.Width ~ ., data = iris, poly = 2), iris)

  expect_silent({
    obj1 <- sparseR(model_matrix = X2, y = y, pre_process = FALSE, max.iter=1e7, family = "coxph")
  })
  expect_silent({
    obj3 <- sparseR(model_matrix = X3, y = y, pre_process = FALSE, max.iter = 1e7, family = "coxph")
  })
})

## Test Sheddon (small) data
data("Sheddon_small")
set.seed(123)

test_that("Different vals of k and poly work, sheddon", {

  srp <- sparseR_prep(formula = ~., data = Z)
  X <- bake(srp, new_data = Z)

  srp2 <- sparseR_prep(survtime~., data = data.frame(survtime = S, Z), family = "coxph")

  sr_cox <- sparseR(survtime~., data = data.frame(survtime = S, Z),
                    family = "coxph", cumulative_k = TRUE)

  plot(sr_cox)
  expect_warning(effect_plot(sr_cox, coef_name = "Age", by = "Sex"))
  apl_cox <- sparseR(model_matrix = X, y = S, gamma = 0,
                     pre_process = FALSE, family = "coxph")
})

