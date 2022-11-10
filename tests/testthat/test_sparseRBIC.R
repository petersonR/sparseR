context("sparseRBIC formula specification")
skip_on_cran() # takes awhile, so only test on GitHub

data(iris)
iris <- iris[1:50,]

# Add another unbalanced factor
iris$Group <- factor(sample(c('A', 'B'), nrow(iris), replace = TRUE))

# Add a nzv variable
iris$NotUseful <- 2

# Add a binary variable
iris$BV <- rbinom(nrow(iris), 1, prob = .5)

# Add an unbalanced binary variable
iris$UBV <- rbinom(nrow(iris), 1, prob = .02)

# Add missing continuous data
iris$Sepal.Length[5] <- NA

# Add missing factor data
iris$Group[2] <- NA

# Add missing binary data
iris$BV[12] <- NA

test_that("Different vals of k and poly work, general formula, RBIC", {
  expect_silent({
    obj1 <- sparseRBIC_step(Sepal.Width ~ ., data = iris, message = FALSE)
  })
  expect_silent({
    obj2 <- sparseRBIC_step(Sepal.Width ~ ., data = iris, k = 2, poly = 2, message = FALSE)
  })
  expect_silent({
    obj3 <- sparseRBIC_step(Sepal.Width ~ ., data = iris, k = 1, poly = 1, message = FALSE)
  })
  expect_silent({
    obj4 <- sparseRBIC_step(Sepal.Width ~ ., data = iris, k = 0, poly = 2, message = FALSE)
  })

  expect_error(sparseRBIC_step(Sepal.Width ~ ., data = iris, k = 1, poly = 0))
  expect_error(sparseRBIC_step(Sepal.Width ~ ., data = iris, k = 1, poly = NULL))
  expect_error(sparseRBIC_step(Sepal.Width ~ ., data = iris, k = 1, poly = 0))
  expect_error(sparseRBIC_step(Sepal.Width ~ ., data = iris, k = NULL, poly = NULL))

})

test_that("Different vals of k and poly work, specific formulae, RBIC", {

  formula <- Sepal.Width ~ Petal.Width + BV + Petal.Length + Group
  expect_silent({
    obj1 <- sparseRBIC_step(formula, data = iris, message = FALSE)
    obj2 <- sparseRBIC_step(formula, data = iris, k = 2, poly = 2, message = FALSE)
    obj3 <- sparseRBIC_step(formula, data = iris, k = 1, poly = 1, message = FALSE)
    obj4 <- sparseRBIC_step(formula, data = iris, k = 0, poly = 2, message = FALSE)
  })

  formula <- Sepal.Width ~ BV + UBV + Group
  expect_silent({
    obj1 <- sparseRBIC_step(formula, data = iris, message = FALSE)
    obj2 <- sparseRBIC_step(formula, data = iris, k = 2, poly = 2, message = FALSE)
    obj3 <- sparseRBIC_step(formula, data = iris, k = 1, poly = 1, message = FALSE)
    obj4 <- sparseRBIC_step(formula, data = iris, k = 0, poly = 2, message = FALSE)
    obj5 <- sparseRBIC_step(formula, data = iris, k = 5, poly = 5, message = FALSE)
  })

})

test_that("Different vals of k and poly work, specific formulae, RAIC", {

  formula <- Sepal.Width ~ Petal.Width + BV + Petal.Length + Group
  expect_silent({
    obj1 <- sparseRBIC_step(formula, data = iris, ic = "RAIC", message = FALSE)
    obj2 <- sparseRBIC_step(formula, data = iris, k = 2, poly = 2, ic = "RAIC", message = FALSE)
    obj3 <- sparseRBIC_step(formula, data = iris, k = 1, poly = 1, ic = "RAIC", message = FALSE)
    obj4 <- sparseRBIC_step(formula, data = iris, k = 0, poly = 2, ic = "RAIC", message = FALSE)
  })

  formula <- Sepal.Width ~ BV + UBV + Group
  expect_silent({
    obj1 <- sparseRBIC_step(formula, data = iris, ic = "RAIC", message = FALSE)
    obj2 <- sparseRBIC_step(formula, data = iris, k = 2, poly = 2, ic = "RAIC", message = FALSE)
    obj3 <- sparseRBIC_step(formula, data = iris, k = 1, poly = 1, ic = "RAIC", message = FALSE)
    obj4 <- sparseRBIC_step(formula, data = iris, k = 0, poly = 2, ic = "RAIC", message = FALSE)
    obj5 <- sparseRBIC_step(formula, data = iris, k = 5, poly = 5, ic = "RAIC", message = FALSE)
  })

})

test_that("Different vals of k and poly work, specific formulae, EBIC", {

  formula <- Sepal.Width ~ Petal.Width + BV + Petal.Length + Group
  expect_silent({
    obj1 <- sparseRBIC_step(formula, data = iris, ic = "EBIC", message = FALSE)
    obj2 <- sparseRBIC_step(formula, data = iris, k = 2, poly = 2, ic = "EBIC", message = FALSE)
    obj3 <- sparseRBIC_step(formula, data = iris, k = 1, poly = 1, ic = "EBIC", message = FALSE)
    obj4 <- sparseRBIC_step(formula, data = iris, k = 0, poly = 2, ic = "EBIC", message = FALSE)
  })

  formula <- Sepal.Width ~ BV + UBV + Group
  expect_silent({
    obj1 <- sparseRBIC_step(formula, data = iris, ic = "EBIC", message = FALSE)
    obj2 <- sparseRBIC_step(formula, data = iris, k = 2, poly = 2, ic = "EBIC", message = FALSE)
    obj3 <- sparseRBIC_step(formula, data = iris, k = 1, poly = 1, ic = "EBIC", message = FALSE)
    obj4 <- sparseRBIC_step(formula, data = iris, k = 0, poly = 2, ic = "EBIC", message = FALSE)
    obj5 <- sparseRBIC_step(formula, data = iris, k = 5, poly = 5, ic = "EBIC", message = FALSE)
  })

})

test_that("Different vals of k and poly work, specific formulae, AIC/BIC", {

  formula <- Sepal.Width ~ Petal.Width + BV + Petal.Length + Group
  expect_silent({
    obj1 <- sparseRBIC_step(formula, data = iris, ic = "BIC", message = FALSE)
    obj2 <- sparseRBIC_step(formula, data = iris, k = 2, poly = 2, ic = "BIC", message = FALSE)
    obj3 <- sparseRBIC_step(formula, data = iris, k = 1, poly = 1, ic = "BIC", message = FALSE)
    obj4 <- sparseRBIC_step(formula, data = iris, k = 0, poly = 2, ic = "BIC", message = FALSE)
  })

  formula <- Sepal.Width ~ BV + UBV + Group
  expect_silent({
    obj1 <- sparseRBIC_step(formula, data = iris, ic = "AIC", message = FALSE)
    obj2 <- sparseRBIC_step(formula, data = iris, k = 2, poly = 2, ic = "AIC", message = FALSE)
    obj3 <- sparseRBIC_step(formula, data = iris, k = 1, poly = 1, ic = "AIC", message = FALSE)
    obj4 <- sparseRBIC_step(formula, data = iris, k = 0, poly = 2, ic = "AIC", message = FALSE)
  })



})

## Test Detrano use-case

data("Detrano")

# Quicken compute time
cleveland <- cleveland[1:100,]

cleveland$thal <- factor(cleveland$thal)
cleveland$case <- 1*(cleveland$num > 0)
cleveland$num <- NULL

# Convert variables into factor variables if necessary!
cleveland$sex <- factor(cleveland$sex)
cleveland$fbs <- factor(cleveland$fbs)
cleveland$exang <- factor(cleveland$exang)

# Simulate missing data
cleveland$thal[2] <- cleveland$thalach[1] <- NA


test_that("Different vals of k and poly work, cleveland, RBIC", {
  expect_silent({
    obj1 <- sparseRBIC_step(formula = case ~ ., data = cleveland, family = "binomial", message = FALSE)
  })
  expect_silent({
    obj2 <- sparseRBIC_step(case ~ ., data = cleveland, k = 2, poly = 2, family = "binomial", message = FALSE)
  })
  expect_silent({
    obj3 <- sparseRBIC_step(case ~ ., data = cleveland, k = 1, poly = 1, family = "binomial", message = FALSE)
  })
  expect_silent({
    obj4 <- sparseRBIC_step(case ~ ., data = cleveland, k = 0, poly = 2, family = "binomial", message = FALSE)
  })
})

test_that("Different vals of k and poly work, specific formulae, cleveland, RBIC", {

  formula <- case ~ trestbps + cp + thalach + thal
  expect_silent({
    obj1 <- sparseRBIC_step(formula, data = cleveland, message = FALSE)
    obj2 <- sparseRBIC_step(formula, data = cleveland, k = 2, poly = 2, message = FALSE)
    obj3 <- sparseRBIC_step(formula, data = cleveland, k = 1, poly = 1, message = FALSE)
    obj4 <- sparseRBIC_step(formula, data = cleveland, k = 0, poly = 2, family = "binomial", message = FALSE)
    obj5 <- sparseRBIC_step(formula, data = cleveland, k = 5, poly = 5, message = FALSE)
  })

  formula <- case ~ thal + cp + chol
  expect_silent({
    obj1 <- sparseRBIC_step(formula, data = cleveland, message = FALSE)
    obj2 <- sparseRBIC_step(formula, data = cleveland, k = 2, poly = 2, message = FALSE)
    obj3 <- sparseRBIC_step(formula, data = cleveland, k = 1, poly = 1, message = FALSE)
    obj4 <- sparseRBIC_step(formula, data = cleveland, k = 0, poly = 2, family = "binomial", message = FALSE)
  })

  formula <- case ~ sex + thal
  expect_silent({
    obj1 <- sparseRBIC_step(formula, data = cleveland, message = FALSE)
    obj2 <- sparseRBIC_step(formula, data = cleveland, k = 2, poly = 2, message = FALSE)
    obj3 <- sparseRBIC_step(formula, data = cleveland, k = 1, poly = 1, message = FALSE)
    obj4 <- sparseRBIC_step(formula, data = cleveland, k = 0, poly = 2, family = "binomial", message = FALSE)
  })

})

test_that("Detrano RBIC functionality", {
  expect_silent(SRL <- sparseRBIC_step(formula = case ~ .,
                                        data = cleveland, message = FALSE))
  expect_silent(MEM <- sparseRBIC_step(formula = case ~ .,
                                        data = cleveland, k = 0,
                                        family = "binomial", message = FALSE))

  expect_output(n <- sparseRBIC_step(formula = case ~ .,
                                        data = cleveland, k = 0, message = FALSE,
                                        family = "binomial", trace = TRUE))

  formula <- case ~ sex + thal
  expect_silent(SRL_b <- sparseRBIC_step(formula, data = cleveland, message = FALSE))
  expect_silent(MEM_b <- sparseRBIC_step(formula, data = cleveland, k = 0, family = "binomial", message = FALSE))

  expect_equal(length(coef(SRL)), length(SRL$fit$coef))
  expect_equal(length(coef(SRL)), length(SRL$fit$coef))

  # Bootstrapping
  expect_silent(f <- sparseRBIC_bootstrap(SRL, B = 3, quiet = TRUE))
  expect_silent(f <- sparseRBIC_bootstrap(MEM, B = 3, quiet = TRUE))
  expect_silent(f <- sparseRBIC_bootstrap(SRL_b, B = 3, quiet = TRUE))
  expect_silent(f <- sparseRBIC_bootstrap(MEM_b, B = 3, quiet = TRUE))

  # Sample splitting
  expect_silent(f <- sparseRBIC_sampsplit(SRL, S = 3, quiet = TRUE))
  expect_silent(f <- sparseRBIC_sampsplit(MEM, S = 3, quiet = TRUE))
  expect_silent(f <- sparseRBIC_sampsplit(SRL_b, S = 3, quiet = TRUE))
  expect_silent(f <- sparseRBIC_sampsplit(MEM_b, S = 3, quiet = TRUE))

})

test_that("Custom ICs work as intended", {
  expect_silent({
    obj1 <- sparseRBIC_step(formula = case ~ ., data = cleveland, family = "binomial", message = FALSE)
  })
  expect_equal(EBIC(obj1), 102.0671, tolerance = .01)
  expect_equal(RBIC(obj1), 90.61904, tolerance = .01)
  expect_equal(RAIC(obj1)[1], 78.50217, tolerance = .01)
})
