context("sparseRBIC formula specification")

data(iris)

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
    obj1 <- sparseRBIC_step(Sepal.Width ~ ., data = iris)
  })
  expect_silent({
    obj2 <- sparseRBIC_step(Sepal.Width ~ ., data = iris, k = 2, poly = 2)
  })
  expect_silent({
    obj3 <- sparseRBIC_step(Sepal.Width ~ ., data = iris, k = 1, poly = 1)
  })
  expect_silent({
    obj4 <- sparseRBIC_step(Sepal.Width ~ ., data = iris, k = 0, poly = 2)
  })

  expect_error(sparseRBIC_step(Sepal.Width ~ ., data = iris, k = 1, poly = 0))
  expect_error(sparseRBIC_step(Sepal.Width ~ ., data = iris, k = 1, poly = NULL))
  expect_error(sparseRBIC_step(Sepal.Width ~ ., data = iris, k = 1, poly = 0))
  expect_error(sparseRBIC_step(Sepal.Width ~ ., data = iris, k = NULL, poly = NULL))

})

test_that("Different vals of k and poly work, specific formulae, RBIC", {

  formula <- Sepal.Width ~ Petal.Width + BV + Petal.Length + Group
  expect_silent({
    obj1 <- sparseRBIC_step(formula, data = iris)
    obj2 <- sparseRBIC_step(formula, data = iris, k = 2, poly = 2)
    obj3 <- sparseRBIC_step(formula, data = iris, k = 1, poly = 1)
    obj4 <- sparseRBIC_step(formula, data = iris, k = 0, poly = 2)
  })

  formula <- Sepal.Width ~ BV + UBV + Group
  expect_silent({
    obj1 <- sparseRBIC_step(formula, data = iris)
    obj2 <- sparseRBIC_step(formula, data = iris, k = 2, poly = 2)
    obj3 <- sparseRBIC_step(formula, data = iris, k = 1, poly = 1)
    obj4 <- sparseRBIC_step(formula, data = iris, k = 0, poly = 2)
    obj5 <- sparseRBIC_step(formula, data = iris, k = 5, poly = 5)
  })

})

test_that("Different vals of k and poly work, specific formulae, RAIC", {

  formula <- Sepal.Width ~ Petal.Width + BV + Petal.Length + Group
  expect_silent({
    obj1 <- sparseRBIC_step(formula, data = iris, ic = "RAIC")
    obj2 <- sparseRBIC_step(formula, data = iris, k = 2, poly = 2, ic = "RAIC")
    obj3 <- sparseRBIC_step(formula, data = iris, k = 1, poly = 1, ic = "RAIC")
    obj4 <- sparseRBIC_step(formula, data = iris, k = 0, poly = 2, ic = "RAIC")
  })

  formula <- Sepal.Width ~ BV + UBV + Group
  expect_silent({
    obj1 <- sparseRBIC_step(formula, data = iris, ic = "RAIC")
    obj2 <- sparseRBIC_step(formula, data = iris, k = 2, poly = 2, ic = "RAIC")
    obj3 <- sparseRBIC_step(formula, data = iris, k = 1, poly = 1, ic = "RAIC")
    obj4 <- sparseRBIC_step(formula, data = iris, k = 0, poly = 2, ic = "RAIC")
    obj5 <- sparseRBIC_step(formula, data = iris, k = 5, poly = 5, ic = "RAIC")
  })

})

## Test Detrano use-case

data("Detrano")
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
    obj1 <- sparseRBIC_step(formula = case ~ ., data = cleveland, family = "binomial")
  })
  expect_silent({
    obj2 <- sparseRBIC_step(case ~ ., data = cleveland, k = 2, poly = 2, family = "binomial")
  })
  expect_silent({
    obj3 <- sparseRBIC_step(case ~ ., data = cleveland, k = 1, poly = 1, family = "binomial")
  })
  expect_silent({
    obj4 <- sparseRBIC_step(case ~ ., data = cleveland, k = 0, poly = 2, family = "binomial")
  })
})

test_that("Different vals of k and poly work, specific formulae, cleveland, RBIC", {

  formula <- case ~ trestbps + cp + thalach + thal
  expect_silent({
    obj1 <- sparseRBIC_step(formula, data = cleveland)
    obj2 <- sparseRBIC_step(formula, data = cleveland, k = 2, poly = 2)
    obj3 <- sparseRBIC_step(formula, data = cleveland, k = 1, poly = 1)
    obj4 <- sparseRBIC_step(formula, data = cleveland, k = 0, poly = 2, family = "binomial")
    obj5 <- sparseRBIC_step(formula, data = cleveland, k = 5, poly = 5)
  })

  formula <- case ~ thal + cp + chol
  expect_silent({
    obj1 <- sparseRBIC_step(formula, data = cleveland)
    obj2 <- sparseRBIC_step(formula, data = cleveland, k = 2, poly = 2)
    obj3 <- sparseRBIC_step(formula, data = cleveland, k = 1, poly = 1)
    obj4 <- sparseRBIC_step(formula, data = cleveland, k = 0, poly = 2, family = "binomial")
  })

  formula <- case ~ sex + thal
  expect_silent({
    obj1 <- sparseRBIC_step(formula, data = cleveland)
    obj2 <- sparseRBIC_step(formula, data = cleveland, k = 2, poly = 2)
    obj3 <- sparseRBIC_step(formula, data = cleveland, k = 1, poly = 1)
    obj4 <- sparseRBIC_step(formula, data = cleveland, k = 0, poly = 2, family = "binomial")
  })

})

test_that("Detrano RBIC functionality", {
  expect_silent(SRL <<- sparseRBIC_step(formula = case ~ ., data = cleveland))
  expect_silent(MEM <<- sparseRBIC_step(formula = case ~ ., data = cleveland, k = 0, family = "binomial"))

  formula <- case ~ sex + thal
  expect_silent(sparseRBIC_step(formula, data = cleveland))
  expect_silent(sparseRBIC_step(formula, data = cleveland, k = 0, family = "binomial"))
})

test_that("Detrano RBIC stepwise predict (coef) functionality", {
  expect_equal(length(coef(SRL)), length(SRL$fit$coef))
  expect_equal(length(coef(SRL)), length(SRL$fit$coef))
})


test_that("Detrano RBIC bootstrap functionality", {
  expect_silent(f <- sparseRBIC_bootstrap(SRL, B = 25, quiet = T))
  expect_silent(f <- sparseRBIC_bootstrap(MEM, B = 25, quiet = T))

  formula <- case ~ sex + thal
  expect_silent(sparseRBIC_step(formula, data = cleveland))
  expect_silent(sparseRBIC_step(formula, data = cleveland, k = 0, family = "binomial"))
})
