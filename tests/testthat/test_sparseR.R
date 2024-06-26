context("sparseR formula specification")
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

set.seed(123)

test_that("Different vals of k and poly work, general formula", {
  expect_silent({
    obj1 <- sparseR(Sepal.Width ~ ., data = iris)
    obj2 <- sparseR(Sepal.Width ~ ., data = iris, k = 2, poly = 2)
    obj3 <- sparseR(Sepal.Width ~ ., data = iris, k = 1, poly = 1)
    obj4 <- sparseR(Sepal.Width ~ ., data = iris, k = 0, poly = 2)
  })

  expect_error(sparseR(Sepal.Width ~ ., data = iris, k = 1, poly = 0))
  expect_error(sparseR(Sepal.Width ~ ., data = iris, k = 1, poly = NULL))
  expect_error(sparseR(Sepal.Width ~ ., data = iris, k = 1, poly = 0))
  expect_error(sparseR(Sepal.Width ~ ., data = iris, k = NULL, poly = NULL))
})

test_that("Different vals of k and poly work, specific formulae", {

  formula <- Sepal.Width ~ Petal.Width + BV + Petal.Length + Group
  expect_silent({
    obj1 <- sparseR(formula, data = iris)
    obj2 <- sparseR(formula, data = iris, k = 2, poly = 2)
    obj3 <- sparseR(formula, data = iris, k = 1, poly = 1)
    obj4 <- sparseR(formula, data = iris, k = 0, poly = 2)
  })

  formula <- Sepal.Width ~ BV + UBV + Group
  expect_silent({
    obj1 <- sparseR(formula, data = iris)
    obj2 <- sparseR(formula, data = iris, k = 2, poly = 2)
    obj3 <- sparseR(formula, data = iris, k = 1, poly = 1)
    obj4 <- sparseR(formula, data = iris, k = 0, poly = 2)
    obj5 <- sparseR(formula, data = iris, k = 5, poly = 5)
  })

  formula <- Sepal.Width ~ BV + UBV
  expect_warning({
    obj1 <- sparseR(formula, data = iris)
    obj2 <- sparseR(formula, data = iris, k = 2, poly = 2)
    obj3 <- sparseR(formula, data = iris, k = 1, poly = 1)
    obj4 <- sparseR(formula, data = iris, k = 0, poly = 2)
  })

})

## Test Detrano use-case

data("Detrano")

# Quicken compute time
cleveland <- cleveland[1:50,]

cleveland$thal <- factor(cleveland$thal)
cleveland$case <- 1*(cleveland$num > 0)
cleveland$num <- NULL

# Convert variables into factor variables if necessary!
cleveland$sex <- factor(cleveland$sex)
cleveland$fbs <- factor(cleveland$fbs)
cleveland$exang <- factor(cleveland$exang)

# Simulate missing data
cleveland$thal[2] <- cleveland$thalach[1] <- NA


test_that("Different vals of k and poly work, cleveland", {
  expect_silent({
    obj1 <- sparseR(case ~ ., data = cleveland)
  })
  expect_silent({
    obj2 <- sparseR(case ~ ., data = cleveland, k = 2, poly = 2)
  })
  expect_silent({
    obj3 <- sparseR(case ~ ., data = cleveland, k = 1, poly = 1)
  })
  expect_silent({
    obj4 <- sparseR(case ~ ., data = cleveland, k = 0, poly = 2)
  })
})

test_that("Different vals of k and poly work, specific formulae, cleveland", {

  formula <- case ~ trestbps + cp + thalach + thal
  expect_silent({
    obj1 <- sparseR(formula, data = cleveland)
    obj2 <- sparseR(formula, data = cleveland, k = 2, poly = 2)
    obj3 <- sparseR(formula, data = cleveland, k = 1, poly = 1)
    obj4 <- sparseR(formula, data = cleveland, k = 0, poly = 2, family = "binomial")
    obj5 <- sparseR(formula, data = cleveland, k = 5, poly = 5)
  })

  formula <- case ~ thal + cp + chol
  expect_silent({
    obj1 <- sparseR(formula, data = cleveland)
    obj2 <- sparseR(formula, data = cleveland, k = 2, poly = 2)
    obj3 <- sparseR(formula, data = cleveland, k = 1, poly = 1)
    obj4 <- sparseR(formula, data = cleveland, k = 0, poly = 2, family = "binomial")
  })

  formula <- case ~ sex + thal
  expect_silent({
    obj1 <- sparseR(formula, data = cleveland)
    obj2 <- sparseR(formula, data = cleveland, k = 2, poly = 2)
    obj3 <- sparseR(formula, data = cleveland, k = 1, poly = 1)
    obj4 <- sparseR(formula, data = cleveland, k = 0, poly = 2, family = "binomial")
  })

})

test_that("Detrano lasso functionality", {
  expect_silent(SRL <- sparseR(formula = case ~ ., data = cleveland))
  expect_silent(APL <- sparseR(formula = case ~ ., data = cleveland, gamma = 0))
  expect_silent(MEM <- sparseR(formula = case ~ ., data = cleveland, k = 0, family = "binomial"))
  expect_silent(SRLp <- sparseR(formula = case ~ ., data = cleveland, poly = 2))

  formula <- case ~ sex + thal
  expect_silent(sparseR(formula, data = cleveland))
  expect_silent(sparseR(formula, data = cleveland, gamma = 0))
  expect_silent(sparseR(formula, data = cleveland, k = 0, family = "binomial"))
  expect_silent(sparseR(formula, data = cleveland, poly = 2))

  expect_equal(nrow(predict(SRL, type = "coef")), nrow(SRL$fit$fit$beta))
  expect_equal(nrow(predict(APL, type = "coef")), nrow(SRL$fit$fit$beta))
  expect_equal(nrow(predict(MEM, type = "coef")), nrow(MEM$fit$fit$beta))
  expect_equal(nrow(predict(SRLp, type = "coef")), nrow(SRLp$fit$fit$beta))

  expect_equal(length(coef(SRL)), nrow(SRL$fit$fit$beta))
  expect_equal(length(coef(APL)), nrow(SRL$fit$fit$beta))
  expect_equal(length(coef(MEM)), nrow(MEM$fit$fit$beta))
  expect_equal(length(coef(SRLp)), nrow(SRLp$fit$fit$beta))

  expect_equal(nrow(predict(SRL,  at = "cv1se", type = "coef")), nrow(SRL$fit$fit$beta))
  expect_equal(nrow(predict(APL,  at = "cv1se", type = "coef")), nrow(SRL$fit$fit$beta))
  expect_equal(nrow(predict(MEM,  at = "cv1se", type = "coef")), nrow(MEM$fit$fit$beta))
  expect_equal(nrow(predict(SRLp, at = "cv1se",  type = "coef")), nrow(SRLp$fit$fit$beta))

  expect_equal(length(coef(SRL, at = "cv1se")), nrow(SRL$fit$fit$beta))
  expect_equal(length(coef(APL, at = "cv1se")), nrow(SRL$fit$fit$beta))
  expect_equal(length(coef(MEM, at = "cv1se")), nrow(MEM$fit$fit$beta))
  expect_equal(length(coef(SRLp, at = "cv1se")), nrow(SRLp$fit$fit$beta))
})


test_that("Detrano MCP functionality", {
  expect_silent(SRL <- sparseR(formula = case ~ .,  penalty = "MCP", data = cleveland))
  expect_silent(APL <- sparseR(formula = case ~ .,  penalty = "MCP", data = cleveland, gamma = 0))
  expect_silent(MEM <- sparseR(formula = case ~ .,  penalty = "MCP", data = cleveland, k = 0))
  expect_silent(SRLp <- sparseR(formula = case ~ ., penalty = "MCP", data = cleveland, poly = 2))

  formula <- case ~ sex + thal
  expect_silent(sparseR(formula, penalty = "MCP", data = cleveland))
  expect_silent(sparseR(formula, penalty = "MCP", data = cleveland, gamma = 0))
  expect_silent(sparseR(formula, penalty = "MCP", data = cleveland, k = 0))
  expect_silent(sparseR(formula, penalty = "MCP", data = cleveland, poly = 2))

  expect_silent(sparseR(formula, penalty = "MCP", ncvgamma = 4, data = cleveland))
  expect_silent(sparseR(formula, penalty = "MCP", ncvgamma = 4, data = cleveland, gamma = 0))
  expect_silent(sparseR(formula, penalty = "MCP", ncvgamma = 4, data = cleveland, k = 0))
  expect_silent(sparseR(formula, penalty = "MCP", ncvgamma = 4, data = cleveland, poly = 2))

  expect_equal(nrow(predict(SRL, type = "coef")), nrow(SRL$fit$fit$beta))
  expect_equal(nrow(predict(APL, type = "coef")), nrow(SRL$fit$fit$beta))
  expect_equal(nrow(predict(MEM, type = "coef")), nrow(MEM$fit$fit$beta))
  expect_equal(nrow(predict(SRLp, type = "coef")), nrow(SRLp$fit$fit$beta))

  expect_equal(length(coef(SRL)), nrow(SRL$fit$fit$beta))
  expect_equal(length(coef(APL)), nrow(SRL$fit$fit$beta))
  expect_equal(length(coef(MEM)), nrow(MEM$fit$fit$beta))
  expect_equal(length(coef(SRLp)), nrow(SRLp$fit$fit$beta))

  expect_equal(nrow(predict(SRL,  at = "cv1se", type = "coef")), nrow(SRL$fit$fit$beta))
  expect_equal(nrow(predict(APL,  at = "cv1se", type = "coef")), nrow(SRL$fit$fit$beta))
  expect_equal(nrow(predict(MEM,  at = "cv1se", type = "coef")), nrow(MEM$fit$fit$beta))
  expect_equal(nrow(predict(SRLp, at = "cv1se",  type = "coef")), nrow(SRLp$fit$fit$beta))

  expect_equal(length(coef(SRL, at = "cv1se")), nrow(SRL$fit$fit$beta))
  expect_equal(length(coef(APL, at = "cv1se")), nrow(SRL$fit$fit$beta))
  expect_equal(length(coef(MEM, at = "cv1se")), nrow(MEM$fit$fit$beta))
  expect_equal(length(coef(SRLp, at = "cv1se")), nrow(SRLp$fit$fit$beta))
})

## Try adding ridging?
test_that("Detrano MCPnet functionality", {
  expect_silent(SRL <- sparseR(formula = case ~ .,  penalty = "MCP", alpha = .2, data = cleveland))
  expect_silent(APL <- sparseR(formula = case ~ .,  penalty = "MCP", alpha = .2, data = cleveland, gamma = 0))
  expect_silent(MEM <- sparseR(formula = case ~ .,  penalty = "MCP", alpha = .2, data = cleveland, k = 0))
  expect_silent(SRLp <- sparseR(formula = case ~ ., penalty = "MCP", alpha = .2, data = cleveland, poly = 2))

  # Print/summary
  expect_output(print(SRL))
  expect_output(print(APL))
  expect_output(print(MEM))
  expect_visible(summary(SRL))
  expect_visible(summary(SRLp))
})

# Test pooling of penalties
test_that("pooling works", {
  expect_silent(SRL <- sparseR(formula = Sepal.Width ~ ., pool = TRUE, data = iris, poly = 2))
  expect_length(val1 <- unique(SRL$results$penalty[SRL$results$Vartype == "Order 1 interaction"]), 1)
  expect_length(val2 <- unique(SRL$results$penalty[SRL$results$Vartype == "Order 2 polynomial"]), 1)
  expect_equal(val1, val2)
})
