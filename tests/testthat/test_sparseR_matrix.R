context("sparseR matrix specification")

data(iris)


X1 <- model.matrix(Sepal.Width ~ . - 1, data = iris)
X2 <- bake(sparseR_prep(Sepal.Width ~ ., data = iris), iris)
X3 <- bake(sparseR_prep(Sepal.Width ~ ., data = iris, poly = 2), iris)

test_that("Different vals of k and poly work (matrix)", {
  expect_silent({
    obj1 <- sparseR(model_matrix = X1, y = iris$Sepal.Width, pre_process = FALSE)
  })
  expect_silent({
    obj2 <- sparseR(model_matrix = X2, y = iris$Sepal.Width, pre_process = FALSE)
  })
  expect_silent({
    obj3 <- sparseR(model_matrix = X3, y = iris$Sepal.Width, pre_process = FALSE, max.iter = 1e7)
  })
})

## Test Detrano use-case

data("Detrano")

# Quicken compute time
cleveland <- cleveland[1:100,]

cleveland$thal <- factor(cleveland$thal)
cleveland$case <- 1*(cleveland$num > 0)

# Convert variables into factor variables if necessary!
cleveland$sex <- factor(cleveland$sex)
cleveland$fbs <- factor(cleveland$fbs)
cleveland$exang <- factor(cleveland$exang)
cleveland$num <- NULL

X1 <- model.matrix(case ~ . - 1, data = cleveland)
X2 <- bake(sparseR_prep(case ~ ., data = cleveland), cleveland)
X3 <- bake(sparseR_prep(case ~ ., data = cleveland, poly = 2), cleveland)
y <- cleveland$case[complete.cases(cleveland)]

test_that("Different vals of k and poly work, cleveland matrix", {
  expect_silent({
    obj1 <- sparseR(model_matrix = X1, y = y, pre_process = FALSE, family = "binomial")
  })
  expect_silent({
    obj2 <- sparseR(model_matrix = X2, y = cleveland$case, pre_process = FALSE)
  })
  expect_silent({
    obj3 <- sparseR(model_matrix = X3, y = cleveland$case, pre_process = FALSE)
  })
})


test_that("Detrano lasso predict (coef) functionality, matrix", {
  expect_silent(SRL <- sparseR(model_matrix = X2, y = cleveland$case, pre_process = FALSE))
  expect_silent(APL <- sparseR(model_matrix = X2, y = cleveland$case, pre_process = FALSE, gamma = 0))
  expect_silent(MEM <- sparseR(model_matrix = X1, y = y, pre_process = FALSE, family = "binomial"))
  expect_silent(SRLp <- sparseR(model_matrix = X3, y = cleveland$case, pre_process = FALSE))

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


test_that("Detrano MCP functionality (matrix)", {
  expect_silent(SRL <- sparseR(model_matrix = X2, y = cleveland$case, pre_process = FALSE, penalty = "MCP"))
  expect_silent(APL <- sparseR(model_matrix = X2, y = cleveland$case, pre_process = FALSE, gamma = 0, penalty = "MCP"))
  expect_silent(MEM <- sparseR(model_matrix = X1, y = y, pre_process = FALSE, family = "binomial", penalty = "MCP"))
  expect_silent(SRLp <- sparseR(model_matrix = X3, y = cleveland$case, pre_process = FALSE, penalty = "MCP"))

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

