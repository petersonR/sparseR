context("sparseR S3 methods")

data(iris)
set.seed(123)

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

## Test print method
test_that("print method works", {

  expect_no_message({
    obj2 <- sparseR(Sepal.Width ~ ., data = iris)
    capture.output(print(obj2))
  })

  expect_output({
    obj1 <- sparseR(Sepal.Width ~ ., data = iris)
    print(obj1)
  })
})


## Test summary method
test_that("summary method works", {

  expect_no_message({
    obj2 <- sparseR(Sepal.Width ~ ., data = iris)
    summary(obj2)
  })

  expect_output({
    obj1 <- sparseR(Sepal.Width ~ ., data = iris)
    print(summary(obj1))
  })
})

## Predict method test
test_that("predict method works", {

  expect_length({
    obj2 <- sparseR(Sepal.Width ~ ., data = iris)
    predict(obj2)
  }, nrow(iris))

  expect_gt({
    obj2 <- sparseR(Sepal.Width ~ ., data = iris)
    cor(predict(obj2), iris$Sepal.Width)
  }, 0.82)

  expect_gt({
    obj2 <- sparseR(Sepal.Width ~ ., data = iris)
    cor(predict(obj2, at = "cv1se"), iris$Sepal.Width)
  }, .78)
})

## Coef method works
test_that("coef method works", {
  expect_equal({
    obj2 <- sparseR(Sepal.Width ~ ., data = iris)
    b <- coef(obj2, at = "cvmin")
    sum(b != 0)
  }, 18)

  expect_equal({
    sum(coef(obj2, at = "cv1se") != 0)
  }, 11)
})

