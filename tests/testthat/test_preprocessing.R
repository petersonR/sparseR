context("Preprocessing Functionality")
skip_on_cran() # takes awhile, so only test on GitHub

data(iris)
iris <- iris[1:50,]

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

test_that("Different vals of k and poly work, general formula", {
  expect_silent({
    obj1 <- sparseR_prep(Sepal.Width ~ ., data = iris)
  })
  expect_silent({
    obj2 <- sparseR_prep(Sepal.Width ~ ., data = iris, k = 2, poly = 2)
  })
  expect_silent({
    obj3 <- sparseR_prep(Sepal.Width ~ ., data = iris, k = 1, poly = 1)
  })
  expect_silent({
    obj4 <- sparseR_prep(Sepal.Width ~ ., data = iris, k = 0, poly = 2)
  })
  expect_silent({
    obj5 <- sparseR_prep(Sepal.Width ~ ., data = iris, k = 5, poly = 5)
  })

  expect_error(sparseR_prep(Sepal.Width ~ ., data = iris, k = 1, poly = 0))
  expect_error(sparseR_prep(Sepal.Width ~ ., data = iris, k = 1, poly = NULL))
  expect_error(sparseR_prep(Sepal.Width ~ ., data = iris, k = 1, poly = 0))
  expect_error(sparseR_prep(Sepal.Width ~ ., data = iris, k = NULL, poly = NULL))

})

test_that("Different vals of k and poly work, specific formulae", {

  formula <- Sepal.Width ~ Petal.Width + BV + Petal.Length + Group
  expect_silent({
    obj1 <- sparseR_prep(formula, data = iris)
    obj2 <- sparseR_prep(formula, data = iris, k = 2, poly = 2)
    obj3 <- sparseR_prep(formula, data = iris, k = 1, poly = 1)
    obj4 <- sparseR_prep(formula, data = iris, k = 0, poly = 2)
    obj5 <- sparseR_prep(formula, data = iris, k = 5, poly = 5)
  })

  formula <- Sepal.Width ~ BV + UBV + Group
  expect_silent({
    obj1 <- sparseR_prep(formula, data = iris)
    obj2 <- sparseR_prep(formula, data = iris, k = 2, poly = 2)
    obj3 <- sparseR_prep(formula, data = iris, k = 1, poly = 1)
    obj4 <- sparseR_prep(formula, data = iris, k = 0, poly = 2)
    obj5 <- sparseR_prep(formula, data = iris, k = 5, poly = 5)
  })

  formula <- Sepal.Width ~ BV + UBV
  expect_warning({
    obj1 <- sparseR_prep(formula, data = iris)
    obj2 <- sparseR_prep(formula, data = iris, k = 2, poly = 2)
    obj3 <- sparseR_prep(formula, data = iris, k = 1, poly = 1)
    obj4 <- sparseR_prep(formula, data = iris, k = 0, poly = 2)
    obj5 <- sparseR_prep(formula, data = iris, k = 5, poly = 5)
  })

  expect_error(sparseR_prep(formula, data = iris, k = 1, poly = 0))
  expect_error(sparseR_prep(formula, data = iris, k = 1, poly = NULL))
  expect_error(sparseR_prep(formula, data = iris, k = .5, poly = 1))
  expect_error(sparseR_prep(formula, data = iris, k = NULL, poly = NULL))

})

## Test Detrano use-case

data("Detrano")
# smaller data set
cleveland <- cleveland[1:50,]

cleveland$thal <- factor(cleveland$thal)
cleveland$case <- 1*(cleveland$num > 0)

# Convert variables into factor variables if necessary!
cleveland$sex <- factor(cleveland$sex)
cleveland$fbs <- factor(cleveland$fbs)
cleveland$exang <- factor(cleveland$exang)

# Simulate missing data
cleveland$thal[2] <- cleveland$thalach[1] <- NA


test_that("Different vals of k and poly work, cleveland", {
  expect_silent({
    obj1 <- sparseR_prep(case ~ ., data = cleveland)
  })
  expect_silent({
    obj2 <- sparseR_prep(case ~ ., data = cleveland, k = 2, poly = 2)
  })
  expect_silent({
    obj3 <- sparseR_prep(case ~ ., data = cleveland, k = 1, poly = 1)
  })
  expect_silent({
    obj4 <- sparseR_prep(case ~ ., data = cleveland, k = 0, poly = 2)
  })


  expect_error(sparseR_prep(case ~ ., data = cleveland, k = 1, poly = 0))
  expect_error(sparseR_prep(case ~ ., data = cleveland, k = 1, poly = NULL))
  expect_error(sparseR_prep(case ~ ., data = cleveland, k = 1, poly = 0))
  expect_error(sparseR_prep(case ~ ., data = cleveland, k = NULL, poly = NULL))

})

test_that("Different vals of k and poly work, specific formulae", {

  formula <- case ~ trestbps + cp + thalach + thal
  expect_silent({
    obj1 <- sparseR_prep(formula, data = cleveland)
    obj2 <- sparseR_prep(formula, data = cleveland, k = 2, poly = 2)
    obj3 <- sparseR_prep(formula, data = cleveland, k = 1, poly = 1)
    obj4 <- sparseR_prep(formula, data = cleveland, k = 0, poly = 2)
    obj5 <- sparseR_prep(formula, data = cleveland, k = 5, poly = 5)
  })

  formula <- case ~ thal + num + chol
  expect_silent({
    obj1 <- sparseR_prep(formula, data = cleveland)
    obj2 <- sparseR_prep(formula, data = cleveland, k = 2, poly = 2)
    obj3 <- sparseR_prep(formula, data = cleveland, k = 1, poly = 1)
    obj4 <- sparseR_prep(formula, data = cleveland, k = 0, poly = 2)
    obj5 <- sparseR_prep(formula, data = cleveland, k = 5, poly = 5)
  })

  formula <- case ~ sex + thal
  expect_silent({
    obj1 <- sparseR_prep(formula, data = cleveland)
    obj2 <- sparseR_prep(formula, data = cleveland, k = 2, poly = 2)
    obj3 <- sparseR_prep(formula, data = cleveland, k = 1, poly = 1)
    obj4 <- sparseR_prep(formula, data = cleveland, k = 0, poly = 2)
    obj5 <- sparseR_prep(formula, data = cleveland, k = 5, poly = 5)
  })

  expect_error(sparseR_prep(formula, data = cleveland, k = 1, poly = 0))
  expect_error(sparseR_prep(formula, data = cleveland, k = 1, poly = NULL))
  expect_error(sparseR_prep(formula, data = cleveland, k = .5, poly = 1))
  expect_error(sparseR_prep(formula, data = cleveland, k = NULL, poly = NULL))

})

test_that("Centering to minimum works", {
  cc <- iris %>%
    dplyr::select(Sepal.Length, Petal.Length, Petal.Width) %>%
    apply(2, min, na.rm = TRUE)

  p1 <- sparseR_prep(Sepal.Width ~ ., iris, k = 0, extra_opts = list(centers = cc))
  c2min <- bake(p1, iris)
  p2 <- sparseR_prep(Sepal.Width ~ ., iris, k = 0, extra_opts = list(center_fn = min))
  c2min2 <- bake(p2, iris)
  expect_identical(c2min2, c2min)

  # testing print + tidy methods for center_to
  expect_output(print(p2$steps[[3]]))
  expect_silent(t1 <- tidy(p2$steps[[3]]))
})


