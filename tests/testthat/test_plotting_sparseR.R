context("Plotting")
data("iris")

set.seed(12321)

obj1 <- sparseR(Sepal.Width ~ ., data = iris)
obj2 <- sparseR(Sepal.Width ~ ., k = 0, data = iris)
obj3 <- sparseR(Sepal.Width ~ ., k = 2, data = iris)

test_that("Plotting sparseR runs without error", {
  expect_silent(plot(obj1))
  expect_silent(plot(obj2))
  expect_silent(plot(obj3))
  expect_silent(plot(obj1, plot_type = "cv"))
  expect_silent(plot(obj1, plot_type = "path"))
})

test_that("Effect plots for sparseR run without error", {
  expect_warning(effect_plot(obj1, coef_name = "Species"))
  expect_silent(effect_plot(obj1, coef_name = "Petal.Width", by = "Species"))
  expect_warning(effect_plot(obj1, coef_name = "Petal.Width"))
  expect_silent(effect_plot(obj1, coef_name = "Petal.Length", by = "Sepal.Length"))

  expect_silent(effect_plot(obj2, coef_name = "Petal.Width"))
  expect_silent(effect_plot(obj2, coef_name = "Species"))
  expect_silent(effect_plot(obj2, "Petal.Width", by = "Species"))

  expect_warning(effect_plot(obj3, coef_name = "Petal.Width"))
  expect_warning(effect_plot(obj3, coef_name = "Species"))
  expect_silent(effect_plot(obj3, "Petal.Width", by = "Species"))

  expect_error(effect_plot(obj2, "Species", by = "Petal.Width"))
})

obj1 <- sparseRBIC_step(Sepal.Width ~ ., data = iris, message = FALSE)
obj2 <- sparseRBIC_step(Sepal.Width ~ ., k = 0, data = iris, message = FALSE)

test_that("Effect plots for sparseRBIC run without error", {
  expect_warning(effect_plot(obj1, coef_name = "Species"))
  expect_silent(effect_plot(obj1, "Petal.Width", by = "Species")) # might need to update to warn
  expect_silent(effect_plot(obj1, coef_name = "Petal.Width"))
  expect_silent(effect_plot(obj1, "Petal.Length", by = "Sepal.Length")) # might need to update to warn

  expect_silent(effect_plot(obj1, "Sepal.Length", by = "Species"))
  expect_warning(effect_plot(obj1, coef_name = "Sepal.Length"))

  expect_silent(effect_plot(obj2, coef_name = "Petal.Width"))
  expect_silent(effect_plot(obj2, coef_name = "Species"))
  expect_silent(effect_plot(obj2, "Petal.Width", by = "Species"))

  expect_error(effect_plot(obj2, "Species", by = "Petal.Width"))
})


