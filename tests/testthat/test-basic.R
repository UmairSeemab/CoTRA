test_that("installed app is present", {
  expect_true(dir.exists(system.file("app", package = "CoTRA")))
  expect_true(file.exists(system.file("app", "app.R", package = "CoTRA")))
})

test_that("dependency checker returns expected structure", {
  deps <- check_cotra_dependencies(quiet = TRUE)
  expect_named(deps, c("ok", "installed", "missing"))
  expect_type(deps$ok, "logical")
})
