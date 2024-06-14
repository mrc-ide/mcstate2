test_that("can select sensible values for progress", {
  withr::with_options(list(mcstate2.progress = TRUE), {
    expect_false(show_progress_bar(FALSE))
    expect_true(show_progress_bar(TRUE))
    expect_true(show_progress_bar(NULL))
  })

  withr::with_options(list(mcstate2.progress = FALSE), {
    expect_false(show_progress_bar(FALSE))
    expect_true(show_progress_bar(TRUE))
    expect_false(show_progress_bar(NULL))
  })

  withr::with_options(list(mcstate2.progress = NULL), {
    expect_false(show_progress_bar(FALSE))
    expect_true(show_progress_bar(TRUE))
    expect_true(show_progress_bar(NULL))
  })
})


test_that("null progress bar does nothing", {
  p <- progress_bar(10, 10, FALSE)(1)
  expect_silent(p(1))
})


test_that("can format overall progress", {
  withr::local_options(cli.num_colors = 1)
  overall <- progress_overall(4, 100, TRUE)
  expect_equal(
    overall(c(100, 30, 10, 0)),
    paste0("[",
           cli::symbol$lower_block_8,
           cli::symbol$lower_block_3,
           cli::symbol$lower_block_1,
           cli::symbol$lower_block_1,
           "]"))
  expect_equal(
    overall(c(100, 30, 10, 100)),
    paste0("[",
           cli::symbol$lower_block_8,
           cli::symbol$lower_block_3,
           cli::symbol$lower_block_1,
           cli::symbol$lower_block_8,
           "]"))
})


test_that("overall progress is empty with one chain", {
  overall <- progress_overall(1, 100, TRUE)
  expect_equal(overall(0), "")
  expect_equal(overall(50), "")
  expect_equal(overall(100), "")
})


test_that("overall progress is empty if disabled", {
  overall <- progress_overall(4, 100, FALSE)
  expect_equal(overall(0), "")
  expect_equal(overall(50), "")
  expect_equal(overall(100), "")
})
