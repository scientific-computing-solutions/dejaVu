context("Testing expandEventCount")

test_that("rejects non-numeric input",
          {
    expect_error(expandEventCount(c(0, 0), c("bad", "good")))
    expect_error(expandEventCount(c("bad", "good"), c(0, 0)))
})


test_that("rejects non-positive times",
          {
    expect_error(expandEventCount(c(0, 0), c(0, 1)))
})


test_that("rejects missing times",
          {
    expect_error(expandEventCount(c(0, 0), c(NA, 1)))
})


test_that("rejects negative event counts",
          {
    expect_error(expandEventCount(c(0, -1), c(1, 1)))
})


test_that("rejects missing event counts",
          {
    expect_error(expandEventCount(c(0, NA), c(1, 1)))
})


test_that("rejects nonconformant times and counts",
          {
    expect_error(expandEventCount(c(0, 1, 2), c(1, 2, 3, 4)))
})


test_that("warns when it builds data",
          {
    expect_warning(expandEventCount(c(0), c(1)))
})


test_that("builds data with correct lengths",
          {
    counts <- seq_len(10) - 1
    times  <- rep.int(1, 10)
    suppressWarnings(result <- expandEventCount(counts, times))
    expect_equal(vapply(result, length, 1), counts)
})


test_that("expands times correctly",
          {
    counts <- seq_len(10) - 1
    times  <- rep.int(2, 10)
    suppressWarnings(result <- expandEventCount(counts, 2))
    suppressWarnings(result.2 <- expandEventCount(counts, times))

    expect_equal(result, result.2)
})


test_that("caps times correctly",
          {
    counts <- seq_len(10)
    times  <- seq_len(length(counts))

    suppressWarnings(result <- expandEventCount(counts, times))
    expect_equal(vapply(result, max, 1), times)
})


test_that("event times are monotonically increasing for each subject",
          {
    counts <- seq_len(10)
    times <- seq_len(length(counts))
    
    suppressWarnings(result <- expandEventCount(counts, times))

    vapply(result, function(x) {all(diff(x) > 0)}, TRUE)

    expect_true(all(vapply(result, function(x) {all(diff(x) > 0)}, TRUE)))
})    
