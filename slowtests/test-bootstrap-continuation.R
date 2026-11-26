library(testthat)
library(rtmle)

test_that("appending subsamples reproduces a single seeded bootstrap run", {
    data("rtmle_object", package = "rtmle", envir = environment())
    boot <- function(x, B, ...) {
        cheap_bootstrap(x, targets = "Treatment effect", time_horizon = 1,
                        B = B, seeds = 431, ...)
    }
    first <- boot(rtmle_object, 2)
    expect_error(boot(first, 1, M = nrow(first$prepared_data) - 1L), "add = FALSE")
    combined <- boot(first, 1)
    all_at_once <- boot(rtmle_object, 3)
    expect_equal(combined$estimate$Cheap_bootstrap, all_at_once$estimate$Cheap_bootstrap)
    expect_equal(combined$estimate$Main_analysis, all_at_once$estimate$Main_analysis)
})
