### test_adding_data.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: okt  6 2026 (08:26) 
## Version: 
## Last-Updated: okt  6 2026 (08:37) 
##           By: Thomas Alexander Gerds
##     Update #: 4
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
library(rtmle)
library(testthat)
library(data.table)
test_that("add_long_data checks",{
    # check that entries without id and date or (start_date, end_date, value) or (date, value) are rejected
    x <- rtmle_init(time_grid = 0:2,
                    name_id = "id",
                    name_outcome = "Y",
                    name_competing = "D",
                    name_censoring = "C")

    expect_error(
        add_long_data(x, timevar_data = list(A = data.table(value = 1:2))),
        "id"
    )
    expect_error(
        add_long_data(x, timevar_data = list(A = data.table(id = 1:2))),
        "accepted formats"
    )
    expect_error(
        add_long_data(x, timevar_data = list(A = data.table(
            id = 1:2, date = 1:2, other = 1:2
        ))),
        "accepted formats"
    )
    expect_error(
        add_long_data(x, outcome_data = data.table(id = 1:2, value = 1:2)),
        "date"
    )

    # check that data are getting sorted according to date (start_date)
    x <- add_long_data(
        x,
        timevar_data = list(
            A = data.table(
                id = c(2, 1, 2, 1),
                date = c(4, 2, 3, 1),
                value = c(4, 2, 3, 1)
            ),
            exposure = data.table(
                id = c(2, 1, 2, 1),
                start_date = c(4, 2, 3, 1),
                end_date = c(5, 3, 4, 2)
            )
        )
    )
    expect_equal(x$long_data$timevar_data$A$id, c(1, 1, 2, 2))
    expect_equal(x$long_data$timevar_data$A$date, c(1, 2, 3, 4))
    expect_equal(x$long_data$timevar_data$exposure$id, c(1, 1, 2, 2))
    expect_equal(
        x$long_data$timevar_data$exposure$start_date,
        c(1, 2, 3, 4)
    )

    # check that in case of multiple outcome, censoring and competing events per subject only the first is kept
    x <- add_long_data(
        x,
        outcome_data = data.table(
            id = c(2, 1, 1),
            date = c(5, 3, 1)
        ),
        censored_data = data.table(
            id = c(2, 1, 1),
            date = c(5, 4, 2)
        ),
        competing_data = data.table(
            id = c(2, 1, 1),
            date = c(5, 6, 3)
        )
    )
    expect_equal(x$long_data$outcome_data$id, c(1, 2))
    expect_equal(x$long_data$outcome_data$date, c(1, 5))
    expect_equal(x$long_data$censored_data$id, c(1, 2))
    expect_equal(x$long_data$censored_data$date, c(2, 5))
    expect_equal(x$long_data$competing_data$id, c(1, 2))
    expect_equal(x$long_data$competing_data$date, c(3, 5))
})

######################################################################
### test_adding_data.R ends here
