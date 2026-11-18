test_that("the continuous estimating equation validates input dimensions", {
    X <- matrix(1, nrow = 5, ncol = 1)
    Y <- rep(0.5, 5)
    offset <- rep(0, 5)

    expect_error(estimating_equation_cpp(X, Y[-1], 0, offset),
                 "Length of Y must match number of rows in X")
    expect_error(estimating_equation_cpp(X, Y, numeric(2), offset),
                 "Length of beta must match number of columns in X")
    expect_error(estimating_equation_cpp(X, Y, 0, offset[-1]),
                 "Length of offset must match number of rows in X")
    expect_error(estimating_equation_cpp(X, Y, 0, offset, rep(1, 4)),
                 "Length of weights_ must match number of rows in X")
})

test_that("the continuous estimating equation accepts vector weights", {
    X <- matrix(1, nrow = 5, ncol = 1)
    Y <- seq(0.1, 0.5, length.out = 5)
    weights <- as.numeric(scale(seq_len(5), center = FALSE))
    fit <- estimating_equation_cpp(X, Y, 0, rep(0, 5), weights)

    expect_type(fit, "double")
    expect_length(fit, 1)
    expect_equal(as.numeric(stats::plogis(fit)), stats::weighted.mean(Y, weights),
                 tolerance = 1e-8)
})

test_that("the continuous estimating equation matches weighted logistic regression", {
    d <- data.frame(Y = rep(c(0.1, 0.4, 0.6, 0.8, 1), 5),
                    L = seq(-1, 1, length.out = 25),
                    offset = sin(seq_len(25)),
                    weights = seq_len(25))
    X <- stats::model.matrix(~ L, d)
    reference <- stats::glm(Y ~ L + offset(offset), data = d,
                            weights = weights, family = stats::quasibinomial())
    fit <- estimating_equation_cpp(X, d$Y, numeric(2), d$offset, d$weights)

    expect_equal(as.numeric(fit), as.numeric(stats::coef(reference)),
                 tolerance = 1e-7)
    score <- crossprod(X, d$weights *
                        (d$Y - stats::plogis(X %*% fit + d$offset)))
    expect_lt(max(abs(score)), 1e-7)
})

test_that("the continuous estimating equation allows IPCW outcomes above one", {
    Y <- c(0, 0, 0.5, 1.5)
    fit <- estimating_equation_cpp(matrix(1, 4, 1), Y, 0, numeric(4))
    expect_equal(as.numeric(stats::plogis(fit)), mean(Y), tolerance = 1e-8)
})

test_that("the continuous estimating equation handles boundary outcomes", {
    X <- matrix(1, 5, 1)
    expect_no_warning(fit_zero <- estimating_equation_cpp(X, numeric(5), 0,
                                                         numeric(5)))
    expect_no_warning(fit_one <- estimating_equation_cpp(X, rep(1, 5), 0,
                                                        numeric(5)))
    expect_lt(stats::plogis(fit_zero)[1], 1e-8)
    expect_gt(stats::plogis(fit_one)[1], 1 - 1e-8)

    # An already solved equation needs no invertible information matrix.
    expect_no_warning(fit <- estimating_equation_cpp(X, rep(1, 5), 0,
                                                    rep(1000, 5)))
    expect_equal(as.numeric(fit), 0)
})

test_that("the continuous estimating equation rejects invalid numeric input", {
    X <- matrix(1, 5, 1)
    Y <- rep(0.5, 5)
    offset <- numeric(5)
    expect_error(estimating_equation_cpp(matrix(numeric(), 0, 1), numeric(),
                                         0, numeric()), "non-empty")
    expect_error(estimating_equation_cpp(matrix(numeric(), 5, 0), Y,
                                         numeric(), offset), "non-empty")
    expect_error(estimating_equation_cpp(replace(X, 1, NA), Y, 0, offset), "finite")
    expect_error(estimating_equation_cpp(X, replace(Y, 1, Inf), 0, offset), "finite")
    expect_error(estimating_equation_cpp(X, Y, NA_real_, offset), "finite")
    expect_error(estimating_equation_cpp(X, Y, 0, replace(offset, 1, NA)), "finite")
    expect_error(estimating_equation_cpp(X, Y, 0, offset, c(-1, 1, 1, 1, 1)),
                 "non-negative")
    expect_error(estimating_equation_cpp(X, Y, 0, offset, rep(0, 5)), "positive")
    expect_error(estimating_equation_cpp(X, Y, 0, offset, rep(Inf, 5)), "finite")
    expect_error(estimating_equation_cpp(X, Y, 0, offset, maxit = 0), "maxit")
    expect_error(estimating_equation_cpp(X, Y, 0, offset, tol = 0), "tol")
    expect_error(estimating_equation_cpp(X, Y, 0, offset, tol = NA_real_), "tol")
})

test_that("the continuous estimating equation safeguards Newton steps", {
    X <- matrix(1, 10, 1)
    Y <- rep(0.2, 10)
    fit <- estimating_equation_cpp(X, Y, 5, numeric(10))
    expect_equal(as.numeric(fit), stats::qlogis(0.2), tolerance = 1e-7)
})

test_that("the continuous estimating equation reports unsolved roots", {
    X <- matrix(1, 5, 1)
    expect_error(estimating_equation_cpp(X, rep(0.2, 5), 0, numeric(5),
                                         maxit = 1), "did not converge")
    expect_error(estimating_equation_cpp(cbind(X, X), rep(0.2, 5), c(0, 0),
                                         numeric(5)), "singular")
    expect_error(estimating_equation_cpp(X, rep(2, 5), 0, numeric(5)),
                 "could not be solved|did not converge")
})

test_that("the continuous estimating equation converges from nonzero coefficients", {
    set.seed(314)
    for (iteration in seq_len(5)) {
        X <- cbind(1, matrix(stats::rnorm(150), 50, 3))
        Y <- stats::runif(50)
        weights <- stats::runif(50, 0.1, 5)
        offset <- stats::rnorm(50)
        reference <- stats::glm.fit(X, Y, weights = weights, offset = offset,
                                    family = stats::quasibinomial(),
                                    control = stats::glm.control(epsilon = 1e-12))
        fit <- estimating_equation_cpp(X, Y, stats::rnorm(4), offset, weights)
        expect_equal(as.numeric(fit), as.numeric(stats::coef(reference)),
                     tolerance = 1e-7)
    }
})
