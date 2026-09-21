### test-cheap-bootstrap.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Apr  9 2025 (10:02) 
## Version: 
## Last-Updated: sep 21 2026 (11:13) 
##           By: Thomas Alexander Gerds
##     Update #: 29
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
library(testthat)
library(data.table)
library(rtmle)
test_that("Cheap bootstrap confidence intervals",{
    tau <- 2
    set.seed(37)
    ld <- simulate_long_data(n = 91,number_visits = 20,beta = list(A_on_Y = -.2,A0_on_Y = -0.3,A0_on_A = 6),register_format = TRUE)
    x <- rtmle_init(time_grid = seq(0,1500,30.45*12),name_id = "id",name_outcome = "Y",name_competing = "Dead",name_censoring = "Censored",censored_label = "censored")
    x <- add_long_data(x,
                    outcome_data=ld$outcome_data,
                    censored_data=ld$censored_data,
                    competing_data=ld$competing_data,
                    timevar_data=ld$timevar_data)
    x <- add_baseline_data(x,data=ld$baseline_data)
    x <- discretize_data(x,start_followup_date = 0)
    x <- regime(x,name = "Always_A",intervention = data.frame("time"=x$intervention_nodes,"A" = factor("1",levels = c("0","1"))),verbose = FALSE)
    x <- regime(x,name = "Never_A",intervention = data.frame("time"=x$intervention_nodes,"A" = factor("0",levels = c("0","1"))),verbose = FALSE)
    x <- prepare_rtmle_data(x)
    x <- target(x,name = "Outcome_risk",estimator = "tmle",regimes = c("Always_A","Never_A"))
    x <- model_formula(x)
    x <- run_rtmle(x,learner = "learn_glmnet",time_horizon = 1,verbose = FALSE)
    x <- run_rtmle(x,learner = "learn_glmnet",time_horizon = 1:tau,verbose = FALSE)
    x <- cheap_bootstrap(x,B = 2,M = 71)
    a = x$estimate$Main_analysis[,.(Time_horizon,Regime,Bootstrap_lower,Bootstrap_upper)]
    b = x$estimate$Cheap_bootstrap$Main_analysis[B == 2][,.(Time_horizon,Regime,Bootstrap_lower,Bootstrap_upper)]
    setkey(a,Time_horizon,Regime)
    setkey(b,Time_horizon,Regime)
    expect_equal(a,b)
})


test_that("cheap bootstrap matches a manual loop with and without replacement", {
    # Use the existing test fixture, including competing events and censoring.
    tau <- 2
    set.seed(37)
    ld <- simulate_long_data(
        n = 91, number_visits = 20,
        beta = list(A_on_Y = -.2, A0_on_Y = -.3, A0_on_A = 6),
        register_format = TRUE
    )
    x <- rtmle_init(time_grid = seq(0,1500,30.45*12),name_id = "id",name_outcome = "Y",name_competing = "Dead",name_censoring = "Censored",censored_label = "censored")
    x <- add_long_data(
        x, outcome_data = ld$outcome_data, censored_data = ld$censored_data,
        competing_data = ld$competing_data, timevar_data = ld$timevar_data
    )
    x <- add_baseline_data(x, data = ld$baseline_data)
    x <- long_to_wide(x,start_followup_date = 0)
    x <- protocol(x, name = "Always_A", intervention = data.frame(
        time = x$intervention_nodes, A = factor("1", levels = c("0", "1"))
    ), verbose = FALSE)
    x <- prepare_rtmle_data(x)
    x <- target(x, name = "Outcome_risk", estimator = "tmle",
                regimes = "Always_A")
    x <- model_formula(x, verbose = FALSE)
    x <- run_rtmle(x, learner = "learn_glmnet", time_horizon = 2)

    N <- nrow(x$prepared_data)
    M <- 71
    B <- 3
    seeds <- 101:103
    theta <- x$estimate$Main_analysis$Estimate

    for (replace in c(FALSE, TRUE)) {
        boot <- numeric(B)
        for (b in seq_len(B)) {
            set.seed(seeds[b])
            inbag <- sample.int(N, M, replace = replace)
            expect_equal(anyDuplicated(inbag) > 0L, replace)

            xb <- copy(x)
            xb$prepared_data <- xb$prepared_data[inbag]
            xb$followup <- xb$followup[inbag]
            xb$prepared_data[, id := seq_len(.N)]
            xb$followup[, id := seq_len(.N)]
            # Recompute protocol matching from the sampled data.
            xb$regimes$Always_A$intervention_match <- NULL
            xb$estimate <- NULL
            xb <- run_rtmle(xb, learner = "learn_glmnet",
                            time_horizon = 2, refit = TRUE)
            boot[b] <- xb$estimate$Main_analysis$Estimate
        }
        expect_true(all(is.finite(c(theta, boot))))

        # Calculate the interval independently of cheap_bootstrap().
        scale <- if (replace) sqrt(M / N) else sqrt(M / (N - M))
        half_width <- qt(.975, df = B) * scale * sqrt(mean((boot - theta)^2))

        result <- suppressMessages(cheap_bootstrap(
            copy(x), time_horizon = 2, B = B, M = M, seeds = seeds,
            replace = replace, add = FALSE
        ))
        draws <- result$estimate$Cheap_bootstrap$Main_analysis
        expect_equal(draws$Bootstrap_estimate[order(draws$B)], boot,
                     tolerance = 1e-6)
        expect_equal(result$estimate$Main_analysis$Bootstrap_lower,
                     pmax(0, theta - half_width), tolerance = 1e-6)
        expect_equal(result$estimate$Main_analysis$Bootstrap_upper,
                     pmin(1, theta + half_width), tolerance = 1e-6)
    }
})



######################################################################
### test-cheap-bootstrap.R ends here
