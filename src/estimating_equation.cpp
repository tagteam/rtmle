// Migrated from contICEIPCW's logistic estimating-equation solver.
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// Evaluate expit without overflowing for large negative linear predictors.
arma::vec expit(const arma::vec& eta) {
    arma::vec mu(eta.n_elem);
    for (arma::uword i = 0; i < eta.n_elem; ++i) {
        if (eta[i] >= 0) {
            mu[i] = 1.0 / (1.0 + std::exp(-eta[i]));
        } else {
            const double exponential = std::exp(eta[i]);
            mu[i] = exponential / (1.0 + exponential);
        }
    }
    return mu;
}

// The logistic objective is concave even for IPCW responses outside [0, 1].
double logistic_objective(const arma::vec& eta, const arma::vec& Y,
                          const arma::vec& weights) {
    double objective = 0;
    for (arma::uword i = 0; i < eta.n_elem; ++i) {
        if (eta[i] >= 0)
            objective += weights[i] * ((Y[i] - 1.0) * eta[i] -
                                       std::log1p(std::exp(-eta[i])));
        else
            objective += weights[i] * (Y[i] * eta[i] -
                                       std::log1p(std::exp(eta[i])));
    }
    return objective;
}

// Solve X' weights (Y - expit(X beta + offset)) = 0. IPCW responses
// can exceed one, so Y is deliberately not restricted to [0, 1].
// [[Rcpp::export]]
arma::vec estimating_equation_cpp(
    const arma::mat& X,
    const arma::vec& Y,
    arma::vec beta,
    arma::vec offset,
    Rcpp::Nullable<arma::vec> weights_ = R_NilValue,
    int maxit = 100,
    double tol = 1e-8,
    bool verbose = false) {
    const arma::uword n = X.n_rows;
    const arma::uword p = X.n_cols;
    if (n == 0 || p == 0)
        Rcpp::stop("X must be a non-empty matrix");
    if (beta.n_elem != p)
        Rcpp::stop("Length of beta must match number of columns in X");
    if (Y.n_elem != n)
        Rcpp::stop("Length of Y must match number of rows in X");
    if (offset.n_elem != n)
        Rcpp::stop("Length of offset must match number of rows in X");
    if (!X.is_finite() || !Y.is_finite() || !beta.is_finite() ||
        !offset.is_finite())
        Rcpp::stop("X, Y, beta, and offset must contain only finite values");
    if (maxit < 1)
        Rcpp::stop("maxit must be positive");
    if (!std::isfinite(tol) || tol <= 0)
        Rcpp::stop("tol must be finite and positive");

    arma::vec weights(n, arma::fill::ones);
    if (weights_.isNotNull()) {
        weights = Rcpp::as<arma::vec>(weights_);
        if (weights.n_elem != n)
            Rcpp::stop("Length of weights_ must match number of rows in X");
        if (!weights.is_finite())
            Rcpp::stop("weights_ must contain only finite values");
        if (arma::any(weights < 0))
            Rcpp::stop("weights_ must be non-negative");
        if (weights.max() == 0)
            Rcpp::stop("weights_ must contain a positive value");
        // A common scale does not change the solution.
        weights /= weights.max();
    }

    for (int iteration = 0; iteration < maxit; ++iteration) {
        Rcpp::checkUserInterrupt();
        const arma::vec eta = X * beta + offset;
        if (!eta.is_finite())
            Rcpp::stop("The linear predictor is not finite");
        const arma::vec mu = expit(eta);
        const arma::vec score = X.t() * (weights % (Y - mu));
        const double score_norm = arma::max(arma::abs(score));
        if (score_norm < tol)
            return beta;

        const arma::vec variance_weights = weights % mu % (1.0 - mu);
        const arma::mat information = X.t() * (X.each_col() % variance_weights);
        arma::vec step;
        if (!arma::solve(step, information, score, arma::solve_opts::no_approx) ||
            !step.is_finite())
            Rcpp::stop("The estimating equation could not be solved: singular information matrix");

        // Newton steps can overshoot even when all coefficients are finite.
        // Halve the step until it improves the concave logistic objective.
        const double objective = logistic_objective(eta, Y, weights);
        const double directional_derivative = arma::dot(score, step);
        bool improved = false;
        double step_factor = 1.0;
        for (int attempt = 0; attempt < 30; ++attempt) {
            const arma::vec candidate = beta + step_factor * step;
            const arma::vec candidate_eta = X * candidate + offset;
            if (candidate.is_finite() && candidate_eta.is_finite()) {
                const double candidate_objective =
                    logistic_objective(candidate_eta, Y, weights);
                // Near a root, improvement can be smaller than floating-point
                // rounding in the summed objective.
                const double rounding_tolerance =
                    1e-12 * (1.0 + std::abs(objective));
                if (std::isfinite(candidate_objective) &&
                    candidate_objective + rounding_tolerance >= objective +
                        1e-4 * step_factor * directional_derivative) {
                    beta = candidate;
                    improved = true;
                    break;
                }
            }
            step_factor *= 0.5;
        }
        if (!improved)
            Rcpp::stop("The estimating equation did not converge: no improving Newton step");
        if (verbose)
            Rcpp::Rcout << "Iteration " << iteration + 1
                        << ": score norm = " << score_norm << '\n';
    }
    // The last permitted step may have solved the equation.
    const arma::vec final_score = X.t() *
        (weights % (Y - expit(X * beta + offset)));
    if (arma::max(arma::abs(final_score)) < tol)
        return beta;
    Rcpp::stop("The estimating equation did not converge within maxit iterations");
}
