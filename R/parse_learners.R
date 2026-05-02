### parse_learners.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Nov  9 2024 (09:55) 
## Version: 
## Last-Updated: apr 10 2026 (11:34) 
##           By: Thomas Alexander Gerds
##     Update #: 125
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
##' Parse learner specifications
##'
##' Normalizes learner specifications for \code{\link{run_rtmle}} and
##' \code{\link{superlearn}}. Available names are used to fill missing
##' \code{fun} elements, and available \code{fun} elements are used to fill
##' missing names. The returned object receives an attribute so repeated parsing
##' leaves it unchanged.
##'
##' @title Parse learner specifications
##' @param learners A learner specification: a character vector of learner names,
##'   a single learner list, or a super-learner list containing \code{folds} and
##'   \code{learners}.
##' @return A normalized learner list.
##' @seealso \code{\link{superlearn}}, \code{\link{run_rtmle}}
##' @examples
##' parse_learners(c("learn_glm","learn_glmnet"))
##' parse_learners(list(folds=5,
##'                learners=
##'                         list("learn_glm",
##'                         "glm100"=list(maxit=200,fun="learn_glm"))))
##' parse_learners(list(folds=10,
##'                learners=list(
##'                         "learn_glmnet",
##'                         "glm"=list(learn_variables="A",fun="glm"),
##'                         list(name="learn_ranger", fun="learn_ranger",num.trees=5))))
##' @export 
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
parse_learners <- function(learners){
    ## Make sure that formatting is obeyed
    ## 1. If learners is a single string
    ## 2. Learners is a vector of strings
    ## 3. Learners does not contain a list of learners, but must have a name and a fun element
    ## 4. Learners contains a list of learners named learners and a folds element for the super learner.
    # learners must be one of the allowed formats
    if (is.character(learners)) {
        # Case 1: single string OR Case 2: vector of strings
        if (!all(nzchar(learners))) {
            stop("All learner names must be non-empty strings.")
        }
    } else if (is.list(learners)) {
        # Case 3: single learner specification
        if (!("learners" %in% names(learners))) {
            if (!all(c("name", "fun") %in% names(learners))) {
                stop("rtmle::parse_learners: Argument 'learners' can either be a list of learners or a single learner which must be a list that contains 'name' and 'fun'.")
            }
        } else {
            # Case 4: super learner specification
            if (!is.list(learners$learners)) {
                stop("'learners' must be a list of learner specifications.")
            }
            if (!is.numeric(learners$folds)) {
                stop("'folds' must be numeric.")
            }
        }
    } else {
        stop("Invalid 'learners' specification.")
    }
    if (is.character(learners) && length(learners) == 1){
        fun <- learners
        if (inherits(try(has_fun <- do.call("inherits",list(x = as.name(fun),"function")),silent = TRUE),
                     "try-error")|| !has_fun){
            stop(paste0("Cannot parse the fun as a function for learner '",learners,"'\nYou may need to load a library first?"))
        }
        # NORMALIZE: use list() instead of NULL so it matches list-form learners
        parsed_learners <- list(name = learners, fun = learners, args = list())
    } else {
        if (is.character(learners)) {
            parsed_learners <- list(
                name = "superlearn",
                learners = lapply(learners, FUN = parse_learners),
                # NORMALIZE: include learner="superlearn" (so list-form and char-vector form align)
                args = list(learner = "superlearn", folds = 10)
            )
        } else {
            stopifnot(is.list(learners))
            if (is.null(learners$learners)){
                #
                # a single learner, not a super learner
                #
                if (length(learners$name)>0){
                    learner_name <- learners$name
                }else{
                    learner_name <- "unnamed"
                }
                if (inherits(learners$fun,"function")){
                    fun <- learners$fun
                } else{
                    if (!is.character(learners$fun)){
                        stop("The learners must specify a learner function as element fun.\n",
                             "Problem with learner: ",
                             paste(utils::capture.output(utils::str(learners)), collapse = "\n"))
                    }
                    fun <- learners$fun
                    if (inherits(try(has_fun <- do.call("inherits",list(x = as.name(fun),"function")),silent = TRUE),
                                 "try-error")|| !has_fun){
                        stop(paste0("Cannot parse the fun as a function for learner: ",
                                    learner_name,
                                    "\nYou may need to load a library first?"))
                    }
                }
                learner_args <- learners
                learner_args$name <- NULL
                learner_args$fun <- NULL

                # NORMALIZE: ensure args is list() (never NULL)
                if (is.null(learner_args)) learner_args <- list()

                parsed_learners <- list(name = learner_name, fun = fun, args = learner_args)
            } else {
                # super learner
                learner_args <- learners

                # make sure that all learners have names
                list_names <- names(learners$learners)
                for (l in 1:length(learners$learners)){
                    if (!is.character(learners$learners[[l]])){
                        if (length(learners$learners[[l]]$name) == 0){
                            if (length(list_names[[l]])>0){
                                learners$learners[[l]]$name <- list_names[[l]]
                            } else {
                                learners$learners[[l]]$name <- paste0("noname_learner_",l)
                            }
                        }
                    }
                }

                learners <- lapply(X = learners$learners, FUN = parse_learners)

                learner_args$learners <- NULL
                learner_args$name <- NULL

                # NORMALIZE: if learner missing, set it
                if (length(learner_args$learner) == 0){
                    learner_args$learner <- "superlearn"
                }

                if (length(learner_args$folds) == 0){
                    warning("parse_learners: Argument `folds' of function super_learn is missing but needed for cross-fitting.\nDefaults to 10-fold cross-validation.")
                    learner_args$folds <- 10
                }
                if (length(learner_args$ensemble_method) == 0){
                    warning("parse_learners: Argument `ensemble_method' of function super_learn is missing. Defaults to `discrete' which implements a discrete superlearner.")
                    learner_args$ensemble_method <- "discrete"
                }

                parsed_learners <- list(
                    name = "superlearn",
                    learners = learners,
                    args = learner_args
                )
            }
        }
    }
    parsed_learners
}
######################################################################
### parse_learners.R ends here
