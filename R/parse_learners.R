### parse_learners.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Nov  9 2024 (09:55) 
## Version: 
## Last-Updated: feb 26 2026 (13:27) 
##           By: Thomas Alexander Gerds
##     Update #: 113
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
##' Parsing list of learners for function superlearn. Available names are 
##' substituted for missing learner_fun elements and available elements learner_fun
##' are substituted for missing names.
##'
##' As a side effect an attribute is given to the object such that when running
##' on its own output the function does nothing.
##' @title Parsing list of learners
##' @param learners List of learners 
##' @return List of learners
##' @seealso \code{\link{superlearn}}, \code{\link{run_rtmle}}
##' @examples
##' parse_learners(c("learn_glm","learn_glmnet"))
##' parse_learners(list(folds=5,
##'                learners=
##'                         list("learn_glm",
##'                         "glm100"=list(maxit=200,learner_fun="learn_glm"))))
##' parse_learners(list(folds=10,
##'                learners=list(
##'                         "learn_glmnet",
##'                         "glm"=list(learn_variables="A",learner_fun="glm"),
##'                         list(name="learn_ranger", learner_fun="learn_ranger",num.trees=5))))
##' @export 
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
parse_learners <- function(learners){
    if (is.character(learners) && length(learners) == 1){
        learner_fun <- learners
        if (inherits(try(has_fun <- do.call("inherits",list(x = as.name(learner_fun),"function")),silent = TRUE),
                     "try-error")|| !has_fun){
            stop(paste0("Cannot parse the learner_fun as a function for learner '",learners,"'\nYou may need to load a library first?"))
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
                if (inherits(learners$learner_fun,"function")){
                    learner_fun <- learners$learner_fun
                } else{
                    if (!is.character(learners$learner_fun)){
                        stop("The learners must specify a learner function as element learner_fun.\n",
                             "Problem with learner: ",
                             paste(utils::capture.output(utils::str(learners)), collapse = "\n"))
                    }
                    learner_fun <- learners$learner_fun
                    if (inherits(try(has_fun <- do.call("inherits",list(x = as.name(learner_fun),"function")),silent = TRUE),
                                 "try-error")|| !has_fun){
                        stop(paste0("Cannot parse the learner_fun as a function for learner: ",
                                    learner_name,
                                    "\nYou may need to load a library first?"))
                    }
                }
                learner_args <- learners
                learner_args$name <- NULL
                learner_args$learner_fun <- NULL

                # NORMALIZE: ensure args is list() (never NULL)
                if (is.null(learner_args)) learner_args <- list()

                parsed_learners <- list(name = learner_name, fun = learner_fun, args = learner_args)
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
                    warning("Argument folds for super learner is missing but needed for cross-fitting.\nDefaults to 10-fold cross-validation.")
                    learner_args$folds <- 10
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
