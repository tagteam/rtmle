### parse_learners.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Nov  9 2024 (09:55) 
## Version: 
## Last-Updated: dec  2 2025 (09:35) 
##           By: Thomas Alexander Gerds
##     Update #: 104
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
        parsed_learners <- list(name = learners,fun = learners,args = NULL)
    }else{
        if (is.character(learners)) {
            parsed_learners <- lapply(learners,FUN = parse_learners)
        } else {
            # FIXME: write a helpful error message
            stopifnot(is.list(learners))
            if (is.null(learners$learners)){
                #
                # a single learner, not a super learner
                # 
                learner_name <- learners$name
                if (inherits(learners$learner_fun,"function")){
                    learner_fun <- learners$learner_fun # as.character(substitute(learners$learner_fun))
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
                # the remaining must be arguments of learner_fun
                learner_args <- learners
                # FIXME: could/should check if arguments of the function
                #        match the arguments character_formula, intervened_data etc.
                #        could also check if the other provided arguments
                #        can be passed (requires that they are named!)
                learner_args$name <- NULL
                learner_args$learner_fun <- NULL
                parsed_learners <- list(name = learner_name,fun = learner_fun,args = learner_args)
            }else{
                # super learner
                # make sure that all learners have names
                learner_args <- learners
                # learners may be a named list
                # and these names should be used if the learner
                # does not have an argument 'name'
                list_names <- names(learners$learners)
                for (l in 1:length(learners$learners)){
                    if (!is.character(learners$learners[[l]])){
                        if (length(learners$learners[[l]]$name) == 0){
                            if (length(list_names[[l]])>0){
                                learners$learners[[l]]$name <- list_names[[l]]
                            }else{
                                learners$learners[[l]]$name <- paste0("noname_learner_",l)
                            }
                        }
                    }
                }
                learners <- lapply(X = learners$learners,FUN = parse_learners)
                learner_args$learners <- NULL
                learner_args$name <- NULL
                if (length(learner_args$folds) == 0){
                    stop("Argument folds for super learner is missing")
                }
                parsed_learners <- list(name = "superlearn",
                                        learners = learners,
                                        args = learner_args)
            }
        }
    }
    parsed_learners
}


######################################################################
### parse_learners.R ends here
