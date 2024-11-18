### parse_learners.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Nov  9 2024 (09:55) 
## Version: 
## Last-Updated: Nov 16 2024 (16:52) 
##           By: Thomas Alexander Gerds
##     Update #: 24
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
##' parse_learners(c("learn_glm","learn_glm"))
##' parse_learners(c(list("learn_glmnet",
##'                  "glm"=list(learn_variables="A")),list("learn_ranger"),
##'                  list(list(learner_fun="learn_ranger",num.trees=5))))
##' @export 
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
parse_learners <- function(learners){
    if (length(attr(learners,"parsed"))>0) return(learners)
    # make sure that all learners have names
    learner_names <- names(learners)
    if (length(learner_names) == 0) {
        names(learners) <- learner_names <- rep("noname_learner",length(learners))
    }
    if (length(missing_names <- which(is.na(names(learners)) | names(learners) == ""))>0){
        tmp_names <- sapply(missing_names,function(x){
            if (is.character(learners[[x]])) {
                learners[[x]]
            } else {
                "noname_learner"
            }
        })
        names(learners)[missing_names] <- tmp_names
    }
    # now checking learner_fun
    for (this_learner in 1:length(learners)){
        ## if no learner function is specified then assume that the name
        ## is a function
        if (is.character(learners[[this_learner]])){
            learners[[this_learner]] <- list(learner_fun = learners[[this_learner]])
        }else{
            if (length(learners[[this_learner]]$learner_fun) == 0){
                learners[[this_learner]]$learner_fun <- names(learners)[[this_learner]]
            }
        }
        # FIXME: could/should check if argumnents of the function
        #        match the arguments character_formula, intervened_data etc.
        #        could also check if the other provided arguments
        #        can be passed (requires that they are named!)
        if (inherits(learners[[this_learner]]$learner_fun,"function")){
            learners[[this_learner]]$learner_fun <- as.character(substitute(learners[[this_learner]]$learner_fun))
        } else{
            if (!is.character(learners[[this_learner]]$learner_fun))
                stop("All learners must specify a learner function either as the name or as element learner_fun.\n",
                     "Problem with learner: ",names(learners)[[this_learner]])
            if (inherits(try(has_fun <- do.call("inherits",list(x = as.name(learners[[this_learner]][["learner_fun"]]),"function")),silent = TRUE),
                         "try-error")|| !has_fun){
                stop(paste0("Cannot parse the learner_fun as a function for learner number: ",
                            this_learner,
                            "\nYou may need to load a library first?"))
            }
        }
    }
    # make sure names are unique
    names(learners) <- make.names(names(learners),unique = TRUE)
    data.table::setattr(learners,"parsed",TRUE)
    learners
}


######################################################################
### parse_learners.R ends here
