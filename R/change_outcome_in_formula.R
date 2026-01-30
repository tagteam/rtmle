### change_outcome_in_formula.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: jan 22 2026 (10:21) 
## Version: 
## Last-Updated: jan 22 2026 (10:34) 
##           By: Thomas Alexander Gerds
##     Update #: 2
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
change_outcome_in_formula <- function(character_formula,outcome_variable){
    if (inherits(character_formula,"formula")){
        fml <- character_formula
    }else{
        if (length(character_formula) != 1){
            stop("Argument character_formula has to be a single character string")
        }
        fml <- stats::formula(character_formula)
    }
    new_fml <- stats::update(fml,paste0(outcome_variable,"~."))
    out <- deparse1(new_fml, collapse = "")
    attr(out,"all_vars") <- all.vars(new_fml)
    out
}


######################################################################
### change_outcome_in_formula.R ends here
