### delete_variables_from_formula.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Nov  3 2024 (14:39) 
## Version: 
## Last-Updated: Nov  3 2024 (14:45) 
##           By: Thomas Alexander Gerds
##     Update #: 3
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
delete_variables_from_formula  <- function(character_formula, delete_vars) {
    allvars <- all.vars(stats::formula(character_formula))
    remaining_vars <- setdiff(allvars[-1],delete_vars)
    paste(allvars[[1]],
          "~",
          paste(remaining_vars, collapse = " + "))
}             
######################################################################
### delete_variables_from_formula.R ends here
