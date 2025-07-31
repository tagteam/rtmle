### delete_variables_from_formula.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Nov  3 2024 (14:39) 
## Version: 
## Last-Updated: Jul 31 2025 (07:34) 
##           By: Thomas Alexander Gerds
##     Update #: 25
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
#
# delete_variables_from_formula("Y~x+Z","Z")
# 
delete_variables_from_formula  <- function(character_formula,
                                           delete_vars) {
    if (inherits(character_formula,"formula")){
        fml <- character_formula
    }else{
        if (length(character_formula) != 1){
            stop("Argument character_formula has to be a single character string")
        }
        fml <- stats::formula(character_formula)
    }
    rhs_terms <- setdiff(attr(stats::terms(fml),"term.labels"),delete_vars)
    if ((nvars <-  length(rhs_terms)) == 0) {
        rhs_terms <- "1"
    }
    out <- deparse1(stats::reformulate(rhs_terms,response = deparse(fml[[2]])),
                    collapse = "")
    attr(out,"number_rhs_variables") <- nvars
    out
}             
######################################################################
### delete_variables_from_formula.R ends here
