### learner_output.R ---
#----------------------------------------------------------------------
## Internal helpers for learner return values
#----------------------------------------------------------------------
parse_learner_output <- function(x,
                                 predicted_values,
                                 diagnostics = NULL,
                                 fit = NULL,
                                 object,
                                 ...){
    if (missing(x)){
        out <- list(predicted_values = as.numeric(predicted_values),
                    diagnostics = diagnostics,
                    fit = fit)
        extras <- list(...)
        if (length(extras)>0){
            out <- c(out,extras)
        }
        if (!missing(object)){
            out$object <- object
        }
        return(out)
    }else{
        if (is.list(x) && "predicted_values"%in%names(x)){
            predicted_values <- x$predicted_values
            diagnostics <- x$diagnostics
            fit <- x$fit
            object <- x$object
            extras <- x[setdiff(names(x),c("predicted_values","diagnostics","fit","object"))]
            out <- c(list(predicted_values = as.numeric(predicted_values),
                          diagnostics = diagnostics,
                          fit = fit),
                     extras)
            if ("object"%in%names(x)){
                out$object <- object
            }
            return(out)
        }
        attrs <- attributes(x)
        extras <- attrs[setdiff(names(attrs),c("names","dim","dimnames","class","diagnostics","fit"))]
        c(list(predicted_values = as.numeric(x),
               diagnostics = attr(x,"diagnostics",exact = TRUE),
               fit = attr(x,"fit",exact = TRUE)),
          extras)
    }
}



######################################################################
### learner_output.R ends here
