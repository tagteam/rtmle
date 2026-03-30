### coef.rtmle.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul  3 2024 (13:48) 
## Version: 
## Last-Updated: mar 25 2026 (17:24) 
##           By: Thomas Alexander Gerds
##     Update #: 60
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
#' @export
#' @method coef rtmle
coef.rtmle <- function(object,time_horizon,...){
    if (!(object$learner$fun[[1]] %chin% c("learn_glmnet","learn_glm"))){
        stop("Can only plot regression coefficients when fitter is either learn_glm or learn_glmnet")
    }
    if (missing(time_horizon)) {
        time_horizon <- min(max(object$times),max(object$run_time_horizons))
    }
    df <- do.call(rbind,lapply(object$intervention_nodes,function(k){
        current_outcome <- paste0(object$names$outcome,"_",k+1)
        time_name <- paste0("time_",k)
        time_block <- object$models[[time_name]]
        do.call(rbind,lapply(names(time_block),function(node_name){
            if (node_name == "outcome") {
                model_names <- paste0("sequence_time_",time_horizon)
            }else{
                model_names <- names(time_block[[node_name]])
            }
            df <- do.call(rbind,lapply(model_names,function(this_outcome){
                if (node_name == "outcome"){
                    outfit <- do.call(rbind,lapply(names(time_block[[node_name]][[current_outcome]]$fit),function(protocol){
                        fit <- time_block[[node_name]][[current_outcome]]$fit[[protocol]][[this_outcome]]
                        if (!is.null(fit)
                            ## && 
                            ## (inherits(fit,"dgCMatrix") || inherits(fit,"data.frame") || inherits(fit,"matrix"))
                            && (is.numeric(coefs <- fit[,1]))){
                            names(coefs) <- rownames(fit)
                            data.table(time = time_name,
                                       protocol = protocol,
                                       node = node_name,
                                       outcome = current_outcome,
                                       terms = names(coefs),
                                       beta = as.numeric(coefs),
                                       stringsAsFactors = FALSE)
                        }
                    }))
                }else{
                    fit <- time_block[[node_name]][[this_outcome]]$fit
                    if (!is.null(fit)
                        ## && (inherits(fit,"dgCMatrix") || inherits(fit,"data.frame"))
                        && (is.numeric(coefs <- fit[,1]))){
                        names(coefs) <- rownames(fit)
                        outfit <- data.table(time = time_name,
                                             protocol = node_name,
                                             node = node_name,
                                             outcome = this_outcome,
                                             terms = names(coefs),
                                             beta = as.numeric(coefs),
                                             stringsAsFactors = FALSE)
                    } else{
                        outfit <- NULL
                    }
                }
            }))
        }))
    }))
    df[]
}
######################################################################
### coef.rtmle.R ends here
