#
# lava model
#
make_regression_model <- function(outcome_variables,parameter_values) {
    m <- lava::lvm()
    for (v in names(outcome_variables)){
        if (outcome_variables[[v]] == "constant"){
            lava::distribution(m, v) <- lava::constant.lvm(value = parameter_values[[paste0("intercept_",v)]])
        }
        if (outcome_variables[[v]] == "binomial"){
            lava::distribution(m, v) <- lava::binomial.lvm("logit")
            lava::intercept(m,v) <- parameter_values[[paste0("intercept_",v)]]
        }
        if (outcome_variables[[v]] == "normal"){
            lava::distribution(m, v) <- lava::normal.lvm()
            if (length(intercept_v <- parameter_values[[paste0("intercept_",v)]]) == 0){
                # if not specified the intercept is zero
                intercept_v <- 0
            }
            lava::intercept(m,v) <- intercept_v
            if (length((var_v = parameter_values[[paste0("var_",v)]]) == 1)){
                lava::variance(m,v) <- (var_v)^2
            }
        }
        if (outcome_variables[[v]] == "lognormal"){
            lava::distribution(m, v) <- lava::lognormal.lvm()
            lava::intercept(m,v) <- parameter_values[[paste0("intercept_",v)]]
            if (length((var_v <- parameter_values[[paste0("var_",v)]]) == 1)){
                lava::variance(m,v) <- (var_v)^2
            }
        }
        if (outcome_variables[[v]] %chin% c("Exponential","Weibull")){
            if (length(parameter_values[[paste0("scale_",v)]]) == 0){
                stop("rtmle::make_regression_model: missing parameter ",paste0("scale_",v))
            }
            if (outcome_variables[[v]] == "Exponential"){shape <- 1} else {shape <- 2}
            lava::distribution(m, v) <- lava::coxWeibull.lvm(
                shape = shape,
                scale = parameter_values[[paste0("scale_",v)]]
            )
        }
        # Regression parameters (if any)
        v_effects <- parameter_values[grep(paste0("^effect_.*_",v,"$"),names(parameter_values),value = TRUE)]
        v_effects <- v_effects[v_effects != 0]
        v_effects <- setNames(v_effects,sub("^effect_","",names(v_effects)))
        v_effects <- setNames(v_effects,sub(paste0("_",v,"$"),"",names(v_effects)))
        if (length(v_effects)>0){
            lava::regression(m) <- formula(paste0(v," ~ ",paste0(sapply(names(v_effects),function(e){paste0("f(",e,",",v_effects[[e]],")")}),collapse = "+")))
        }
    }
    m
}
