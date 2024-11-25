### test-rtmle-versus-ltmle.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Nov 16 2024 (17:04) 
## Version: 
## Last-Updated: Nov 25 2024 (08:06) 
##           By: Thomas Alexander Gerds
##     Update #: 31
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
library(testthat)
library(ltmle)
library(rtmle)
library(data.table)

test_that("single time point compare rtmle with ltmle",{
    set.seed(8)
    rexpit <- function(x) rbinom(n=length(x), size=1, prob=plogis(x))
    # Single time point Example
    n <- 1000
    W <- rnorm(n)
    A <- rexpit(-1 + 2 * W)
    Y <- rexpit(W + A)
    data <- data.frame(W, A, Y)
    result1 <- ltmle(data, Anodes="A", Ynodes="Y", abar=1,estimate.time = FALSE,gbounds = c(0,1))
    summary(result1)
    x <- rtmle_init(intervals = 1,name_id = "id",name_outcome = "Y",name_competing = NULL,name_censoring = NULL)
    rdata <- cbind(id = 1:NROW(data),data)
    setDT(rdata)
    setnames(rdata,c("id","W","A_0","Y_1"))
    x$prepared_data <- rdata
    x$names$name_baseline_covariates <- "W"
    x$names$name_time_covariates <- "A"
    protocol(x) <- list(name = "A",treatment_variables = "A",intervention = 1)
    target(x) <- list(name = "Outcome_risk",strategy = "additive",estimator = "tmle",protocols = "A")
    x <- run_rtmle(x,refit = TRUE)
    expect_equal(result1$fit$g[["A"]], as.matrix(x$models[[1]][["A_0"]]$fit))
    ## all.equal(result1$fit$Q[["Y"]], as.matrix(x$models[[1]][["Y_1"]]$fit))
    expect_equal(result1$estimates[["tmle"]],x$estimate$Outcome_risk$A$Estimate,tolerance = 0.0001)
    expect_equal(result1$IC$tmle,x$IC$Outcome_risk$A$time_horizon_1,tolerance = 0.0001)
})

test_that("longitudinal data compare rtmle with ltmle",{
    set.seed(17)
    ld <- simulate_long_data(n = 1291,number_visits = 20,beta = list(A_on_Y = -.2,A0_on_Y = -0.3,A0_on_A = 6),register_format = TRUE)
    x <- rtmle_init(intervals = 3,name_id = "id",name_outcome = "Y",name_competing = "Dead",name_censoring = "Censored",censored_label = "censored")
    x$long_data <- ld[c("outcome_data","censored_data","competing_data","timevar_data")]
    add_baseline_data(x) <- ld$baseline_data[,start_followup_date:=0]
    x <- long_to_wide(x,intervals = seq(0,2000,30.45*12))
    protocol(x) <- list(name = "Always_A",treatment_variables = "A",intervention = 1)
    prepare_data(x) <- list()
    target(x) <- list(name = "Outcome_risk",strategy = "additive",estimator = "tmle",protocols = "Always_A")
    x <- run_rtmle(x,learner = "learn_glm",time_horizon = 1:3)
    ## summary(x)
    ldata <- copy(x$prepared_data)
    ldata[,id := NULL]
    ldata[,rtmle_predicted_outcome := NULL]
    ldata[,start_followup_date := NULL]
    detQ <- function(data, current.node, nodes, called.from.estimate.g){
        death.index <- grep("Dead_",names(data))
        if(length(death.index)==0){
            ## message("No death/terminal event node found")
            return(NULL)
        }
        hist.death.index <- death.index[death.index < current.node]
        if(length(hist.death.index)==0)
            return(NULL)
        else{
            is.deterministic <- Reduce("+",lapply(data[,hist.death.index,drop=FALSE],
                                                  function(dd){x=dd;x[is.na(dd)] <- 0;x}))>=1
            # should be unnecessary to exclude those who readily
            # have a missing value for death, but it does not hurt either
            is.deterministic[is.na(is.deterministic)] <- FALSE
            list(is.deterministic=is.deterministic, Q.value=0)
        }
    }
    suppressWarnings(suppressMessages(y1 <- ltmle(data = ldata[,1:6,with = FALSE],Anodes = c("A_0"),Cnodes = c("Censored_1"),Ynodes = c("Y_1"),Lnodes = c("L_0"),gform = c(A_0 = "A_0 ~ sex + age + L_0",Censored_1 = "Censored_1 ~ sex + age + L_0 + A_0"),Qform = c(Y_1 = "Q.kplus1 ~ sex + age + L_0 + A_0"),survivalOutcome = TRUE,deterministic.Q.function = detQ,abar = 1,estimate.time = FALSE)))
    expect_equal(y1$estimates[["tmle"]],x$estimate$Outcome_risk$Always_A$Estimate[[1]])
    expect_equal(y1$IC[["tmle"]],x$IC$Outcome_risk$Always_A[[1]])
    suppressWarnings(suppressMessages(y2 <- ltmle(data = ldata[,1:11,with = FALSE],Anodes = c("A_0","A_1"),Cnodes = c("Censored_1","Censored_2"),Ynodes = c("Y_1","Y_2"),Lnodes = c("L_0","Dead_1","L_1"),gform = c(A_0 = "A_0 ~ sex + age + L_0",Censored_1 = "Censored_1 ~ sex + age + L_0 + A_0",A_1 = "A_1 ~ sex + age + L_0 + A_0",Censored_2 = "Censored_2 ~ sex + age + L_0 + A_0 + L_1 + A_1"),Qform = c(Y_1 = "Q.kplus1 ~ sex + age + L_0 + A_0",Y_2 = "Q.kplus1 ~ sex + age + L_0 + A_0 + L_1 + A_1"),survivalOutcome = TRUE,deterministic.Q.function = detQ,abar = rep(1,2),estimate.time = FALSE)))
    # fitted g-models
    for (g in c("A_0","Censored_1","A_1","Censored_2")){
        a <- y2$fit$g[[g]]
        b <- as.matrix(x$models[["Always_A"]][[g]]$fit)
        b <- b[match(rownames(b),rownames(a)),]
        expect_equal(a,b)
    }
    expect_equal(y2$estimates[["tmle"]],x$estimate$Outcome_risk$Always_A$Estimate[[2]])
    expect_equal(y2$IC[["tmle"]],x$IC$Outcome_risk$Always_A[[2]])
    suppressWarnings(suppressMessages(y3 <- ltmle(data = ldata,Anodes = c("A_0","A_1","A_2"),Cnodes = c("Censored_1","Censored_2","Censored_3"),Ynodes = c("Y_1","Y_2","Y_3"),Lnodes = c("L_0","Dead_1","L_1","Dead_2","L_2"),gform = c(A_0 = "A_0 ~ sex + age + L_0",Censored_1 = "Censored_1 ~ sex + age + L_0 + A_0",A_1 = "A_1 ~ sex + age + L_0 + A_0",Censored_2 = "Censored_2 ~ sex + age + L_0 + A_0 + L_1 + A_1",A_2 = "A_2 ~ sex + age + L_0 + A_0 + L_1 + A_1",Censored_3 = "Censored_3 ~ sex + age + L_0 + A_0 + L_1 + A_1 + L_2 + A_2"),Qform = c(Y_1 = "Q.kplus1 ~ sex + age + L_0 + A_0",Y_2 = "Q.kplus1 ~ sex + age + L_0 + A_0 + L_1 + A_1",Y_3 = "Q.kplus1 ~ sex + age + L_0 + A_0 + L_1 + A_1 + L_2 + A_2"),survivalOutcome = TRUE,deterministic.Q.function = detQ,abar = rep(1,3),estimate.time = FALSE)))
    expect_equal(y3$estimates[["tmle"]],x$estimate$Outcome_risk$Always_A$Estimate[[3]])
    expect_equal(y3$IC[["tmle"]],x$IC$Outcome_risk$Always_A[[3]])
})




######################################################################
### test-rtmle-versus-ltmle.R ends here
