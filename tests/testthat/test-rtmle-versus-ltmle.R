### test-rtmle-versus-ltmle.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Nov 16 2024 (17:04) 
## Version: 
## Last-Updated: Jul 30 2025 (06:31) 
##           By: Thomas Alexander Gerds
##     Update #: 70
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
    suppressMessages(result1 <- ltmle(data, Anodes="A", Ynodes="Y", abar=1,estimate.time = FALSE,gbounds = c(0,1)))
    ## summary(result1)
    x <- rtmle_init(intervals = 1,name_id = "id",name_outcome = "Y",name_competing = NULL,name_censoring = NULL,censored_label = "censored")
    rdata <- cbind(id = 1:NROW(data),data)
    setDT(rdata)
    setnames(rdata,c("id","W","A_0","Y_1"))
    x$prepared_data <- rdata
    x$names$name_baseline_covariates <- "W"
    x$names$name_time_covariates <- "A"
    suppressMessages(x <- protocol(x,name = "A",treatment_variables = "A",intervention = 1))
    x <- target(x,name = "Outcome_risk",estimator = "tmle",protocols = "A")
    x <- model_formula(x)
    x <- run_rtmle(x,refit = TRUE,verbose = FALSE)
    expect_equal(result1$fit$g[["A"]], as.matrix(x$models[["A_0"]]$fit))
    ## all.equal(result1$fit$Q[["Y"]], as.matrix(x$models[["Y_1"]]$fit))
    expect_equal(result1$estimates[["tmle"]],x$estimate$Main_analysis$Estimate,tolerance = 0.0001)
    expect_equal(result1$IC$tmle,x$IC$Outcome_risk$A$time_horizon_1,tolerance = 0.0001)
})

test_that("longitudinal data compare rtmle with ltmle",{
    set.seed(17)
    ld <- simulate_long_data(n = 1291,number_visits = 20,beta = list(A_on_Y = -.2,A_0_on_Y = -0.3,A_0_on_A = 6),register_format = TRUE)
    x <- rtmle_init(intervals = 3,name_id = "id",name_outcome = "Y",name_competing = "Dead",name_censoring = "Censored",censored_label = "censored")
    x <- add_long_data(x,outcome_data=ld$outcome_data,censored_data=ld$censored_data,competing_data=ld$competing_data,timevar_data=ld$timevar_data)
    x <- add_baseline_data(x,data=ld$baseline_data)
    suppressMessages(x <- long_to_wide(x,intervals = seq(0,2000,30.45*12)))
    suppressMessages(x <- protocol(x,name = "Always_A",treatment_variables = "A",intervention = 1))
    x <- prepare_data(x)
    x <- target(x,name = "Outcome_risk",estimator = "tmle",protocols = "Always_A")
    x <- model_formula(x)
    suppressWarnings(x <- run_rtmle(x,learner = "learn_glm",time_horizon = 1:3,verbose = FALSE))
    ## summary(x)
    ldata <- copy(x$prepared_data)
    ldata[,id := NULL]
    ldata[,rtmle_predicted_outcome := NULL]
    ldata[,c("A_0","A_1","A_2") := lapply(.SD,function(a){1*(a == 1)}),.SDcols = c("A_0","A_1","A_2")]
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
    expect_equal(y1$estimates[["tmle"]],x$estimate$Main_analysis$Estimate[[1]])
    expect_equal(y1$IC[["tmle"]],x$IC$Outcome_risk$Always_A[[1]])
    suppressWarnings(suppressMessages(y2 <- ltmle(data = ldata[,1:11,with = FALSE],Anodes = c("A_0","A_1"),Cnodes = c("Censored_1","Censored_2"),Ynodes = c("Y_1","Y_2"),Lnodes = c("L_0","Dead_1","L_1"),gform = c(A_0 = "A_0 ~ sex + age + L_0",Censored_1 = "Censored_1 ~ sex + age + L_0 + A_0",A_1 = "A_1 ~ sex + age + L_0 + A_0",Censored_2 = "Censored_2 ~ sex + age + L_0 + A_0 + L_1 + A_1"),Qform = c(Y_1 = "Q.kplus1 ~ sex + age + L_0 + A_0",Y_2 = "Q.kplus1 ~ sex + age + L_0 + A_0 + L_1 + A_1"),survivalOutcome = TRUE,deterministic.Q.function = detQ,abar = rep(1,2),estimate.time = FALSE)))
    # fitted g-models
    for (g in c("A_0","Censored_1","A_1","Censored_2")){
        a <- y2$fit$g[[g]]
        b <- as.matrix(x$models[[g]]$fit)
        rownames(b) <- sub("01","0",rownames(b))
        rownames(b) <- sub("11","1",rownames(b))
        b <- b[match(rownames(a),rownames(b)),]
        expect_equal(a,b)
    }
    expect_equal(y2$estimates[["tmle"]],x$estimate$Main_analysis$Estimate[[2]])
    expect_equal(y2$IC[["tmle"]],x$IC$Outcome_risk$Always_A[[2]])
    suppressWarnings(suppressMessages(y3 <- ltmle(data = ldata,Anodes = c("A_0","A_1","A_2"),Cnodes = c("Censored_1","Censored_2","Censored_3"),Ynodes = c("Y_1","Y_2","Y_3"),Lnodes = c("L_0","Dead_1","L_1","Dead_2","L_2"),gform = c(A_0 = "A_0 ~ sex + age + L_0",Censored_1 = "Censored_1 ~ sex + age + L_0 + A_0",A_1 = "A_1 ~ sex + age + L_0 + A_0",Censored_2 = "Censored_2 ~ sex + age + L_0 + A_0 + L_1 + A_1",A_2 = "A_2 ~ sex + age + L_0 + A_0 + L_1 + A_1",Censored_3 = "Censored_3 ~ sex + age + L_0 + A_0 + L_1 + A_1 + L_2 + A_2"),Qform = c(Y_1 = "Q.kplus1 ~ sex + age + L_0 + A_0",Y_2 = "Q.kplus1 ~ sex + age + L_0 + A_0 + L_1 + A_1",Y_3 = "Q.kplus1 ~ sex + age + L_0 + A_0 + L_1 + A_1 + L_2 + A_2"),survivalOutcome = TRUE,deterministic.Q.function = detQ,abar = rep(1,3),estimate.time = FALSE)))
    expect_equal(y3$estimates[["tmle"]],x$estimate$Main_analysis$Estimate[[3]])
    expect_equal(y3$IC[["tmle"]],x$IC$Outcome_risk$Always_A[[3]])
})
if (FALSE){
    test_that("longitudinal data compare rtmle with ltmle two-dimensional treatment",{
        rexpit <- function (x) rbinom(n = length(x), size = 1, prob = plogis(x))
        set.seed(9)
        n <- 1000
        age <- rbinom(n, 1, 0.5)
        gender <- rbinom(n, 1, 0.5)
        L_0 <- 2*age - 3*gender + rnorm(n)
        A_0 <- rexpit(age + gender)
        B_0 <- rexpit(age + 0.5*gender)
        Y_1 <- rexpit(-2+plogis(-1 + age - gender + L_0 - 0.5*B_0 - A_0 + rnorm(n)))
        A_1 <- rexpit(age + gender)
        B_1 <- rexpit(age + 0.5*gender - A_1)
        Y_2 <- rexpit(-2+plogis(-1 + age - gender + L_1 - 0.5*B_1 - A_0 - 0.8*A_1 + rnorm(n)))
        L_1 <- 2*age - 3*gender + 2*A_0 + rnorm(n)
        A_2 <- rexpit(age + gender + A_1 - L_0 - B_1)
        B_2 <- rexpit(age + 0.5*gender - A_1 - A_2)
        L_2 <- 2*age - 3*gender + 2*A_1 + rnorm(n)
        Y_3 <- rexpit(-2+plogis(-1 + age - gender + L_2 - 0.5*B_2 - A_1 - 0.8*A_2 + rnorm(n)))
        # once outcome has occurred it stays
        Y_2[Y_1 == 1] <- 1
        Y_3[Y_2 == 1] <- 1
        data <- data.frame(age, gender, L_0, A_0, B_0, Y_1, A_1, B_1, L_1, Y_2, A_2, B_2, L_2, Y_3)
        result <- ltmle(data,Anodes=c("A_0","B_0","A_1","B_1","A_2","B_2"),Lnodes=c("L_0", "L_1", "L_2"),Ynodes=grep("^Y", names(data),value = TRUE),survivalOutcome = TRUE,estimate.time = FALSE,abar=list(c(1,0,1,0,1,0),c(0,1,0,1,0,1)))
        summary(result)
        x <- rtmle_init(intervals = 3,name_id = "id",name_outcome = "Y",name_competing = NULL,name_censoring = NULL)
        rdata <- cbind(id = 1:NROW(data),data)
        setDT(rdata)
        ## setnames(rdata,c("id","age","gender","A_0","B_0","L_1","A_1","B_1","Y_1","L_2","A_2","B_2","Y_2"))
        x <- add_baseline_data(x,rdata[,.(id,age,gender)])
        x <- add_wide_data(x,outcome_data = rdata[,.(id,Y_1,Y_2,Y_3)],timevar_data = list(A = rdata[,.(id,A_0,A_1,A_2)],B = rdata[,.(id,B_0,B_1,B_2)],L = rdata[,.(id,L_0,L_1,L_2)]))
        x <- prepare_data(x)
        x <- protocol(x,name = "AnotB",treatment_variables = c("A","B"),intervention = c(1,0))
        x <- protocol(x,name = "BnotA",treatment_variables = c("A","B"),intervention = c(0,1))
        x <- target(x,name = "Outcome_risk",estimator = "tmle",protocols = c("AnotB","BnotA"))
        x <- model_formula(x,exclusion_rules = make_exclusion_rule(3,"A","B",recursive = FALSE),inclusion_rules = make_exclusion_rule(3,"B","A",recursive = FALSE))
        ## x <- model_formula(x)
        x <- run_rtmle(x,refit = TRUE,times = 3)
        # fitted g-models
        for (g in c("A_0","B_0","A_1","B_1")){
            a <- result$fit$g[[1]][[g]]
            b <- as.matrix(x$models[[g]]$fit)
            rownames(b) <- sub("01","0",rownames(b))
            rownames(b) <- sub("11","1",rownames(b))
            b <- b[match(rownames(a),rownames(b)),]
            expect_equal(a,b)
        }
        x$estimate
    })
}

######################################################################
### test-rtmle-versus-ltmle.R ends here
