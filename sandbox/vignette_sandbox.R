## devtools::install()
library(rtmle)
library(data.table)
devtools::load_all("~/Rpackages/rtmle/")

diab_data <- fread("./data/diabetes_sim_data.csv")
diab_data[event == "MACE" & MACE == 1][, .(id, date)]
diab_data[event == "death",.(id,date)]
diab_data[order(id,date), .SD[.N], by = id][event%in%c("visit","dropout"),.(id,date)]
diab_data[last == 1 &event %in% c("visit","dropout"),.(id,date)]
## diab_data[event%in%c("baseline","MACE","visit"),.(id, date = time, HbA1c)]
tv_covs <- c("HbA1c","GLP1","SGLT2","DPP4")
    tv <- lapply(tv_covs, function(tv){
	if (tv %in% c("GLP1","SGLT2","DPP4")){
	    d <- diab_data[,c("id","time",tv),with = FALSE]
	    setnames(d,c("id","start_exposure","value"))
	    d[,end_exposure := shift(x = start_exposure,n = 1,type = "lead"),by = id]
	    # when treatment starts at the last time point in a time-series
	    # then the duration is unknown
	    d <- d[!is.na(end_exposure)]
	    d <- d[value != 0]
	    d[,value := NULL]
	}else{
	    d <- diab_data[event%in%c("baseline","MACE","visit"),c("id","time",tv),with = FALSE]
	    setnames(d,c("id","date","value"))
	}
	d[]
    })   
names(tv) <- tv_covs
outcome_data <- diab_data[event == "MACE", .(date = min(time)), by = id]


## Get first MACE event time for each id
x <- rtmle_init(time_grid = seq(0, 5*6, by = 6),
		name_id = "id",
		name_competing = "death",
                name_censoring = "censoring",
		name_outcome = "MACE")
x <- add_baseline_data(x,data=diab_data[time == 0, .(id, sex, age)])
x <- add_long_data(x, outcome_data = outcome_data)
x <- add_long_data(x, HbA1c = tv[["HbA1c"]])
x <- add_long_data(x, SGLT2 = tv[["SGLT2"]])

x <- long_to_wide(x,
                  HbA1c = "last",
                  SGLT2 = list(exposed = "exposed",
                               exposed_3month = function(data, grid, name, id) exposure_time(data, grid, name, id, threshold = 3)))

x$data$timevar_data

## FIx this bug:
tv[["SGLT2"]][id == 5]
x$data$timevar_data$SGLT2[id == 5]

tv[["SGLT2"]][id == 1]
x$data$timevar_data$SGLT2[id == 1]


x1 <- add_long_data(x, timevar_data = tv[1])
x2 <- add_long_data(x, HbA1c = tv[["HbA1c"]])
identical(x1, x2)

x <- add_long_data(x, outcome_data = outcome_data)
x <- add_long_data(x, timevar_data = tv[1], outcome_data = outcome_data, censored_data = NULL, competing_data = NULL)
x <- add_long_data(x, timevar_data = tv[1], outcome_data = outcome_data, censored_data = NULL, competing_data = NULL)
y <- add_long_data(x, timevar_data = tv[1], outcome_data = outcome_data, censored_data = NULL, competing_data = NULL)




x <- add_long_data(x, outcome_data = outcome_data, timevar_data = NULL, censored_data = NULL, competing_data = NULL)
x <- add_long_data(x, timevar_data = tv[1], outcome_data = outcome_data, censored_data = NULL, competing_data = NULL)
x <- add_long_data(x, timevar_data = tv[2], outcome_data = outcome_data, censored_data = NULL, competing_data = NULL)
x <- add_long_data(x, timevar_data = tv[3], outcome_data = outcome_data, censored_data = NULL, competing_data = NULL)

x <- long_to_wide(x,breaks = time_grid,start_followup_date = 0)

## Fix something about date/time internal naming.

add_long_data(x, outcome_data = outcome_data, timevar_data = NULL, censored_data = NULL)

diab_data[event == "MACE", .(id, time)]


