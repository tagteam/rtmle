map_grid <- function(grid,
                     data,
                     name,
                     rollforward,
                     values=c(1,0),
                     fun_aggregate = NULL,
                     fill=NA,
                     value_is_factor=FALSE,
                     id){
    value = NULL
    if (length(data)==0) return(NULL)
    setkeyv(grid,c(id,"date"))
    setkeyv(data,c(id,"date"))
    data = copy(data)
    if (length(values) == 2){
        data[,value:=values[[1]]]
    } else{
        stopifnot("value"%in%names(data))
    }
    grid <- data[grid,roll=rollforward]
    if (length(values) == 2){
        # missing value means no event in this interval
        grid[is.na(grid$value),value:=values[[2]]]
    }else{
        # missing value means no measurement can be rolled forward
        # but nothing is done about it here
    }
    setkeyv(grid,c(id,"interval"))
    # note: need do.call because otherwise fun_aggregate is not
    #       interpreted correctly
    wide <- do.call(data.table::dcast,
                    list(grid,
                         stats::formula(paste0(id,"~interval")),
                         value.var="value",
                         sep="_",
                         fun.aggregate = fun_aggregate,
                         fill=fill))
    if (value_is_factor) {
        # this is for ltmle censored/uncensored
        for (cc in names(wide)[-1]){
            set(wide,j=cc,value=factor(wide[[cc]],levels=values))
        }
    }
    grid[,value:=NULL]
    data[,value:=NULL]
    # dcast assigns numeric column names when value.var
    # has length one
    setnames(wide,c(id,paste0(name,"_",names(wide)[-1])))
    setkeyv(wide,id)
    wide[]
}
