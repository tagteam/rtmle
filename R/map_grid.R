map_grid <- function(grid,
                     data,
                     name,
                     rollforward,
                     values=c(1,0),
                     fun_aggregate = NULL,
                     fill=NA,
                     id){
    value = NULL
    stopifnot(rollforward[[1]]>0)
    if (length(data)==0) return(NULL)
    setkeyv(grid,c(id,"date"))
    setkeyv(data,c(id,"date"))
    data = copy(data)
    if (length(values) == 2){
        data[,value:=values[[1]]]
    } else{
        stopifnot("value"%in%names(data))
    }
    # when rollforward = Inf then the last value is carried forward
    # with equidistant time grids rollforward can be the length of the interval
    
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
    wide <- fast_cast(grid,
                      id = id,
                      value_col = "value",
                      fill = fill,
                      fun_aggregate = fun_aggregate)
    setkey(wide,id)
    ## print(all.equal(wide,wide1))
    grid[,value:=NULL]
    data[,value:=NULL]
    # dcast assigns numeric column names when value.var
    # has length one
    setnames(wide,c(id,paste0(name,"_",names(wide)[-1])))
    setkeyv(wide,id)
    wide[]
}
