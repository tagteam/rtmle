map_grid <- function(grid,
                     data,
                     name,
                     rollforward,
                     values=c(1,0),
                     fun.aggregate = NULL,
                     fill=NA,
                     X.factor=FALSE,
                     id){
    X = NULL
    if (length(data)==0) return(NULL)
    setkeyv(grid,c(id,"date"))
    setkeyv(data,c(id,"date"))
    data = copy(data)
    data[,X:=values[[1]]]
    grid <- data[grid,roll=rollforward]
    # missing value means no event in this interval
    grid[is.na(grid$X),X:=values[[2]]]
    setkeyv(grid,c(id,"interval"))
    # note: need do.call because otherwise fun.aggregate is not
    #       interpreted correctly
    wide <- do.call(data.table::dcast,list(grid,
                               stats::formula(paste0(id,"~interval")),
                               value.var="X",
                               sep="_",
                               fun.aggregate = fun.aggregate,
                               fill=fill))
    if (X.factor) {
        # this is for ltmle censored/uncensored
        for (cc in names(wide)[-1]){
            set(wide,j=cc,value=factor(wide[[cc]],levels=values))
        }
    }
    grid[,X:=NULL]
    data[,X:=NULL]
    # dcast assigns numeric column names when value.var
    # has length one
    setnames(wide,c(id,paste0(name,"_",names(wide)[-1])))
    setkeyv(wide,id)
    wide[]
}
