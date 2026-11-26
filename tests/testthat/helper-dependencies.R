# Give a copy of a package function private dependency bindings for one test.
with_test_dependencies <- function(fun, ...) {
    environment(fun) <- list2env(list(...), parent = environment(fun))
    fun
}

without_package <- function(fun, unavailable) {
    with_test_dependencies(fun, requireNamespace = function(package, ...) {
        if (identical(package, unavailable)) return(FALSE)
        base::requireNamespace(package, ...)
    })
}
