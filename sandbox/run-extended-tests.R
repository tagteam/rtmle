run_extended_file <- function(path, file) {
    pkgload::load_all(path, quiet = TRUE, export_all = FALSE)
    setwd(tempdir())

    if (basename(dirname(file)) == "slowtests") {
        # test_file() evaluates in the file's directory.
        test_file <- basename(file)
        if (!file.copy(file, test_file)) stop("Could not copy ", file)
        Sys.setenv(NOT_CRAN = "true")
        testthat::local_edition(3)
        testthat::test_file(test_file, reporter = "summary", stop_on_failure = TRUE)
        return(invisible(NULL))
    }

    extension <- tolower(tools::file_ext(file))
    if (extension == "org") {
        lines <- readLines(file, warn = FALSE)
        starts <- grep("^[[:space:]]*#\\+BEGIN_SRC[[:space:]]+R([[:space:]]|$)",
                       lines, ignore.case = TRUE)
        ends <- grep("^[[:space:]]*#\\+END_SRC[[:space:]]*$", lines, ignore.case = TRUE)
        blocks <- lapply(starts, function(start) {
            end <- ends[ends > start][1L]
            if (is.na(end)) stop("Unclosed R block at line ", start)
            lines[start + seq_len(end - start - 1L)]
        })
        file <- "vignette.R"
        writeLines(as.character(unlist(blocks)), file)
    } else if (extension == "rmd") {
        file <- knitr::purl(file, output = "vignette.R", quiet = TRUE, documentation = 0)
    }
    source(file, local = globalenv(), print.eval = TRUE)
    invisible(NULL)
}

run_extended_tests <- function(path = ".") {
    path <- normalizePath(path, mustWork = TRUE)
    slowtests <- list.files(file.path(path, "slowtests"), "\\.[rR]$", full.names = TRUE)
    vignettes <- list.files(file.path(path, "draft-vignettes"), "\\.(org|Rmd|R)$",
                            full.names = TRUE, ignore.case = TRUE)
    if (!length(slowtests) || !length(vignettes)) {
        stop("Run from the rtmle source checkout or supply its path.")
    }

    failed <- character()
    for (file in c(slowtests, vignettes)) {
        label <- file.path(basename(dirname(file)), basename(file))
        message("\nRunning ", label)
        tryCatch(
            callr::r(run_extended_file, args = list(path, file), show = TRUE,
                     user_profile = FALSE),
            error = function(error) {
                message(conditionMessage(error))
                failed <<- c(failed, label)
            }
        )
    }
    if (length(failed)) {
        stop("Extended checks failed in: ", paste(failed, collapse = ", "), call. = FALSE)
    }
    message("\nAll extended checks passed.")
    invisible(NULL)
}

if (sys.nframe() == 0L) run_extended_tests()
