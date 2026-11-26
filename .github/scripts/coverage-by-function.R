write_function_coverage_summary <- function(coverage, file) {
    line_coverage <- covr::tally_coverage(coverage, by = "line")
    named <- !is.na(line_coverage$functions) & nzchar(line_coverage$functions)
    line_coverage <- line_coverage[named, ]

    filenames <- unique(line_coverage$filename)
    sources <- setNames(lapply(filenames, readLines, warn = FALSE), filenames)
    groups <- split(line_coverage, paste(line_coverage$filename,
                                         line_coverage$functions, sep = "\r"))

    sections <- lapply(groups, function(function_lines) {
        filename <- function_lines$filename[1L]
        line_numbers <- seq.int(min(function_lines$line), max(function_lines$line))
        hits <- function_lines$value[match(line_numbers, function_lines$line)]
        missed <- !is.na(hits) & hits == 0

        label <- sprintf("%s (%s) — %d of %d lines missed",
                         function_lines$functions[1L], filename,
                         sum(missed), nrow(function_lines))
        code <- sprintf("%s %4d  %s", ifelse(missed, "-", " "), line_numbers,
                        sources[[filename]][line_numbers])

        c(
            "<details>", "",
            paste0("<summary>", htmltools::htmlEscape(label), "</summary>"), "",
            "```diff", code, "```", "",
            "</details>", ""
        )
    })

    cat(
        "### Coverage by function", "",
        "Expand a function to see its source. Lines marked with `-` contain code that neither test suite executed.", "",
        unlist(sections, use.names = FALSE),
        file = file, append = TRUE, sep = "\n"
    )
}
