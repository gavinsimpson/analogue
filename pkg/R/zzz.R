.First.lib <- function(lib, pkg)  {
    cat(paste("This is analogue",
              utils::packageDescription("analogue")$Version,
              "\n"))
    invisible()
}
