.First.lib <- function(lib, pkg)  {
    packageStartupMessage("This is analogue ",
                          utils::packageDescription("analogue",
                                                    field="Version"),
                          appendLF = TRUE)
}
