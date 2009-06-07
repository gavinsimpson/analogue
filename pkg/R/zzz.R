.First.lib <- function(lib, pkg)  {
    library.dynam("analogue", pkg, lib)
    packageStartupMessage("This is analogue ",
                          utils::packageDescription("analogue",
                                                    field="Version"),
                          appendLF = TRUE)
}
