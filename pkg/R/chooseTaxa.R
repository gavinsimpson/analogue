chooseTaxa <- function(object, ...) {
    UseMethod("chooseTaxa")
}

chooseTaxa.default <- function(object, n.occ = 1, max.abun = 0,
                               type = c("AND","OR"), ...) {
    if(missing(type))
        type <- "AND"
    type <- match.arg(type)
    occ.want <- colSums(object > 0) >= n.occ
    abun.want <- apply(object, 2, max) >= max.abun
    want <- if(isTRUE(all.equal(type, "AND"))) {
        occ.want & abun.want
    } else {
        occ.want & abun.want
    }
    return(object[,want])
}
