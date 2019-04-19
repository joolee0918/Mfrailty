# ========================================================= 
# Extract offset name
# =========================================================

untangle.offset <- function(tt) {
    spc <- attr(tt, "offset")
    if (length(spc) == 0) 
        vars = character(0)
    facs <- attr(tt, "factors")
    fname <- dimnames(facs)
    vars = (fname[[1]])[spc]
    return(vars)
}
