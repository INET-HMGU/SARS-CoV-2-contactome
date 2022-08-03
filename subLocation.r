# subcellular location proportion calculation

subLocaPtoM <- function(x) {
    dt_raw <- str_split(x$Subcellular.main.location, pattern = ", ")
    dt <- data.frame(subcell = tolower(unlist(dt_raw)))
    dt$subcell[is.na(dt$subcell)] <- "not detected"
    dt2 <- merge(sub_major, dt, by.x = "subcellular", by.y = "subcell")
    dt2_t <- prop.table(table(dt2$major))
    return(as.data.frame(dt2_t))
}

subLocaCtoM <- function(x) {
    dt_raw <- str_split(x$Subcellular.main.location, pattern = ", ")
    dt <- data.frame(subcell = tolower(unlist(dt_raw)))
    dt$subcell[is.na(dt$subcell)] <- "not detected"
    dt2 <- merge(sub_major, dt, by.x = "subcellular", by.y = "subcell")
    dt2_t <- table(dt2$major)
    return(as.data.frame(dt2_t))
}
