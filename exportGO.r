export <- function(go, checklist) {
    dfa1 <- list()
    for (i in checklist) {
        if(length(go[[i]]) == 19) {
            dfa1[i] <- go[i]
            dfa1[[i]]$cluster <- i
        } else {
            next
        }
    }
    return(dfa1)
}
