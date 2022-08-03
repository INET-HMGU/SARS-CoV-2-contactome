funcEnrich <- function(ocg, organism, querylist) {
    go <- list()
    df <- data.frame("name" = names(ocg$clustsizes), "clustsize" = ocg$clustsizes)
    for (i in querylist) {
        query_list <- ocg[[4]][ocg[[4]][, 2] == i, "node"]
        goquery <- gost(query = query_list, organism = organism, correction_method = "bonferroni", evcodes = TRUE)
        goquery$result$inCommunity <- paste0( goquery$meta$query_metadata$queries$query_1, collapse = ",")
        goquery$result$annotatedInCommunity <- paste0(goquery$meta$genes_metadata$query$query_1$ensgs, collapse = ",")
        target <- query_list[query_list %in% binary$node]
        goquery$result$viralTarget <- paste0(target, collapse = ",")
        go[[i]] <- goquery$result
    }
    return(go)
}
