# from particular node lists to their 1st interactors and generate subnetwork
combineNetwork <- function(network, node) {
    gwas_random_g <- make_ego_graph(network, nodes = node, order = 1, mode = "all")
    gwas_random_list_df <- lapply(gwas_random_g, as_data_frame)
    gwas_random_df <- do.call(rbind, gwas_random_list_df)
    gwas_random_g_merge <- graph_from_data_frame(gwas_random_df, directed = FALSE)
    gwas_random_final <- simplify(induced_subgraph(network, names(V(gwas_random_g_merge))), remove.loops = FALSE)
    return(gwas_random_final)
}
