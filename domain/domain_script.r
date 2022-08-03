COL_OK_HIGH <- "#a63603"
COL_OK_LOW <- "#fdae6b"
COL_KO <- "#fee6ce"
COL_BG <- "lightgrey"

faux_to_domain_intx <- function(intx_data, domains_per_protein)
{
  tmp <- apply(intx_data,1,
               function(x)
               {
                 tmp <- domains_per_protein[[ x[["human_uniprot"]] ]]
                 if (length(tmp) > 0)
                 {
                   t(sapply(tmp, function(y) { c(x, y) }))
                 }
               })
  domain_intx <- as.data.frame(do.call("rbind", tmp[ !sapply(tmp, is.null)]), stringsAsFactors = F)
  colnames(domain_intx)[[6]] <- "human_pfam"
  domain_intx$intx_domain_ID <- sprintf("%s_%s", domain_intx$virus, domain_intx$human_pfam)
  return (domain_intx)
}

f_analyze_domain_intx <- function(intx_data, domains_per_protein)
{
  domains_per_human_prot <- domains_per_protein[ unique(intx_data$human_uniprot) ]
  domain_intx <- faux_to_domain_intx(intx_data, domains_per_protein)
  # get intx that do not have a domain
  no_domain_intx <- intx_data$intx_ID[! intx_data$intx_ID %in% domain_intx$intx_ID ]
  # table w/ intx per domain-intx
  protein_intx_per_domain_intx <- sapply(split(domain_intx$intx_ID, domain_intx$intx_domain_ID), unique)
  protein_intx_per_domain_intx <- protein_intx_per_domain_intx[ order(sapply(protein_intx_per_domain_intx, length), decreasing = T)]
  # count intx that have an overlapping domain-virus intx; those that do not overlap; and those w/o domains
  w_overlapping_domain_intx <- sapply(split(domain_intx, domain_intx$intx_ID),
                                      function(x)
                                      {
                                        this_intx_ID <- x[["intx_ID"]][[1]]
                                        this_domain_intx <- x[["intx_domain_ID"]]
                                        tmp <- sapply(this_domain_intx,
                                                      function(y)
                                                      {
                                                        sum((y == domain_intx$intx_domain_ID) & (domain_intx$intx_ID != this_intx_ID))
                                                      })
                                        (sum(tmp) > 0)
                                      })
  counts <- c("overlapping" = sum(w_overlapping_domain_intx),
              "unique" = sum(!w_overlapping_domain_intx),
              "no domain" = length(no_domain_intx))
  to_return <- list("domain_intx" = domain_intx,
                    "counts" = counts,
                    "protein_intx_per_domain_intx" = protein_intx_per_domain_intx,
                    "domains_per_human_prot" = domains_per_human_prot)
  return (to_return)
}

faux_piechart_redundancy_summary <- function(counts)
{
  labels <- sprintf("%s\n(%.0f%%, n=%s)", names(counts), 100*counts/sum(counts), counts)
  pie(counts, col=c(COL_OK_HIGH, COL_KO, COL_BG), clockwise = T, border="white", labels=labels)
  title(main=sprintf("Recurrence of virus-human PPIs\n(one datapoint per virus-human PPI)\nn_PPI=%s", sum(counts)), cex.main=1, line=-1)
}

faux_plot_real_vs_random <- function(real_counts, rnd_counts, pvalue, XLAB)
{
  hist_rnd <- hist(rnd_counts, breaks = (0:(max(rnd_counts)+1)), plot = F, right = F) #T, include.lowest = F)
  density_rnd <- density(rnd_counts)
  MAXY <- max(density_rnd$y)*1.2
  MAXX <- max(c(density_rnd$x, real_counts))*1.2
  plot(hist_rnd, freq=F, las=1, col="grey90", border="white", ylim=c(0, MAXY), xlim=c(0,MAXX), xlab=XLAB, xaxt='n', main="")
  xaxis0 <- par("xaxp")[[1]]
  xaxis1 <- par("xaxp")[[2]]
  xaxis_bins <- par("xaxp")[[3]]
  xaxis_ticks <- seq(from = xaxis0, to = xaxis1, length.out = xaxis_bins+1)
  axis(1, tick=T, at=xaxis_ticks+0.5, labels=xaxis_ticks, xpd=T)
  lines(density_rnd$x+0.5, density_rnd$y, lwd=2, col=COL_BG) #"grey")
  arrows(x0 = real_counts+0.5, x1 = real_counts+0.5, y1 = 0, y0=MAXY/10, col=COL_OK_HIGH, lwd=2, length=0.1)
  text(x = real_counts+0.5, y=MAXY/10, real_counts, pos = 3, col=COL_OK_HIGH, font = 2)
  
  abline(v=median(rnd_counts)+0.5, lty=2)
  text(x = median(rnd_counts)+0.5, y=MAXY, pos = 4, sprintf("median(rnd)=%d", median(rnd_counts)), col=COL_BG, xpd=T, font=2)
  title(main=sprintf("HuSCI vs randomized networks\nn_rnd=%s; p=%.3f", length(rnd_counts), pvalue), cex.main=1)
}

f_plot_result <- function(to_return, OUTPUT_PDF)
{
  pdf(OUTPUT_PDF, width=10, height=3.5)
  par("mfrow"=c(1,3))
  par_orig <- par("mar")
  par("mar"=c(3,4,3,3))
  faux_piechart_redundancy_summary(to_return$real_res$counts)
  par("mar" = par_orig)
  faux_plot_real_vs_random(to_return$real_res$counts[[1]], to_return$rnd_res$rnd_counts$overlapping, to_return$rnd_res$pvalue, "Number of PPIs w/ redundant virus-domain interactions")

  ## Analysis per virus protein: get type of PPIs (partner w/ redundant domain, unique domain, or no domain) per virus protein
  data <- to_return$real_res$domain_intx
  per_virus_protein <- split(data, data$virus)
  overlapping_res <- list()
  for (this_virus_protein in names(per_virus_protein))
  {
    this_domain_intx <- per_virus_protein[[this_virus_protein]]
    per_PPI <- split(this_domain_intx, this_domain_intx$intx_ID)
    overlapping_PPI <- sapply(per_PPI, 
                              function(x)
                              {
                                tmp <- sapply(per_PPI, 
                                              function(y)
                                              { # check other PPIs of this virus: do they have the same domain?
                                                sum(x[["human_pfam"]] %in% y[["human_pfam"]]) > 0
                                              } )
                                sum(tmp) > 1 # we always have at least (self identify)
                              })
    overlapping_res[[this_virus_protein]] <- overlapping_PPI
  }
  N_overlapping_per_virus_prot <- sapply(overlapping_res, sum)
  N_intx_per_virus_prot <- sapply(overlapping_res, length)
  # count also those intx w/ human proteins that do not have a domain
  N_intx_per_virus_total <- sapply(split(to_return$intx_data$human, to_return$intx_data$virus), length)
  all_virus_proteins <- names(N_intx_per_virus_total)
  to_show <- rbind("overlapping"=N_overlapping_per_virus_prot[all_virus_proteins],
                   "unique"=(N_intx_per_virus_prot-N_overlapping_per_virus_prot)[all_virus_proteins],
                   "all" = N_intx_per_virus_total[all_virus_proteins])
  colnames(to_show) <- all_virus_proteins
  to_show[ is.na(to_show)] <- 0
  to_show <- rbind(to_show, "no_domain"=to_show["all",] - to_show["unique",] - to_show["overlapping",])
  to_show <- to_show[c("overlapping","unique","no_domain"),]
  MAXY <- max(apply(to_show, 2, sum))*1.2
  b <- barplot(to_show, beside=F, border=NA, col=c(COL_OK_HIGH, COL_KO, COL_BG), las=2, ylim=c(0,MAXY), ylab="Number of PPIs")
  title(main="Recurrence of PPIs per virus protein", cex.main=1)
  grid(nx=NA, ny=NULL)
  legend("topright", fill=rev(c(COL_OK_HIGH, COL_KO, COL_BG)), border=NA, bty='n', rev(c("overlapping", "unique", "no domain")))

  dev.off()
}

faux_generate_randomized_networks <- function(intx_data, N_SWITCHES_FACTOR=10, NUM_RANDOM=1000, R_object="rnd_networks.R_object")
{
  if (!is.null(R_object) && file.exists(R_object))
  {
    load(R_object)
  }else
  {
    library("BiRewire")

    # save uniprot mappings
    uniprot_mappings <- intx_data$human_uniprot
    names(uniprot_mappings) <- intx_data$human
    uniprot_mappings <- uniprot_mappings[ !duplicated(names(uniprot_mappings))]
    # generate adjacency matrix w/ PPI info
    adj_matrix <- tapply(rep(T, nrow(intx_data)), intx_data[,c("virus", "human")], as.logical, default = F)
    n_edge_switches <- sum(adj_matrix)*N_SWITCHES_FACTOR

    all_rnd_pairs <- list()
    for (i_r in 1:NUM_RANDOM)
    {
      rnd_network <- birewire.rewire.bipartite(adj_matrix, max.iter = n_edge_switches, exact = T)
      rnd_pairs <- as.data.frame(as.table(rnd_network), stringsAsFactors=F)
      rnd_pairs <- rnd_pairs[ rnd_pairs$Freq,]
      rnd_pairs$Freq <- NULL
      colnames(rnd_pairs) <- c("virus", "human")
      rnd_pairs$group <- "rnd"
      rnd_pairs$human_uniprot <- uniprot_mappings[rnd_pairs$human]
      rnd_pairs$intx_ID <- sprintf("%s_%s", rnd_pairs$virus, rnd_pairs$human)
      all_rnd_pairs[[i_r]] <- rnd_pairs
    }
    if (!is.null(R_object))
    {
      save(all_rnd_pairs, file=R_object)
    }
  }
  return (all_rnd_pairs)
}

f_analyze_domain_intx_vs_random <- function(intx_data, domain_intx, domains_per_protein, dataset_ID, NUM_RANDOM=1000)
{
  # generate randomized networks
  rnd_intx <- faux_generate_randomized_networks(intx_data, NUM_RANDOM=NUM_RANDOM)
  # process all randomized networks
  rnd_res <- lapply(rnd_intx, function(x) { f_analyze_domain_intx(x, domains_per_protein) })
  # compare results of real vs random networks
  real_res <- domain_intx$counts[[1]]
  rnd_counts <- as.data.frame(t(sapply(rnd_res, function(x) { x[["counts"]]})), stringsAsFactors = F)
  pvalue <- mean(real_res <= rnd_counts$overlapping)
  to_return <- list("rnd_intx" = rnd_intx,
                    "rnd_counts" = rnd_counts,
                    "rnd_res" = rnd_res,
                    "pvalue" = pvalue)
  return (to_return)
}

f_process_dataset <- function(intx_data, domains_per_protein, dataset_ID, OUTPUT_PDF=NULL)
{
  R_object <- sprintf("tmp_files/%s_results.R_object", dataset_ID)
  if (file.exists(R_object))
  {
    load(R_object)
  }else
  {
    domain_intx <- f_analyze_domain_intx(intx_data, domains_per_protein)
    rnd_res <- f_analyze_domain_intx_vs_random(intx_data, domain_intx, domains_per_protein, dataset_ID)
    to_return <- list("intx_data" = intx_data,
                      "real_res" = domain_intx,
                      "rnd_res" = rnd_res)
    save(to_return, file=R_object)
  }
  if (!is.null(OUTPUT_PDF))
  {
    f_plot_result(to_return, OUTPUT_PDF)
  }
  return (to_return)
}

faux_evaluate_significance <- function(domain_intx, all_prots, all_pfams)
{
  res_fisher <- list()
  for (this_human_pfam in all_pfams)
  {
    for (this_virus_prot in all_prots)
    {
      P1 <- sum((domain_intx$virus == this_virus_prot) & (domain_intx$human_pfam == this_human_pfam)) 
      N1 <- sum((domain_intx$virus == this_virus_prot) & (domain_intx$human_pfam != this_human_pfam)) 
      P2 <- sum((domain_intx$virus != this_virus_prot) & (domain_intx$human_pfam == this_human_pfam)) 
      N2 <- sum((domain_intx$virus != this_virus_prot) & (domain_intx$human_pfam != this_human_pfam)) 
      this_fisher <- fisher.test(matrix(c(P1, N1, P2, N2), nrow=2))
      res_fisher[[sprintf("%s_%s", this_virus_prot, this_human_pfam)]] <- c(this_virus_prot, this_human_pfam, P1, N1, P2, N2, this_fisher$p.value)
    }
  }
  df_summary <- as.data.frame(do.call("rbind", res_fisher), stringsAsFactors=F)
  colnames(df_summary) <- c("virus", "human_pfam", "P1", "N1", "P2", "N2", "pvalue")
  for (this_f in c("P1", "N1", "P2", "N2", "pvalue"))
  {
    df_summary[[this_f]] <- as.numeric(df_summary[[this_f]])
  }
  return (df_summary)
}

## read pfam data: get domains per protein and ID_to_name
pfam_data <- read.table("input/9606.tsv", comment.char = "", skip=3, header=F, sep = "\t", colClasses = "character")
pfam_data <- pfam_data[ pfam_data$V8 == "Domain",]
domains_per_protein <- sapply(split(pfam_data$V6, pfam_data$V1), unique)
domain_info <- pfam_data$V7
names(domain_info) <- pfam_data$V6
domain_info <- domain_info[ !duplicated(names(domain_info)) ]

## read virus-human interactions 
all_intx <- read.csv("input/binary_edge_1126.csv", header=T, colClasses = "character", sep = "\t")
all_intx$intx_ID <- sprintf("%s_%s", all_intx$virus, all_intx$human)

## process dataset: get recurrent PPIs in real vs randomized networks
HuSCI_res <- f_process_dataset(all_intx, domains_per_protein, "HuSCI", "HuSCI_domain_analysis.pdf")

## evaluate significance of protein-domain associations
all_pfams <- unique(HuSCI_res$real_res$domain_intx$human_pfam)
all_prots <- unique(HuSCI_res$intx_data$virus)
df_summary <- faux_evaluate_significance(HuSCI_res$real_res$domain_intx, all_prots, all_pfams)
df_summary$PPIs <- ifelse(row.names(df_summary) %in% names(HuSCI_res$real_res$protein_intx_per_domain_intx),
                          sapply(HuSCI_res$real_res$protein_intx_per_domain_intx[ row.names(df_summary) ], paste, collapse=","),
                          "")
save(df_summary, file="tmp_files/df_summary.R_object")
df_summary$domain_name <- domain_info[ df_summary$human_pfam ]
write.table(df_summary, file = "husci_virus_domain_associations.csv", quote=FALSE, sep = "\t", col.names = T, row.names = F)
selected <- df_summary[ (df_summary$P1 > 1) & (df_summary$pvalue < 0.05),]
dim(selected)
write.table(selected, file = "husci_virus_domain_associations.filtered.csv", quote=FALSE, sep = "\t", col.names = T, row.names = F)

## evaluate significant associations in random networks
rnd_df_summary_file <- "tmp_files/rnd_df_summary.R_object"
if (file.exists(rnd_df_summary_file))
{
  load(rnd_df_summary_file)
}else
{
  rnd_df_summary <- list()
  NUM_RANDOM <- length(HuSCI_res$rnd_res$rnd_res)
  for (i_r in 1:NUM_RANDOM)
  {
    rnd_df_summary[[i_r]] <- faux_evaluate_significance(HuSCI_res$rnd_res$rnd_res[[i_r]]$domain_intx, all_prots, all_pfams)
    if (i_r %% 10 == 0) { cat(sprintf("Evaluating associations in randomized networks: %s/%s (%s)\n", i_r, NUM_RANDOM, date()))}
  }
  save(rnd_df_summary, file=rnd_df_summary_file)
}

## require at least 2 interactions (P1 >= 2; to avoid spurious associations) and p-value < 0.05
pdf("number_of_significant_virus_domain_associations.pdf", width=4, height=4)
PVALUE_CUTOFF <- 0.05
MIN_INTX <- 2
n_sign_real <- sum((df_summary$pvalue < PVALUE_CUTOFF) & df_summary$P1 >= MIN_INTX)
n_sign_rnd <- sapply(rnd_df_summary, function(x) { sum((x[["pvalue"]] < PVALUE_CUTOFF) & (x[["P1"]] >= MIN_INTX)) })
pvalue <- sum(n_sign_real <= n_sign_rnd)/length(n_sign_rnd)
faux_plot_real_vs_random(n_sign_real, n_sign_rnd, pvalue, "Number of significant virus-domain associations")
dev.off()

