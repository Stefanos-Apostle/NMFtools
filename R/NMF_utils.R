#' @title NMF Manhattan
#'
#' @description Function to create NMF Factor manhattan plot
#' @param NMF_model NMF object
#' @param nmf_factor Numerical position or column name of individual factor
#' @param subset_w option to select a subset of the w matrix
#' @return A ggplot formatted manhattan plot
#' @examples
#' plot <- NMF_manhattan(NMF_model, nmf_factor = 21)
#' @export
NMF_manhattan <- function(NMF_model, nmf_factor,  subset_w = NA) {

  # subset w functionality not built in yet

  nmf <- data.frame(SNP = names(NMF_model$w[, nmf_factor]),
                    contr = NMF_model$w[, nmf_factor],
                    chrom = NMF_model$w_metadata$chrom,
                    pos = NMF_model$w_metadata$pos)
  nmf$chrom <- factor(nmf$chrom, levels = c(as.character(1:22), "X"))

  plot <- ggplot(data = nmf, aes(x = as.numeric(pos), y = contr)) +
    geom_point(size = .01, pch = ".") +
    facet_wrap(~chrom, switch = "x", nrow = 1, scales = "free_x") +
    xlab("Genomic Position") +
    #ylim(c(.0000006, max(nmf$contr))) +
    ylim(c(.0000006, max(NMF_model$w))) +
    theme_classic() +
    ggtitle(label = paste("nmf", nmf_factor)) +
    theme(text = element_text(face = "bold"), axis.line.y = element_line(linewidth = 1), axis.ticks.y = element_line(size = 1, lineend = "square"),
          plot.title = element_text(hjust = 0.5),axis.title.y = element_blank(), axis.line.x = element_blank(),
          axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    coord_cartesian()
  #scale_y_continuous(limits = c(0, max(nmf$contr)), expand = c(0,0))
  #scale_y_continuous(limits = c(0, max(SNPs_df$w)), breaks = seq(0, max(SNPs_df$w), by = .00005), expand = c(0,0))

  return(plot)
}

#' @title Multiple Ridges
#'
#' @description Function to create skyline tracks of NMF factors
#' @param NMF_model NMF object
#' @param nmf_factors Numerical positions or column names of individual factor
#' @param window Number of base pairs that the function will smooth the manhattan skyline over
#' @param chromosome Which chromosome to plot skyline over
#' @param custom_range Array of length 2 with start and stop positions along chromosome.
#' @return A dataframe of window with corresponding max contribution in the window. To be used for plotting with plot_multiple_ridges()
#' @examples
#' # Calculating max contributions in factors 21, 54, and 108 at the NNAT locus (GRCh37) in 0.1 Mb windows
#' multiple_ridges_df <- multiple_ridges(NMF_model, nmf_factor = c(21, 54, 108), window = 100000, chromosome = 20, custome_range = c(36149617,36152092))
#' @export
multiple_ridges <- function(NMF_model, nmf_factors, window = 10000, chromosome, custom_range = NA) {

  ridge_df <- data.frame()
  for (f in nmf_factors) {

    nmf <- data.frame(SNP = names(NMF_model$w[, f]),
                      contr = NMF_model$w[, f],
                      chrom = NMF_model$w_metadata$chrom,
                      pos = NMF_model$w_metadata$pos)
    nmf$chrom <- factor(nmf$chrom, levels = c(as.character(1:22), "X"))

    nmf_c <- nmf[which(nmf$chrom == chromosome), ]
    nmf_reg <- nmf_c[which(nmf_c$pos >= custom_range[1] & nmf_c$pos <= custom_range[2]),]


    # window = 100000
    if (is.na(custom_range[1])) {
      chunks = seq(from = min(nmf_reg$pos), to = max(nmf_reg$pos), by = window)
    }else if (length(custom_range) == 2 & is.numeric(custom_range)){
      chunks = seq(from = custom_range[1], to = custom_range[2], by = window)
    }else{
      stop("custom_range must be a vector of length 2 and be of class 'numeric'")
    }


    max_df <- data.frame()
    for (i in chunks) {
      chunk_df <- nmf_reg[which(nmf_reg$pos >= i & nmf_reg$pos < (i +window)), ]

      if(nrow(chunk_df) != 0) {
        tdf <- data.frame(refsnp_id = i,
                          chr_name = chromosome,
                          chrom_start = i +(window/2),
                          chrom_end = i +(window/2),
                          w = max(chunk_df$contr))
        max_df <- rbind(max_df, tdf)
      }
    }
    max_df$nmf_factor <- as.character(f)
    ridge_df <- rbind(ridge_df, max_df)
  }
  return(ridge_df)
}

#' @title Plot Multiple Ridges
#'
#' @description Function to create skyline tracks of NMF factors
#' @param mult_ridge_df Dataframe of chromosome window with associated max factor contributions. Output of multiple_ridges()
#' @param squish Geom_ridgeline parameter to move stacked plots closer together by squishing the y axis
#' @param lwd Plot line width
#' @param alpha Transparency of plot color fill
#' @param add_chr Option to add ideogram at the bottom of the plot to show where this region is.
#' @return A plot showing the skyline between factors across a specific region
#' @examples
#' # Calculating max contributions in factors 21, 54, and 108 at the NNAT locus (GRCh37) in 0.1 Mb windows
#' multiple_ridges_df <- multiple_ridges(NMF_model, nmf_factor = c(21, 54, 108), window = 100000, chromosome = 20, custome_range = c(36149617,36152092))
#' plot_mutliple_ridges(multiple_ridges_df)
#' @export
plot_multiple_ridges <- function(mult_ridge_df, squish = 10000, lwd = 1, alpha = 1, add_chr = F) {
  num_fs <- length(unique(mult_ridge_df$nmf_factor))
  mult_ridge_df$peak_id <- unlist(lapply(X = 1:num_fs, FUN = function(X){rep(X, nrow(mult_ridge_df)/num_fs)}))

  plot <- ggplot(data = mult_ridge_df, aes(x = chrom_start, y = peak_id/squish, height = w, group = peak_id, fill = factor(peak_id))) +
    geom_ridgeline(size = lwd, alpha = alpha) +
    scale_fill_discrete(labels = unique(mult_ridge_df$nmf_factor)) +
    labs(fill = "NMF Factor") +
    theme_classic()+
    ylab("Contribution") +
    xlab("Chr Position") +
    theme(axis.text.y = element_blank(), plot.title = element_text(hjust = 0.5), text = element_text(face = "bold")) +
    ggtitle(paste("Chr", unique(mult_ridge_df$chr_name), ": ", min(mult_ridge_df$chrom_start), "-", max(mult_ridge_df$chrom_start), sep = ""))

  if (add_chr == T) {
    p.ideo <- Ideogram(genome = "hg19",subchr =  "chr8", zoom.region = c(min(mult_ridge_df$chrom_start), max(mult_ridge_df$chrom_start)), alpha = 0.2)
    plot <- ggarrange(plot, p.ideo@ggplot +rremove("y.text"), ncol = 1, heights = c(5,1))
  }

  return(plot)
}

#' @title Test Perplexity
#'
#' @description Function to plot dimenstional reduction at different perplexities
#' @param NMF_model NMF object
#' @param perplexity Perplexity value passed to Rtsne()
#' @param return_df Option to return plots or data frame of dimensional reduction
#' @return A plot of dimensional reduction at the specified perplexity or a dataframe if return_df=T
#' @examples
#' plot_p40 <- test_perplexity(NMF_model, perplexity = 40)
#' @export
test_perplexity <- function(NMF_model, perplexity, return_df = F) {
  ht <- NMF$h
  colnames(ht) <- c(1:ncol(ht))
  h_tsne_test <- Rtsne(t(ht), perplexity = perplexity, check_duplicates = F)
  layout <- h_tsne_test$Y
  colnames(layout) <- c("tSNE1", "tSNE2")

  if (return_df) {
    return(as.data.frame(layout))
  }else{
    plot <- ggplot(as.data.frame(layout), aes(tSNE1, tSNE2)) +
      geom_point(size = 0.2) +
      theme_classic() +
      theme(aspect.ratio = 1)

    return(plot)
  }
}

#' @title sNN HClust Clusters
#'
#' @description Function to calc hclust of sNN clusters, as a way to compare shared nearest neighbors cluters.
#' @param sNN_clusters
#' @param trait_dists
#' @param clusters_to_calc
#' @param method
#' @return Hclust object ordering the sNN clusters
#' @examples
#' sNN_hc_clusters(...)
#' @export
sNN_hc_clusters <- function(sNN_clusters, trait_dists, clusters_to_calc = "All", method = "euclidean"){
  names(sNN_clusters) <- c(1:length(sNN_clusters))
  cluster_dists <- data.frame()
  for (i in c(0:max(sNN_clusters))) {

    ts <- which(sNN_clusters == i)
    tsn <- names(sNN_clusters)[ts]

    td <- trait_dists[as.numeric(tsn), ]
    tdc <- colMeans(td)

    cluster_dists <- rbind(cluster_dists, tdc)

    #nrow(cluster_dists) == length(unique(sNN_clusters))
  }

  if (length(clusters_to_calc) == 1) {
    cds <- dist(cluster_dists, method = method)
  }else{
    cds <- dist(cluster_dists[as.numeric(clusters_to_calc), ], method = method)
  }
  chc <- hclust(cds)
  return(chc)

}

#' @title Count Factors
#'
#' @description Function to calulcate the number of factors with non zero weight for each phenotype
#' @param dimred_df Dataframe of dimensional reduction of phenotype matrix
#' @param h_matrix Phenotypic matrix (NMF_model$h)
#' @return Dimred_df with new column of count values
#' @examples
#' count_factors(dimred_df, h_matrix = NMF_model$h)
#' @export
count_factors <- function(dimred_df, h_matrix) {
  pls <- c()
  for (i in c(1:ncol(h_matrix))) {
    pl <- length(which(h_matrix[, i] != 0))
    pls <- c(pls, pl)
  }
  names(pls) <- colnames(h_matrix)
  dimred_df$num_factors <- pls
  return(dimred_df)
}

#' @title Plot Trait Keyword
#'
#' @description Function to highlight traits related to an array of keywords on a dimensional reduction map.
#' @param dimred_df Dataframe of dimensional reduction of phenotype matrix
#' @param keyword Array of character strings to use for trait searching
#' @param h_matrix Phenotypic matrix (NMF_model$h)
#' @param dimred Which type of dimensional reduction to plot
#' @param size Values to pass to ggplot to change size of poiints
#' @param clusters Cluster assignments for each trait (can use sNN clustering output)
#' @return Ggplot of dimensional reduction highlighting the traits related to supplied keywords
#' @examples
#' plot_kw(dimred_df, keyword = c("BMI", "Type 2 diabetes", "diet"), h_matrix = NMF_model$h)
#' @export
plot_kw <- function(dimred_df, keyword, h_matrix,dimred = "tSNE", size = "num_factors", clusters = F) {

  dimred_df <- count_factors(dimred_df, h_matrix)

  select <- rep(NA, nrow(dimred_df))
  for (k in keyword) {
    traits <- grep(k, dimred_df$description, ignore.case = T)
    #print(k)
    #rint(traits)
    select[traits] <- k
  }
  dimred_df$select <- select

  if (clusters != F) {
    # string of cluster attributes for each trait
    dimred_df$clusters <- clusters
  }

  if (dimred == "tSNE") {
    plot <- ggplot(dimred_df, aes(x = TSNE1, y = TSNE2, size = num_factors, color = select, alpha = ifelse(is.na(select), 0.2, 1)))
  }else if (dimred == "UMAP") {
    plot <- ggplot(dimred_df, aes(x = UMAP1, y = UMAP2, size = num_factors, color = select, alpha = ifelse(is.na(select), 0.2, 1)))
  }else{
    stop("dimred must be one of the following; c('tSNE', 'UMAP'")
  }

  if (size == "num_factors") {
    plot <- plot +
      geom_point()
  }else if (is.numeric(size)) {
    plot <- plot +
      geom_point(size = size)
  }else{
    stop("Size must be default or a single number")
  }

  plot <- plot +
    #geom_point() +
    #ggtitle(keyword) +
    theme_classic() +
    labs(color = "Keyword") +
    scale_radius(range = c(0.5, 5)) +
    scale_alpha(guide = "none") +
    theme(aspect.ratio = 1)
  return(plot)
}


#' @title Plot Trait Keyword
#'
#' @description Function to update dimred_df highlighting traits related to keyword()
#' @param dimred_df Dataframe of dimensional reduction of phenotype matrix
#' @param keyword Array of character strings to use for trait searching
#' @param clusters Cluster assignments for each trait (can use sNN clustering output)
#' @param h_matrix Phenotypic matrix (NMF_model$h)
#' @param exact If false traits contianing keyword are selected, else if exact = T only traits that match keyword are selected.
#' @return Updated dimred_df with new column to identify if traits realted to keyword.
#' @examples
#' dimred_kw(dimred_df, keyword = c("BMI", "Type 2 diabetes", "diet"), h_matrix = NMF_model$h)
#' @export
dimred_kw <- function(dimred_df, keyword, clusters, h_matrix, exact = F) {
  "function to update dimred_df for a set of keywords"

  dimred_df <- count_factors(dimred_df, h_matrix)
  dimred_df$clusters <- clusters

  select <- rep(NA, nrow(dimred_df))
  for (k in keyword) {
    if (exact == F) {
      traits <- grep(k, dimred_df$description, ignore.case = T)
      select[traits] <- dimred_df$clusters[traits]
    }else if (exact == T) {
      traits <- which(tolower(k) == tolower(dimred_df$description))
      select[traits] <- dimred_df$clusters[traits]
    }
  }
  dimred_df$select <- select

  dimred_df <- dimred_df[-which(is.na(dimred_df$select)), ]
  return(dimred_df)
}

#' @title Hclust on the sNN object to create a dendrogram of sNN clusters.
#'
#' @description Function to update dimred_df highlighting traits related to keyword()
#' @param dimred_df Dataframe of dimensional reduction of phenotype matrix
#' @param keyword Array of character strings to use for trait searching
#' @param clusters Cluster assignments for each trait (can use sNN clustering output)
#' @param sNN shared nearest neighbors object.
#' @param return_plotlist If false traits one combined plot, else if T individual plots are returned as a list.
#' @return Updated dimred_df with new column to identify if traits realted to keyword.
#' @examples
#' cluster_dotplot(dimred_df, keyword = c("BMI", "Type 2 diabetes", "diet"), clusters = clusters, sNN = sNN_object)
#' @export
cluster_dotplot <- function(dimred_df, keyword, clusters, sNN, return_plotlist = F) {

  dimred_df <- dimred_kw(dimred_df, keyword, clusters, h_matrix)
  clusters_plotted <- unique(clusters[as.numeric(rownames(dimred_df))])

  ## dendrogram

  hc <- sNN_hc_clusters(sNN_clusters = clusters, trait_dists = sNN$dist, clusters_to_calc = clusters_plotted)
  dend_plot <- ggdendrogram(hc, rotate = F, labels = F) + scale_y_reverse()


  ## dotplot
  plot <- ggplot(dimred_df, aes(y = factor(description), x = factor(select, levels = as.numeric(hc$labels[hc$order])), size = num_factors)) +
    geom_point() +
    theme_classic() +
    xlab("Cluster") +
    ylab("") +
    labs(size = "No. Factors") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5), axis.text.y = element_text(size = 5))

  if (return_plotlist == F) {
    comb_plot <- egg::ggarrange(plot, dend_plot, ncol = 1)
    return(comb_plot)
  }else if (return_plotlist == T) {
    plot_list <- list(plot, dend_plot)
    return(plot_list)
  }
}

#' @title Highlight clusters.
#'
#' @description Function to plot the clusters from the keywords to accompany the dotplot/hclut figure.
#' @param dimred_df Dataframe of dimensional reduction of phenotype matrix
#' @param keyword Array of character strings to use for trait searching
#' @param clusters Cluster assignments for each trait (can use sNN clustering output)
#' @param h_matrix Phenotypic matrix (NMF_model$h)
#' @param exact If false traits contianing keyword are selected, else if exact = T only traits that match keyword are selected.
#' @return Updated dimred_df with new column to identify if traits realted to keyword.
#' @examples
#' highlight_clusters(dimred_df, keyword = c("BMI", "Type 2 diabetes", "diet"), h_matrix = NMF$h, clusters = clusters, sNN = sNN_object)
#' @export
highlight_clusters <- function(dimred_df, keyword, h_matrix, clusters, exact = F) {
  "Function to plot the clusters from the keywords to accompany the dotplot/hclut figure"
  dimred_df$clusters <- clusters
  plot <- plot_kw(dimred_df, keyword, h_matrix)

  dimred_df2 <- dimred_kw(dimred_df, keyword, clusters, h_matrix, exact = exact)
  clusters_plotted <- unique(clusters[as.numeric(rownames(dimred_df2))])

  dimred_df$clusters[-which(dimred_df$clusters %in% dimred_df2$clusters)] <- NA

  plot <- plot +
    geom_point(aes(color = factor(dimred_df$clusters), alpha = ifelse(is.na(dimred_df$clusters), 0.5, 1)), size = 0.5) +
    labs(color = "Cluster")
  plot$layers <- plot$layers[2]
  return(plot)
}

#' @title sNN Clusters Panel
#'
#' @description Function to plot a panel of figures including highlighted traits, dendrogram of clusters, and dimensional reduction ggpoint of there.
#' @param dimred_df Dataframe of dimensional reduction of phenotype matrix
#' @param keyword Array of character strings to use for trait searching
#' @param clusters Cluster assignments for each trait (can use sNN clustering output)
#' @param h_matrix Phenotypic matrix (NMF_model$h)
#' @param sNN shared nearest neighbors object.
#' @return Updated dimred_df with new column to identify if traits realted to keyword.
#' @examples
#' sNN_clusters_plots(dimred_df, keyword = c("BMI", "Type 2 diabetes", "diet"), h_matrix = NMF$h, clusters = clusters, sNN = sNN_object)
#' @export
sNN_clusters_plots <- function(dimred_df, keyword, clusters, sNN, h_matrix) {

  dotplot <- cluster_dotplot(dimred_df, keyword, clusters, sNN, return_plotlist = T)
  dotplot[[1]]
  dotplot[[2]]
  tsne_plot <- highlight_clusters(dimred_df, keyword, h_matrix, clusters)

  plot <- egg::ggarrange(dotplot[[1]], dotplot[[2]], tsne_plot, ncol = 2, byrow = F)
  return(plot)
}

#' @title Plot keyword UMAP
#'
#' @description Function to plot traits related to keyword on UMAP of the h matrix
#' @param h_umap UMAP output of h_matrix
#' @param keyword Array of character strings to use for trait searching
#' @return Ggpoint plot of UMAP and highlighted traits
#' @examples
#' plot_kw_umap(h_umap, keyword = c("BMI", "Type 2 diabetes", "diet"))
#' @export
plot_kw_umap <- function(h_umap, keyword) {

  select <- rep(NA, nrow(h_umap))
  for (k in keyword) {
    traits <- grep(k, h_umap$description, ignore.case = T)
    select[traits] <- k
  }
  h_umap$select <- select

  plot <- ggplot(h_tsne, aes(x = 1, y = 2, size = num_factors, color = select, alpha = ifelse(is.na(select), 0.2, 1))) +
    geom_point() +
    #ggtitle(keyword) +
    theme_classic() +
    labs(color = "Keyword") +
    scale_radius(range = c(0.5, 5)) +
    scale_alpha(guide = "none") +
    theme(aspect.ratio = 1)
  return(plot)
}

#' @title Phenotype Matrix Heatmap
#'
#' @description Function to plot traits related to keyword on UMAP of the h matrix
#' @param NMF NMF object
#' @param keyword Array of character strings to use for trait searching
#' @param nmf_hc Hclust object on the phenotypic matrix
#' @param out_dir Path to figure output directory
#' @param exact If false traits contianing keyword are selected, else if exact = T only traits that match keyword are selected.
#' @param custom_filename If false default name by keywords is used, if character string provided this will be used for file prefix
#' @param column Which column of the phenoitypic metadata to use for searching trait names
#' @return Pheatmap of traits related to keyword columns ordered by hclust
#' @examples
#' pheno_hm(NMF, keyword = c("BMI", "Type 2 diabetes", "diet"), nmf_hc, out_dir = "./path_to_output/")
#' @export
pheno_hm <- function(NMF, keyword, nmf_hc, out_dir, exact = F, custom_filename = F, column = "description") {

  if (length(term) > 1) {
    name = paste(term, collapse = "_")
  }else{
    name = term
  }

  for (t in term) {
    if (exact == F) {
      ct <- grep(t, NMF$h_metadata[,column], ignore.case = T)
    }else if (exact == T) {
      ct <- which(tolower(NMF$h_metadata[,column]) == tolower(t))
    }

    tdf <- data.frame(col_num = ct,
                      Description = NMF$h_metadata[ct,column],
                      Search_Term = t)
    if (t == term[1]) {
      term_df <- tdf
    }else{
      term_df <- rbind(term_df, tdf)
    }
  }
  #rownames(term_df) <- paste("trait", term_df$col_num, sep = "")

  h_col <- c()
  for (kw_row in term_df$col_num) {
    #kw_row <- which(NMF$h_metadata$description == descr)
    x <- NMF$h_metadata[kw_row,1:5]
    tn <- paste(x[-which(x == "")], collapse = "-")
    h_col <- c(h_col, grep(tn, colnames(NMF$h)))
  }
  term_df$h_col <- h_col
  rownames(term_df) <- paste("trait", term_df$h_col, sep = "")


  #trait_descr_df <- data.frame(Description = manifest$description[term_df$col_num])
  #rownames(trait_descr_df) <- paste("trait", term_df$col_num, sep = "")
  h_matrix <- NMF$h
  colnames(h_matrix) <- paste("trait", c(1:ncol(h_matrix)), sep = "")

  if (length(term_df$h_col) == 0) {
    stop(paste(term, "not found in manifest files"))
  }else if (length(term_df$h_col) == 1) {
    h_df <- data.frame(trait = h_matrix[, term_df$h_col],
                       factor = as.character(c(1:150)))
  }else{
    h_df <- as.data.frame(h_matrix[, term_df$h_col])
    h_df$factor <- as.character(c(1:150))
  }
  rownames(h_df) <- h_df$factor

  #melt_data <- melt(h_df, id.vars = c("factor"))
  if (custom_filename == F) {
    filename <- paste(out_dir, name, "_heatmap.pdf", sep = "")
  }else{
    filename = paste(out_dir, custom_filename, "_heatmap.pdf", sep = "")
  }

  pdf(filename, width = 35, height = 15)
  h_df <- h_df[, -which(colnames(h_df) == "factor")]
  rg <- max(abs(h_df))
  if (exact == F) {
    b_plot <- pheatmap(t(h_df), cellwidth = 10, cellheight = 10, border_color = NA, cluster_cols = nmf_hc,
                       cluster_rows = F,
                       treeheight_col = 200, fontsize_col = 8,
                       color = colorRampPalette(c("blue", "white", "red"))(100),
                       breaks = seq(-.005, .005, length.out = 100),
                       main = ifelse(custom_filename == F, name, custom_filename),
                       annotation_row = term_df[, -which(colnames(term_df) %in% c("col_num", "h_col"))])#, labels_row = term_df$Description)
  }else if (exact == T) {
    b_plot <- pheatmap(t(h_df), cellwidth = 10, cellheight = 10, border_color = NA, cluster_cols = nmf_hc,
                       cluster_rows = F,
                       treeheight_col = 200, fontsize_col = 8, labels_col = c(1:150),
                       color = colorRampPalette(c("blue", "white", "red"))(100),
                       #annotation_row = term_df[, -which(colnames(term_df) == "col_num")],
                       breaks = seq((-1*rg), rg, length.out = 100),
                       main = ifelse(custom_filename == F, name, custom_filename))#, labels_row = term_df$Description)
  }

  print(b_plot)
  dev.off()

}

#' @title NMF Factor Udler Spider Plot
#'
#' @description Function to plot spider plots to compare with Udler bNMF paper.
#' @param factor NMF factor
#' @param h_df Data.frame of the NMF h_matrix
#' @param maxs Max values for each of these traits
#' @param mins Min values for each of these traits
#' @param custom_title Option to set a specific title name for the plot
#' @return spider plot for a specific nmf factor
#' @examples
#' factor_spider_plot(factor = 21, h_df = NMF$h, maxs, mins)
#' @export
factor_spider_plot <- function(factor, h_df, maxs, mins, custom_title = F) {
  final_df <- data.frame()
  for (k in factor){
    factor_df <- h_df[k,]

    HOMA_B <- factor_df$trait29/factor_df$trait15
    TG <- factor_df$trait28
    HDL <- factor_df$trait17
    WC <- factor_df$trait4372
    BMI <- factor_df$trait4109
    FastedIns <- factor_df$trait29
    HOMA_IR <- TG/HDL
    ProIns <- factor_df$trait4475

    values <- c(HOMA_B, ProIns, HOMA_IR, FastedIns, BMI, WC, HDL, TG)
    values[which(values > unique(maxs))] <- unique(maxs)
    values[which(values < unique(mins))] <- unique(mins)

    spider_df <- data.frame(row.names = c("HOMA_B", "ProIns", "HOMA_IR", "FastedIns", "BMI", "WC", "HDL", "TG"),
                            maxs, mins,
                            values)
    colnames(spider_df) <- c("max", "min", k)
    spider_df[,3][which(is.na(spider_df[,3]))] <- 0

    if (ncol(final_df) == 0) {
      final_df <- spider_df
    }else{
      final_df <- cbind(final_df, spider_df[, 3])
    }
  }
  colnames(final_df) <- c("min", "max", factor)


  if (custom_title == F) {
    plot_title = paste("Factor ", paste(factor, collapse = ", "), sep = "")
  }else{
    plot_title = custom_title
  }
  plot <- radarchart(as.data.frame(t(final_df)), title = plot_title, plty = 1, seg = 2, pcol = 1)
  return(plot)
}

#' @title Gene of interest manhattan plot
#'
#' @description Function to plot manhattans of genes of interest for a specific factor, faceted by nearest gene
#' @param top_bmi_genes Genes to plot (Need to change this object name)
#' @param NMF NMF object
#' @param factor NMF Factor to plot
#' @return GGplot of contribution manhattans by genes for a factor
#' @examples
#' goi_manhattan(top_bmi_genes = c("ADIPOQ", "FTO", "LEP"), NMF = NMF_object)
#' @export
goi_manhattan <- function(top_bmi_genes, NMF, factor) {
  top_bmi_w <- as.matrix(NMF$w[which(NMF$w_metadata$nearest_genes %in% top_bmi_genes), ])
  top_bmi_meta <- data.frame(snp = rownames(top_bmi_w),
                             nearest_gene = NMF$w_metadata$nearest_genes[which(NMF$w_metadata$nearest_genes %in% top_bmi_genes)],
                             pos = NMF$w_metadata$pos[which(NMF$w_metadata$nearest_genes %in% top_bmi_genes)])

  top_bmi_meta$contr = top_bmi_w[match(rownames(top_bmi_w), top_bmi_meta$snp), factor]

  plot <- ggplot(top_bmi_meta, aes(x = pos, y = contr)) +
    geom_point(size = 0.5) +
    facet_wrap(~nearest_gene, scales = "free_x", nrow = 1) +
    theme_classic() +
    #ylim(c(0, max(top_bmi_w))) +
    theme(text = element_text(face = "bold"), strip.text = element_text(size = 5, face = "bold", angle = 90), axis.line.y = element_line(linewidth = 1), axis.ticks.y = element_line(size = 1, lineend = "square"),
          plot.title = element_text(hjust = 0.5),axis.title.y = element_blank(), axis.line.x = element_blank(),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none",
          strip.background = element_rect( fill="white", color="white"
          )) +
    coord_cartesian()

  return(plot)
}

#' @title Plot Trait Weights
#'
#' @description Function to plot a trait's weights across all factors
#' @param NMF NMF object
#' @param trait Specific name of a trait
#' @return GGplot of a trait's weight across all factors
#' @examples
#' trait_factor_weights(NMF = NMF_object, trait = "BMI")
#' @export
trait_factor_weights <- function(NMF, trait) {


  ## need to finish this function
  ct <- which(tolower(NMF$h_metadata$description) == tolower(trait))
  Description = NMF$h_metadata$description[ct]

  x <- NMF$h_metadata[ct,1:5]
  tn <- paste(x[-which(x == "")], collapse = "-")
  h_col <- grep(tn, colnames(NMF$h))
  trait_num <- paste("trait", h_col, sep = "")

  trait_col <-  which(colnames(NMF$h) == paste("phecode", trait_num, "both_sexes.rds_pos", "both_sexes", sep = "-"))
  weights <- NMF$h[, trait_col]


  rank_df <- data.frame(factor = names(weights),
                        weight = weights,
                        order = rank(weights, ties.method = "random"))
  plot <- ggplot(rank_df, aes(x = order,y = weight, label = ifelse(abs(weight) > 0.005, as.character(rank_df$factor), ""))) +
    geom_text_repel(box.padding = 0.9, point.padding = 1, max.overlaps = Inf, min.segment.length = 0) +
    geom_point() +
    ylab("Factor Weight") +
    xlab("Order") +
    ggtitle(paste("'", trait, "' Factor Weights", sep = "")) +
    theme_classic() +
    theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5), text = element_text(face = "bold"))
  return(plot)

}

#' @title Get top Factor Genes
#'
#' @description Function to find the genes associated with the top SNPs in a factor
#' @param NMF_model NMF object
#' @param factor Specific nmf factor
#' @param cutoff Contribution cutoff to find genes assocatiated to SNPs above this cutoff
#' @return list of unique gene names
#' @examples
#' factor21_genes <- get_nmf_genes(NMF_model = NMF_object, factor = 21)
#' @export
get_nmf_genes <- function(NMF_model, factor, cutoff = .00005) {

  nmf <- NMF_model$w[, factor]
  sig_nmf <- names(nmf[which(nmf > cutoff)])

  #ensembl <- useEnsembl("snp",dataset = "hsapiens_snp",GRCh = "37")
  #bm_df <- getBM(attributes=c("refsnp_id",
  #                 "associated_gene"),
  #  filters ="snp_filter", values =sig_nmf, mart = ensembl, uniqueRows=TRUE)


  nmf_genes <- unique(NMF$w_metadata$nearest_genes[which(NMF$w_metadata$rsid %in% sig_nmf)])

  utts <- c()
  for (u in nmf_genes) {
    if (length(grep(" and ", u)) > 0) {
      k <- unlist(strsplit(u, " and "))
    } else if (length(grep(",", u)) > 0) {
      k <- unlist(strsplit(u, ","))
    } else{
      k <- u
    }
    utts <- c(utts, k)
  }

  r <- which(utts == "" | is.na(utts) | utts == "Intergenic")
  if (length(r) > 0) {utts <- utts[-r]}

  return(utts)
}

#' @title Fix Duplicates
#'
#' @description Function to fix the problem of merging creating replicates, need to find a proper way to merge that doesnt cause this
#' @param go_df Gene ontology results dataframe
#' @return Gene ontology data frame with duplicates removed
#' @examples
#' go_df <- fix_dups(go_df)
#' @export
fix_dups <- function(go_df) {

  u <- unique(go_df$term_name)
  for (k in u) {
    rs <- which(go_df$term_name == k)
    if (length(rs) > 1) {
      sub <- go_df[rs[1], ]
      go_df <- go_df[-rs, ]
      go_df <- rbind(go_df, sub)
    }
  }
  return(go_df)
}

#' @title Gene Ontology Enrichment
#'
#' @description Function to calc GO enrichment for each factor of the NMF model
#' @param NMF_model NMF_object
#' @param cutoff Contribution cutoff for variants to be used in GO enrichment
#' @return Gene ontology data frame
#' @examples
#' go_df <- calc_go_w(NMF_model = NMF_object)
#' @export
calc_go_w <- function(NMF_model, cutoff = 0.0001) {
  pb <- txtProgressBar(1, ncol(NMF_model$w), style = 3)
  #ensembl <- useEnsembl("snp",dataset = "hsapiens_snp",GRCh = "37")
  go_df <- data.frame()
  for (i in c(1:ncol(NMF$w))) {
    setTxtProgressBar(pb, i)
    nmf_genes <- get_nmf_genes(NMF_model = NMF, factor = i, cutoff = cutoff)

    if (length(nmf_genes) > 0) {
      nmf_go <- gost(nmf_genes, organism = "hsapiens", significant = F)

      if (is.null(nmf_go) == F) {
        tdf <- nmf_go$result[, c("term_name", "p_value")]
        colnames(tdf)[2] <- i

        if (nrow(go_df) == 0) {
          go_df <- tdf
        } else{
          go_df <- merge(go_df, tdf, by = "term_name", all = T)
          go_df <- fix_dups(go_df)
        }
      }
    }
  }
  return(go_df)
}

#' @title Format GO Dataframe
#'
#' @description Function to format the output go_df
#' @param go_df Gene ontology dataframe, output of calc_go_w()
#' @param cutoff Significance cutoff for GO enrichments
#' @return Gene ontology data frame with only significant enrichments
#' @examples
#' go_df <- format_go_df(go_df)
#' @export
format_go_df <- function(go_df, cutoff = 0.05){

  # only keeping the terms that are significant in atleast one of the nmf factors
  keep <- c()
  for (r in c(1:nrow(go_df))) {
    tl <- which(go_df[r,] > -log(cutoff))
    if (length(tl) > 0) {
      keep <- c(keep, r)
    }
  }
  go_df <- go_df[keep, ]
  rownames(go_df) <- go_df$term_name
  go_df$term_name <- NULL
  go_df[is.na(go_df)] <- 1

  # calc_go_w skips factors that do not have any assocaited genes based on cutoff to reduce time in connecting to biomart servers, so this function formats the go_df to contain these skipped factors
  missing <- which(!(c(1:150) %in% colnames(go_df)))
  for (m in missing) {
    go_df <- cbind(go_df, data.frame(m = 1))
    colnames(go_df)[which(colnames(go_df) == "m")] <- m
  }

  # log transformation and final prep of column order for pheatmap
  go_df <- -log(go_df)
  go_df <- go_df[, order(as.numeric(colnames(go_df)))]

  return(go_df)
}

#' @title Plot GO Enrichment
#'
#' @description Function to plot a heatmap of GO enrichments
#' @param NMF NMF object
#' @param filename Name of file to save heatmap to
#' @param hclust_obj Option to provide an hclust object to cluster NMF factors
#' @param cutoff Variant contribution cutoff for GO enrichment calc
#' @return Heatmap of GO enrichments across NMF factors
#' @examples
#' plot_go_heatmap(NMF = NMF_object, filename = "/path_to_outoput/GO_heatmap.pdf")
#' @export
plot_go_heatmap <- function(NMF, filename, hclust_obj = F, cutoff = 0.0001) {

  go_df <- calc_go_w(NMF, cutoff = cutoff)
  #go_df_save <- go_df
  go_df <- format_go_df(go_df)

  pdf(filename, width = 25, height = 15)
  pheatmap(go_df, border_color = NA, cluster_cols = nmf_hc, cluster_rows = T,
           treeheight_col = 200, fontsize_col = 8, color = colorRampPalette(c("white", "red"))(100),
           breaks = seq(-log(0.05), -log(0.001), length.out = 100)
  )
  dev.off()

  return(go_df)
}

#' @title Plot Only Significant Enrichments
#'
#' @description Function to Plot unique significant GO enrichments for specific NMF factors
#' @param go_df Formatted GO enrichment dataaframe
#' @param factors Array of NMF factor integer names to plot
#' @param hclust_obj Option to provide an hclust object to cluster NMF factors
#' @param cutoff Variant contribution cutoff for GO enrichment calc
#' @return Heatmap of significant GO enrichments across specified NMF factors
#' @examples
#' only_sig(go_df, factors = c(21, 87, 109))
#' @export
only_sig <- function(go_df, factors) {


  go_df <- go_df[, factors]
  keep <- c()
  for (r in c(1:nrow(go_df))) {
    tl <- which(go_df[r,] > -log(0.05))
    if (length(tl) > 0) {
      keep <- c(keep, r)
    }
  }
  go_df <- go_df[keep, ]


  pheatmap(go_df, border_color = NA, cluster_rows = T, cluster_cols = F, show_rownames = F,
           cellwidth = 7, cellheight = .1, treeheight_row = 0,
           treeheight_col = 200, fontsize_col = 8, color = colorRampPalette(c("white", "red"))(100),
           breaks = seq(-log(0.05), -log(0.0000001), length.out = 100), clustering_distance_rows = "euclidean",
  )
}

#' @title Get Unique GO Enrichments
#'
#' @description Function to return a list of GO terms that are specifically significant in that factor only
#' @param go_df Formatted GO enrichment dataaframe
#' @param factor Name of factor that we want specific enrichments for
#' @param pval_cutoff Significance cutoff for GO terms
#' @return List of GO terms that are specifically significant in that factor only
#' @examples
#' get_unique_GOterms(go_df, factor = 21)
#' @export
get_unique_GOterms <- function(go_df, factor, pval_cutoff = 0.05) {
  if (ncol(go_df) > 2) {
    oth_max <- apply(X=go_df[, which(colnames(go_df) != factor)], MARGIN=1, FUN=max)
  }else{
    oth_max <- go_df[, which(colnames(go_df) != factor)]
  }
  uniq_sig <- which(go_df[,which(colnames(go_df) == factor)] > -log(pval_cutoff) & oth_max < -log(pval_cutoff))
  uns <- rownames(go_df)[uniq_sig]
  return(uns)
}

#' @title Factor Rank List
#'
#' @description Function to created a ranked list for gsea
#' @param NMF_model NMF object
#' @param factor Name of factor that we want to create a ranked list for
#' @return Names list of gene and associated contribution
#' @examples
#' RNK_list <- factor_ranks(NMF_model = NMF_object, factor = 21)
#' @export
factor_ranks <- function(NMF_model, factor) {
  factor_contr <- NMF$w[, factor]
  nearest_genes <- NMF$w_metadata$nearest_genes

  rank_list <- c()
  for (nrg in unique(nearest_genes)) {
    ror <- which(nearest_genes == nrg)
    mg <- mean(factor_contr[ror])
    names(mg) <- nrg
    rank_list <- c(rank_list, mg)
  }
  rank_list <- rank_list[order(abs(rank_list), decreasing = T)]
  return(rank_list)
}

#' @title Calc GSEA
#'
#' @description Function to calculate GSEA for each factor (this is running GO for some reason, where did I use this one at?)
#' @param NMF_model NMF object
#' @param cutoff Variant contribution cutoff to be used for database curreation
#' @return Names list of gene and associated contribution
#' @examples
#' RNK_list <- factor_ranks(NMF_model = NMF_object, factor = 21)
#' @export
calc_fgsea_w <- function(NMF_model, cutoff = 0.0001) {
  pb <- txtProgressBar(1, ncol(NMF_model$w), style = 3)
  #ensembl <- useEnsembl("snp",dataset = "hsapiens_snp",GRCh = "37")
  go_df <- data.frame()
  for (i in c(1:ncol(NMF$w))) {
    setTxtProgressBar(pb, i)
    nmf_genes <- get_nmf_genes(NMF_model = NMF, factor = i, cutoff = cutoff)

    if (length(nmf_genes) > 0) {
      nmf_gsea <- fgsea()
      nmf_gsea <- gost(nmf_genes, organism = "hsapiens", significant = F)

      if (is.null(nmf_go) == F) {
        tdf <- nmf_go$result[, c("term_name", "p_value")]
        colnames(tdf)[2] <- i

        if (nrow(go_df) == 0) {
          go_df <- tdf
        } else{
          go_df <- merge(go_df, tdf, by = "term_name", all = T)
          go_df <- fix_dups(go_df)
        }
      }
    }
  }
  return(go_df)
}

#' @title Create UCSC Tracks
#'
#' @description Function to format a table in UCSC compatible format
#' @param file filename
#' @param NMF_model NMF_object
#' @param factor Specific factor name to format table for
#' @param name Track name
#' @param description Track description
#' @param window Base pair length of window size to create contribution skyline
#' @return Wriets table to file
#' @examples
#' create_ucsc_track(file = "Factor21.bedGraph", NMF_model = NMF_object, factor = 21, name = "Factor 21", description = "Contribution skyline for factor 21 in a 5000bp window.")
#' @export
create_ucsc_track <- function(file, NMF_model, factor, name, description, window = 5000) {

  ucsc_df <- data.frame()
  for (chr in c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X')) {
    map <- multiple_ridges(nmf_factors = factor, NMF_model, window, chromosome = chr)
    tdf <- data.frame(chr = paste("chr", map$chr_name, sep = ""),
                      start = map$chrom_start-(window/2),
                      end = map$chrom_start + (window/2),
                      value = map$w)
    ucsc_df <- rbind(ucsc_df, tdf)
  }

  header_df <- data.frame(head = c(paste("browser position chr", unique(map$chr_name), ":", min(ucsc_df$start), "-", max(ucsc_df$end), sep = ""),
                                   "browser hide all",
                                   "browser pack refGene encodeRegions",
                                   "browser full altGraph",
                                   paste('track type=bedGraph name="', name, '" description="', description, '" visibility=full color=200,100,0 altColor=0,100,200 priority=20', sep='')))


  write.table(x = header_df, file = file, row.names = F, quote = F, col.names = F)
  write.table(x = ucsc_df, file = file, row.names = F, quote = F, col.names = F, append = T)

}

#' @title Write UCSC Tracks
#'
#' @description Function to write UCSC tracks of NMF factors
#' @param out_dir Path to output directory
#' @param NMF_model NMF_object
#' @param window Base pair length of window size to create contribution skyline
#' @return Writes tracks for each of the NMF factors
#' @examples
#' write_tracks(out_dir = "/path_to_output/", NMF_model = NMF_object)
#' @export
write_tracks <- function(out_dir, NMF_model, window = 5000) {

  for (f in c(1:ncol(NMF_model$w))) {
    print(f)
    out_file <- paste(out_dir, "/NMF_factor_", f, "_skyline.bedGraph", sep = "")
    name = paste("NMF_factor_", f, sep = "")
    description = paste("Contribution skyline for factor ",f, " with smoothing window ", window, sep = "")
    create_ucsc_track(out_file, NMF_model,factor = f, name, description, window)
  }
}

#' @title Tissue Enrichment by GSEA
#'
#' @description Function to calculate tissue enrichment by gsea
#' @param NMF_model NMF_object
#' @param closest path to bedtools closest output
#' @param narrowpeaks path to scATAC seqs narrowpeaks file
#' @param num_top_peaks Number of top peaks to use for tissue specificity from tissue peak calling
#' @return Dataframe of tissue enrichments for each factor
#' @examples
#' calc_gsea_per_tissue(NMF = NMF_object, closes = "path_to_closest_file", narrowpeaks = "path_to_narrowpeaks_file")
#' @export
tissue_enrichment_df <- calc_gsea_per_tissue <- function(NMF, closest, narrowpeaks, num_top_peaks = 100) {
  cumul_gsea <- c()
  for (f in c(1:ncol(NMF@w))) {

    factor <- NMF@w[, f]
    names(factor) <- closest$num_id[match(names(factor), closest$rsID)]

    fna <- which(is.na(names(factor)))
    if (length(fna) >0){factor <- factor[-fna]}
    #factor <- tapply(factor, names(factor), mean)
    factor <- factor[order(factor, decreasing = T)]

    ## 'gene set' of the top x tissue peaks
    set_df <- narrowpeaks[order(narrowpeaks$log_pValue, decreasing = T), ]
    set <- as.character(set_df$num_id[1:num_top_peaks])


    # enrichment of tissue in factor
    #gsea <- fgsea::calcGseaStat(stats = factor, selectedStats = match(set, names(factor)), gseaParam = 0)
    #cumul_gsea <- c(cumul_gsea, gsea)

    res <- fgsea::fgseaMultilevel(pathways = list(factor = set), stats = factor)

    cumul_gsea <- c(cumul_gsea, res$padj)
  }
  return(cumul_gsea)
}

#' @title Split peak name
#'
#' @description Function to clean up narrowpeak file names
#' @param narrowpeaks Narrowpeaks dataframe
#' @param col_name Name to pull from metadata
#' @return DUpdated dataframe with formatted names
#' @examples
#' narrowpeaks <- split_peak_name(narrowpeaks)
#' @export
split_peak_name <- function(narrowpeaks, col_name = "name") {
  new_names <- unlist(lapply(X = narrowpeaks[[col_name]], FUN = function(X){strsplit(X, "bed_peak_")[[1]][2]}))
  narrowpeaks$peak_id <- new_names
  return(narrowpeaks)
}

#' @title Enrichment Calculation
#'
#' @description Function to calculate overepresentation enrichment
#' @param peaks dataframe where one column is a list of peak names and the other column is a c("yes", "no") of if it a top peak
#' @param databse top X peak names of the tissue/cell type
#' @return p values for enrichment
#' @examples
#' r <- enrichment_calc(peaks, database)
#' @export
enrichment_calc <- function(peaks, database) {
  '
  peaks = dataframe where one column is a list of peak names and the other column is a c("yes", "no") of if it a top peak
  databse = top X peak names of the tissue/cell type

  m n
  k p

  '

  query <- peaks$name[which(peaks$top == "yes")]
  nonquery <- peaks$name[which(peaks$top == "no")]

  m <- length(which(query %in% database)) # query peaks in database peaks
  n <- length(which(nonquery %in% database)) # non query peaks in the database peaks
  k <- length(query) - m # query peaks not in the database peaks
  p <- length(nonquery) - n # non query peaks not in the database peas

  enrich_table <- matrix(c(m, n, k, p),
                         nrow = 2,
                         dimnames = list(topPeaks = c("yes", "no"),
                                         SetPeaks = c("in", "out")))

  res <- fisher.test(enrich_table, alternative = "greater")
  return(res$p.value)
}

#' @title Factor Enrichments
#'
#' @description Function to calculate overepresentation enrichment for each factor in NMF model
#' @param NMF NMF object
#' @param filtered_closest dataframe of bedtools closest object
#' @param narrowpeaks dataframe of scATAC seq peak calling object
#' @param num_closest Number of closest peaks top use for enrichment calculation
#' @param num_peaks Number of top peaks to use for atac specificity and enrichment calc
#' @param factos Option to select only a specific set of NMF factors to calculate enrichment for
#' @return Data frame of factors by cell type or tissue enrichment
#' @examples
#' enrichment_df <- factor_enrichments(NMF = NMF_object, filtered_closest, narrowpeaks)
#' @export
factor_enrichments <- function(NMF, filtered_closest, narrowpeaks, num_closest = 500, num_peaks = 1000, factors = NA) {
  enrich_pval <-c()
  w_matrix <- as.matrix(NMF$w)
  if (is.na(factors[1])) {
    factors <- c(1:ncol(w_matrix))
  }else {
    factors <- factors
  }
  for (f in factors) {
    factor <- w_matrix[,f]
    filtered_closest$factor_contr <- factor[match(filtered_closest$rsID, names(factor))]
    top_closest <- filtered_closest[order(filtered_closest$factor_contr, decreasing = T), ]
    top_closest <- top_closest[1:num_closest, ]

    peaks <- data.frame(name = closest_peaks$peak,
                        top = ifelse(closest_peaks$peak %in% top_closest$peak, "yes", "no"))

    top_peaks <- narrowpeaks[order(narrowpeaks$log_pValue, decreasing = T), ]
    tissue_db <- top_peaks$peak[c(1:num_peaks)]

    pval <- enrichment_calc(peaks = peaks, database = tissue_db)
    enrich_pval <- c(enrich_pval, pval)

  }
  return(enrich_pval)
}

#' @title Fix Pval Float
#'
#' @description Function to replace pval = 0 with smallest r value
#' @param enrichment_df Output of enrichment_calc
#' @return Updated enrichment_df
#' @examples
#' enrichment_df <- fix_pval_float(enirchment_df)
#' @export
fix_pval_float <- function(enrichment_df) {
  t_enrichdf <- enrichment_df
  for (c in c(1:ncol(t_enrichdf))) {
    vals <- t_enrichdf[, c]
    vals[which(vals == 0)] <- .Machine$double.xmin
    t_enrichdf[, c] <- vals
  }
  return(t_enrichdf)
}

#' @title Format Factor Rank
#'
#' @description Function to format factor rank
#' @param factor_rank factor_rank
#' @return Updated factor rank
#' @examples
#' factor_rank <- format_factor_rank(factor_rank)
#' @export
format_factor_rank <- function(factor_rank) {
  split_genes <- c()
  for (i in c(1:length(factor_rank))) {
    n <- names(factor_rank)[i]
    lsc <- grep(",", n)
    if (length(lsc) > 0) {
      ns <- unlist(strsplit(n, ","))
      nl <- rep(factor_rank[i], length(ns))
      names(nl) <- ns
      split_genes <- c(split_genes, nl)
    }
  }
  new_rank <- factor_rank[-grep(",", names(factor_rank))]
  new_rank <- c(new_rank, split_genes)
  new_rank <- new_rank[order(abs(new_rank), decreasing = T)]
  return(new_rank)
}
















