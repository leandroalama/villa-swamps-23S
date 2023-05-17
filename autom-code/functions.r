load_packages <- function() {
  packages <- c("tidyverse", "phyloseq", "DESeq2", "dendextend", "viridis", "ampvis2", "data.table")
  lapply(packages, require, character.only = TRUE)
}

load_sample_info_tab <- function(file_path_md) {
  sample_info_tab <- read.csv(file_path_md, header = TRUE, row.names = 1, check.names = FALSE, sep = ",")
  sample_info_tab$color <- as.character(sample_info_tab$COLOUR)
  return(sample_info_tab)
}

load_count_tab <- function(file_path_ct) {
    count_tab <- as.matrix(read.table(file_path_ct, header = T, 
                           row.names = 1, check.names = F, sep = "\t"))
    return(count_tab)
}

load_tax_tab <- function(file_path_tt) {
    tax_tab <- as.matrix(read.table(file_path_tt, header=T,
                                    row.names=1, check.names=F, sep="\t"))
    return(tax_tab)
}

automate_deseq <- function(count_tab, sample_info_tab){
    
  deseq_counts <- DESeqDataSetFromMatrix(count_tab, colData = sample_info_tab, design = ~1)
  deseq_counts <- estimateSizeFactors(deseq_counts, type = "poscounts")
  deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)
  vst_trans_count_tab <- assay(deseq_counts_vst)
  return(vst_trans_count_tab)
}



create_dendrogram <- function(data, sample_info_tab) {
  euc_dist <- dist(t(data))
  euc_clust <- hclust(euc_dist, method="ward.D2")
  euc_dend <- as.dendrogram(euc_clust, hang=0.1)
  dend_cols <- as.character(sample_info_tab$color[order.dendrogram(euc_dend)])
  labels_colors(euc_dend) <- dend_cols
  plot(euc_dend, main="Dendogram", ylab="VST Euc. dist.")
}

create_phyloseq <- function(count_tab, sample_info_tab) {
  vst_count_phy <- otu_table(count_tab, taxa_are_rows=T)
  sample_info_tab_phy <- sample_data(sample_info_tab)
  vst_physeq <- phyloseq(vst_count_phy, sample_info_tab_phy)
  return(vst_physeq)
}

get_pcoa_eigenvals <- function(physeq) {
  vst_pcoa <- ordinate(physeq, method="MDS", distance="euclidean")
  eigen_vals <- vst_pcoa$values$Eigenvalues
  return(eigen_vals)
}

perform_pcoa <- function(count_tab, sample_info) {
  # Crear objeto phyloseq
  count_phy <- otu_table(count_tab, taxa_are_rows=T)
  sample_info_phy <- sample_data(sample_info)
  physeq <- phyloseq(count_phy, sample_info_phy)
  
  # Realizar PCoA
  pcoa <- ordinate(physeq, method="MDS", distance="euclidean")
  eigen_vals <- pcoa$values$Eigenvalues 
  
  # Crear gráfico
  plot_ordination(physeq, pcoa, color="TYPE") + 
    geom_point(size=1) + labs(col="TYPE") + 
    geom_text(aes(label=rownames(sample_info_phy), hjust=0.3, vjust=-0.4)) + 
    coord_fixed(sqrt(eigen_vals[2]/eigen_vals[1])) + ggtitle("PCoA") + 
    scale_color_manual(values=unique(sample_info$COLOUR[order(sample_info$COLOUR)])) + 
    theme(legend.text = element_text(size=10))
}

calculate_richness <- function(count_tab, tax_tab, sample_info_tab, x_var="TYPE", color_var="TYPE", measures=c("Chao1", "Shannon", "Simpson")) {
  
  # Crear objeto phyloseq
  count_tab_phy <- otu_table(count_tab, taxa_are_rows=T)
  tax_tab_phy <- tax_table(tax_tab)
  sample_info_tab_phy <- sample_data(sample_info_tab)
  ASV_physeq <- phyloseq(count_tab_phy, tax_tab_phy, sample_info_tab_phy)
  
  # Calcular riqueza
  richness <- plot_richness(ASV_physeq, x = x_var , color = color_var, measures = measures) + geom_boxplot()
  
  return(richness)
}


create_stack_plot <- function(count_tab, tax_tab, sample_info_tab, taxrank="class", threshold=5) {
  # create phyloseq object
  count_tab_phy <- otu_table(count_tab, taxa_are_rows=T)
  tax_tab_phy <- tax_table(tax_tab)
  sample_info_tab_phy <- sample_data(sample_info_tab)
  ASV_physeq <- phyloseq(count_tab_phy, tax_tab_phy, sample_info_tab_phy)
  
  # calculate proportions by taxonomic group
  clases_counts_tab <- otu_table(tax_glom(ASV_physeq, taxrank=taxrank))
  clases_tax_vec <- as.vector(tax_table(tax_glom(ASV_physeq, taxrank=taxrank))[,3])
  rownames(clases_counts_tab) <- as.vector(clases_tax_vec)
  unclassified2_tax_counts <- colSums(count_tab) - colSums(clases_counts_tab)
  clases_and_unidentified_counts_tab <- rbind(clases_counts_tab, "Unclassified"=unclassified2_tax_counts)
  clases_taxa_counts_tab <- clases_and_unidentified_counts_tab
  clases_taxa_proportions_tab <- apply(clases_taxa_counts_tab, 2, function(x) x/sum(x)*100)
  
  # filter by threshold and create "Other" category
  temp_filt_clases_taxa_proportions_tab <- data.frame(clases_taxa_proportions_tab[apply(clases_taxa_proportions_tab, 1, max) > threshold, ])
  filtered2_proportions <- colSums(clases_taxa_proportions_tab) - colSums(temp_filt_clases_taxa_proportions_tab) 
  filt_clases_taxa_proportions_tab <- rbind(temp_filt_clases_taxa_proportions_tab, "Other"=filtered2_proportions)
  clases_stack <- filt_clases_taxa_proportions_tab
  clases_stack$Major_Taxa <- row.names(clases_stack)
  clases_stack.g <- gather(clases_stack, Sample, Proportion, -Major_Taxa)
  
  # create plot
  datos_stack<-data.frame("Sample"=row.names(sample_info_tab), "char"=sample_info_tab$TYPE, "color"=sample_info_tab$color, stringsAsFactors=F) 
  stack_plot <- ggplot(clases_stack.g, aes(x=Sample, y=Proportion, fill=Major_Taxa)) +
    geom_bar(width=0.6, stat="identity") +
    theme_bw() +
    theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1), legend.title=element_blank()) +
    labs(x="Sample", y="% of 16S rRNA gene copies recovered", title="All samples") +
    facet_grid(~1, scales = 'free_x', space = 'free_x')
  
  return(stack_plot)
}

heatmap_function <- function(count_tab, tax_tab, sample_info_tab,
                         group_by, facet_by, tax_aggregate, tax_add,
                         plot_colorscale, plot_legendbreaks){
  
  #Convertir datos en phyloseq object
  count_tab_phy <- otu_table(count_tab, taxa_are_rows=T)
  tax_tab_phy <- tax_table(tax_tab)
  sample_info_tab_phy <- sample_data(sample_info_tab)
  ASV_physeq <- phyloseq(count_tab_phy, tax_tab_phy, sample_info_tab_phy)
  
  #Crear otutable y metadataX
  otutable <- data.frame(OTU = rownames(phyloseq::otu_table(ASV_physeq)@.Data),
                         phyloseq::otu_table(ASV_physeq)@.Data,
                         phyloseq::tax_table(ASV_physeq)@.Data,
                         check.names = FALSE)
  
  metadataX <- data.frame(phyloseq::sample_data(ASV_physeq), 
                           check.names = FALSE)
  
  metadataX <-setDT(metadataX, keep.rownames = TRUE)[]
  
  #Cargar datos en ampvis2
  av2 <- amp_load(otutable, metadataX)
  
  #Crear heatmap con ampvis2
  heat_map_plot <- amp_heatmap(av2, 
                               group_by = group_by, 
                               facet_by = facet_by, 
                               plot_values = FALSE,
                               tax_show = 10,
                               tax_aggregate = tax_aggregate,
                               tax_add = tax_add,
                               plot_colorscale = plot_colorscale,
                               plot_legendbreaks = plot_legendbreaks)
  
  return(heat_map_plot)
}








