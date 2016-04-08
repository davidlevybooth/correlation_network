### Weighted correlation network analysis with igraph
### D. Levy-Booth 04/08/2016

############################################################

### From Wikipedia: Weighted correlation network analysis, also known
### as weighted gene co-expression network analysis (WGCNA), is a widely
### used data mining method especially for studying biological networks
### based on pairwise correlations between variables.
### https://en.wikipedia.org/wiki/Weighted_correlation_network_analysis 

### Only positive correlations are used in this analysis to study
### co-occurance and co-abundance. Exclusion is much more difficult to parse
### with a sparse dataset such as OTU abundance. False-discovery-rate (FDR)
### correction is used to only plot significant correlations. 

### In addition, environmental niches can be explored by mapping factor
### "preference" as a mean value for each OTU to the network (Robador et al. 2016)
### http://www.nature.com/ismej/journal/v10/n4/full/ismej2015157a.html

### Analysis Steps: 
### 1. Optional: Compile OTU, Taxonomy and Sample metadata
### 2. Optional: Filter OTU data, e.g., removal of low-abundance OTUs
### 3. Create and filter correlation matrix
### 4. Create correlation network
### 5. Optional: Calculate OTU means for environmental factors
### 6. Apply OTU means as node attributes
### 7. Plot graphs with OTU variable means as node colors 
### 8. Optional: Test node coordinant dissimiarity against OTU means
### 9. Optional: Export network to ploting software (e.g., Cytoscape)
### 10, Optional: plot OTU variable means against variable distributions in samples 
############################################################

### Load libraries 
require(igraph)
require(phyloseq)
require(Hmisc)
require(ecodist)

### 1. Optional: Compile OTU, Taxonomy and Sample metadata
### Load OTU BIOM table into phyloseq package
biom_file = "biom_table.json.biom"
phyloseq_object = import_biom(biom_file, parse_function=parse_taxonomy_greengenes)

### Add metadata to phyloseq
MAP=import_qiime_sample_data("metadata_table.txt")

### Create and edit phyloseq object
phyloseq_object <- merge_phyloseq(phyloseq_object, MAP)

### 2. Optional: Filter OTU data, e.g., removal of low-abundance OTUs
### Phyloseq not required if OTU table already in proper form 
phyloseq_object_pruned <- prune_samples(sample_data(phyloseq_object)$Subset_variable=="attribute", phyloseq_object) # To subset samples
phyloseq_object_merged <- merge_samples(phyloseq_object_pruned, sample_data(phyloseq_object_pruned)$Merge_variable) # To merge samples
phyloseq_object_filtered <- filter_taxa(phyloseq_object_merged, function(x) sum(x) > 100, TRUE) # To remove OTUS with < 100 counts
OTU_labels <- tax_table(phyloseq_object_filtered)[,"Taxonomic_Level"] # Create label vector at desired taxonomic level 
phyloseq_object_transformed = transform_sample_counts(phyloseq_object_filtered, function(x) x/sum(x)) # calcualte relative abundance
otu_table <- as.matrix(otu_table(phyloseq_object_transformed)) #extract just OTUs

#Pull out metadata 
metadata <- sample_data(phyloseq_object)

### 3. Create and filter correlation matrix
otu_cor <- rcorr(otu_table) # calculate correlation table # Suggest keeping default Pearson correlations
otu_cor_r <- otu_cor$r
cor_len <- length(otu_cor_r) - nrow(otu_cor_r) # Remove self-correlations from n 
adusted_p <- p.adjust(otu_cor$P, method="fdr") #FDR-correct p-values

# Option 1. To select correlations by FDR-corrected p-values 
otu_cor_r[ adust_p > .05 ] <- 0
## OR ##
# Option 2. To select correlations by r-threshold (here > 0.8)
rna_cor[ rna_cor < .8 ] <- 0

# Remove self-correlations
diag(rna_cor_r) <- 0


### 4. Create correlation network
otu_graph <- graph.adjacency(otu_cor_r, weighted=TRUE, mode="lower")


### 5. Optional: Calculate OTU means for environmental factors
### 5a. "plotmax" Function to calculate the plot in which each OTU is at maximum 
### OTU table rownames must be sample names 
plotmax <- function(matrix) {
  plotmaxs <- rep(0, ncol(matrix))
  for (i in 1:ncol(matrix)) {
    row_index = which.max( matrix[,i] )
    row_name = rownames(matrix[row_index])
    plotmaxs[i] = row_name
  }
  return(plotmaxs)
}

plot_max <- plotmax(otu_table)


### 5b. "var_weights" function to calcualte the "environmental preference" of each
### OTU. variable must be in the form of an atomic vector: i.e., metadata$variable
var_weights <- function (otu_table, variable) {
  otu_var = variable * otu_table # 1
  otu_sum <- colSums(otu_table) # 2
  var_sum <- colSums(otu_var) # 3
  var_weight <- var_sum/otu_sum # 4
 return(var_weight)
}

variable_weight <- var_weights(otu_table, metadata$variable)

### 6. Apply OTU means as node attributes
V(otu_graph)$variable <- variable_weight

### 7. Plot graphs with OTU variable means as node colors 
### Set color palette for variable 
V(otu_graph)$color <- 1:(length(variable_weight)) 
otu_graph$palette <- scales::dscale(V(otu_graph)$variable %>% cut(9), sequential_pal)

### Plot that shit
plot(otu_graph, layout=layout.fruchterman.reingold, vertex.size=6, vertex.label=NA)


### 8. Optional: Test node coordinate dissimiarity against OTU means
### Calculate node X,Y coordinates (default: Fruchterman-Reingold layout)
graph_layout <- layout.spring(otu_graph)

mantel(dist(graph_layout) ~ dist(variable_weight), nperm=9999)

### 9. Optional: Export network to ploting software (e.g., Cytoscape)
write.graph(otu_graph, file="filename.graphml", format="graphml")


### 10. Optional: plot OTU variable means against variable distributions in samples
variable <- metadata$variable

### Function to plot density distributions 
dist_plot <- function (actual, weight, title) {
  
  red_AD_Col <- rgb(1,0,0,0.2)
  red_WD_Col <- rgb(0,0,1,0.2)
  
  red_AD <- density(actual)
  red_WD <- density(weight)
  xlim <- range(red_AD$x)
  
  par(mar=c(5, 4, 4, 5) + 0.1)
  plot(red_AD, xlim = xlim, col="red", ylab = paste("Sample", title, "Density",sep=" "), xlab = title, main="")
  polygon(red_AD, border = red_AD_Col, density = -1, col = red_AD_Col)
  par(new=TRUE)
  plot(red_WD,xlim = xlim,col="blue", axes = FALSE, bty = "n", xlab = "", ylab = "", main=paste("Mean OTU", title, "Distribution", sep=" "))
  polygon(red_WD, border = red_WD_Col, density = -1, col = red_WD_Col)
  axis(4)
  mtext(paste("Mean", title, "Density", sep=" "), side=4, line=3)
  legend('topleft',c(paste("Sample", title, sep=" "),paste("Mean OTU", title, sep=" ")),
         fill = c(red_AD_Col, red_WD_Col), 
         col= c("red", "blue"),
         bty = 'n',
         border = NA)
  
}

#Plot it baby!
dist_plot(variable, variable_weight, "Variable Name")
