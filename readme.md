The purpose of this document is to show the correlation matrix of pathway genes from real data. The correlation is displayed as heatmaps. Genes are ordered by hierarchical clustering using pearson correlation distance and complete linkage (default in hclust).

*  BRCA data from [BRB-ArrayTools](https://brb.nci.nih.gov/BRB-ArrayTools/) built-in dataset
*  Run KEGG pathway analysis BRCA1 vs BRCA2 using the default options
*  15 KEGG pathways were found

Output

*  PathwayClassComparison folder - BRB-ArrayTools Output
*  15 svg files - heatmap of correlation matrix of 15 pathways
*  [1 svg file](hsa04110_hsa05213.svg) - heatmap of cor matrix from combining the 1st and 2nd pathways
*  [1 svg file](hsa04110_hsa05213_hsa05219_hsa05223.svg) - heatmap of cor matrix from combining 4 pathways (102 genes)

```R
# x <- read.delim("clipboard", as.is=T)
# y[[1]] <- x[1:cumsum(n)[1], 4]
# y[[2]] <- x[(cumsum(n)[1]+1):cumsum(n)[2], 4]
# for(i in 3:15) y[[i]] <- x[(cumsum(n)[i-1]+1):cumsum(n)[i], 4]
# pathwayid <- y
# names(pathwayid) <- scan("clipboard", "")
# dump("pathwayid", file = "pathwayid.Rdmped")

library("ComplexHeatmap")
# I can add a color bar at the top of the heatmap (cf heatmaply and gplots)
# I can control the clustering
# Better use of space compared to gplots::heatmap.2()

source("pathwayid.Rdmped")
# expdesign <- read.delim("EXPDESIGN.txt", as.is = TRUE)
geneid <- read.delim("GENEID.txt", as.is = TRUE)
lograt <- read.delim("LOGRAT.txt", as.is = TRUE, header = FALSE)
for(i in seq_along(pathwayid)) {
  ind <- match(pathwayid[[i]], geneid[, 1])
  mcor <- cor(t(lograt[ind, ]))
  colnames(mcor) <- rownames(mcor) <- pathwayid[[i]]
  ord <- hclust(as.dist(1-mcor))$order
  svg(paste0(names(pathwayid)[i], ".svg"))
  Heatmap( mcor[ord, ord], cluster_rows = FALSE, cluster_columns = FALSE,
     row_names_gp = gpar(fontsize = 6),
     column_names_gp = gpar(fontsize = 6),
     column_title = paste0(names(pathwayid)[i], " (g=", length(pathwayid[[i]]), ")"), name = "value")  
  dev.off()
}
```

Put two gene sets together to understand inter-geneset correlation
```R
i <- 1
ind1 <- match(pathwayid[[i]], geneid[, 1])
mcor <- cor(t(lograt[ind1, ]))
ord1 <- hclust(as.dist(1-mcor))$order
i <- 2
ind2 <- match(pathwayid[[i]], geneid[, 1])
mcor <- cor(t(lograt[ind2, ]))
ord2 <- hclust(as.dist(1-mcor))$order

mcor <- cor(t(lograt[c(ind1[ord1], ind2[ord2]), ]))
colnames(mcor) <- rownames(mcor) <- c(pathwayid[[1]][ord1], pathwayid[[2]][ord2])

ha_column = HeatmapAnnotation(df = data.frame(KEGG = c(rep("hsa04110", 40), rep("hsa05213", 24))),
    col = list(KEGG = c("hsa04110" =  "seagreen", "hsa05213" = "darkorange")))
ht1 <- Heatmap( mcor, cluster_rows = FALSE, cluster_columns = FALSE, name = "value", top_annotation = ha_column, row_names_gp = gpar(fontsize = 6), column_names_gp = gpar(fontsize = 6))
svg("hsa04110_hsa05213.svg")
draw(ht1, annotation_legend_side = "bottom")
dev.off()
```

Combine 4 gene sets
```R
ind <- vector("list", 4)
ord <- vector("list", 4)
for(i in 1:4) {
  ind[[i]] <- match(pathwayid[[i]], geneid[, 1])
  mcor <- cor(t(lograt[ind[[i]], ]))
  ord[[i]] <- hclust(as.dist(1-mcor))$order
}

mcor <- cor(t(lograt[unlist(mapply(function(x,y) x[y], ind, ord)), ]))
colnames(mcor) <- rownames(mcor) <- unlist(mapply(function(x, y) x[y], pathwayid[1:4], ord))

ha_column = HeatmapAnnotation(df = data.frame(KEGG = rep(names(pathwayid[1:4]), sapply(pathwayid[1:4], length))),
    col = list(KEGG = c("hsa04110" =  "seagreen", "hsa05213" = "darkorange", "hsa05219" = "yellow", "hsa05223" = "pink")))
ht1 <- Heatmap(mcor, cluster_rows = FALSE, cluster_columns = FALSE,
               name = "value",
               top_annotation = ha_column,
               row_names_gp = gpar(fontsize = 3),
               column_names_gp = gpar(fontsize = 3))
svg("hsa04110_hsa05213_hsa05219_hsa05223.svg")
draw(ht1, annotation_legend_side = "bottom")
dev.off()
```
