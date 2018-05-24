The purpose of this document is to show the correlation matrix of pathway genes from real data. The correlation is displayed as heatmaps. Genes are ordered by hierarchical clustering using pearson correlation distance and complete linkage (default in hclust).

*  BRCA data from BRB-ArrayTools built-in dataset
*  Run KEGG pathway analysis BRCA1 vs BRCA2
*  15 KEGG pathways were found using the default options.

Output

*  PathwayClassComparison folder - BRB-ArrayTools Output
*  15 HTML files - interactive heatmap of correlation matrix of 15 pathways

```R
# x <- read.delim("clipboard", as.is=T)
# y[[1]] <- x[1:cumsum(n)[1], 4]
# y[[2]] <- x[(cumsum(n)[1]+1):cumsum(n)[2], 4]
# for(i in 3:15) y[[i]] <- x[(cumsum(n)[i-1]+1):cumsum(n)[i], 4]
# pathwayid <- y
# names(pathwayid) <- scan("clipboard", "")
# dump("pathwayid", file = "pathwayid.Rdmped")

source("pathwayid.Rdmped")
library(heatmaply)
# expdesign <- read.delim("EXPDESIGN.txt", as.is = TRUE)
geneid <- read.delim("GENEID.txt", as.is = TRUE)
lograt <- read.delim("LOGRAT.txt", as.is = TRUE, header = FALSE)
for(i in seq_along(pathwayid)) {
  ind <- match(pathwayid[[i]], geneid[, 1])
  mcor <- cor(t(lograt[ind, ]))
  colnames(mcor) <- rownames(mcor) <- pathwayid[[i]]
  heatmaply(mcor, margins = c(60, 80, 40, 0),
            k_col = 1, k_row = 1,
            colors = BrBG,
            limits = c(-1,1),
            distfun = "pearson",
      main = paste0(names(pathwayid)[i], " (g=", length(pathwayid[[i]]), ")"),
      file = paste0(names(pathwayid)[i], ".html"))
}
```
