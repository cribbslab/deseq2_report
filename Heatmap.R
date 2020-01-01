

dds<- DESeqDataSetFromMatrix(countData=df_mRNA,
                             colData=meta_data,
                             design= ~treatment)

# make the DeSeqDataSet
dds <- DESeq(dds)

# use the log transform on the data set
rld <- rlog(dds, blind=F)
topVarianceGenes <- head(order(rowVars(assay(rld)), decreasing=T),100)
matrix <- assay(rld)
[ topVarGenes ]
matrix <- matrix - rowMeans(matrix)

# select the 'contrast' you want
annotation_data <- as.data.frame(colData(rld)[c("ConditionA","ConditionB")])
pheatmap(matrix, annotation_col=annotation_data)

plotUpDownSigGenes <- function(results, colNums, rld, title) {
  
  # make the lists
  upgenes <- rownames(head(results[ order( results$log2FoldChange ), ], n=20))
  downgenes <- rownames(head(results[ order( -results$log2FoldChange ), ], n=20))
  
  # this gives us the rows we want
  rows <- match(upgenes, row.names(rld))
  mat <- assay(rld)[rows,colNums]
  mat <- mat - rowMeans(mat)
  
  # the labels are hard coded at the moment :(
  df <- as.data.frame(colData(rld)[c("labelA","labelB")])
  pheatmap(mat, fontsize=5, annotation_col=df, main=paste(title,"top 20 up genes"))
  
  # this gives us the rows we want
  rows <- match(downgenes, row.names(rld))
  mat <- assay(rld)[rows,colNums]
  mat <- mat - rowMeans(mat)
  
  df <- as.data.frame(colData(rld)[c("labelA","labelB")])
  pheatmap(mat, fontsize=5, annotation_col=df, main=paste(title,"top 20 down genes"))
}


contrastDEGenes <- subset(results(dds, contrast=c("A","B")), padj < 0.05)

# this part is kind of funky
# the function needs to know which columns
# correspond to the samples (to pull from rld)
aCols <- c(1,2,3)
bCols <- c(4,5,6)

# get the log transforms again
rld <- rlog(dds, blind=F)

# call
plotUpDownSigGenes(
  contrastDEGenes,
  c(aCols, bCols),
  rld,
  "Title for the plot"
)