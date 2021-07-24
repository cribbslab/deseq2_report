

ensembl_to_symbol <- function(dataframe, ensembl_column){
  
  dataframe_tmp <- dataframe %>% 
    select(ensembl_column)
  data <- unlist(dataframe_tmp)
  data <- as.vector(data)
  annots <-  AnnotationDbi::select(org.Hs.eg.db, keys=data,
                                   columns=c("SYMBOL","GENENAME"), keytype = "ENSEMBL")

  result <- merge(dataframe, annots, by.x=ensembl_column, by.y="ENSEMBL")
  return(result)
  }


setClass(Class="filter_genes_return",
         representation(
           res="data.frame",
           sig="data.frame"
         )
)

filter_genes <- function(result, name, species='human'){
  
  test <- as.data.frame(result)
  
  data <- as.vector(rownames(test))
  
  if(species == 'human'){
    annots <-  AnnotationDbi::select(org.Hs.eg.db, keys=data,
                                     columns="SYMBOL", keytype = "ENSEMBL")
  }else{
    annots <-  AnnotationDbi::select(org.Mm.eg.db, keys=data,
                                     columns="SYMBOL", keytype = "ENSEMBL")
  }

  result <- merge(test, annots, by.x="row.names", by.y="ENSEMBL")
  res <- result %>% 
    dplyr::select(log2FoldChange, SYMBOL, baseMean, padj, Row.names) %>% 
    na.omit()
  
  sig <- res %>% 
    dplyr::filter(log2FoldChange > 1 | log2FoldChange < -1) %>% 
    dplyr::filter(padj < 0.05)
  dir.create(file.path("results"), showWarnings = FALSE)
  sig_name = paste("results/", name,"_sig.csv", sep="")
  sig_name_tsv = paste("results/", name,"_sig.tsv", sep="")
  res_name = paste("results/",name,"_res.csv", sep="")
  res_name_tsv = paste("results/", name,"_res.tsv", sep="")
  write_csv(sig, sig_name)
  write_csv(res, res_name)
  write_tsv(sig, sig_name_tsv)
  write_tsv(res, res_name_tsv)
  return(new("filter_genes_return",
             res=res,
             sig=sig))
}


filter_genes_single <- function(result, name){
  
  test <- as.data.frame(result)
  
  data <- as.vector(rownames(test))
  annots <-  AnnotationDbi::select(org.Hs.eg.db, keys=data,
                                   columns="SYMBOL", keytype = "ENSEMBL")
  
  result <- merge(test, annots, by.x="row.names", by.y="ENSEMBL")
  res <- result %>% 
    dplyr::select(log2FoldChange, SYMBOL, baseMean, padj, Row.names) %>% 
    na.omit()
  
  sig <- res %>% 
    filter(log2FoldChange > 2 | log2FoldChange < -2)
  
  sig_name = paste("results/", name,"_sig.csv", sep="")
  sig_name_tsv = paste("results/", name,"_sig.tsv", sep="")
  res_name = paste("results/",name,"_res.csv", sep="")
  res_name_tsv = paste("results/", name,"_res.tsv", sep="")
  write_csv(sig, sig_name)
  write_csv(res, res_name)
  write_tsv(sig, sig_name_tsv)
  write_tsv(res, res_name_tsv)
  return(list("sig"= sig, "res"= res))
}

run_deseq2_full <- function(df_mRNA, meta_data, model){
  
  
  dds<- DESeqDataSetFromMatrix(countData=df_mRNA,
                               colData=meta_data,
                               design=as.formula(model)) 
  
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  
  dds <- DESeq(dds, parallel=FALSE)
  
  return(dds)
}





setClass(Class="DESeq2_return",
         representation(
           res="DESeqResults",
           dds="DESeqDataSet"
         )
)


run_deseq2 <- function(df_mRNA, meta_data, control="untreated", test="treated", value, model){
  
  df_mRNA = df_mRNA[,rownames(meta_data)]
  
  
  dds<- DESeqDataSetFromMatrix(countData=df_mRNA,
                               colData=meta_data,
                               design=as.formula(model))
  
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  
  dds <- DESeq(dds)
  
  res <- results(dds, contrast = c(value, test, control))
  
  return(new("DESeq2_return",
             res=res,
             dds=dds))
}

setClass(Class="DESeq2_lrt_return",
         representation(
           res="DESeqResults",
           dds="DESeqDataSet"
         )
)


run_deseq2_lrt <- function(df_mRNA, meta_data, control="untreated", test="treated", value, model,
                           reduced){
  
  df_mRNA = df_mRNA[,rownames(meta_data)]
  
  
  dds<- DESeqDataSetFromMatrix(countData=df_mRNA,
                               colData=meta_data,
                               design=as.formula(model))
  
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  
  dds <- DESeq(dds, test="LRT", reduced =as.formula(reduced))
  
  res <- results(dds, contrast = c(value, test,control))
  
  return(new("DESeq2_lrt_return",
             res=res,
             dds=dds))
}



plot_volcano <- function(res, species="human"){
  
  test <- as.data.frame(res)
  
  data <- as.vector(rownames(test))
  
  if(species == 'human'){
  annots <-  AnnotationDbi::select(org.Hs.eg.db, keys=data,
                                   columns="SYMBOL", keytype = "ENSEMBL")
  }else{
  annots <-  AnnotationDbi::select(org.Mm.eg.db, keys=data,
                                     columns="SYMBOL", keytype = "ENSEMBL")
  }
  result <- merge(test, annots, by.x="row.names", by.y="ENSEMBL")
  res <- result %>% 
    dplyr::select(log2FoldChange, SYMBOL, baseMean, padj, Row.names) %>% 
    na.omit()
  
  
  mutateddf <- mutate(res, sig=ifelse(res$padj<0.01, "padj<0.01", "Not Sig")) #Will have different colors depending on significance
  input <- cbind(gene=rownames(mutateddf), mutateddf )
  input <- input %>% 
    arrange(input$padj)
  
  symbol_data <- head(input, 30)

  #convert the rownames to a column
  volc = ggplot(input, aes(log2FoldChange, -log10(padj))) + #volcanoplot with log2Foldchange versus pvalue
    geom_point(aes(col=sig)) + #add points colored by significance
    geom_point(data=symbol_data, aes(log2FoldChange, -log10(padj)), colour="red") +
    ggtitle("Volcano") #e.g. 'Volcanoplot DESeq2'
  
  #setEPS()
  #postscript("MUG_volcano.eps")
  volcano <- volc+geom_text_repel(data=symbol_data, aes(label=`SYMBOL`)) + scale_colour_Publication() + theme_bw()#adding text for the genes
  return(volcano)
}




plotPCA34 <- function (object, ...) 
{
  .local <- function (object, intgroup = "condition", ntop = 500, 
                      returnData = FALSE) 
  {
    rv <- rowVars(assay(object))
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                       length(rv)))]
    pca <- prcomp(t(assay(object)[select, ]))
    percentVar <- pca$sdev^2/sum(pca$sdev^2)
    if (!all(intgroup %in% names(colData(object)))) {
      stop("the argument 'intgroup' should specify columns of colData(dds)")
    }
    intgroup.df <- as.data.frame(colData(object)[, intgroup, 
                                                 drop = FALSE])
    group <- if (length(intgroup) > 1) {
      factor(apply(intgroup.df, 1, paste, collapse = " : "))
    }
    else {
      colData(object)[[intgroup]]
    }
    d <- data.frame(PC3 = pca$x[, 3], PC4 = pca$x[, 4], group = group, 
                    intgroup.df, name = colnames(object))
    if (returnData) {
      attr(d, "percentVar") <- percentVar[1:2]
      return(d)
    }
    ggplot(data = d, aes_string(x = "PC3", y = "PC4", color = "group")) + 
      geom_point(size = 3) + xlab(paste0("PC3: ", round(percentVar[1] * 
                                                          100), "% variance")) + ylab(paste0("PC4: ", round(percentVar[2] * 
                                                                                                              100), "% variance")) + coord_fixed()
  }
  .local(object, ...)
}

theme_Publication <- function(base_size=14, base_family="arial") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.title = element_text(face="italic"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}


pca <- function(vsd, gene_prop = 1, batch = NA, 
                       remove_batch_effect = FALSE){
  
  ##Using PCAtools 
  
  ###Number of genes to use in PCA
  if(gene_prop < 0 | gene_prop > 1){stop("gene_prop must be between 0 and 1")}
  
  selectgenes <- floor(gene_prop * nrow(vsd))
  
  
  
  ##To correct for batch effect: 
  if(is.na(batch) & remove_batch_effect){stop("Batch column must be provided to remove batch effect")}
  if(remove_batch_effect){
    assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd[[batch]])
  }
  
  
  ##Run PCA       
  pca_data <- PCAtools::pca(mat = assay(vsd), 
                            metadata = colData(vsd), 
                            removeVar = (1-gene_prop)) 
  return(pca_data)
}

pca_biplot <- function(pca_data, vsd, gene_prop = 1, design, batch = NA, title, 
                       remove_batch_effect = FALSE){
  
  selectgenes <- floor(gene_prop * nrow(vsd))
  
  ##plot PC1 vs PC2
  Biplot <- biplot(pca_data, 
                   showLoadings = FALSE, 
                   lab = NULL,
                   #lab = pca_data$metadata$rowid,
                   #drawConnectors = TRUE,
                   #boxedLabels = TRUE,
                   #labhjust = 2,
                   #labvjust = 4,
                   colby = design, 
                   #shape = NULL,
                   shape = if(!is.na(batch)){batch},
                   #encircle = TRUE,
                   #ellipse = TRUE,
                   #ellipseConf = 0.95,
                   #ellipseFill = TRUE,
                   #ellipseAlpha = 1/4,
                   #ellipseLineSize = 1.0,
                   #showLoadingsNames = TRUE,
                   #boxedLoadingsNames = TRUE,
                   #ntopLoadings = 10,
                   #sizeLoadingsNames = 2,
                   #colLoadingsNames = "black",
                   #fillBoxedLoadings = "white",
                   #drawConnectorsLoadings = TRUE,
                   #alphaLoadingsArrow = 1,
                   #lengthLoadingsArrowsFactor = 1,
                   #legendPosition = 'right', 
                   legendLabSize = 10, 
                   legendIconSize = 4.0) %>%
    
    ggpar( title = title,
           caption = paste0("PCA using: ", selectgenes, " genes") ,
           legend = "right",
           legend.title = list(color = design, shape = batch),
           palette = c("red4", "blue4", "green4", "orange", "black", "brown", "violet"),
           shapekey = c(20, 4, 6),
           font.title = c("bold", "brown", 18),
           font.subtitle = c("bold", "dark grey", 12),
           font.caption = c("bold.italic", "royal blue", 10),
           font.legend = c("bold", "black", 8),
           font.x = c("bold", "black", "11"), 
           font.y = c("bold", "black", "11"),
           font.tickslab = c("bold", "black", "8"),
           font.family = "Arial",
           x.text.angle = 0,
           ggtheme = theme_pubr()
    )
  #Bp <<- Biplot
  #return(Bp)
  ifelse(remove_batch_effect,
         Biplot$labels$subtitle <- paste0(Biplot$labels$subtitle, ", batch corrected"),
         FALSE)
  return(Biplot)  
  #ggsave(paste0(plots_dir, title, ifelse(remove_batch_effect, "_batch_corrected_", "_"), "PCA_all_samples.png"))
}

Plot_pca_loadings <- function(vsd, 
                              gene_prop = 1, 
                              nPCs = 6,
                              title = Plot_title,
                              batch = NA,
                              remove_batch_effect = FALSE,
                              results_dir = "results_pca/"
) {
  
  dir.create(file.path(results_dir), showWarnings = FALSE)
  ##Using PCAtools 
  
  ###Number of genes to use in PCA
  if(gene_prop < 0 | gene_prop > 1){stop("gene_prop must be between 0 and 1")}
  
  selectgenes <- floor(gene_prop * nrow(vsd))
  
  ##To correct for batch effect: 
  if(is.na(batch) & remove_batch_effect){stop("Batch column must be provided to remove batch effect")}
  if(remove_batch_effect){
    assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd[[batch]])
  }
  
  
  
  
  ##Run PCA       
  pca_data <- PCAtools::pca(mat = assay(vsd), 
                            metadata = colData(vsd), 
                            removeVar = (1-gene_prop)) 
  
  ##Create a data.frame of selected PC loadings
  Loads_PCs <- PCAtools::getLoadings(pca_data, 1) %>% 
    as.data.frame() %>%
    tibble::rownames_to_column("ensembl_gene_id")
  
  for (i in 2:nPCs) {
    Loads_PCs <- Loads_PCs %>% dplyr::mutate(getLoadings(pca_data, i))
    
  }
  
  ##save in results folder
  write.table(Loads_PCs,
              paste0(results_dir, title, ifelse(remove_batch_effect, "_batch_corrected_", "_"),
                     "PCA_loadings.tsv"),
              col.name=TRUE,
              sep="\t",
              na = "NA",
              row.names=FALSE,
              quote=FALSE)
  
  ##Make a loadings plot for selected loadings
  PLoadings <- plotloadings(pca_data,
                            components = getComponents(pca_data, 1:nPCs),
                            rangeRetain = 0.1,
                            labSize = 3.0,
                            #title = 'Loadings plot',
                            #subtitle = 'PC1, PC2',
                            #caption = 'Top 10% variables',
                            shape = 23,
                            shapeSizeRange = 1,1,
                            col = c("dark red", "white", "dark blue"),
                            colMidpoint = 0,
                            
                            drawConnectors = TRUE,
                            lengthConnectors = unit(0.005, "npc"),
                            positionConnectors = "right",
                            #labvjust = 0.5, 
                            
  ) %>%
    
    
    ggpar( title = title,
           subtitle = "All samples",
           caption = "PC Loadings - Top 10% variables: ",
           legend = "right",
           font.title = c("bold", "brown", 18),
           font.subtitle = c("bold", "dark grey", 14),
           font.caption = c("bold.italic", "royal blue", 10),
           font.legend = c("bold", "black", 8),
           font.x = c("bold", "black", "11"), 
           font.y = c("bold", "black", "11"),
           font.tickslab = c("bold", "black", "8"),
           font.family = "Comic Sans MS",
           x.text.angle = 0,
           
           
           
           ggtheme = theme_pubr()
           
    )
  
  
  return(PLoadings)
  ##ggsave(paste0(plots_dir, title, ifelse(remove_batch_effect, "_batch_corrected", FALSE), "_PCA_all_samples.png"))
  
}



GSEA_plots <-  function(pathways, title, results_dir){
  
  ##read in the annotated results table
  results_annotated <- results_dir
  
  ##create ranked gene list for fgsea
  gseaInput <-  results_annotated %>% 
    dplyr::filter(!is.na(ENTREZID), !is.na(log2FoldChange)) %>% 
    arrange(log2FoldChange)
  
  ranks <- pull(gseaInput,log2FoldChange)
  names(ranks) <- gseaInput$SYMBOL
  
  
  ###choose pathway databases to load. 
  
  ##REACTOME GENE SETS
  HS_CP_REACTOME <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
  HS_CP_REACTOME <- split(x = HS_CP_REACTOME$gene_symbol, f = HS_CP_REACTOME$gs_name)
  
  
  ##HALLMARKS GENE SETS
  HS_HALLMARK <- msigdbr(species = "Homo sapiens", category = "H")
  HS_HALLMARK <- split(x = HS_HALLMARK$gene_symbol, f = HS_HALLMARK$gs_name)
  
  ##KEGG GENE SETS
  HS_CP_KEGG <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")
  HS_CP_KEGG <- split(x = HS_CP_KEGG$gene_symbol, f = HS_CP_KEGG$gs_name)
  
  
  ##Choose gene sets list
  #pathways <- HS_HALLMARK
  pathway_dbs <- list("HS_HALLMARK" = HS_HALLMARK, "HS_CP_KEGG" = HS_CP_KEGG, "HS_CP_REACTOME" = HS_CP_REACTOME)
  a <-c("HS_HALLMARK", "HS_CP_KEGG", "HS_CP_REACTOME")
  b <- a[str_which(a, str_to_upper(pathways))]
  pathways <- pathway_dbs[[b]]
  ##GSEA analysis using the "ranks" table against the pathways in the database.
  
  fgseaRes <- fgsea(pathways, ranks, minSize=15, maxSize = 500, nperm=1000)
  
  
  
  ## Show results in a nice table
  
  fgseaResTidy <- fgseaRes %>%
    as_tibble() %>%
    arrange(desc(NES)) %>%
    dplyr::filter(padj <= 0.05) %>%
    dplyr::mutate(Pathway = str_replace_all(pathway, "_", " "),
                  Adjusted_p_value = padj, 
                  Normalised_Enrichment_Score = NES,
                  Size = size,
                  Key_genes = leadingEdge) %>%
    dplyr::select(Pathway, Normalised_Enrichment_Score, Size, Adjusted_p_value) 
  #%>% dplyr::filter(abs(Normalised_Enrichment_Score) > 2.3)
  
  ##save results table
  #return(fgseaResTidy)
  write.table(fgseaResTidy,
              file = paste0("results_gsea/", title, "fgsea_", b ,"_pathways.tsv"),
              col.name=TRUE,
              sep="\t",
              na = "NA",
              append = FALSE,
              row.names=FALSE,
              quote=FALSE)
  
  
  
  if(dim(fgseaResTidy)[1] != 0){
    FGSEA_Results <- fgseaResTidy %>% kbl(col.names = c("Pathway", "Normalised Enrichment Score", "Size", "Adjusted p-value"), 
                                          caption = "L363 Car-resistant cells, Top Differentially Expressed Pathways",
                                          digits = 2,
                                          escape = F, 
                                          align = "c"
    ) %>% 
      kable_paper(html_font = "Comic Sans MS",
                  full_width = T, bootstrap_options = c("striped", "bordered", "hover")) %>%
      column_spec(2, color = "white", background = ifelse(fgseaResTidy$Normalised_Enrichment_Score < 0, "red", "green")) %>% footnote("Using Kable")
    
    
    
    ###Plots enrichment scores
    
    FGSEA_plot <- (ggbarplot(fgseaResTidy,
                             x = "Pathway",
                             y = "Normalised_Enrichment_Score",
                             fill = "Adjusted_p_value",
                             #fill = "Size",
                             #palette = "jco",
                             sort.val = "desc",
                             label = FALSE,
                             xlab = FALSE,
                             size = 1,
                             width = 1,
    )
    + scale_x_discrete(labels = function(x) str_wrap(x, 15))
    
    + scale_fill_steps(
      low = "#132B43",
      high = "#56B1F7",
      space = "Lab",
      na.value = "grey50",
      guide = "coloursteps",
      aesthetics = "fill",
      name = "P-value",
      n.breaks = 5
    )) %>%
      
      ggpar(#legend.title = list(color = "P-values"),
        rotate = TRUE,
        title = title,
        submain = paste0(control, " vs ", test),
        caption = paste0("GSEA with ", pathways, " gene sets from MSigDB"),
        ggtheme = theme_pubr(),
        #ggtheme = theme(legend.direction = "vertical"),
        xlab = "",
        ylab = "Normalized Enrichment Score",
        font.title = c("bold", "brown", 18),
        font.subtitle = c("bold", "dark grey", 14),
        font.caption = c("bold.italic", "royal blue", 12),
        font.legend = c("bold", "black", 6),
        font.x = c("bold", "black", "11"),
        #x.text.angle = 90,
        font.y = c("bold", "black", "11"),
        font.tickslab = c("bold", "black", "8"),
        font.family = "Calibri",
        
      )
    
    
    #ggsave(plots_dir, title, "_", pathways, "_FGSEA_plot.png")
    GSEA_return <- list(FGSEA_plot)
    FGSEA_plot
  }
  ##kable table of results
  
}


Annotate_genes_results <- function(res){
  
  if(species=='human'){
  res$GENENAME <- mapIds(org.Hs.eg.db,
                         keys=row.names(res),
                         column="GENENAME",
                         keytype="ENSEMBL",
                         multiVals="first")
  res$SYMBOL <- mapIds(org.Hs.eg.db,
                       keys=row.names(res),
                       column="SYMBOL",
                       keytype="ENSEMBL",
                       multiVals="first")
  
  res$ENTREZID <- mapIds(org.Hs.eg.db,
                         #EnsDb.Hsapiens.v86,
                         keys=row.names(res),
                         column="ENTREZID",
                         keytype="ENSEMBL",
                         multiVals="first")
  return(res)
  }else{
  res$GENENAME <- mapIds(org.Mm.eg.db,
                         keys=row.names(res),
                         column="GENENAME",
                         keytype="ENSEMBL",
                         multiVals="first")
  res$SYMBOL <- mapIds(org.Mm.eg.db,
                       keys=row.names(res),
                       column="SYMBOL",
                       keytype="ENSEMBL",
                       multiVals="first")
  
  res$ENTREZID <- mapIds(org.Mm.eg.db,
                         #EnsDb.Hsapiens.v86,
                         keys=row.names(res),
                         column="ENTREZID",
                         keytype="ENSEMBL",
                         multiVals="first")
  return(res)
  
  }
  
  }


gsea_rnk <- function(results_annotated, results_dir, test="test", control="control"){

gseaInput <-  results_annotated %>% 
  dplyr::filter(!is.na(ENTREZID), !is.na(log2FoldChange)) %>% 
  arrange(log2FoldChange)


gsea_input_rnk <- dplyr::select(gseaInput, ENTREZID, log2FoldChange)

write.table(gsea_input_rnk, results_dir,col.name=FALSE,sep="\t",row.names=FALSE,quote=FALSE)

}


GSEA_plots <-  function(pathways, title, results_dir, control="control", test="test"){
  
  ##read in the annotated results table
  results_annotated <- results_dir
  
  ##create ranked gene list for fgsea
  gseaInput <-  results_annotated %>% 
    dplyr::filter(!is.na(ENTREZID), !is.na(log2FoldChange)) %>% 
    arrange(log2FoldChange)

  ranks <- pull(gseaInput,log2FoldChange)
  names(ranks) <- gseaInput$SYMBOL
  
  
  ###choose pathway databases to load. 
  
  ##REACTOME GENE SETS
  HS_CP_REACTOME <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
  HS_CP_REACTOME <- split(x = HS_CP_REACTOME$gene_symbol, f = HS_CP_REACTOME$gs_name)
  
  
  ##HALLMARKS GENE SETS
  HS_HALLMARK <- msigdbr(species = "Homo sapiens", category = "H")
  HS_HALLMARK <- split(x = HS_HALLMARK$gene_symbol, f = HS_HALLMARK$gs_name)
  
  ##KEGG GENE SETS
  HS_CP_KEGG <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")
  HS_CP_KEGG <- split(x = HS_CP_KEGG$gene_symbol, f = HS_CP_KEGG$gs_name)
  
  
  ##Choose gene sets list
  #pathways <- HS_HALLMARK
  pathway_dbs <- list("HS_HALLMARK" = HS_HALLMARK, "HS_CP_KEGG" = HS_CP_KEGG, "HS_CP_REACTOME" = HS_CP_REACTOME)
  a <-c("HS_HALLMARK", "HS_CP_KEGG", "HS_CP_REACTOME")
  b <- a[str_which(a, str_to_upper(pathways))]
  pathways <- pathway_dbs[[b]]
  ##GSEA analysis using the "ranks" table against the pathways in the database.
  
  fgseaRes <- fgsea(pathways, ranks, minSize=15, maxSize = 500, nperm=1000)
  
  
  
  ## Show results in a nice table
  
  fgseaResTidy <- fgseaRes %>%
    as_tibble() %>%
    arrange(desc(NES)) %>%
    dplyr::filter(padj <= 0.05) %>%
    dplyr::mutate(Pathway = str_replace_all(pathway, "_", " "),
                  Adjusted_p_value = padj, 
                  Normalised_Enrichment_Score = NES,
                  Size = size,
                  Key_genes = leadingEdge) %>%
    dplyr::select(Pathway, Normalised_Enrichment_Score, Size, Adjusted_p_value) 
  #%>% dplyr::filter(abs(Normalised_Enrichment_Score) > 2.3)
  
  ##save results table
  #return(fgseaResTidy)
  write.table(fgseaResTidy,
              file = paste0("results_gsea/", title, "_fgsea_", b ,"_pathways.tsv"),
              col.name=TRUE,
              sep="\t",
              na = "NA",
              append = FALSE,
              row.names=FALSE,
              quote=FALSE)
  
  
  if(dim(fgseaResTidy)[1] == 0){
    return("no data")
  }
  
  ##kable table of results
  FGSEA_Results <- fgseaResTidy %>% kbl(col.names = c("Pathway", "Normalised Enrichment Score", "Size", "Adjusted p-value"),
                                        digits = 2,
                                        escape = F, 
                                        align = "c"
  ) %>% 
    kable_paper(html_font = "Comic Sans MS",
                full_width = T, bootstrap_options = c("striped", "bordered", "hover")) %>%
    column_spec(2, color = "white", background = ifelse(fgseaResTidy$Normalised_Enrichment_Score < 0, "red", "green")) %>% footnote("Using Kable")
  
  
  
  ###Plots enrichment scores
  
  FGSEA_plot <- (ggbarplot(fgseaResTidy,
                           x = "Pathway",
                           y = "Normalised_Enrichment_Score",
                           fill = "Adjusted_p_value",
                           #fill = "Size",
                           #palette = "jco",
                           sort.val = "desc",
                           label = FALSE,
                           xlab = FALSE,
                           size = 1,
                           width = 1,
  )
  + scale_x_discrete(labels = function(x) str_wrap(x, 15))
  
  + scale_fill_steps(
    low = "#132B43",
    high = "#56B1F7",
    space = "Lab",
    na.value = "grey50",
    guide = "coloursteps",
    aesthetics = "fill",
    name = "P-value",
    n.breaks = 5
  )) %>%
    
    ggpar(#legend.title = list(color = "P-values"),
      rotate = TRUE,
      title = title,
      submain = paste0(control, " vs ", test),
      caption = paste0("GSEA with ", pathways, " gene sets from MSigDB"),
      ggtheme = theme_pubr(),
      #ggtheme = theme(legend.direction = "vertical"),
      xlab = "",
      ylab = "Normalized Enrichment Score",
      font.title = c("bold", "brown", 18),
      font.subtitle = c("bold", "dark grey", 14),
      font.caption = c("bold.italic", "royal blue", 12),
      font.legend = c("bold", "black", 6),
      font.x = c("bold", "black", "11"),
      #x.text.angle = 90,
      font.y = c("bold", "black", "11"),
      font.tickslab = c("bold", "black", "8"),
      font.family = "Arial",
      
    )
  
  FGSEA_plot
  #ggsave(plots_dir, title, "_", pathways, "_FGSEA_plot.png")
  GSEA_return <- list(FGSEA_Results, FGSEA_plot)
  return(GSEA_return)
  
}
