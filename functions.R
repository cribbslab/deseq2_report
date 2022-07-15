

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
  }else if(species == 'mouse'){
    annots <-  AnnotationDbi::select(org.Mm.eg.db, keys=data,
                                     columns="SYMBOL", keytype = "ENSEMBL")
  }else if(species == 'macaque'){
    annots <-  AnnotationDbi::select(org.Mmu.eg.db, keys=data,
                                     columns="SYMBOL", keytype = "ENSEMBL")
  }else if(species == 'rabbit'){
    ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL",
                      dataset="ocuniculus_gene_ensembl", 
                      host="uswest.ensembl.org")
    annots <- as.data.frame(getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'), mart=ensembl))
    colnames(annots) <- c('ENSEMBL','SYMBOL')
  }else if(species == 'pig'){
    ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL",
                      dataset="sscrofa_gene_ensembl", 
                      host="uswest.ensembl.org")
    annots <- as.data.frame(getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'), mart=ensembl))
    colnames(annots) <- c('ENSEMBL','SYMBOL')
  }else{
    print('please supply a valid species within the config.yml file')
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
  
  res_names <- resultsNames(dds)
  
  resLFC <- lfcShrink(dds, coef=res_names[2], type="apeglm")
  
  return(new("DESeq2_return",
             res=res,
             dds=dds))
}

setClass(Class="DESeq2_lrt_return",
         representation(
           resLFC="DESeqResults",
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
  }else if(species == 'mouse'){
    annots <-  AnnotationDbi::select(org.Mm.eg.db, keys=data,
                                     columns="SYMBOL", keytype = "ENSEMBL")
  }else if(species == 'macaque'){
    annots <-  AnnotationDbi::select(org.Mmu.eg.db, keys=data,
                                     columns="SYMBOL", keytype = "ENSEMBL")
  }else if(species == 'rabbit'){
    ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL",
                      dataset="ocuniculus_gene_ensembl", 
                      host="uswest.ensembl.org")
    annots <- as.data.frame(getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'), mart=ensembl))
    colnames(annots) <- c('ENSEMBL','SYMBOL')
  }else if(species == 'pig'){
    ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL",
                      dataset="sscrofa_gene_ensembl", 
                      host="uswest.ensembl.org")
    annots <- as.data.frame(getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'), mart=ensembl))
    colnames(annots) <- c('ENSEMBL','SYMBOL')
  }else{
    print('please supply a valid species within the config.yml file')
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


Annotate_genes_results <- function(res, species='human'){
  
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
  }else if(species=='mouse'){
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
  
  }else if(species=='macaque'){
    res$GENENAME <- mapIds(org.Mmu.eg.db,
                           keys=row.names(res),
                           column="GENENAME",
                           keytype="ENSEMBL",
                           multiVals="first")
    res$SYMBOL <- mapIds(org.Mmu.eg.db,
                         keys=row.names(res),
                         column="SYMBOL",
                         keytype="ENSEMBL",
                         multiVals="first")
    
    res$ENTREZID <- mapIds(org.Mmu.eg.db,
                           #EnsDb.Hsapiens.v86,
                           keys=row.names(res),
                           column="ENTREZID",
                           keytype="ENSEMBL",
                           multiVals="first")
    return(res)
    
  }else if(species=='pig'){

    res$SYMBOL <- mapIds(org.Ss.eg.db,
                         keys=row.names(res),
                         column="SYMBOL",
                         keytype="ENSEMBL",
                         multiVals="first")
    
    res$ENTREZID <- mapIds(org.Mmu.eg.db,
                           #EnsDb.Hsapiens.v86,
                           keys=row.names(res),
                           column="ENTREZID",
                           keytype="ENSEMBL",
                           multiVals="first")
    return(res)
    
  }else if(species=='rabbit'){
    res$GENENAME <- mapIds(org.Mmu.eg.db,
                           keys=row.names(res),
                           column="GENENAME",
                           keytype="ENSEMBL",
                           multiVals="first")
    res$SYMBOL <- mapIds(org.Mmu.eg.db,
                         keys=row.names(res),
                         column="SYMBOL",
                         keytype="ENSEMBL",
                         multiVals="first")
    
    res$ENTREZID <- mapIds(org.Mmu.eg.db,
                           #EnsDb.Hsapiens.v86,
                           keys=row.names(res),
                           column="ENTREZID",
                           keytype="ENSEMBL",
                           multiVals="first")
    return(res)
    
  }else{
    print('please supply a valid species within the config.yml file')
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



' Function to visualise enrichment results using a barplot
#'
#' \code{xEnrichBarplot} is supposed to visualise enrichment results using a barplot. It returns an object of class "ggplot".
#'
#' @param eTerm an object of class "eTerm"
#' @param top_num the number of the top terms (sorted according to FDR or adjusted p-values). If it is 'auto', only the significant terms (see below FDR.cutoff) will be displayed
#' @param displayBy which statistics will be used for displaying. It can be "fc" for enrichment fold change (by default), "adjp" or "fdr" for adjusted p value (or FDR), "pvalue" for p value, "zscore" for enrichment z-score
#' @param FDR.cutoff FDR cutoff used to declare the significant terms. By default, it is set to 0.05. This option only works when setting top_num (see above) is 'auto'
#' @param bar.label logical to indicate whether to label each bar with FDR. By default, it sets to true for bar labelling
#' @param bar.label.size an integer specifying the bar labelling text size. By default, it sets to 3
#' @param bar.color either NULL or fill color names ('lightyellow-orange' by default)
#' @param bar.width bar width. By default, 80% of the resolution of the data
#' @param wrap.width a positive integer specifying wrap width of name
#' @param font.family the font family for texts
#' @param signature logical to indicate whether the signature is assigned to the plot caption. By default, it sets TRUE showing which function is used to draw this graph
#' @return an object of class "ggplot"
#' @note none
#' @export
#' @seealso \code{\link{xEnricherGenes}}, \code{\link{xEnricherSNPs}}, \code{\link{xEnrichViewer}}
#' @include xEnrichBarplot.r
#' @examples
#' \dontrun{
#' # Load the XGR package and specify the location of built-in data
#' library(XGR)
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata/"
#' 
#' # 1) load eQTL mapping results: cis-eQTLs significantly induced by IFN
#' cis <- xRDataLoader(RData.customised='JKscience_TS2A', RData.location=RData.location)
#' ind <- which(cis$IFN_t > 0 & cis$IFN_fdr < 0.05)
#' df_cis <- cis[ind, c('variant','Symbol','IFN_t','IFN_fdr')]
#' data <- df_cis$variant
#' 
#' # 2) Enrichment analysis using Experimental Factor Ontology (EFO)
#' # Considering LD SNPs and respecting ontology tree
#' eTerm <- xEnricherSNPs(data, ontology="EF", include.LD="EUR", LD.r2=0.8, ontology.algorithm="lea", RData.location=RData.location)
#'
#' # 3) Barplot of enrichment results
#' bp <- xEnrichBarplot(eTerm, top_num="auto", displayBy="fc")
#' #pdf(file="enrichment_barplot.pdf", height=6, width=12, compress=TRUE)
#' print(bp)
#' #dev.off()
#' }
#' 
#' # 4) use font family (Arial)
#' \dontrun{
#' BiocManager::install("extrafont")
#' library(extrafont)
#' font_import()
#' fonttable()
#' ## creating PDF files with fonts
#' library(extrafont)
#' loadfonts()
#' bp <- xEnrichBarplot(eTerm, top_num="auto", displayBy="fc", font.family="Arial Black")
#' pdf(file="enrichment_barplot_fonts.pdf", height=6, width=12, family="Arial Black")
#' print(bp)
#' dev.off()
#' }

xEnrichBarplot <- function(eTerm, top_num=10, displayBy=c("fc","adjp","fdr","zscore","pvalue"), FDR.cutoff=0.05, bar.label=TRUE, bar.label.size=3, bar.color='lightyellow-orange', bar.width=0.8, wrap.width=NULL, font.family="sans", signature=TRUE) 
{
  
  displayBy <- match.arg(displayBy)
  
  if(is.null(eTerm)){
    warnings("There is no enrichment in the 'eTerm' object.\n")
    return(NULL)
  }
  
  ## when 'auto', will keep the significant terms
  df <- xEnrichViewer(eTerm, top_num="all")
  if(top_num=='auto'){
    top_num <- sum(df$adjp<FDR.cutoff)
    if(top_num <= 1){
      top_num <- 10
    }
  }
  df <- xEnrichViewer(eTerm, top_num=top_num, sortBy="adjp")
  
  ## text wrap
  if(!is.null(wrap.width)){
    width <- as.integer(wrap.width)
    res_list <- lapply(df$name, function(x){
      x <- gsub('_', ' ', x)
      y <- strwrap(x, width=width)
      if(length(y)>1){
        paste0(y[1], '...')
      }else{
        y
      }
    })
    df$name <- unlist(res_list)
  }
  
  name <- height <- direction <- hjust <- NULL
  adjp <- zscore <- pvalue <- fc <- NULL
  
  ##########
  ## consider the direction of z-score
  df <- df %>% dplyr::mutate(direction=ifelse(zscore>0,1,-1))
  ##########	
  
  if(displayBy=='adjp' | displayBy=='fdr'){
    df <- df %>% dplyr::arrange(direction, desc(adjp), zscore) %>% dplyr::mutate(height=-1*log10(adjp)) %>% dplyr::mutate(hjust=1)
    df$name <- factor(df$name, levels=df$name)
    ####
    if(length(df$height[!is.infinite(df$height)])==0){
      df$height <- 10
    }else{
      df$height[is.infinite(df$height)] <- max(df$height[!is.infinite(df$height)])
    }
    ####
    p <- ggplot(df, aes(x=name, y=height))
    p <- p + ylab(expression(paste("Enrichment significance: ", -log[10]("FDR"))))
    
  }else if(displayBy=='fc'){
    df <- df %>% dplyr::arrange(direction, fc, desc(adjp)) %>% dplyr::mutate(height=log2(fc)) %>% dplyr::mutate(hjust=ifelse(height>=0,1,0))
    df$name <- factor(df$name, levels=df$name)
    p <- ggplot(df, aes(x=name, y=height))
    p <- p + ylab(expression(paste("Enrichment changes: ", log[2]("FC"))))
    
  }else if(displayBy=='pvalue'){
    df <- df %>% dplyr::arrange(direction, desc(pvalue), zscore) %>%  dplyr::mutate(height=-1*log10(pvalue)) %>% dplyr::mutate(hjust=1)
    df$name <- factor(df$name, levels=df$name)
    ####
    if(length(df$height[!is.infinite(df$height)])==0){
      df$height <- 10
    }else{
      df$height[is.infinite(df$height)] <- max(df$height[!is.infinite(df$height)])
    }
    ####
    p <- ggplot(df, aes(x=name, y=height))
    p <- p + ylab(expression(paste("Enrichment significance: ", -log[10]("p-value"))))
    
  }else if(displayBy=='zscore'){
    df <- df %>% dplyr::arrange(direction, zscore, desc(adjp)) %>% dplyr::mutate(height=zscore) %>% dplyr::mutate(hjust=ifelse(height>=0,1,0))
    df$name <- factor(df$name, levels=df$name)
    p <- ggplot(df, aes(x=name, y=height))
    p <- p + ylab("Enrichment z-scores")
    
  }
  
  if(is.null(bar.color)){
    bp <- p + geom_col(color='grey80',fill='transparent', width=bar.width)
  }else{
    bar.color <- unlist(strsplit(bar.color, "-"))
    if(length(bar.color)==2){
      bp <- p + geom_col(aes(fill=height), width=bar.width) + scale_fill_gradient(low=bar.color[1],high=bar.color[2]) 
    }else{
      bp <- p + geom_col(color='grey80',fill='transparent', width=bar.width)
    }
  }
  #bp <- p + geom_col(aes(fill=height)) + scale_fill_gradient2(low="cyan", mid="grey", high="orange", midpoint=0)
  bp <- bp + theme_bw() + theme(legend.position="none",axis.title.y=element_blank(), axis.text.y=element_text(size=10,color="black"), axis.title.x=element_text(size=12,color="black")) + coord_flip()
  
  bp <- bp + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  if(bar.label){
    ## get text label
    to_scientific_notation <- function(x) {
      res <- format(x, scientific=T)
      res <- sub("\\+0?", "", res)
      sub("-0?", "-", res)
    }
    label <- to_scientific_notation(df$adjp)
    df$label <- paste('FDR', as.character(label), sep='=')
    
    ## hjust==1
    df_sub <- df %>% dplyr::filter(hjust==1)
    bp <- bp + geom_text(data=df_sub, aes(label=label),hjust=1,size=bar.label.size,family=font.family)
    ## hjust==0
    df_sub <- df %>% dplyr::filter(hjust==0)
    bp <- bp + geom_text(data=df_sub, aes(label=label),hjust=0,size=bar.label.size,family=font.family)
  }
  
  ## caption
  if(signature){
    caption <- paste("Created by xEnrichBarplot from XGR version", utils ::packageVersion("XGR"))
    bp <- bp + labs(caption=caption) + theme(plot.caption=element_text(hjust=1,face='bold.italic',size=8,colour='#002147'))
  }
  
  ## change font family to 'Arial'
  bp <- bp + theme(text=element_text(family=font.family))
  
  ## put arrows on x-axis
  bp <- bp + theme(axis.line.x=element_line(arrow=arrow(angle=30,length=unit(0.25,"cm"), type="open")))
  
  ## x-axis (actually y-axis) position
  bp <- bp + scale_y_continuous(position="top")
  
  return(bp)
}
