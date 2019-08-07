#                          _        ____
#  _ __ ___  ___ ___ _ __ | |_ ___ |  _ \
# | '__/ _ \/ __/ _ \ '_ \| __/ _ \| |_) |
# | | |  __/ (_|  __/ |_) | || (_) |  _ <
# |_|  \___|\___\___| .__/ \__\___/|_| \_\
#                   |_|
#
# August 2019 receptoR v 1.3
## Last update: 2019-08-07, Derek Toms
## functions.R


#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$
desat = function(cols, sat=0.5) {
    X <- diag(c(1, sat, 1)) %*% rgb2hsv(col2rgb(cols))
    hsv(X[1,], X[2,], X[3,])
}

#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$
#$# Data processing 31 July 2018
##  Uploaded count table processing, similar to above but without downloading CEL files
processDataUpload = function(finished_table, read_counts, datasetID, userComments, gpl){
    ## categories need to have R-safe names so we'll store the user input as a new column
    finished_table$category.labels <- finished_table$category
    safeNames <- make.names(levels(factor(finished_table$category)),unique=TRUE)
    finished_table$category <- factor(finished_table$category, levels = levels(factor(finished_table$category)), labels=safeNames)
    ## colours
    catCol <- factor(brewer.pal(9,"Set1")[1:length(levels(factor(finished_table$category)))])
    finished_table$colours <- finished_table$category
    finished_table$colours <- factor(finished_table$colours, levels=levels(factor(finished_table$category)), labels=catCol)
    
    ## timestamp
    timeStamp <- strftime(Sys.time(),"%Y%m%d-%H%M")
          
withProgress(
   message = "Processing user count data", value = 0,
   {
# Convert ------------------------------------------------------------
        incProgress(0.1, message = "Converting to expression set")
        sample_names <- as.character(finished_table$samples)
        final_counts <- tbl_df(read_counts) %>% dplyr::select(sample_names)
        uploaded_features <- as.character(read_counts[,1])
        matrix_counts <- as.matrix(final_counts)
        rownames(matrix_counts) <- make.names(uploaded_features, unique=TRUE)   ##### This is currently a problem with non-unique row names: solved, but not the best method
        all_eset = ExpressionSet(matrix_counts)
        all_pData<-data.frame(tissue=finished_table$category)
        rownames(all_pData)<-finished_table$samples 
        
        ## Change of object name
        all_eset_final<-all_eset
        pData(all_eset_final)<-all_pData

# Differentially Expressed Gene (DEG) Analysis ------------------------------------------------------------
        incProgress(0.1, message = "Determining differential gene expression")

        eset = all_eset_final

        tissue = as.factor(pData(eset)$tissue)
        design = model.matrix(~0 + tissue)
        colnames(design) = levels(tissue)
        
        incProgress(0.1, message = "Fitting model")
        fit = lmFit(eset, design)
        matrices<-t(combn(levels(tissue),2))
        contrasts <- paste(matrices[,1],matrices[,2],sep='-')
        contrast_matrix = makeContrasts(contrasts=contrasts, levels = design)

        fit2 = contrasts.fit(fit, contrast_matrix)
        efit = eBayes(fit2)
        
        tfit = treat(fit2, lfc = 1)

        results = decideTests(efit)
        results_lfc = decideTests(tfit)

        coefs = colnames(contrast_matrix)
        sig_genes = lapply(coefs, function(x) {
            topTable(efit, coef = x, number = Inf, p.value = 0.01, sort.by = "p")
        })

        sig_genes_lfc = lapply(coefs, function(x) {
            topTreat(tfit, coef = x, number = Inf, p.value = 0.01, sort.by = "p")
        })
        
        names(sig_genes) = coefs
        names(sig_genes_lfc) = coefs

        sapply(sig_genes, nrow)
        sapply(sig_genes_lfc, nrow)

        # get list of DEG
        de_choices = names(sig_genes_lfc)
        # set groups
        categories = levels(tissue)

# Save user-generated experiments -----------------------------------   
        incProgress(0.1, message = "Saving processed data")
        mylist = list("uploaded_features"=uploaded_features, "eset"=eset, "de_choices"=de_choices, "sig_genes_lfc"=sig_genes_lfc, "categories"=categories,"timeStamp"=timeStamp,"species"=gpl)
        incProgress(0.1, message = "Success!")
    })
    return(mylist)
}
loadUserDatasets <- function() {
    data <- read.csv(system.file("extdata","userDB.csv",package="receptoR"))
    return(data)
 }

 # end processing
#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$

 # Get de genes from topTable output that match a genelist -------------------------------------

 get_de_genes = function(gene_list, de_choice, sig_genes_lfc) {
       sig_genes_lfc[[de_choice]] %>% tibble::rownames_to_column("Symbol") %>% filter(Symbol %in% gene_list)
   }

 # Get summary expression table ----------------------------------------------------------------

get_expression_summary = function(eset, gene_list) {
   
    geneData <- get_gene_data(eset, gene_list)
    # cat(file=stderr(),"Is the data null?\n",!is.null(geneData),"\n")
    if(!is.null(geneData)){
    geneData <- geneData %>% 
    group_by(Symbol, tissue) %>% 
    summarise(expression = mean(expression)) %>% 
    spread(tissue, expression)
    return(geneData)
    } else {
        return(NULL)
    }
  
 }

gene_heatmap = function(eset, subset_probes, anno_col = 'tissue', probe_level = FALSE, gsm_show = TRUE, ...) {

    # cat(file=stderr(), "incoming ", length(subset_probes), " probes\n")
    mat = exprs(eset) %>% as.data.frame() %>% tibble::rownames_to_column("Symbol") %>% filter(Symbol %in% subset_probes)
    row_labs=mat$Symbol
    mat = data.matrix(mat[,-1])
    rownames(mat) = row_labs
    mat = mat[!(apply(mat,1,function(y) all(y==0))),]
    row_labs = rownames(mat)
    mat = aggregate(mat, list(genes = rownames(mat)), mean)
    rownames(mat) = mat$genes
    mat = as.matrix(mat[,-1])
    row_labs=rownames(mat)
    anno_df = pData(eset)[anno_col]
    # another fix for Tibble deprecating row names
    anno_df <- data.frame(anno_df)
    rownames(anno_df) <- colnames(mat)
    # end fix
    vars = as.character(unique(anno_df[[anno_col]]))
    anno_colours = brewer.pal(9, "Set1")[1:length(vars)]
    names(anno_colours) = vars
    anno_colours = list(anno_colours)
    names(anno_colours) = anno_col 
    pheatmap(mat, show_colnames = gsm_show, annotation_col = anno_df, annotation_colors = anno_colours, labels_row = row_labs,...)
  
 }


 # Get expression data from eset for given genes -----------------------------------------------

 get_gene_data = function(eset, gene_list) {
   ph = pData(eset) %>% tibble::rownames_to_column("Sample")
   ph$Sample <- colnames(exprs(eset)) # 2018-04-17 fix missing row names
       keepCols <- c("Symbol")
       geneData <- exprs(eset) %>% as.data.frame() %>% tibble::rownames_to_column("Symbol") %>% filter(Symbol %in% gene_list) %>% gather(Sample, expression, -one_of(keepCols)) %>% left_join(ph, by = "Sample")
       # cat(file=stderr(),"%%%%%% gene list (",species,"):\n",gene_list,"\n%%%%%% length of expression data frame:\n", length(rownames(geneData)),"\n")
       if (length(rownames(geneData))>2){
       return(geneData)} else {return(NULL)}
 }

 # Create a top N gene list by abs LFC ---------------------------------------------------------

 top_genes = function(de_data, top = 25) {
   de_data %>% 
     group_by(Symbol) %>% 
     summarise(meanLFC = mean(logFC)) %>% 
     top_n(top, abs(meanLFC))
 }


 # Make a boxplot faceted by gene --------------------------------------------------------------

 by_gene_boxplot = function(gene_data, tissues = categories) {
   gene_data %>% 
     filter(tissue %in% tissues) %>%
     ggplot(aes(x = tissue, y = expression)) +
       geom_boxplot(aes(fill = tissue)) +
       facet_wrap(~Symbol) +
       theme_bw() + theme(axis.text.x = element_blank()) +
       scale_fill_brewer(palette = "Set1")
 }
 
 by_gene_violplot = function(gene_data, tissues = categories) {
    gene_data %>% 
      filter(tissue %in% tissues) %>%
      ggplot(aes(x = tissue, y = expression)) +
        geom_violin(aes(fill = tissue)) +
        facet_wrap(~Symbol) +
        theme_bw() + theme(axis.text.x = element_blank()) +
        scale_fill_brewer(palette = "Set1")
  }

 # Make a box plot for the expression of given genes, one plot per tissue -----------------------

 overall_expression_boxplot = function(gene_data, tissues = categories) {

   tissue_cols = brewer.pal(3, "Set1")[1:length(tissues)]
   names(tissue_cols) = tissues

   gene_plots = gene_data %>% 
     filter(tissue %in% tissues) %>% 
     group_by(tissue, Symbol) %>% mutate(meanExp = mean(expression)) %>% 
     group_by(tissue) %>% 
     do(plot = mutate(., Symbol = factor(Symbol, levels = unique(Symbol[order(meanExp, decreasing = TRUE)]))) %>% 
       ggplot(aes(x = Symbol, y = expression)) +
         geom_boxplot(fill = tissue_cols[unique(.$tissue)]) +
         theme_bw() +
         theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
         labs(title = unique(.$tissue), x = "")
     )
   plot_grid(plotlist = gene_plots$plot, ncol = 1)    
    
 }


 # Perform sPLS-DA on selected genes -----------------------------------------------------------

 get_plsda = function(eset, genes, probe) {
     
         exp = exprs(eset) %>% as.data.frame() %>% tibble::rownames_to_column("Symbol") %>% filter(Symbol %in% genes)
         row_labs=exp$Symbol
         exp = data.matrix(exp[,-1])
         rownames(exp) = row_labs
         exp = exp[!(apply(exp,1,function(y) all(y==0))),]

     
     tissue = factor(pData(eset)$tissue)
     tissue_grps = pData(eset)$tissue
    
     if (probe == FALSE) {
       exp = aggregate(exp, list(genes = rownames(exp)), mean)
       rownames(exp) = exp$genes
       exp = as.matrix(exp[,-1])
     }
    
     exp = t(exp)
         return(list(
           result = splsda(exp, tissue, ncomp = 2),
           tissue_grps = tissue_grps,
           varNames = colnames(exp)
         ))
 } 