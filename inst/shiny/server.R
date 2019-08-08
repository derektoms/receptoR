#                          _        ____
#  _ __ ___  ___ ___ _ __ | |_ ___ |  _ \
# | '__/ _ \/ __/ _ \ '_ \| __/ _ \| |_) |
# | | |  __/ (_|  __/ |_) | || (_) |  _ <
# |_|  \___|\___\___| .__/ \__\___/|_| \_\
#                   |_|
#
# August 2019 receptoR v 1.3
## Last update: 2019-08-07, Derek Toms
## server.R



########################################
#$#$#$#$#$#$    HEADER     $#$#$#$#$#$#$
########################################

# App structural packages:

#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @import reshape2
#' @import RColorBrewer
#' @import readr
#' @import stringr
#' @import shiny
#' @import shinythemes
#' @import shinyjs
#' @import dbplyr
#' @import DT
#' @import writexl

#' @import BiocManager
# Bioinformatics packages installed BiocManager:

#' @import Biobase
#' @import limma
#' @import annotate
#' @import pheatmap
#' @import mixOmics
#' @import cowplot
#' @import affy

source("functions.R")
global <- reactiveValues (DatasetTable = loadUserDatasets())
options(shiny.maxRequestSize=30*1024^2) ## 30 MB max file upload
########################################
#$#$#$#$#$#$#    SERVER    #$#$#$#$#$#$#
########################################

server <- function(input, output, session) {
    
#   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-
#  / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \
# `-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'

#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$
## This is the database search begins
#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$

# Quick link from the main page
observeEvent(input$linkSearch, {
  updateNavbarPage(session, "receptorMain", selected="searchPanel")
})


# Set up colour environment
#_,.-'~'-.,__,.-'~'-.,__,.-'~'-.,__,.-'~'-.,__,.-'~'-.,_
  catCol <- RColorBrewer::brewer.pal(3, "Set1")
  rowCol <-desat(catCol)
  userID <- NULL
 
  ## Set up reactive table to store experimental samples
  userSamples <- reactiveValues()
  userSamples$df <- data.frame()
    
  # 2019-07-31 Upload user data
  # Upload read count table
  #_,.-'~'-.,__,.-'~'-.,__,.-'~'-.,__,.-'~'-.,__,.-'~'-.,_
  userEset <- reactive({
      inFile <- input$eset_upload
      if (is.null(inFile))
          return(NULL)
      df<- read.csv(inFile$datapath,header=TRUE,sep=",")
      return(df)
  })
  
  output$upload_table <- DT::renderDataTable({
      df <- userEset()
      datatable(df, options=list(
          searching=TRUE, 
          paging=TRUE,
          scrollX=TRUE, 
          scrollY='25vh',
          scrollCollapse=TRUE,
          fixedHeader=TRUE,
          autoWidth=FALSE))
      })
  
  observeEvent(input$uploadButton, {
      showModal(modalDialog(title="Select your data to upload for analysis","Please ensure your data looks right before proceeding. A properly formatted count table should be free of missing values, and not be normalised (i.e. raw counts). Gene symbols are the required identifier and should be based on human (HGNC) and mouse (MGI) nomenclature.", radioButtons("speciesSelection", "Choose species for gene symbol annotation:", choices = c("Mouse" = "mouse", "Human" = "human")),
          fileInput('eset_upload','Choose file to upload', accept = c('text/csv','text/comma-separated-values','.csv')),
          DT::dataTableOutput("upload_table"),
          easyClose = TRUE,
          footer = tagList(
              actionButton("uploaded","Upload read table"))))
  })
  
  # 'uploaded file' flag
  eset_is_uploaded = FALSE
  
  observeEvent(input$uploaded, {
      eset_is_uploaded <<- TRUE
      removeModal()
      uploadSamples <- userEset()
      tableRows <- ncol(uploadSamples)
      userSamples$df <<- data.frame(samples = colnames(uploadSamples), category = rep("Not yet assigned", tableRows), features = rep(nrow(uploadSamples), tableRows), description = rep("User uploaded samples",tableRows))
      updateTabsetPanel(session = session, inputId = "searchpanel", selected = "2")  ## jump to 'Assign' tab
  })

  observeEvent(input$clear_upload, {
      eset_is_uploaded <<- FALSE
      userSamples$df <<- data.frame()
      removeModal()
  })

# Assign categories to each sample (GSM)
#_,.-'~'-.,__,.-'~'-.,__,.-'~'-.,__,.-'~'-.,__,.-'~'-.,_
  
  observeEvent(input$assignButton, {
        userSamples$df[,"category"] <<- as.character(userSamples$df[,"category"])
        userSamples$df[input$gsm_table_rows_selected,"category"] <<- input$selection
        userSamples$df[,"category"] <<- as.factor(userSamples$df[,"category"])
  })      
   
  output$gsm_table <- DT::renderDataTable({
      if(input$assignButton == 0){
         return (datatable(userSamples$df, extensions = 'Buttons', options=list(
               dom = 'Bfrtip',
               buttons = list(list(extend = 'colvis')),
               searching=TRUE, 
               paging=FALSE,
               scrollX=TRUE, 
               scrollY='60vh', 
               scrollCollapse=TRUE,
               fixedHeader=TRUE,
               autoWidth=TRUE,
               columnDefs=list(list(
             targets = "_all",
             render = JS(
                 "function(data, type, row, meta) {",
                     "return type === 'display' && typeof data === 'string' && data.length > 100 ?",
                     "'<span title=\"' + data + '\">' + data.substr(0, 100) + '...</span>' : data;",
                     "}")
                     )))))
      } else {
         return (datatable(userSamples$df, extensions = 'Buttons', options=list(
             dom = 'Bfrtip',
             buttons = list(list(extend = 'colvis')),
             searching=TRUE, 
             paging=FALSE,
             scrollX=TRUE, 
             scrollY='60vh', 
             scrollCollapse=TRUE,
             fixedHeader=TRUE,
             autoWidth=TRUE,
             columnDefs=list(list(
             targets = "_all",
             render = JS(
                 "function(data, type, row, meta) {",
                     "return type === 'display' && typeof data === 'string' && data.length > 100 ?",
                     "'<span title=\"' + data + '\">' + data.substr(0, 100) + '...</span>' : data;",
                     "}")
                     )))) %>%
                     formatStyle('category', target="row", backgroundColor=styleEqual(c(input$cat1, input$cat2, input$cat3), c(rowCol[1], rowCol[2], rowCol[3]))))
      }
  })
                 
  proxy.gsm = DT::dataTableProxy('gsm_table')
 
  observeEvent(input$assignButton,{
      proxy.gsm %>% selectRows(NULL)
  }) 
  
  
  ## UI output

    output$categorySelect <- renderUI(
      fluidRow(
        column(12,
               selectizeInput("selection", "Select a Category",
                           c("category1" <- {input$cat1},
                             "category2" <- {input$cat2},
                             "category3" <- {input$cat3},
                             "category4" <- "Not included")
                             # , options = list(create=TRUE, plugins = list("remove_button")))  ### <- "remove_button" isn't what I thought it was. I would also like the "create" option but I will need to link this to the table as cat1-3 are linked (otherwise new variables are not coloured or sent along for processing)
        )
      )     ### 2018-12-10 I'd like to have a button to add category 3
    )
    )  


# Finished table, to ultimately lead to CEL download
#_,.-'~'-.,__,.-'~'-.,__,.-'~'-.,__,.-'~'-.,__,.-'~'-.,_
userSamples$finishedtable <- NULL

observeEvent(input$assignButton, {
    userSamples$finishedtable <<- dplyr::filter(userSamples$df, category %in% c(input$cat1, input$cat2, input$cat3))
})
 
  output$finishedtable <- DT::renderDataTable({
      if(!is.null(userSamples$finishedtable)){
      datatable(userSamples$finishedtable,
      options=list(
          searching=FALSE, 
          paging=FALSE,
          scrollX=TRUE, 
          scrollY='60vh', 
          scrollCollapse=TRUE,
          fixedHeader=TRUE,
          autoWidth=TRUE,
          columnDefs=list(list(
          targets = "_all",
          render = JS(
              "function(data, type, row, meta) {",
                  "return type === 'display' && typeof data === 'string' && data.length > 100 ?",
                  "'<span title=\"' + data + '\">' + data.substr(0, 100) + '...</span>' : data;",
                  "}")
          )))) %>%
      formatStyle('category',target="row",
      backgroundColor=styleEqual(c(input$cat1,input$cat2,input$cat3),c(rowCol[1],rowCol[2],rowCol[3]))
  )}})

rv <- reactiveValues(download_flag = 0)


  output$report <- downloadHandler(
      filename = function(){paste(input$downloadId,"processed_data.rda",sep="_")},
      content = function(file){
          save(mylist,file = file)
          rv$download_flag <- rv$download_flag + 1
      })



# Modal confirming download, and processing function
#_,.-'~'-.,__,.-'~'-.,__,.-'~'-.,__,.-'~'-.,__,.-'~'-.,_     
observeEvent(input$downloadRDA, {
            shinyjs::disable("downloadRDA")
            mylist<<-processDataUpload(userSamples$finishedtable, userEset(), input$downloadId, input$comments, input$speciesSelection)
            showModal(modalDialog(title="Your dataset was successfully processed!","Next, save the processed data file to your local computer using the \"Download\" button. This file can then be uploaded from the menu on the \"Load Expression Datasets\" tab.",
            easyClose = TRUE,
            footer = tagList(
                actionButton("processedOK","OK"),downloadButton("report","Download"))))# modal
        })


# Reset button, modal confirmation
#_,.-'~'-.,__,.-'~'-.,__,.-'~'-.,__,.-'~'-.,__,.-'~'-.,_
  observeEvent(input$linkReset, {
      showModal(modalDialog(title="Important! Are you sure you want to reset everything?","All data and categorized samples will be lost. This can not be undone.",
      footer = tagList(
          modalButton("Cancel"),
          actionButton("buttonReset","Yes, reset."))))# modal
      # confirm reset (all categories, sample search, gone)
      observeEvent(input$buttonReset, {
          shinyjs::enable("gplSelection")
          userSamples$df <<- userSamples$df[0,]
          reset("searchText")
          reset("cat1")
          reset("cat2")
          reset("cat3")
          reset("downloadId")
          replaceData(proxy.search, NULL)
          replaceData(proxy.gsm, NULL)
          userSamples$finishedtable <<- NULL
          removeModal()
          updateTabsetPanel(session = session, inputId = "searchpanel", selected = "1")
        })

    
  })
  


#   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-
#  / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \
# `-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'

#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$
## This is where the analysis begins
#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$

# Quick link from the modal
observeEvent(input$processedOK, {
  updateNavbarPage(session, "receptorMain", selected="expressionPanel")
})

# Conditional nav tabs
hideTab(inputId = "receptorMain", target = "Gene-level Expression")
hideTab(inputId = "receptorMain", target = "Sample-level Expression")

observeEvent(input$user_data,{
    if(input$user_data!="none"){
        showTab(inputId = "receptorMain", target = "Gene-level Expression")
        showTab(inputId = "receptorMain", target = "Sample-level Expression")
    } else {
        hideTab(inputId = "receptorMain", target = "Gene-level Expression")
        hideTab(inputId = "receptorMain", target = "Sample-level Expression")

    }
})

# Load dataset
#_,.-'~'-.,__,.-'~'-.,__,.-'~'-.,__,.-'~'-.,__,.-'~'-.,_
output$loadUserExperiments = renderUI({
    selectizeInput(inputId="user_data",label="Select an experiment for analysis",choices=c("none"="none","upload processed data file (.rda)"="upload",paste(global$DatasetTable$desc)),selected="none")
})


userProcessed <- reactive({
    inFile2 <- input$processed_upload
    if (is.null(inFile2))
        return(NULL)
    load(inFile2$datapath,.GlobalEnv)
    return(mylist)
})

observeEvent(input$user_data,{
    id <- NULL
    datasetToLoad <- NULL
   
   if(input$user_data=="none"){
        assign(x = "uploaded_features", value = NULL, envir = .GlobalEnv) ## ^ this set will be used for uploaded data
        assign(x = "eset", value = NULL, envir = .GlobalEnv)
        assign(x = "de_choices", value = NULL, envir = .GlobalEnv)
        assign(x = "sig_genes_lfc", value = NULL, envir = .GlobalEnv)

    } else if (input$user_data == "upload"){
            
        showModal(modalDialog(title="Select your data to upload for visualization",
            fileInput('processed_upload','Choose file to upload', accept = c('.rda')),
            easyClose = TRUE, footer = tagList(actionButton("uploadAnalysis","Analyse"),modalButton("Cancel"))))
            
        } else {
        id <- global$DatasetTable$userID[which(global$DatasetTable$desc == input$user_data)]
        assign(x = "species", value = global$DatasetTable$species[which(global$DatasetTable$desc == input$user_data)], envir = .GlobalEnv)
        assign(x = "uploaded_features", value = NULL, envir = .GlobalEnv) ## I think this might have to be here for non-uploaded datasets
        datasetToLoad <- system.file("extdata",paste("app_data_", id, ".rda", sep=''),package="receptoR")
        withProgress(message="Loading dataset",value=0.2,{
            load(datasetToLoad,envir=.GlobalEnv)
            incProgress(0.3, message = "Loading contrasts")
            # cat(file=stderr(),"Groups are the issue:\n",categories,"\nDE choices:\n",de_choices,"\n")
            updateCheckboxGroupInput(session, "tissues", choices = categories, selected = categories)
            updateCheckboxGroupInput(session, "pls_tissues", choices = categories, selected = categories)
            updateCheckboxGroupInput(session, "de", choices = de_choices, selected = de_choices[1])

            incProgress(0.3, message ="Loading genelists")
            updateCheckboxGroupInput(session, "genelist", label = NULL, choices = names(gene_lists[[species]]), selected = NULL, inline = FALSE)
            incProgress(0.2, message = "Loading gene names")
            updateSelectInput(session, "gene", choices = make.names(uploaded_features))
        })
        
    }
    
})

observeEvent(input$uploadAnalysis, {
    
    withProgress(message="Loading dataset",value=0.2,{
        incProgress(0.3, message = "Loading contrasts")
        assign(x = "uploaded_features", value = userProcessed()$uploaded_features, envir = .GlobalEnv)
        assign(x = "eset", value = userProcessed()$eset, envir = .GlobalEnv)
        assign(x = "de_choices", value = userProcessed()$de_choices, envir = .GlobalEnv)
        assign(x = "sig_genes_lfc", value = userProcessed()$sig_genes_lfc, envir = .GlobalEnv)
        assign(x = "categories", value = userProcessed()$categories, envir = .GlobalEnv)
        assign(x = "timeStamp", value = userProcessed()$timeStamp, envir = .GlobalEnv)
        assign(x = "species", value = userProcessed()$species, envir = .GlobalEnv)
        updateCheckboxGroupInput(session, "tissues", choices = categories, selected = categories)
        updateCheckboxGroupInput(session, "pls_tissues", choices = categories, selected = categories)
        updateCheckboxGroupInput(session, "de", choices = de_choices, selected = de_choices[1])

        incProgress(0.3, message ="Loading genelists")
        updateCheckboxGroupInput(session, "genelist", label = NULL, choices = names(gene_lists[[species]]), selected = NULL, inline = FALSE)
        incProgress(0.2, message = "Loading gene names")
        updateSelectInput(session, "gene", choices = make.names(uploaded_features))
        removeModal()
    })
})

# Download DEG
#_,.-'~'-.,__,.-'~'-.,__,.-'~'-.,__,.-'~'-.,__,.-'~'-.,_

rvDEG <- reactiveValues(download_flag = 0)

  # proxy.finishedtable = dataTableProxy('finishedtable')
  output$reportDEG <- downloadHandler(
      filename = paste(input$user_data,"DEG_report.xlsx",sep="_"),
      # filename = paste(input$user_data,"DEG_report.csv",sep="_"),
      content = function(file){
          write_xlsx(sig_genes_lfc, path=file)
          rvDEG$download_flag <- rvDEG$download_flag + 1
      })
      
# Load genes tab
#_,.-'~'-.,__,.-'~'-.,__,.-'~'-.,__,.-'~'-.,__,.-'~'-.,_

  geneList = reactive({
    if (is.null(input$genelist) && is.null(input$gene)) {
      return(NULL)
    }
    
    genes = c()

    if (!is.null(input$genelist)) {
      for (gene in input$genelist) {
        genes = c(genes, gene_lists[[paste(species)]][[gene]])
      }
    }

    if (!is.null(input$gene)) {
      genes = c(genes, input$gene)
    }
    
    return(unname(genes))
  })
  
 summary_gene_data = reactive({
     validate(
       need(input$user_data!="none","No dataset selected. Please select an experiment for analysis."),
       need(geneList(), "No genes selected. Please select receptor type(s) to analyse.")
     )
   get_expression_summary(eset, geneList())
 })


# Gene outputs
#_,.-'~'-.,__,.-'~'-.,__,.-'~'-.,__,.-'~'-.,__,.-'~'-.,_
  
  output$genes = DT::renderDataTable({
      validate(
        need(input$user_data!="none","No dataset selected. Please select an experiment for analysis in 'Load Expression Data'."),
        need(geneList(), "No genes selected. Please select receptor type(s) to analyse in 'Load Expression Data'."),
        need(!is.null(summary_gene_data()), "Gene lists don't match features in your data. Check the input box below, you may have other identifiers that are not gene symbols.")
      )
    
     summary_gene_data() %>% datatable() %>% 
      formatRound(2:4)
    
  })
  
  # single gene plot
 output$singleGenePlot = renderPlot({
     validate(
       need(input$user_data!="none","No dataset selected."),
       need(geneList(), "No genes selected."),
       need(input$genes_rows_selected >= 1, "Please select one or more genes from the 'Average Expression' table to inspect expression by tissue type.")
     )
    
    rows = as.integer(input$genes_rows_selected)
    genes_to_plot = summary_gene_data()$Symbol[rows]
    
    gene_data = get_gene_data(eset, genes_to_plot)
    density <- gene_data %>% group_by(tissue,Symbol) %>% summarise(count=n())
    # cat(file=stderr(),"the density is ", as.matrix(density),"\n")
    if (any(density$count > 3)) {by_gene_violplot(gene_data,tissues=categories)}
        else {by_gene_boxplot(gene_data,tissues=categories)}
    
    
  })

# Expression tab
#_,.-'~'-.,__,.-'~'-.,__,.-'~'-.,__,.-'~'-.,__,.-'~'-.,_
  observe({
    toggle("de_choices", anim = TRUE, condition = input$de_state )
  })
  
  genesToPlot = reactive({
    validate(
      need(input$user_data!="none","No dataset selected. Please select an experiment for analysis in 'Load Expression Data'."),
      need(geneList(), "No genes selected. Please select receptor type(s) to analyse in 'Load Expression Data'.")
    )

    genes = geneList()
    
    if(input$de_state) {
      selected_de = input$de
      de_lists = lapply(selected_de, function(x) {as.character(get_de_genes(genes, x, sig_genes_lfc)$Symbol) })
      genes = Reduce(union, de_lists)
    } 
   
    return(genes) 
  }) 

  
  
# Heatmap plot
#_,.-'~'-.,__,.-'~'-.,__,.-'~'-.,__,.-'~'-.,__,.-'~'-.,_
  output$expressionPlot = renderPlot({
      validate(
          need(input$user_data!="none","No dataset selected. Please select an experiment for analysis in 'Load Expression Data'."),
          need(geneList(), "No genes selected. Please select receptor type(s) to analyse in 'Load Expression Data'."),
          need(input$tissues, "No tissues selected. Please choose at least one tissue to plot receptor heatmap."),
          need(length(genesToPlot())>10, if(input$de_state){paste("Based on the genes selected in 'Load Data', ", length(genesToPlot())," genes were differentially expressed in these comparisons (",paste(input$de, collapse = ", "), "); try choosing additional comparisons, or unselecting the 'Show differentially expressed only' option in the side menu.",sep="")}else{paste("No genes to plot as a heatmap (minimum = 10). Try including more receptor types in 'Load Data'.")})
    )
       
    selected_tissues = input$tissues
    sub_eset = eset[, eset$tissue %in% selected_tissues]
    genes = genesToPlot()
    
    # cat(file=stderr(), "Preparing heatmap:\n Tissues:", paste(input$tissues, collapse = ", "), "\n gene list: ",paste(genesToPlot(),collapse=", "),"\n genes: ", paste(genes,collapse=", "),"\n")
    
    gene_heatmap(sub_eset, genes, scale = "row",
                  probe_level = FALSE,
                  gsm_show = input$hm_gsm,
                  show_rownames = input$hm_rownames,
                  cluster_rows = input$hm_row_cluster,
                  cluster_cols = input$hm_col_cluster,
                  border_color = NA)
    })
  
  output$heatmap_ui = renderUI({
    plotOutput("expressionPlot", height = input$hm_height, width = input$hm_width)
  })

# Overall expression
#_,.-'~'-.,__,.-'~'-.,__,.-'~'-.,__,.-'~'-.,__,.-'~'-.,_
  output$overallPlot = renderPlot({
    validate(
        need(input$user_data!="none","No dataset selected. Please select an experiment for analysis in 'Load Expression Data'."),
        need(geneList(), "No genes selected. Please select receptor type(s) to analyse in 'Load Expression Data'."),
    need(input$tissues, "No tissues selected. Please choose at least one tissue to plot receptor heatmap."),
    need(length(genesToPlot())>0, if(input$de_state){paste("Based on the genes selected in 'Load Data', ", length(genesToPlot())," genes were differentially expressed in these comparisons (",paste(input$de, collapse = ", "), "); try choosing additional comparisons, or unselecting the 'Show differentially expressed only' option in the side menu.",sep="")}else{paste("No genes to plot. Try including more receptor types in 'Load Data'.")}) 
    )
    
    gene_data = get_gene_data(eset, genesToPlot())
    overall_expression_boxplot(gene_data, tissues = input$tissues)
    
  })
  

# By gene boxplots ----------------------------------------------------------------------------

  output$byGenePlot = renderPlot({
      validate(
          need(input$user_data!="none","No dataset selected. Please select an experiment for analysis in 'Load Expression Data'."),
          need(geneList(), "No genes selected. Please select receptor type(s) to analyse in 'Load Expression Data'."),
      need(input$tissues, "No tissues selected. Please choose at least one tissue to plot receptor heatmap."),
      need(length(genesToPlot())>0, if(input$de_state){paste("Based on the genes selected in 'Load Data', ", length(genesToPlot())," genes were differentially expressed in these comparisons (",paste(input$de, collapse = ", "), "); try choosing additional comparisons, or unselecting the 'Show differentially expressed only' option in the side menu.",sep="")}else{paste("No genes to plot. Try including more receptor types in 'Load Data'.")}) 
      )    
    gene_data = get_gene_data(eset, genesToPlot())
    by_gene_boxplot(gene_data, tissues = input$tissues)
  })
  


  plsdaData = reactive({
    
    selected_tissues = input$pls_tissues
    if(length(selected_tissues) < 2) {
      return(NULL)
    }
    
    
    sub_eset = eset[, eset$tissue %in% selected_tissues]
    genes = genesToPlot()
    if(is.null(uploaded_features)){genes = gene2probe(genesToPlot(), mapped_probes)}
    
    if(length(genes) < 10) {
        # cat(file=stderr(),"genes are: ",genes,"\n")
      return(NULL)
    }
    
    # probe = input$pls_probe
    
    get_plsda(sub_eset, genes, probe = FALSE) 
    
  })

# PCA plot ----------------------------------------------------------------------------
  output$indPlot = renderPlot({
    validate(
       need(input$user_data!="none","No dataset selected. Please select an experiment for analysis in 'Load Expression Data'."),
       need(geneList(), "No genes selected. Please select receptor type(s) to analyse in 'Load Expression Data'."),
       need(length(input$pls_tissues) >= 2, "Please select at least two tissues for a PLS-DA plot."),
       need(!is.null(plsdaData()$result), if(input$de_state){paste("Based on the genes selected in 'Load Data', ", length(genesToPlot())," genes were differentially expressed in these tissues (",paste(input$tissues, collapse = ", "), "); try unselecting that option in the side menu.",sep="")}else{paste("No genes to plot. Try including more receptor types in 'Load Data'.")})
    )
    
    plotIndiv(plsdaData()$result, ind.names = FALSE, group = factor(plsdaData()$tissue_grps), pch = 16, 
              col.per.group = brewer.pal(3, "Set1")[1:length(input$pls_tissues)], legend = TRUE, cex = 2, ellipse=TRUE, title="Plot of individual samples",style="graphics")
  })

# Correlation Circle plot ----------------------------------------------------------------------------  
  output$varPlot = renderPlot({
      validate(
         need(input$user_data!="none","No dataset selected. Please select an experiment for analysis in 'Load Expression Data'."),
         need(geneList(), "No genes selected. Please select receptor type(s) to analyse in 'Load Expression Data'."),
         need(length(input$pls_tissues) >= 2, "Please select at least two tissues for a Correlation circle plot."),
         need(!is.null(plsdaData()$result), if(input$de_state){paste("Based on the genes selected in 'Load Data', ", length(genesToPlot())," genes were differentially expressed in these tissues (",paste(input$tissues, collapse = ", "), "); try unselecting that option in the side menu.",sep="")}else{paste("No genes to plot. Try including more receptor types in 'Load Data'.")})
      )
      comp = as.integer(input$pls_ncomp)
    plotVar(plsdaData()$result, var.names = list(plsdaData()$varNames), comp.select=comp, cex = 1, overlap=FALSE, col="grey",title="Correlation circle between genes and discriminant components", style="graphics")
    
  })

  output$numGenesUI = renderUI({
    numericInput("pls_num_genes", "Select number of genes to show contributions for", 
                 value = 25, min = 1, max = length(geneList()), step = 1)
  })
  
# Loadings plot ----------------------------------------------------------------------------
  output$contribPlot = renderPlot({
      validate(
         need(input$user_data!="none","No dataset selected. Please select an experiment for analysis in 'Load Expression Data'."),
         need(geneList(), "No genes selected. Please select receptor type(s) to analyse in 'Load Expression Data'."),
         need(length(input$pls_tissues) >= 2, "Please select at least two tissues for a Loadings plot."),
         need(!is.null(plsdaData()$result), if(input$de_state){paste("Based on the genes selected in 'Load Data', ", length(genesToPlot())," genes were differentially expressed in these tissues (",paste(input$tissues, collapse = ", "), "); try unselecting that option in the side menu.",sep="")}else{paste("No genes to plot. Try including more receptor types in 'Load Data'.")})
      )

    
    grps = plsdaData()$result$names$colnames$Y
    ndisplay = input$pls_num_genes
    comp = as.integer(input$pls_ncomp)
    plotLoadings(plsdaData()$result, name.var = plsdaData()$varNames, ndisplay = ndisplay, comp = comp, contrib='max', method='mean',legend.color = catCol[1:length(grps)],title=paste("Weight of the top ", ndisplay, " genes contributing to discriminant component ", comp, sep=""),size.title=1)
     
  })
  
#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$  
  ## Kill shinyApp when session closes
  session$onSessionEnded(stopApp)

}
