#                          _        ____
#  _ __ ___  ___ ___ _ __ | |_ ___ |  _ \
# | '__/ _ \/ __/ _ \ '_ \| __/ _ \| |_) |
# | | |  __/ (_|  __/ |_) | || (_) |  _ <
# |_|  \___|\___\___| .__/ \__\___/|_| \_\
#                   |_|
#
# August 2019 receptoR v 1.3
## Last update: 2019-08-07, Derek Toms
## ui.R

########################################
#$#$#$#$#$#$    HEADER     $#$#$#$#$#$#$
########################################

library(shiny)
library(shinythemes)
library(shinyjs)


########################################
#$#$#$#$#$#$  Javascript   $#$#$#$#$#$#$
########################################


# this codes for the 'enter/return' key as an action button
jscode <- '
$(function() {
  var $els = $("[data-proxy-click]");
    $.each(
    $els,
    function(idx, el) {
      var $el = $(el);
      var $proxy = $("#" + $el.data("proxyClick"));
      $el.keydown(function (e) {
        if (e.keyCode == 13) {
          $proxy.click();
        }
      });
    }
  );
});
'

jscode2 <- '
$(function() {
    $(#receptorMain li a[data-value="Gene-level Expression"]).hide();
    $(#receptorMain li a[data-value="Sample-level Expression"]).hide();
    
});
'


########################################
#$#$#$#$#$#$#$    UI     $#$#$#$#$#$#$#$
########################################

ui <- fluidPage(
tags$head(tags$script(HTML(jscode))),
tags$head(tags$script(HTML(jscode2))),
tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "www/receptor.CSS")),
tags$head(tags$link(rel = "stylesheet", href = "https://use.fontawesome.com/releases/v5.6.3/css/all.css",  integrity="sha384-UHRtZLI+pbxtHCWp1t77Bi1L4ZtiqrqD80Kn4Z8NTSRyMA2Fd33n5dQ8lWUE00s/", crossorigin="anonymous")),
shinyjs::useShinyjs(),
navbarPage("receptoR",
    id = "receptorMain",
    theme = shinytheme("spacelab"),

# Start page  ------------------------------------------------------------------------------

    tabPanel("Start here",
       value ="startPanel",
       h3("Welcome to receptoR!"),
       hr(),
       sidebarLayout(
           sidebarPanel(
               # h4("An automated hypothesis generation software to identify cellular signaling pathways from transcriptomics data"),
               p("This software allows you to browse and analyze transcriptomics data. This is based on the idea that each cell type expresses a particular suite of cellular receptors that drive its behaviour."),
               tags$ol(tags$li("A cell transcribes mRNA that will be translated into functional receptor proteins."),tags$li("Isolated RNA from this cell is converted to labeled cDNA, which can be analysed by microarray or next-generation sequencing to measure expression of specific transcripts."),tags$li("Each array represents a snapshot of a specific transcriptome. Thousands of these have been digitized and made publicly available."),tags$li("By mining this data, we can predict which receptors are expressed by our cells or tissues of interest to direct bioengineering strategies.")),
               hr(),
               
               #div
               p("To begin using receptor, ",
               actionLink("linkSearch","upload expression data")),
               hr(),
               p("code created by Derek Toms, Qing Yun Tong and Matthew Workentine"),
               p("Copyright (C) 2019, code licensed under GPLv3")
               #/div
               ),
           mainPanel(
               img(src="www/overview.png",width="100%")
               )
               
        )),

#  Upload data  ------------------------------------------------------------------------------

    tabPanel("Upload Transcriptome Data",
       value = "searchPanel",
       includeCSS("www/receptor.CSS"),
       h3("Upload and categorize expression data"), ## change
       hr(),
       sidebarLayout(
       sidebarPanel(
           # style = "position:fixed;width:30%",
           conditionalPanel(condition="input.searchpanel==1",
           h4("Upload Expression Data"),
           br(),
           actionButton("uploadButton", "Upload data to analyze")),
           
           conditionalPanel(condition="input.searchpanel==2",
           h4("Define the categories to assign each sample for comparison."),
           p("Each sample of interest should be assigned to a category. In this way, experimental comparisons can be performed to determine differential expression between categories. A minimum of two and a maximum of three categories should be defined. If you are only interested in a single sample type it is recommended that this is compared to a 'background' sample to identify enriched receptor genes."),

           tags$div(class="inputWithIcon", textInput("cat1", label=NULL, placeholder="Category 1 (e.g., pancreatic endocrine cells)"), tags$span(style="color:#E41A1C",icon("circle",class="fa-2x"))),

           tags$div(class="inputWithIcon",textInput("cat2", label=NULL, placeholder="Category 2 (e.g., photoreceptors)"), tags$span(style="color:#377EB8",icon("circle",class="fa-2x"))),
           
           tags$div(class="inputWithIcon",textInput("cat3", label=NULL, placeholder="Category 3 (optional)"),
           tags$span(style="color:#4DAF4A",icon("circle",class="fa-2x"))),
                      
           hr(),
           h4("Highlight samples, then click to Assign them to the specificed category."),
           p("Using the table at right and the drop down menu below, click on samples and \'Assign\' them to different categories. Samples can be filtered using the search bar."),
           fluidRow(column(8,uiOutput("categorySelect")),
           column(4,actionButton("assignButton", "Assign")))
           ),
           
           conditionalPanel(condition="input.searchpanel==3",
               h4("Thank you for using receptoR!"), # this should be slightly more informative, along the lines of Process data
               p(" Please enter your name and any comments/bugs/questions/requests in the box below, then click the \'Download and Process\' button to retrieve the raw files from the NCBI server and process them based on their assigned categories."), # move this off of this page; separate "Help/Comments" button
               textAreaInput("comments","Comments",width="100%",height="100px",resize="vertical"), # ditto
               textInput("downloadId","Download ID"), # needs to be mandatory (or filled by default)
               actionButton("downloadRDA","Process")),
           hr(),
               # Help banner on the bottom -------------------------
           # h4("Help me!"),
           p("Click ",actionLink("linkReset","here "),"to start again.")
       ),
       
       mainPanel(
        tabsetPanel(id = "searchpanel",
        tabPanel("Upload", value=1,
            h4("Transcriptome data"), # return uploaded table
            p("Begin by uploading transcriptome data as a read count table. We recommend reads from at least eight samples to improve statistical power. In the next step, these samples will be assigned to one of three categories to determine differential expression between sample types.")
        ),
        # Assign samples to categories ------------------------------------------------------
        tabPanel("Assign", value=2,
            h4("Assign individual samples (columns) to categories of your choosing"),
            DT::dataTableOutput("gsm_table")
        ),
        # This will be where the CEL files are downloaded (confirmation, etc) ------------
        tabPanel("Process", value=3,
        h4("Please confirm samples are properly categorized before proceeding"),
        p("Expression samples annotated:"), DT::dataTableOutput("finishedtable")        
        ))
        )
        
    )),
    
    # Load Gene Expression Data tab -------------------------------------
    tabPanel("Load Expression Datasets",
        value="expressionPanel",
        h3("Pick from previously analyzed experiments to perform analyses"),
        hr(),
        sidebarLayout(
        sidebarPanel(
            h4("Load Experiment"),
            uiOutput("loadUserExperiments"),
            hr(),
            checkboxGroupInput("genelist", "Select a receptor type to analyze", 
                  choices = NULL),
            br(),
            selectInput("gene", "Select additional non-receptor coding gene(s) to include in the analysis.", choices = NULL, multiple = TRUE),
            helpText("Search by gene symbol; availability of a given gene is based on microarray probe annotations."),
            downloadButton("reportDEG","Download differential gene expression analysis"),
            helpText("This Microsoft Excel file (.XSLX) contains all differentially expressed genes among these tissues. It can be further used for downstream analyses including functional enrichment analysis.")
        ),
        mainPanel(
           fluidRow(                
               column(6, h4("Average Gene-by-gene Expression"), DT::dataTableOutput("genes")),            
               column(6, h4("Gene Violin Plot"), plotOutput("singleGenePlot")))
        )
        )
    ),
    
    # Magnitude expression tab ------------------------------------------------------------------------------
    
    tabPanel("Gene-level Expression",
        h3("Compare genes based on absolute expression and differential expression between experimental groups"),
        hr(),
        sidebarLayout(
        sidebarPanel(
            h4("Gene expression"),
            checkboxGroupInput("tissues", label = "Select tissues to include",
            choices = NULL, selected = NULL),
            br(),
            checkboxInput("de_state", label = "Show differentially expressed only", value = TRUE),
            checkboxGroupInput("de", label = "Choose comparison(s) to show", choices = NULL, selected = NULL),
            br(),
            conditionalPanel(condition="input.absexpanel==1",
                h5("Heatmap parameters"),
                checkboxInput("hm_gsm", "Show GSM (column names)", value = TRUE),
                checkboxInput("hm_rownames", "Show gene symbols (row names)", value = TRUE),
                checkboxInput("hm_col_cluster", "Cluster columns", value = TRUE),
                checkboxInput("hm_row_cluster", "Cluster rows", value = TRUE),
                numericInput("hm_width", "Plot width (px)", value = 900, min = 100, max = 2400, step = 10),
                numericInput("hm_height", "Plot height (px)", value = 1200, min = 100, max = 2400, step = 10))
        ),
        mainPanel(
            tabsetPanel(type = "tabs",
            tabPanel("Heatmap", value=1, h4("Cluster analysis and a heatmap representation of gene expression."), p("Genes with similar expression patterns will cluster as rows, while individual microarrays cluster as columns."), uiOutput("heatmap_ui")),
            tabPanel("Summary boxplots", h4("Boxplots of expression data by tissue."), plotOutput("overallPlot", height = 600)),
            tabPanel("By-gene boxplots", h4("Boxplots of expression data by gene."), plotOutput("byGenePlot", height = 600)),
            id = "absexpanel")
        )
        )
    ),

    # Mixomics tab ---------------------------------------------
    tabPanel("Sample-level Expression",
        h3("Compare trends in samples based on experimental groups"),
        hr(),
        sidebarLayout(
        sidebarPanel(
            h4("Sample expression"),
            checkboxGroupInput("pls_tissues", label = "Select tissues to inclued",
            choices = NULL, selected = NULL),
            # checkboxInput("pls_probe", "Perform PLS-DA at probe level", value = FALSE),
            br(),
            h4("Gene contribution plot"),
            uiOutput("numGenesUI"),
            radioButtons("pls_ncomp", "Select component for gene contribution plot", choices = c(1,2)),
            br()

        ),
        mainPanel(
            tabsetPanel(type = "tabs",
            tabPanel("Discriminant Analysis", h4("Sparse Partial Least Squares Discriminant Analysis (sPLS-DA)  selects genes that are informative about a specific group. "), plotOutput("indPlot", height = 800)),
            
            tabPanel("Component loadings plot", h4("Gene contribution to each principle component."), p("The longer the bar (in either direction) the more that gene contributes to that component."), plotOutput("contribPlot", height = 800)),
            
            tabPanel("Circle variance", h4("Circle variance projections onto tissue."), p("Strongly correlated genes are projected in the same direction from the origin; the greater the distance the stronger the association."), plotOutput("varPlot", height = 800))
        ),
        position = c("right","left"),
        fluid = TRUE
        )
        )
    )
)
)