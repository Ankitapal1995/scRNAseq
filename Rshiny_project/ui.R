library(shiny)
library(shinyjs)  # For custom JavaScript actions

ui <- fluidPage(
  useShinyjs(),  # Enable shinyjs
  tags$head(tags$style(HTML("
    #loading {
      position: fixed;
      top: 50%;
      left: 50%;
      transform: translate(-50%, -50%);
      font-size: 24px;
      color: #333;
      background-color: rgba(255, 255, 255, 0.9);
      padding: 20px;
      border-radius: 10px;
      display: none;
      z-index: 9999;
    }
  "))),
  
  div(id = "loading", "Please wait, your job is running..."),
  
  titlePanel("Single-cell RNA-seq Analysis Pipeline"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("barcodes", "Choose Barcodes File (barcodes.tsv or barcodes.tsv.gz)", accept = c(".tsv", ".tsv.gz")),
      fileInput("features", "Choose Features File (genes.tsv or genes.tsv.gz)", accept = c(".tsv", ".tsv.gz")),
      fileInput("matrix", "Choose Matrix File (matrix.mtx or matrix.mtx.gz)", accept = c(".mtx", ".mtx.gz")),
      numericInput("minFeatures", "Minimum Features (genes per cell):", 200),
      numericInput("maxFeatures", "Maximum Features (genes per cell):", 2500),
      numericInput("maxMT", "Maximum Mitochondrial Percent:", 5),
      actionButton("runPipeline", "Run Pipeline")
    ),
    
    mainPanel(
      h4("QC Plots:"),
      plotOutput("qcPlot"),
      h4("UMAP Plot:"),
      plotOutput("umapPlot"),
      h4("Cluster Markers:"),
      DT::dataTableOutput("markersTable")
    )
  )
)
