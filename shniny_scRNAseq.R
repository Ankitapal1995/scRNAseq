library(shiny)
library(shinydashboard)
library(Seurat)
library(plotly)

ui <- dashboardPage(
  dashboardHeader(title = "scRNA-seq Analysis"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Data Input", tabName = "input", icon = icon("file-upload")),
      menuItem("Preprocessing", tabName = "preprocessing", icon = icon("cogs")),
      menuItem("Visualization", tabName = "visualization", icon = icon("chart-bar"))
    )
  ),
  dashboardBody(
    tabItems(
      tabItem(tabName = "input",
              fileInput("file", "Choose an input file", accept = c(".csv", ".rds")),
              actionButton("load_data", "Load Data")
      ),
      tabItem(tabName = "preprocessing",
              actionButton("qc", "Run Quality Control"),
              actionButton("run_pca", "Run PCA")
      ),
      tabItem(tabName = "visualization",
              plotlyOutput("umap_plot"),
              plotlyOutput("pca_plot")
      )
    )
  )
)

server <- function(input, output, session) {
  
  # Reactive expression to load the Seurat object
  seurat_obj <- reactiveVal(NULL)
  
  # Load data when the user clicks the 'Load Data' button
  observeEvent(input$load_data, {
    req(input$file)
    
    # Assuming input file is an RDS file containing a Seurat object
    seurat_data <- readRDS(input$file$datapath)
    seurat_obj(seurat_data)
    
    showNotification("Data Loaded Successfully!")
  })
  
  # Quality control (QC) processing
  observeEvent(input$qc, {
    req(seurat_obj())
    # Run basic QC steps (e.g., filter low-quality cells)
    pbmc <- seurat_obj()
    pbmc <- SCTransform(pbmc)  # Normalize and scale the data
    seurat_obj(pbmc)
    
    showNotification("Quality Control Complete!")
  })
  
  # Run PCA when the user clicks the 'Run PCA' button
  observeEvent(input$run_pca, {
    req(seurat_obj())
    pbmc <- seurat_obj()
    
    # Run PCA
    pbmc <- RunPCA(pbmc)
    seurat_obj(pbmc)
    
    showNotification("PCA Complete!")
  })
  
  # Plot UMAP
  output$umap_plot <- renderPlotly({
    req(seurat_obj())
    pbmc <- seurat_obj()
    
    # Run UMAP for visualization
    pbmc <- RunUMAP(pbmc, dims = 1:10)
    
    # Plot UMAP using Plotly
    umap_data <- Embeddings(pbmc, "umap")
    umap_plot <- plot_ly(x = umap_data[,1], y = umap_data[,2], type = 'scatter', mode = 'markers')
    umap_plot
  })
  
  # Plot PCA
  output$pca_plot <- renderPlotly({
    req(seurat_obj())
    pbmc <- seurat_obj()
    
    # Plot PCA using ggplot2 and then convert to plotly
    pca_data <- Embeddings(pbmc, "pca")
    pca_plot <- ggplot(as.data.frame(pca_data), aes(x = PC_1, y = PC_2)) + geom_point()
    ggplotly(pca_plot)
  })
}

shinyApp(ui, server)
