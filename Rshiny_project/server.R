library(shiny)
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(DT)

server <- function(input, output, session) {
  
  # Reactive value to store Seurat object
  scRNAseqData <- reactiveVal(NULL)
  
  # Validate if all required files are uploaded
  validate_files <- function() {
    req(input$barcodes, input$features, input$matrix)  # Ensure all files are uploaded
    return(TRUE)
  }
  
  observeEvent(input$runPipeline, {
    # Validate that all files are uploaded
    if (!validate_files()) {
      showNotification("Please upload all required files!", type = "error")
      return(NULL)
    }
    
    # Show loading message
    runjs("$('#loading').show();")
    
    # Read the files into the required formats
    barcodes <- read.delim(input$barcodes$datapath, header = FALSE)
    features <- read.delim(input$features$datapath, header = FALSE)
    matrix <- Matrix::readMM(input$matrix$datapath)
    
    # Set rownames and colnames for the matrix
    rownames(matrix) <- features$V1
    colnames(matrix) <- barcodes$V1
    
    # Create Seurat object
    pbmc <- CreateSeuratObject(counts = matrix, project = "PBMC3k", min.cells = 3, min.features = 200)
    
    # Calculate mitochondrial percentage
    pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
    
    # Visualize QC metrics
    output$qcPlot <- renderPlot({
      VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
    })
    
    # Filter cells
    pbmc <- subset(pbmc, subset = nFeature_RNA > input$minFeatures & nFeature_RNA < input$maxFeatures & percent.mt < input$maxMT)
    
    # Normalize data
    pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
    
    # Identify highly variable features
    pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
    
    # Scale data
    pbmc <- ScaleData(pbmc, vars.to.regress = c("percent.mt"))
    
    # Perform PCA
    pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
    
    # Visualize PCA (Optional)
    output$pcaPlot <- renderPlot({
      ElbowPlot(pbmc)
    })
    
    # Find neighbors
    pbmc <- FindNeighbors(pbmc, dims = 1:10)
    
    # Cluster cells
    pbmc <- FindClusters(pbmc, resolution = 0.5)
    
    # Run UMAP
    pbmc <- RunUMAP(pbmc, dims = 1:10)
    
    # Visualize clusters
    output$umapPlot <- renderPlot({
      DimPlot(pbmc, reduction = "umap", label = TRUE)
    })
    
    # Find markers for each cluster
    cluster.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
    
    # Display cluster markers
    output$markersTable <- DT::renderDT({
      datatable(cluster.markers)
    })
    
    # Save Seurat object into reactive value
    scRNAseqData(pbmc)
    
    # Hide loading message after completion
    runjs("$('#loading').hide();")
  })
}
